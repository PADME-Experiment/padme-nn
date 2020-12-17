#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSpectrum.h"
#include "TPDF.h"
#include <typeinfo>

#define NPOINTS 1024
#define BinsPerNs 1024./400.
#define NsPerBin 400./1024.
#define NoEventsToPrint 0//10*3//10 events times 3 analyses

char name[256];
Int_t fProcessing = 3;

void DigitalProcessingRRC(Double_t *uin, Double_t *uout,int npoints, Double_t timebin);

void GraphDrawing(TGraph *graph[], TGraph *procgraph[], Int_t NoEvents);

Int_t CalcChaTime(Double_t fSamples[], Int_t fProcessing, TH1D* hTSpec, TH1D* hPreTSpec, TCanvas* cTSpec[NoEventsToPrint], Int_t& TSpecNo, Int_t EventNo, std::vector<Double_t> *Time, Bool_t PrintFlag, Int_t& TSpecPrintNo);

int main(Int_t argc, char **argv){
  /*  TFile *readfile = TFile::Open("Data.tmp/GeneratedSignals.root","update");
  if (!readfile) {
    std::cout<<"No file named GeneratedSignals.root found"<<std::endl;
    return 0;
    }*/

  TFile *readfile = TFile::Open("Data.tmp/test/200913TestGeneratedSignalsLandauAmp200k.root","update");
  if (!readfile) {
    std::cout<<"No file named 200913TestGeneratedSignalsLandauAmp200k.root found in Data.tmp/test"<<std::endl;
    return 0;
  }
  
  TTree *readtree = (TTree*)readfile->Get("SigTree");
  //readfile->GetObject("SigTree",readtree);

  TFile *writefile = TFile::Open("Data.tmp/200913ReconstructedSignalsLandauAmp200k.root","recreate");
  if (!writefile) {
    std::cout<<"No file named 200913ReconstructedSignalsLandauAmp200k.root found in Data.tmp/test"<<std::endl;
    return 0;
  }

  TTree *writetree = new TTree("RecoTree","Reconstructed Signal Tree");
 
  Int_t TotEvents=readtree->GetEntries();

  Bool_t PrintFlag=1;
  
  Double_t sigfinal[1024];//Final signal - what the DAQ would register, so the sum of all individual hits
  Double_t  sigProc[NPOINTS];

  Int_t EventIndex;//Index of event
  Int_t HitsPerEvent;//Number of hits per event
  Int_t tbins[NPOINTS];//Time in bins
  
  std::vector<Int_t> *HitTimesBins = 0;//Hit arrival time in bins
  std::vector<Double_t> *HitAmp = 0;//Hit amplitude

  //read branches written in SignalGeneratorToy.cpp
  readtree->SetBranchAddress("SigFinal",sigfinal);
  readtree->SetBranchAddress("EventNo",&EventIndex);
  readtree->SetBranchAddress("HitsPerEvent",&HitsPerEvent);
  readtree->SetBranchAddress("HitTimeBins",&HitTimesBins);//Hit arrival time in bins
  readtree->SetBranchAddress("HitAmplitude",&HitAmp);
  readtree->SetBranchAddress("TBins",&tbins);//Runs from 1 to 1024

  std::vector<Double_t> *tUnProcRecoBins =new std::vector<Double_t>; //Reconstructed hit times from TSpectrum on raw data
  std::vector<Double_t> *tProcRecoBins= new std::vector<Double_t>; //Reconstructed hit times from TSpectrum on RRC processed data
  std::vector<Double_t> *tDerivRecoBins= new std::vector<Double_t>; //Reconstructed hit times from TSpectrum on derivative of raw data

  Int_t NoRecoUnProc=0;//Number of hits reconstructed by TSpectrum without digital processing
  Int_t NoRecoProc=0;//Number of hits reconstructed by TSpectrum from digitally processed signals
  Int_t NoRecoDeriv=0;//Number of hits reconstructed by TSpectrum from digitally processed signals

  //TBranch *b = readtree->GetBranch("TRecoUnProc");
  //readtree->GetListOfBranches()->Remove(b);
  // TLeaf* l= readtree->GetLeaf("TRecoUnProc");
  // readtree->GetListOfLeaves()->Remove((TObject *)l);
  // readtree->Write();
  
  //create new branch for reconstructed time
  writetree->Branch("TRecoUnProc","vector<Double_t>",&tUnProcRecoBins);
  writetree->Branch("TRecoProc","vector<Double_t>",&tProcRecoBins);
  writetree->Branch("TRecoDeriv","vector<Double_t>",&tDerivRecoBins);
  writetree->Branch("NoRecoUnProc",&NoRecoUnProc,"NoRecoUnProc/I");
  writetree->Branch("NoRecoProc",&NoRecoProc,"NoRecoProc/I");
  writetree->Branch("NoRecoDeriv",&NoRecoDeriv,"NoRecoDeriv/I");
  
  readtree->Print();
  std::cout<<"TotEvents = "<<TotEvents<<std::endl;

  Int_t TotalTSpecUnProcNo=0;//total number of hits in all event found by TSpectrum **without** RRC Processing
  Int_t TotalTSpecProcNo=0;//total number of hits in all event found by TSpectrum with RRC Processing
  Int_t TotalTSpecDeriv=0;//total number of hits in all event found by TSpectrum on the derivative of the signals
  
  TH1D *hRaw = new TH1D("hRaw","hRaw", 1024,0.,1024.);//Signal before digital processing
  TH1D *hTSpecRecoUnProc = new TH1D("hTSpecRecoUnProc","hTSpecRecoUnProc", 1024,0.,1024.);//Signal after TSpectrum (will have markers of TSpectrum-found peaks)
  TH1D *hTSpecRecoProc = new TH1D("hTSpecRecoProc","hTSpecRecoProc", 1024,0.,1024.);//Signal after TSpectrum && digital processing (will have markers of TSpectrum-found peaks)
  TH1D *hTSpecDeriv = new TH1D("hTSpecDeriv","hTSpecDeriv", 1024,0.,1024.);//Signal after TSpectrum && digital processing (will have markers of TSpectrum-found peaks)

  TCanvas *cTSpectrum[NoEventsToPrint];
 
  TGraph  *grapharr[NoEventsToPrint];
  TGraph  *procgrapharr[NoEventsToPrint];
  
  TCanvas *cFinal= new TCanvas(name,name,0,0,800,800);//TSpectrum results
  cFinal->Print("TSpectrumResults.pdf[");//opens the file but doesn't print anything

  Int_t printno=0;//counter for events to print
  Int_t TSpecPrintNo=0;

  std::vector<Int_t> TruthHitTimeBinDiffs;//Difference in hit arrival time in bins
  Int_t TruthMinHitTimeBinDiffs;//Minimum difference in hit arrival time for each event, in bins
  
  for(Int_t eventno=0;eventno<TotEvents;eventno++){//loop over total events
    readtree->GetEntry(eventno);
    DigitalProcessingRRC(sigfinal,sigProc,1024,0.4);
    if(printno<NoEventsToPrint){
      grapharr[printno] = new TGraph();
      if(fProcessing>=2) procgrapharr[printno] = new TGraph();
      for(int ii=0;ii<NPOINTS;ii++){
	grapharr[printno]->SetPoint(ii,tbins[ii],sigfinal[ii]);
	if(fProcessing>=2) procgrapharr[printno]->SetPoint(ii,tbins[ii],sigProc[ii]);
      }
      printno++;
    }//end of graphing

    //time vectors are cleared every event so that they only ever contain the no. hits found from a single event
    tProcRecoBins->clear();
    tUnProcRecoBins->clear();
    tDerivRecoBins->clear();

    TruthHitTimeBinDiffs.clear();
    for(int ii=0;ii<HitsPerEvent;ii++){
      if(ii<HitsPerEvent-1){
	TruthHitTimeBinDiffs.push_back(HitTimesBins->at(ii+1)-HitTimesBins->at(ii));
      }
    }
    TruthMinHitTimeBinDiffs=*min_element(TruthHitTimeBinDiffs.begin(), TruthHitTimeBinDiffs.end()); 
    //    if(printno<NoEventsToPrint&&TruthMinHitTimeBinDiffs<9*BinsPerNs&&HitsPerEvent>1){
    //std::cout<<"Event no "<<eventno<<" no hits "<<HitsPerEvent<<" min diff "<<TruthMinHitTimeBinDiffs<<std::endl;
    //PrintFlag=1;//I just want to print events with hits that are close together
    //}
    //    else PrintFlag=0;
    
    /*    NoRecoUnProc = CalcChaTime(sigfinal,1, hTSpecRecoUnProc, hRaw, cTSpectrum, TotalTSpecUnProcNo,eventno,tUnProcRecoBins,PrintFlag, TSpecPrintNo);
    std::sort(tUnProcRecoBins->begin(),tUnProcRecoBins->end());
    */
    
    if(fProcessing>=2){
      NoRecoProc = CalcChaTime(sigfinal,2, hTSpecRecoProc, hRaw, cTSpectrum, TotalTSpecProcNo,eventno,tProcRecoBins,PrintFlag, TSpecPrintNo);
      std::sort(tProcRecoBins->begin(),tProcRecoBins->end());
      }

    if(fProcessing>=3){
      NoRecoDeriv = CalcChaTime(sigfinal,3, hTSpecDeriv, hRaw, cTSpectrum, TotalTSpecDeriv,eventno,tDerivRecoBins,PrintFlag, TSpecPrintNo); //array of amps for TSpec, processing type (develop, RRC, deriv), histogram for TSpectrum, histogram for before tspectrum, canvas for tspectrum, no. hits found by tspectrum, EventNo, vector of tspec times
      std::sort(tDerivRecoBins->begin(),tDerivRecoBins->end());
    }

    if(fProcessing==1&&eventno%100==0) std::cout<<"Event "<<eventno<<" HitsPerEvent "<<HitsPerEvent<<" NoRecoUnProc "<<NoRecoUnProc<<std::endl;
    if(fProcessing==2&&eventno%100==0) std::cout<<"Event "<<eventno<<" HitsPerEvent "<<HitsPerEvent<<" NoRecoUnProc "<<NoRecoUnProc<<" NoRecoProc "<<NoRecoProc<<std::endl;
    if(fProcessing==3&&eventno%100==0) std::cout<<"Event "<<eventno<<" HitsPerEvent "<<HitsPerEvent<<" NoRecoUnProc "<<NoRecoUnProc<<" NoRecoProc "<<NoRecoProc<<" NoRecoDeriv "<<NoRecoDeriv<<std::endl;

    writetree->Fill();
  }//end of event loop
  
  cFinal->Print("TSpectrumResults.pdf]");//closes the file
  GraphDrawing(grapharr,procgrapharr,printno);
  
  writetree->Print();
  writetree->Write(nullptr,TObject::kOverwrite);
  delete readfile;
  delete writefile;
}//end of main

void DigitalProcessingRRC(Double_t *uin, Double_t *uout,int npoints, Double_t timebin) { //Beth, implemented from Venelin's idea 06/2019
  // Double_t R1=1300;//ohms
  // Double_t R2=100.; //ohms
  // Double_t C=0.015e-9; //nF
  
  //simulating RRC circuit. Circuit diagram in presentation from Beth 19/06/2019 
  //approximating voltage output to:
  //dQ/dt=(uin(t)/R2)-((R1+R2)/CR1R2)*Q(t)
  //uout(t)=uin(t)-Q(t)/C


  Double_t R1=650;//ohms
  Double_t R2=100; //ohms
  Double_t C=0.03; //nF

  
  //Calculating the output pulse:
  Double_t integr=0;
  
  static Double_t ic[1024];

  Double_t bin_width=timebin;

  integr=0;
  ic[0]= uin[0]/R2;
  
  for(int i=1;i<npoints;i++) {
    integr+=ic[i-1]*bin_width; //integr = intgrated charge = charge of this bin + charge of previous bin + bin before...
    ic[i]= uin[i]/R2 - ((R1+R2)/(C*R1*R2))* integr; //ic = current through capacitor = dQ/dt
    uout[i] = uin[i] - integr/(C);
  }
}

void GraphDrawing(TGraph *graph[], TGraph *procgraph[], Int_t NoEvents){
  TCanvas *canvasarr[NoEventsToPrint];
  TMultiGraph *gfinal;
  for(int ii=0;ii<NoEvents;ii++){
    gfinal = new TMultiGraph(); 

    sprintf(name, "cToy%d", ii);
    canvasarr[ii]=new TCanvas(name,name,0,0,800,800);
    canvasarr[ii]->cd();

    if(fProcessing>=2) procgraph[ii]->SetLineColor(kBlue);
    gfinal->Add(graph[ii],"l");
    if(fProcessing>=2) gfinal->Add(procgraph[ii],"l");
    gfinal->Draw("same");
    gStyle->SetOptStat(0);
    if(ii<NoEvents-1) canvasarr[ii]->Print("SampleSignals.pdf(");
    else canvasarr[ii]->Print("SampleSignals.pdf)");
  }
}

Int_t CalcChaTime(Double_t fSamples[], Int_t fProcessing, TH1D* hTSpec, TH1D* hPreTSpec, TCanvas* cTSpec[NoEventsToPrint], Int_t& TSpecNo, Int_t EventNo, std::vector<Double_t> *Time, Bool_t PrintFlag, Int_t& TSpecPrintNo){
  hTSpec->Reset();
  hPreTSpec->Reset();

  if(fProcessing==1)  sprintf(name,"hRaw%d",EventNo);
  if(fProcessing==2)  sprintf(name,"hRaw->RRCPRocessing%d",EventNo);
  if(fProcessing==3)  sprintf(name,"hRaw->Deriv%d",EventNo);
  hPreTSpec->SetNameTitle(name,name);
  
  Double_t SigTau;//decay time of signal in ns
  Double_t fPeakSearchWidth;
  Double_t ThresholdFraction;
  if(fProcessing!=3)  ThresholdFraction=1./5;
  
  Int_t derivpoints = 10;//Number of bins before over which to take the derivative
  /*if(fProcessing==1){
    SigTau=20;
    fPeakSearchWidth=SigTau*log(1./ThresholdFraction)*BinsPerNs;//fixes sigma so that TSpectrum won't reconstruct one hit as being more than one
    }*/

  //  fPeakSearchWidth=0.3*BinsPerNs;
  if(fProcessing==3) fPeakSearchWidth=5;//6;//0.3;// 5 = working new toy//*BinsPerNs;
  if(fProcessing==2)  fPeakSearchWidth=13;
  if(fProcessing==1) fPeakSearchWidth=6;//6;//6 = develop, 20 next guess
  
  Double_t fTimeBin=0.4;
  Int_t npeaks = 10;
  Int_t nfound=0;
 
  //currently looking for peaks with TSpectrum to obtain multi hit times
  Double_t AbsSamRec[1024];
  Double_t AbsSamRecDP[1024];
  Double_t AbsSamRecDeriv[1024];
  
  for(UShort_t i=0;i<1024;i++){
    AbsSamRec[i] = (fSamples[i]);///4096.*1000.; //in mV positivo
    if(i>=derivpoints)    AbsSamRecDeriv[i]=(AbsSamRec[i]-AbsSamRec[i-derivpoints]);
    else if(i<derivpoints) AbsSamRecDeriv[i]=0;
    if(fProcessing==3)
      if(AbsSamRecDeriv[i]<=0) AbsSamRecDeriv[i] = 0;//to prevent TSpectrum being confused by negative values, set negative values of the derivative to 0
      else if (AbsSamRecDeriv[i]>0) AbsSamRecDeriv[i] = AbsSamRecDeriv[i];
  }

  if(fProcessing==2){
    DigitalProcessingRRC(AbsSamRec,AbsSamRecDP,1024,fTimeBin);
  }
  for(int i=0;i<1025;i++){
    if(i<1000){
      if(fProcessing==1) hTSpec->SetBinContent(i,AbsSamRec[i]);
      else if(fProcessing==2) hTSpec->SetBinContent(i,AbsSamRecDP[i]);
      else if(fProcessing==3) hTSpec->SetBinContent(i,AbsSamRecDeriv[i]);
    }
    else hTSpec->SetBinContent(i,0);
  }
  
  for(int i=0;i<1025;i++){
    if(i<1000) hPreTSpec->SetBinContent(i,AbsSamRec[i]);
    else hPreTSpec->SetBinContent(i,0);
  }

  //cTSpec->cd();
  Double_t VMax = hTSpec->GetMaximum(); //mV
  Double_t VMin = hTSpec->GetMinimum(); //mV

  Double_t VMax2 = hPreTSpec->GetMaximum(); //mV
  Double_t VMin2 = hPreTSpec->GetMinimum(); //mV

  if(fProcessing==3) ThresholdFraction=3./VMax;
  
  TSpectrum *s = new TSpectrum(2*npeaks);//npeaks max number of peaks

  //  if(fProcessing==1) ThresholdFraction=VMax/10;

  if(VMax>3) //any signal with amplitude < 3 is considered noise (given that the noise ~U[-1.5,1.5])
    nfound = s->Search(hTSpec,fPeakSearchWidth,"nobackground, nodraw",ThresholdFraction);
  
  //std::cout<<"TSpectrum finds "<<nfound<<" peaks"<<std::endl;
  if(TSpecPrintNo<NoEventsToPrint&&PrintFlag==1){
    if (fProcessing==1)  sprintf(name,"cTSpectrum%d",EventNo);
    else if (fProcessing==2)  sprintf(name,"cTSpectrumProc%d",EventNo);
    else if (fProcessing==3)  sprintf(name,"cTSpectrumDeriv%d",EventNo);
    cTSpec[EventNo]= new TCanvas(name,name,0,0,800,800);//TSpectrum results
    cTSpec[EventNo]->cd();
    hPreTSpec->Draw();
    hTSpec->Draw("same");
    hPreTSpec->SetLineColor(kBlack);
    TLegend* legend;
    if(fProcessing==2) legend = new TLegend(0.5,0.65,0.9,0.75);//(x1,y1,x2,y2)
    else if(fProcessing==3) legend = new TLegend(0.55,0.8,0.9,0.9);
    else if(fProcessing>3) legend = new TLegend(0.55,0.65,0.9,0.75);
    gStyle->SetLegendTextSize(0.03);
    gStyle->SetOptStat(0);
    if(fProcessing>=2) legend->AddEntry(hPreTSpec,"Raw data");
    if(fProcessing==2) legend->AddEntry(hTSpec,"RRC Processed Data");
    if(fProcessing==3) legend->AddEntry(hTSpec,"Derivative of Data");
    legend->Draw("same");
    cTSpec[EventNo]->Print("TSpectrumResults.pdf");
    TSpecPrintNo++;
  }

  for(Int_t ll=0;ll<nfound;ll++){ //peak loop per channel
    TSpecNo++;
    Double_t xp = (s->GetPositionX())[ll];
    Double_t yp = (s->GetPositionY())[ll];
    //    std::cout<<"TSpectrum peak "<<ll<<" xp "<<xp<<" yp "<<yp<<std::endl;
    Time->push_back(xp);
  }

  //  std::cout<<"nfound "<<nfound<<" Time->size() "<<Time->size()<<std::endl;
  
  return Int_t(Time->size());
}
