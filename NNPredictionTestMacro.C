#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRatioPlot.h"
#include "TLine.h"
 
#define NPOINTS 1024
#define NsPerBin 400./1024.
#define BinsPerNs 1024./400.
#define MaxTimeWindow 25//00. //hits are distributed uniformly within a 200ns window
 
char name[256];

TFile* OpenTruthFile(Int_t RunTypeFlag);
void OpenNNFile(Int_t RunTypeFlag, ifstream& NNFile);

void NNPredictionTestMacro(){

  Int_t RunTypeFlag=14;//0 = General Test (NoisyTestHits), 1 = Well Separated Signals, 2 = Exactly Two General Signals, 3 = Exactly Two Well Separated Signals, 4 = Exactly 3 General Signals, 5 - Exactly 3 Well Separated Signals, 6 = One or Two General Signals, 7 = Landau Distributed Standard Signals, 8 = Signals with 5ns Rise Time, 9 = Signals with 10ns Rise Time, 10 = Signals with 7ns Rise Time && Fixed 50mV Amplitude, 11 = Signals with mean 5 hits per event, 12 = Signals with mean 3 hits per event analysed with 5HitPerEventModel, 13 = Signals with exponential rise and fall, mean 3 hits per event, 14 = standard signals mutliplied by uniformly distributed random "calibration constant"~U[1,3] to represent the different gains of differen PVeto channels

  Bool_t UseDerivFlag=0;
  
  std::cout<<"RunType is "<<RunTypeFlag<<std::endl;
  
  TFile *truthfile;
  truthfile=OpenTruthFile(RunTypeFlag);
  if(truthfile==nullptr) return;
  
  TTree *tree;
  truthfile->GetObject("SigTree",tree);

  //Open NN prediction .txt file
  std::ifstream NNFile;
  OpenNNFile(RunTypeFlag,NNFile);
  
  if(!NNFile) return;
  
  float value;
  std::vector<Int_t> NNHitNumbers;
  while(!NNFile.eof()){
    NNFile>>value;
    NNHitNumbers.push_back(value);
  }

  std::ifstream StandardNNFile;
  /* if(RunTypeFlag==14){
    StandardNNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/RandomCalibration/RandomCalibrationResults/PVetoMCNNTestNoisyHitsWithoutRandomCalibration.txt");//Standard/StandardResults/PVetoMCNNTestNoisyHitsLandauAmp200k.txt");
    if(!StandardNNFile){
      std::cout<<"File /Users/bethlong/Google Drive/PVetoNN/PVetoMCData/Standard/StandardResults/PVetoMCNNTestNoisyHitsLandauAmp200k.txt does not exist."<<std::endl;
      return;
    }
    }*/
  
  float StandardValue;
  std::vector<Int_t> StandardNNHitNumbers;
  while(!StandardNNFile.eof()){
    StandardNNFile>>StandardValue;
    StandardNNHitNumbers.push_back(StandardValue);
    }
  
  //Open derivative reconstruction .root file
  if(UseDerivFlag==1){
    TFile *derivfile = TFile::Open("Data/200913ReconstructedSignalsLandauAmp200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named 200913ReconstructedSignalsLandauAmp200k.root found in test directory"<<std::endl;
      return;
    }
    tree->AddFriend("RecoTree",derivfile);
  }  
  
  //Check if NN & truth files are compatible in size
  if(NNHitNumbers.size()!=tree->GetEntries()+1){//When the NNFile.txt is being written, there's an extra blank line added to the end, therefore it should have TotEvents+1 lines
    std::cout<<"NNFile and truthfile are mismatched"<<std::endl;
    std::cout<<"NNHitNumbers.size() = "<<NNHitNumbers.size()<<" tree->GetEntries() = "<<tree->GetEntries()<<std::endl;
    return;
  }
  
  Long64_t TotEvents=tree->GetEntries();
  std::cout<<"TotEvents: "<<TotEvents/1000<<"k"<<std::endl;
  Int_t EventNumbers;//Index of event

  Int_t TruthHitsPerEvent;//Number of hits per event in MC "truth" data
  std::vector<Double_t> *TruthHitTimes=new std::vector<Double_t>;//Hit arrival time in ns in truth file

  Int_t DerivHitsPerEvent;//Number of hits per event in derivative reconstruction data
  Int_t ProcHitsPerEvent;//Number of hits per event in RRC reconstruction data

  //Read branch info from .root files into vectors
  tree->SetBranchAddress("HitsPerEvent",&TruthHitsPerEvent);  
  tree->SetBranchAddress("HitTimeNs",&TruthHitTimes);
  tree->SetBranchAddress("NoRecoProc",&ProcHitsPerEvent);
  if(UseDerivFlag==1) tree->SetBranchAddress("NoRecoDeriv",&DerivHitsPerEvent);  


  std::vector<Double_t> TruthHitTimeDiffs;//Difference in hit arrival time in ns
  Double_t TruthMinHitTimeDiffs;//Minimum difference in hit arrival time for each event, in ns

  //Set up analysis histograms
  TH1F *hProcHitNumbers = new TH1F("hProcHitNumbers","Hits per event",15,0,14);//number of hits per event from neural network
  TH1F *hDerivHitNumbers;
  if(UseDerivFlag==1)  hDerivHitNumbers = new TH1F("hDerivHitNumbers","Hits per event",15,0,14);//number of hits per event from neural network
  TH1F *hNNHitNumbers = new TH1F("hNNHitNumbers","Hits per event",15,0,14);//number of hits per event from neural network
  TH1F *hTrueHitNumbers = new TH1F("hTrueHitNumbers","Hits per event",15,0,14);//number of hits per event from simulation
  TH1F *hStandardNNHitNumbers;
  //  if(RunTypeFlag==14) hStandardNNHitNumbers = new TH1F("hStandardNNHitNumbers","Hits per event",15,0,14);//number of hits per event from neural network in standard case (as opposed to case with random calibration constant)
  
  TH1F *hProcCorrectEventsVsTrueTime = new TH1F("hProcCorrectEventsVsTrueTime","Minimum time difference between hits per correctly reconstructed event",MaxTimeWindow,0,MaxTimeWindow);//Number of events where derivative reconstructs correct no. hits
  TH1F *hDerivCorrectEventsVsTrueTime;
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime = new TH1F("hDerivCorrectEventsVsTrueTime","Minimum time difference between hits per correctly reconstructed event",MaxTimeWindow,0,MaxTimeWindow);//Number of events where derivative reconstructs correct no. hits
  TH1F *hNNCorrectEventsVsTrueTime = new TH1F("hNNCorrectEventsVsTrueTime","Minimum time difference between hits per correctly reconstructed event",MaxTimeWindow,0,MaxTimeWindow);//Number of events where NN reconstructs correct no. hits
  TH1F *hTrueEventsVsTrueTime = new TH1F("hTrueHitsVsTrueTime","Minimum time difference between hits per event",MaxTimeWindow,0,MaxTimeWindow);//Time spectrum of simulated events
  TH1F *hStandardNNCorrectEventsVsTrueTime;
  //  if(RunTypeFlag==14) hStandardNNCorrectEventsVsTrueTime = new TH1F("hStandardNNCorrectEventsVsTrueTime","Minimum time differencqe between hits per correctly reconstructed event",MaxTimeWindow,0,MaxTimeWindow);//Number of events where NN reconstructs correct no. hits
  
  Int_t NNWrongNo=0;
  Int_t StandardNNWrongNo=0;
  Int_t ProcWrongNo=0;
  Int_t DerivWrongNo=0;
  
  for(Int_t EventIndex=0;EventIndex<TotEvents;EventIndex++){
    TruthHitTimeDiffs.clear();
    tree->GetEntry(EventIndex);
    if(EventIndex%100==0)  std::cout<<"Event: "<<EventIndex<<" Truth hits: "<<TruthHitsPerEvent<<" "<<" NN Predicted hits: "<<NNHitNumbers[EventIndex]<<" RRC Predicted hits: "<<ProcHitsPerEvent<<" Derivative Predicted hits: "<<DerivHitsPerEvent<<std::endl;
    
    hProcHitNumbers->Fill(ProcHitsPerEvent);
    if(UseDerivFlag==1)    hDerivHitNumbers->Fill(DerivHitsPerEvent);
    hNNHitNumbers->Fill(NNHitNumbers[EventIndex]);
    hTrueHitNumbers->Fill(TruthHitsPerEvent);
    //    if(RunTypeFlag==14) hStandardNNHitNumbers->Fill(StandardNNHitNumbers[EventIndex]);
    
    if(TruthHitsPerEvent-NNHitNumbers[EventIndex]!=0) NNWrongNo++;
    //    if(RunTypeFlag==14&&TruthHitsPerEvent-StandardNNHitNumbers[EventIndex]!=0) StandardNNWrongNo++;
    if(TruthHitsPerEvent-ProcHitsPerEvent!=0) ProcWrongNo++;
    if(UseDerivFlag==1&&TruthHitsPerEvent-DerivHitsPerEvent!=0) DerivWrongNo++;

    if(TruthHitsPerEvent<2) continue;//if there are fewer than two hits in an event, there aren't multiple hits for which to evaluate the time difference
    for(int ii=0;ii<TruthHitsPerEvent;ii++){
      if(ii<TruthHitsPerEvent-1) TruthHitTimeDiffs.push_back(TruthHitTimes->at(ii+1)-TruthHitTimes->at(ii));
      TruthMinHitTimeDiffs=*min_element(TruthHitTimeDiffs.begin(), TruthHitTimeDiffs.end());
    }
    hTrueEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
    if(TruthHitsPerEvent==NNHitNumbers[EventIndex])   hNNCorrectEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
    if(TruthHitsPerEvent==ProcHitsPerEvent)   hProcCorrectEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
    if(UseDerivFlag==1&&TruthHitsPerEvent==DerivHitsPerEvent)   hDerivCorrectEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
    //    if(RunTypeFlag==14&&TruthHitsPerEvent==StandardNNHitNumbers[EventIndex])   hStandardNNCorrectEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
  }
  
  TLegend* NNLegend = new TLegend(0.6,0.8,0.9,0.9);
  //gStyle->SetLegendTextSize(0.03);
  NNLegend->AddEntry(hTrueHitNumbers,"ToyMC Truth");
  if(RunTypeFlag==14) NNLegend->AddEntry(hNNHitNumbers,"Neural Network with Random Multiplication");
  else NNLegend->AddEntry(hNNHitNumbers,"Neural Network");
  

  TLegend* StandardNNLegend;
  /*  if(RunTypeFlag==14){
    StandardNNLegend = new TLegend(0.6,0.8,0.9,0.9);
    //gStyle->SetLegendTextSize(0.03);
    StandardNNLegend->AddEntry(hTrueHitNumbers,"ToyMC Truth");
    StandardNNLegend->AddEntry(hStandardNNHitNumbers,"Neural Network without Random Multiplication");
    }*/
  
  TLegend* ProcLegend = new TLegend(0.6,0.8,0.9,0.9);
  //gStyle->SetLegendTextSize(0.03);
  ProcLegend->AddEntry(hTrueHitNumbers,"ToyMC Truth");
  ProcLegend->AddEntry(hProcHitNumbers,"RRC + TSpectrum");

  TLegend* DerivLegend;
  if(UseDerivFlag==1){
    DerivLegend = new TLegend(0.6,0.8,0.9,0.9);
    //gStyle->SetLegendTextSize(0.03);
    DerivLegend->AddEntry(hTrueHitNumbers,"ToyMC Truth");
    DerivLegend->AddEntry(hDerivHitNumbers,"Derivative + TSpectrum");
  }

  TLegend* AllLegend = new TLegend(0.5,0.75,0.9,0.9);
  //  gStyle->SetLegendTextSize(0.03);
  AllLegend->AddEntry(hTrueHitNumbers,"ToyMC Truth");
  //  AllLegend->AddEntry(hProcHitNumbers,"RRC + TSpectrum");
  if(UseDerivFlag==1)  AllLegend->AddEntry(hDerivHitNumbers,"Derivative + TSpectrum");
  if(RunTypeFlag==14){
    //    AllLegend->AddEntry(hStandardNNHitNumbers,"Neural Network without Multiplication");
    AllLegend->AddEntry(hNNHitNumbers,"Neural Network with Random Multiplication");
  }
  else   AllLegend->AddEntry(hNNHitNumbers,"Neural Network");

  TLegend* RatioLegend = new TLegend(0.6,0.8,0.9,0.9);
  if(UseDerivFlag==1)  RatioLegend->AddEntry(hDerivHitNumbers,"Derivative + TSpectrum");
  if(RunTypeFlag==14){
    //    RatioLegend->AddEntry(hStandardNNHitNumbers,"Neural Network without Multiplication");
    RatioLegend->AddEntry(hNNHitNumbers,"Neural Network with Random Multiplication");
  }
  else   RatioLegend->AddEntry(hNNHitNumbers,"Neural Network");
  
  TCanvas *cNNHitNumbers = new TCanvas("cNNHitNumbers","Truth vs NN predicted hits per event",0,0,800,800);
  TRatioPlot *NNEfficiencyRatio = new TRatioPlot(hNNHitNumbers,hTrueHitNumbers);
  NNEfficiencyRatio->Draw();
  NNEfficiencyRatio->GetUpperPad()->cd();
  hTrueHitNumbers->Draw("AH");
  hNNHitNumbers->Draw("sames");
  NNLegend->Draw("sames");
  hNNHitNumbers->SetLineColor(kRed);
  TGraph *NNHitRatio = NNEfficiencyRatio->GetLowerRefGraph();

  TCanvas *cStandardNNHitNumbers;
  TRatioPlot *StandardNNEfficiencyRatio;
  TGraph *StandardNNHitRatio;
  
  /*if(RunTypeFlag==14){
    cStandardNNHitNumbers = new TCanvas("cStandardNNHitNumbers","Truth vs NN predicted hits per event without Multiplication",0,0,800,800);
    StandardNNEfficiencyRatio = new TRatioPlot(hStandardNNHitNumbers,hTrueHitNumbers);
    StandardNNEfficiencyRatio->Draw();
    
    StandardNNEfficiencyRatio->GetUpperPad()->cd();
    hStandardNNHitNumbers->Draw("AH");
    hTrueHitNumbers->Draw("sames");
    StandardNNLegend->Draw("sames");
    hStandardNNHitNumbers->SetLineColor(kBlack);
    hStandardNNHitNumbers->SetLineColor(kBlack);
    StandardNNHitRatio = StandardNNEfficiencyRatio->GetLowerRefGraph();
    }*/

  /*TCanvas *cProcHitNumbers = new TCanvas("cProcHitNumbers","Truth vs derivative predicted hits per event",0,0,800,800);
  TRatioPlot *ProcEfficiencyRatio = new TRatioPlot(hProcHitNumbers,hTrueHitNumbers);
  ProcEfficiencyRatio->Draw();

  ProcEfficiencyRatio->GetUpperPad()->cd();
  hProcHitNumbers->Draw("AH");
  hTrueHitNumbers->Draw("sames");
  ProcLegend->Draw("sames");
  
  ProcEfficiencyRatio->GetLowerPad()->cd();
  hProcHitNumbers->SetLineColor(kRed);
  */
  TGraph *DerivHitRatio;
  if(UseDerivFlag==1){
    TCanvas *cDerivHitNumbers = new TCanvas("cDerivHitNumbers","Truth vs derivative predicted hits per event",0,0,800,800);
    TRatioPlot *DerivEfficiencyRatio = new TRatioPlot(hDerivHitNumbers,hTrueHitNumbers);
    DerivEfficiencyRatio->Draw();

    DerivEfficiencyRatio->GetUpperPad()->cd();
    hDerivHitNumbers->Draw("AH");
    hTrueHitNumbers->Draw("sames");
    DerivLegend->Draw("sames");

    DerivHitRatio = DerivEfficiencyRatio->GetLowerRefGraph();
  }

  TCanvas *cHitRatioCanvas = new TCanvas("cHitRatioCanvas","Deriv/NN efficiency plots",0,0,800,800);
  NNHitRatio->Draw("AP");
  if(UseDerivFlag==1)  DerivHitRatio->Draw("Psame");
  //  if(RunTypeFlag==14) StandardNNHitRatio->Draw("Psame");
  NNHitRatio->SetLineColor(kRed);
  RatioLegend->Draw("same");
  if(UseDerivFlag==1)  DerivHitRatio->GetXaxis()->SetRangeUser(0,12);
  NNHitRatio->GetXaxis()->SetRangeUser(0,12);
  NNHitRatio->GetYaxis()->SetRangeUser(0,1.3);
  //  NNHitRatio->SetTitle("No. Reconstructed Hits/No. Truth Hits vs No. Truth Hits");
  
  TCanvas *cNNCorrectEventTimes = new TCanvas("cNNCorrectEventTimes","Truth/NN predicted events per time diff for correct reconstruction",0,0,800,800);
  TRatioPlot *NNCorrectTimeRatio = new TRatioPlot(hNNCorrectEventsVsTrueTime,hTrueEventsVsTrueTime);
  NNCorrectTimeRatio->Draw();
  NNCorrectTimeRatio->GetUpperPad()->cd();
  hTrueEventsVsTrueTime->Draw("AH");
  hNNCorrectEventsVsTrueTime->Draw("sames");
  NNLegend->Draw("sames");
  hNNCorrectEventsVsTrueTime->SetLineColor(kRed);
  NNCorrectTimeRatio->GetLowerPad()->cd();
  TLine *line = new TLine(0,0.95,MaxTimeWindow,0.95);
  line->SetLineColor(kRed);
  line->Draw("same");
  TGraph *NNTimeRatio = NNCorrectTimeRatio->GetLowerRefGraph();
  
  TCanvas *cStandardNNCorrectEventTimes;
  TRatioPlot *StandardNNCorrectTimeRatio;
  TGraph *StandardNNTimeRatio;
  /*  if(RunTypeFlag==14){
    cStandardNNCorrectEventTimes = new TCanvas("cStandardNNCorrectEventTimes","Truth/NN predicted events per time diff for correct reconstruction",0,0,800,800);
    StandardNNCorrectTimeRatio = new TRatioPlot(hStandardNNCorrectEventsVsTrueTime,hTrueEventsVsTrueTime);
    StandardNNCorrectTimeRatio->Draw();
    StandardNNCorrectTimeRatio->GetUpperPad()->cd();
    hTrueEventsVsTrueTime->Draw("AH");
    hStandardNNCorrectEventsVsTrueTime->Draw("sames");
    StandardNNLegend->Draw("sames");
    hStandardNNCorrectEventsVsTrueTime->SetLineColor(kBlack);
    StandardNNCorrectTimeRatio->GetLowerPad()->cd();
    line->SetLineColor(kRed);
    line->Draw("same");
    StandardNNTimeRatio = StandardNNCorrectTimeRatio->GetLowerRefGraph();
    }*/

  /*
  TCanvas *cProcCorrectEventTimes = new TCanvas("cProcCorrectEventTimes","Truth/RRC predicted events per time diff for correct reconstruction",0,0,800,800);
  TRatioPlot *ProcCorrectTimeRatio = new TRatioPlot(hProcCorrectEventsVsTrueTime,hTrueEventsVsTrueTime);
  ProcCorrectTimeRatio->Draw();

  ProcCorrectTimeRatio->GetUpperPad()->cd();
  hTrueEventsVsTrueTime->Draw("AH");
  hProcCorrectEventsVsTrueTime->Draw("sames");
  ProcLegend->Draw();//"sames");
  
  ProcCorrectTimeRatio->GetLowerPad()->cd();
  line->Draw("same");
  */
  
  TGraph *DerivTimeRatio;
  if(UseDerivFlag==1){
    TCanvas *cDerivCorrectEventTimes = new TCanvas("cDerivCorrectEventTimes","Truth/Derivative predicted events per time diff for correct reconstruction",0,0,800,800);
    TRatioPlot *DerivCorrectTimeRatio = new TRatioPlot(hDerivCorrectEventsVsTrueTime,hTrueEventsVsTrueTime);
    DerivCorrectTimeRatio->Draw();
    
    DerivCorrectTimeRatio->GetUpperPad()->cd();
    hTrueEventsVsTrueTime->Draw("AH");
    hDerivCorrectEventsVsTrueTime->Draw("sames");
    DerivLegend->Draw();//"sames");
    
    DerivCorrectTimeRatio->GetLowerPad()->cd();
    line->Draw("same"); 

    DerivTimeRatio = DerivCorrectTimeRatio->GetLowerRefGraph();
  }
  
  TCanvas *cTimeRatioCanvas = new TCanvas("cTimeRatioCanvas","Deriv/NN efficiency plots",0,0,800,800);
  NNTimeRatio->Draw("AP");
  //  if(RunTypeFlag==14)  StandardNNTimeRatio->Draw("Psame");
  //  StandardNNTimeRatio->Draw("AP");
  if(UseDerivFlag==1)  DerivTimeRatio->Draw("Psame");
  NNTimeRatio->SetLineColor(kRed);
  RatioLegend->Draw("same");
  NNTimeRatio->GetXaxis()->SetRangeUser(0,MaxTimeWindow);
  NNTimeRatio->GetYaxis()->SetRangeUser(0,1);
  if(UseDerivFlag==1)  DerivTimeRatio->GetYaxis()->SetRangeUser(0,1);
  //  NNTimeRatio->SetTitle("No. Reconstructed Hits/No. Truth Hits vs Min Time Diff");
  line->Draw();
  line->SetLineColor(kBlue);

  TCanvas *cStandardTimeRatioCanvas;
  /*  if(RunTypeFlag==14){
    cStandardNNHitNumbers = new TCanvas("cStandardTimeRatioCanvas","Deriv/NN efficiency plots",0,0,800,800);
    StandardNNTimeRatio->Draw("AP");
    //  StandardNNTimeRatio->Draw("AP");
    if(UseDerivFlag==1)  DerivTimeRatio->Draw("Psame");
    StandardNNTimeRatio->SetLineColor(kBlack);
    RatioLegend->Draw("same");
    StandardNNTimeRatio->GetXaxis()->SetRangeUser(0,MaxTimeWindow);
    StandardNNTimeRatio->GetYaxis()->SetRangeUser(0,1);
    if(UseDerivFlag==1)  DerivTimeRatio->GetYaxis()->SetRangeUser(0,1);
    //  NNTimeRatio->SetTitle("No. Reconstructed Hits/No. Truth Hits vs Min Time Diff");
    line->Draw();
    line->SetLineColor(kBlue);
    }*/
  TCanvas *cAllHitNumbers = new TCanvas("cAllHitNumbers","Truth vs RRC vs derivative vs CNN predicted hits per event",0,0,800,800);
  //  hProcHitNumbers->Draw();
  if(UseDerivFlag==1)  {
    hDerivHitNumbers->Draw();//"sames");
    hTrueHitNumbers->Draw("sames");
  }
  else hTrueHitNumbers->Draw();
  AllLegend->Draw("sames");
  hNNHitNumbers->Draw("sames");
  /*  if(RunTypeFlag==14){
    hStandardNNHitNumbers->Draw("sames");
    hStandardNNHitNumbers->SetFillColorAlpha(kBlack,0.35);
    }*/
  //  hProcHitNumbers->SetLineColor(kMagenta);
  if(UseDerivFlag==1)  hDerivHitNumbers->SetLineColor(kBlack);
  hTrueHitNumbers->SetFillColorAlpha(kBlue,0.25);
  if(UseDerivFlag==1)  hDerivHitNumbers->SetFillColorAlpha(kBlack,0.35);
  hNNHitNumbers->SetFillColorAlpha(kRed,0.25);

  TCanvas *cAllCorrectEventTimes = new TCanvas("cAllCorrectEventTimes","Truth/RRC/NN/Derivative predicted events per time diff for correct reconstruction",0,0,800,800);
  hTrueEventsVsTrueTime->Draw();
  hNNCorrectEventsVsTrueTime->Draw("sames");
  //  hProcCorrectEventsVsTrueTime->Draw("sames");
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime->Draw("sames");
  AllLegend->Draw("sames");
  //  hProcCorrectEventsVsTrueTime->SetLineColor(kMagenta);
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime->SetLineColor(kBlack);
  hTrueEventsVsTrueTime->SetFillColorAlpha(kBlue,0.25);
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime->SetFillColorAlpha(kBlack,0.35);
  hNNCorrectEventsVsTrueTime->SetFillColorAlpha(kRed,0.25);
  /*  if(RunTypeFlag==14){
    hStandardNNCorrectEventsVsTrueTime->SetFillColorAlpha(kBlack,0.35);
    hStandardNNCorrectEventsVsTrueTime->Draw("sames");
    }*/
  gStyle->SetOptStat(0);
  
  std::cout<<"Number of wrongly reconstructed events from NN: "<<NNWrongNo<<" out of "<<TotEvents<<" from derivative "<<DerivWrongNo<<" out of "<<TotEvents<<std::endl;
  std::cout<<"NN Efficiency: "<<(TotEvents-NNWrongNo)*100/(1.*TotEvents)<<"%, Derivative Efficiency: "<<(TotEvents-DerivWrongNo)*100/(1.*TotEvents)<<"%"<<std::endl;

}

TFile* OpenTruthFile(Int_t RunTypeFlag){
  //Open truth .root file
  if(RunTypeFlag==0){
    TFile* truthfile = TFile::Open("Data/test/TestGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named TestGeneratedSignals.root found in directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
  
  else if(RunTypeFlag==1){
    TFile* truthfile = TFile::Open("Data/WellSeparatedGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named WellSeparatedGeneratedSignals.root found"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
  
  else if(RunTypeFlag==2){
    TFile* truthfile = TFile::Open("Data/TwoTestGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named TwoTestGeneratedSignals.root found"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
   
  else if(RunTypeFlag==3){
    TFile* truthfile = TFile::Open("Data/TwoWellSeparatedGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named TwoWellSeparatedGeneratedSignals.root found"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
  
  else if(RunTypeFlag==4){
    TFile* truthfile = TFile::Open("Data/ThreeTestGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named ThreeTestGeneratedSignals.root found"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
   
  else if(RunTypeFlag==5){
    TFile* truthfile = TFile::Open("Data/ThreeWellSeparatedGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named ThreeWellSeparatedGeneratedSignals.root found"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }

  else if(RunTypeFlag==6){
    TFile* truthfile = TFile::Open("Data/OneTwoTestGeneratedSignals.root","read");
    if (!truthfile) {
      std::cout<<"No file named OneTwoTestGeneratedSignals.root found"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }

  else if(RunTypeFlag==7||RunTypeFlag==12||RunTypeFlag==14){
    //    TFile* truthfile= TFile::Open("Data/test/201215TestGeneratedSignalsLandauAmp200kRandomCalibration0point7to1point3.root","read");
    TFile* truthfile= TFile::Open("Data/test/200913TestGeneratedSignalsLandauAmp200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named 2012133TestGeneratedSignalsLandauAmp200kRandomCalibration.root found in test directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }

  else if(RunTypeFlag==8){
    TFile* truthfile = TFile::Open("Data/test/TestGeneratedSignalsLandauAmp5nsrise200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named TestGeneratedSignalsLandauAmp5nsrise200k.root found in test directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }

  else if(RunTypeFlag==9){
    TFile* truthfile = TFile::Open("Data/test/TestGeneratedSignalsLandauAmp10nsrise200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named TestGeneratedSignalsLandauAmp10nsrise200k.root found in test directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }

  else if(RunTypeFlag==10){
    TFile* truthfile = TFile::Open("Data/test/TestGeneratedSignalsFixedAmp200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named TestGeneratedSignalsFixedAmp200k.root found in test directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }

  else if(RunTypeFlag==11){
    TFile* truthfile = TFile::Open("Data/test/201003TestGeneratedSignals5HitsPerEvent200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named 201003TestGeneratedSignals5HitsPerEvent200k.root found in test directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
  
  else if(RunTypeFlag==13){
    TFile* truthfile = TFile::Open("Data/test/201126TestGeneratedSignalsLandauAmp200kExpRise.root","read");
    if (!truthfile) {
      std::cout<<"No file named 201126TestGeneratedSignalsLandauAmp200kExpRise.root found in test directory"<<std::endl;
      return nullptr;
    }
    else return truthfile;
  }
  else{
    std::cout<<"No truth file to open"<<std::endl;
    return nullptr;
  }
}

void OpenNNFile(Int_t RunTypeFlag, ifstream& NNFile){
  if(RunTypeFlag==0){
    NNFile.open("Data/PVetoMCNNGeneralTestHits.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTestNoisyHits.txt found"<<std::endl;
      return;
    }
  }

  else if(RunTypeFlag==4){
    NNFile.open("Data/PVetoMCNNThreeTestNoisyHits.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNThreeTestNoisyHits.txt found"<<std::endl;
      return;
    }
  }

  else if(RunTypeFlag==6){
    NNFile.open("Data/PVetoMCNNOneTwoTestNoisyHits.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNOneTwoTestNoisyHits.txt found"<<std::endl;
      return;
    }
  }

  else if(RunTypeFlag==7){
    NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/Standard/StandardResults/Model2Sept200913PVetoMCNNTestNoisyHitsLandauAmp200k.txt");
    if (!NNFile) {
      std::cout<<"No file named Model2Sept200913PVetoMCNNTestNoisyHitsLandauAmp200k.txt found"<<std::endl;
      return;
    }
  }
  
  else if(RunTypeFlag==8){
    NNFile.open("Data/PVetoMCNNTestNoisyHitsLandauAmp5nsrise200k.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTestNoisyHitsLandauAmp5nsrise200k.txt found"<<std::endl;
      return;
    }
  }
  
  else if(RunTypeFlag==9){
    NNFile.open("Data/PVetoMCNNTestNoisyHitsLandauAmp10nsrise200k.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTestNoisyHitsLandauAmp10nsrise200k.txt found"<<std::endl;
      return;
    }
  }
  
  else if(RunTypeFlag==10){
    NNFile.open("Data/PVetoMCNNTestNoisyHitsFixedAmp200k.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTestNoisyHitsFixedAmp200k.txt found"<<std::endl;
      return;
    }
  }

 else if(RunTypeFlag==11){
    NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/5HitsPerEvent/5HitsPerEventResults/PVetoMCNNTest5HitsPerEvent200k.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTest5HitsPerEvent200k.txt found"<<std::endl;
      return;
    }
  }

 else if(RunTypeFlag==12){
    NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/5HitsPerEvent/3HitsPerEventResults/PVetoMCNNTest3HitsPerEventWith5HitsPerEventModel.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTest3HitsPerEventWith5HitsPerEventModel.txt found"<<std::endl;
      return;
    }
  }

 else if(RunTypeFlag==13){
    NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/ExpRise/ExpRiseResults/PVetoMCNNTestLandauAmpExpRise200k.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTest3HitsPerEventWith5HitsPerEventModel.txt found"<<std::endl;
      return;
    }
  }

 else if(RunTypeFlag==14){
   //NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/RandomCalibration/RandomCalibrationResults/PVetoMCNNTestNoisyHitsRandomCalibration0point7to1point3_200k.txt");//standard signals, multiplied by random calibration constant in python script
   NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/RandomCalibration/RandomCalibrationResults/PVetoMCNNTestNoisyHitsStandardSampleRandomCalibration0point7to1point3_200k.txt");//standard signals, multiplied by random calibration constant in python script
       //    NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/RandomCalibration/RandomCalibrationResults/PVetoMCNNTestNoisyHitsWithoutRandomCalibration200k.txt");
    if (!NNFile) {
      std::cout<<"No file named PVetoMCNNTestNoisyHitsRandomCalibration1to3_200k.txt found"<<std::endl;
      return;
    }
  }

 else{
   std::cout<<"No reco file to open"<<std::endl;
   return nullptr;
 }

}
