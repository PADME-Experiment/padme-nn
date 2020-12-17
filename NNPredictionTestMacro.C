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

  Int_t RunTypeFlag=7;//0 = General Test (NoisyTestHits), 1 = Well Separated Signals, 2 = Exactly Two General Signals, 3 = Exactly Two Well Separated Signals, 4 = Exactly 3 General Signals, 5 - Exactly 3 Well Separated Signals, 6 = One or Two General Signals, 7 = Landau Distributed Standard Signals, 8 = Signals with 5ns Rise Time, 9 = Signals with 10ns Rise Time, 10 = Signals with 7ns Rise Time && Fixed 50mV Amplitude, 11 = Signals with mean 5 hits per event, 12 = Signals with mean 3 hits per event analysed with 5HitPerEventModel, 13 = Signals with exponential rise and fall, mean 3 hits per event, 14 = standard signals mutliplied by uniformly distributed random "calibration constant"~U[1,3] to represent the different gains of differen PVeto channels

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

  //Read branch info from .root files into vectors
  tree->SetBranchAddress("HitsPerEvent",&TruthHitsPerEvent);  
  tree->SetBranchAddress("HitTimeNs",&TruthHitTimes);
  if(UseDerivFlag==1) tree->SetBranchAddress("NoRecoDeriv",&DerivHitsPerEvent);  


  std::vector<Double_t> TruthHitTimeDiffs;//Difference in hit arrival time in ns
  Double_t TruthMinHitTimeDiffs;//Minimum difference in hit arrival time for each event, in ns

  //Set up analysis histograms
  TH1F *hDerivHitNumbers;
  if(UseDerivFlag==1)  hDerivHitNumbers = new TH1F("hDerivHitNumbers","Hits per event",15,0,14);//number of hits per event from neural network
  TH1F *hNNHitNumbers = new TH1F("hNNHitNumbers","Hits per event",15,0,14);//number of hits per event from neural network
  TH1F *hTrueHitNumbers = new TH1F("hTrueHitNumbers","Hits per event",15,0,14);//number of hits per event from simulation
  
  TH1F *hDerivCorrectEventsVsTrueTime;
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime = new TH1F("hDerivCorrectEventsVsTrueTime","Minimum time difference between hits per correctly reconstructed event",MaxTimeWindow,0,MaxTimeWindow);//Number of events where derivative reconstructs correct no. hits
  TH1F *hNNCorrectEventsVsTrueTime = new TH1F("hNNCorrectEventsVsTrueTime","Minimum time difference between hits per correctly reconstructed event",MaxTimeWindow,0,MaxTimeWindow);//Number of events where NN reconstructs correct no. hits
  TH1F *hTrueEventsVsTrueTime = new TH1F("hTrueHitsVsTrueTime","Minimum time difference between hits per event",MaxTimeWindow,0,MaxTimeWindow);//Time spectrum of simulated events
  
  Int_t NNWrongNo=0;
  Int_t DerivWrongNo=0;
  
  for(Int_t EventIndex=0;EventIndex<TotEvents;EventIndex++){
    TruthHitTimeDiffs.clear();
    tree->GetEntry(EventIndex);
    if(EventIndex%100==0)  std::cout<<"Event: "<<EventIndex<<" Truth hits: "<<TruthHitsPerEvent<<" "<<" NN Predicted hits: "<<NNHitNumbers[EventIndex]<<" Derivative Predicted hits: "<<DerivHitsPerEvent<<std::endl;
    
    if(UseDerivFlag==1)    hDerivHitNumbers->Fill(DerivHitsPerEvent);
    hNNHitNumbers->Fill(NNHitNumbers[EventIndex]);
    hTrueHitNumbers->Fill(TruthHitsPerEvent);
    
    if(TruthHitsPerEvent-NNHitNumbers[EventIndex]!=0) NNWrongNo++;
    if(UseDerivFlag==1&&TruthHitsPerEvent-DerivHitsPerEvent!=0) DerivWrongNo++;

    if(TruthHitsPerEvent<2) continue;//if there are fewer than two hits in an event, there aren't multiple hits for which to evaluate the time difference
    for(int ii=0;ii<TruthHitsPerEvent;ii++){
      if(ii<TruthHitsPerEvent-1) TruthHitTimeDiffs.push_back(TruthHitTimes->at(ii+1)-TruthHitTimes->at(ii));
      TruthMinHitTimeDiffs=*min_element(TruthHitTimeDiffs.begin(), TruthHitTimeDiffs.end());
    }
    hTrueEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
    if(TruthHitsPerEvent==NNHitNumbers[EventIndex])   hNNCorrectEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
    if(UseDerivFlag==1&&TruthHitsPerEvent==DerivHitsPerEvent)   hDerivCorrectEventsVsTrueTime->Fill(TruthMinHitTimeDiffs);
  }
  
  TLegend* NNLegend = new TLegend(0.6,0.8,0.9,0.9);
  //gStyle->SetLegendTextSize(0.03);
  NNLegend->AddEntry(hTrueHitNumbers,"ToyMC Truth");
  if(RunTypeFlag==14) NNLegend->AddEntry(hNNHitNumbers,"Neural Network with Random Multiplication");
  else NNLegend->AddEntry(hNNHitNumbers,"Neural Network");
  
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
  if(UseDerivFlag==1)  AllLegend->AddEntry(hDerivHitNumbers,"Derivative + TSpectrum");
  AllLegend->AddEntry(hNNHitNumbers,"Neural Network");

  TLegend* RatioLegend = new TLegend(0.6,0.8,0.9,0.9);
  if(UseDerivFlag==1)  RatioLegend->AddEntry(hDerivHitNumbers,"Derivative + TSpectrum");
  RatioLegend->AddEntry(hNNHitNumbers,"Neural Network");
  
  TCanvas *cNNHitNumbers = new TCanvas("cNNHitNumbers","Truth vs NN predicted hits per event",0,0,800,800);
  TRatioPlot *NNEfficiencyRatio = new TRatioPlot(hNNHitNumbers,hTrueHitNumbers);
  NNEfficiencyRatio->Draw();
  NNEfficiencyRatio->GetUpperPad()->cd();
  hTrueHitNumbers->Draw("AH");
  hNNHitNumbers->Draw("sames");
  NNLegend->Draw("sames");
  hNNHitNumbers->SetLineColor(kRed);
  TGraph *NNHitRatio = NNEfficiencyRatio->GetLowerRefGraph();

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
  if(UseDerivFlag==1)  DerivTimeRatio->Draw("Psame");
  NNTimeRatio->SetLineColor(kRed);
  RatioLegend->Draw("same");
  NNTimeRatio->GetXaxis()->SetRangeUser(0,MaxTimeWindow);
  NNTimeRatio->GetYaxis()->SetRangeUser(0,1);
  if(UseDerivFlag==1)  DerivTimeRatio->GetYaxis()->SetRangeUser(0,1);
  //  NNTimeRatio->SetTitle("No. Reconstructed Hits/No. Truth Hits vs Min Time Diff");
  line->Draw();
  line->SetLineColor(kBlue);

  TCanvas *cAllHitNumbers = new TCanvas("cAllHitNumbers","Truth vs RRC vs derivative vs CNN predicted hits per event",0,0,800,800);
  if(UseDerivFlag==1)  {
    hDerivHitNumbers->Draw();//"sames");
    hTrueHitNumbers->Draw("sames");
  }
  else hTrueHitNumbers->Draw();
  AllLegend->Draw("sames");
  hNNHitNumbers->Draw("sames");
  if(UseDerivFlag==1)  hDerivHitNumbers->SetLineColor(kBlack);
  hTrueHitNumbers->SetFillColorAlpha(kBlue,0.25);
  if(UseDerivFlag==1)  hDerivHitNumbers->SetFillColorAlpha(kBlack,0.35);
  hNNHitNumbers->SetFillColorAlpha(kRed,0.25);

  TCanvas *cAllCorrectEventTimes = new TCanvas("cAllCorrectEventTimes","Truth/RRC/NN/Derivative predicted events per time diff for correct reconstruction",0,0,800,800);
  hTrueEventsVsTrueTime->Draw();
  hNNCorrectEventsVsTrueTime->Draw("sames");
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime->Draw("sames");
  AllLegend->Draw("sames");
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime->SetLineColor(kBlack);
  hTrueEventsVsTrueTime->SetFillColorAlpha(kBlue,0.25);
  if(UseDerivFlag==1)  hDerivCorrectEventsVsTrueTime->SetFillColorAlpha(kBlack,0.35);
  hNNCorrectEventsVsTrueTime->SetFillColorAlpha(kRed,0.25);
  gStyle->SetOptStat(0);
  
  std::cout<<"Number of wrongly reconstructed events from NN: "<<NNWrongNo<<" out of "<<TotEvents<<" from derivative "<<DerivWrongNo<<" out of "<<TotEvents<<std::endl;
  std::cout<<"NN Efficiency: "<<(TotEvents-NNWrongNo)*100/(1.*TotEvents)<<"%, Derivative Efficiency: "<<(TotEvents-DerivWrongNo)*100/(1.*TotEvents)<<"%"<<std::endl;

}

TFile* OpenTruthFile(Int_t RunTypeFlag){
  //Open truth .root file
  if(RunTypeFlag==7){
    //    TFile* truthfile= TFile::Open("/Users/bethlong/Documents/PADME/SignalToyMC/Data/test/201215TestGeneratedSignalsLandauAmp200kRandomCalibration0point7to1point3.root","read");
    TFile* truthfile= TFile::Open("/Users/bethlong/Documents/PADME/SignalToyMC/Data/test/200913TestGeneratedSignalsLandauAmp200k.root","read");
    if (!truthfile) {
      std::cout<<"No file named 2012133TestGeneratedSignalsLandauAmp200kRandomCalibration.root found in test directory"<<std::endl;
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
 if(RunTypeFlag==7){
    NNFile.open("/Users/bethlong/Google Drive/PVetoNN/PVetoMCData/Standard/StandardResults/Model2Sept200913PVetoMCNNTestNoisyHitsLandauAmp200k.txt");
    if (!NNFile) {
      std::cout<<"No file named Model2Sept200913PVetoMCNNTestNoisyHitsLandauAmp200k.txt found"<<std::endl;
      return;
    }
  }
 else{
   std::cout<<"No reco file to open"<<std::endl;
   return nullptr;
 }

}
