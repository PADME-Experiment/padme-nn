#include "TMath.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <algorithm>
#include <fstream>

#define NPOINTS 1024
#define NsPerBin 400./1024.
#define BinsPerNs 1024./400.
#define rise 7//ns
#define tau 20//ns
#define UseCalibrationFlag 0//multiply hit amplitude by random constant to represent unequal gain of veto channels
#define VetoChannelTot 90//total number of veto channels

char name[256];

void u_expo(Double_t *t, Double_t *uin){//creates array of 1024 points of linearly rising and then exponentially faling signal after time = risetime
  //Time is in ns!
  for(int i=0;i<NPOINTS;i++) {
    // Step function at 10 ns -> 1V 
    if(t[i] < rise ) {
      uin[i] =t[i]/rise;
      //      std::cout<<"i " <<i<<" t "<<t[i]<<" uin "<<uin[i]<<" case1"<<std::endl;
    } else { 
      uin[i] = exp(-(t[i] - rise)/tau);
      //      std::cout<<"i " <<i<<" t "<<t[i]<<" uin "<<uin[i]<<" case2"<<std::endl;
    }
  }
}

void u_scint(Double_t *t, Double_t *uin){//creates array of 10000 points of exponentially faling signal *after* time = offset with decay constant tau=20e-9s. This exponential doesn't have a non-zero rise time
//Time is in ns!
  for(int i=0;i<NPOINTS;i++) {
    // Step function at 10 ns -> 1V 
    uin[i] =-1*exp(-(t[i]/rise))+exp(-(t[i]/tau));
    //      std::cout<<"i " <<i<<" t "<<t[i]<<" uin "<<uin[i]<<std::endl;
  }
}

void u_gauss(Double_t *t, Double_t *uin){//creates array of 1024 points of a gaussian signal
  double sigma = 0.7;
  float mu = 20.;
  
  for(int i=0;i<NPOINTS;i++){
    t[i]=i;
    uin[i]=TMath::Gaus(i, mu, sigma, 1);//returns the value of a normalalised gaussian distribution at point i
  }
}

int main(Int_t argc, char **argv){

  int eventtot; //number of events to run over

  //input argument order is: No. events, sigfile, hitfile, rootfile
  
  if(argc>1){
    eventtot=atoi(argv[1]);
    std::cout<<"Setting event tot to "<<argv[1]<<std::endl;
  }
  else{
    eventtot=10e3;//100e3;//argc is number of parameters given to code in command line. If no number of events is given, default to running 100k events.
    std::cout<<"Eventtot defaults to "<<eventtot<<std::endl;
  }

  std::ofstream sigfile;
  if(argc>2){
    sigfile.open(argv[2],std::ofstream::out);//set file name
    std::cout<<"Setting signal file name to "<<argv[2]<<std::endl;
  }
  else{
   sigfile.open("PVetoMCNoisySignals.txt", std::ofstream::out);
   std::cout<<"Signal filename defaults to PVetoMCNoisySignals.txt"<<std::endl;
  }

  std::ofstream hitfile;
  if(argc>3){
    hitfile.open(argv[3],std::ofstream::out);//set file name
    std::cout<<"Setting hit file name to "<<argv[3]<<std::endl;
  }
  else{
   hitfile.open("PVetoMCNoisyHits.txt", std::ofstream::out);
   std::cout<<"Hit filename defaults to PVetoMCNoisyHits.txt"<<std::endl;
  }

  TFile *fileout;
  if(argc>4){
   fileout = TFile::Open(argv[4],"RECREATE");
   std::cout<<"Setting ROOT file name to "<<argv[4]<<std::endl;
  }
  else{
    fileout=TFile::Open("GeneratedSignals.root","RECREATE");
    std::cout<<"ROOT filenmae defaults to GeneratedSignals.root"<<std::endl;
  }

  TTree *tree = new TTree("SigTree","Generated Signal Tree");
  
  Int_t tbins[NPOINTS];//samples
  Double_t twindow[NPOINTS];//holds sample times in ns
  Double_t time_bin_width=400./NPOINTS; //ns - 1 signal window = 400ns long  
  for(int i=0;i<NPOINTS;i++) {
    twindow[i] = i*time_bin_width;
  }

  Bool_t NoiseFlag=1;//NoiseFlag=0 doesn't add noise, NoiseFlag=1 adds noise uniformly distributed between [-0.1,0.1], +/-1% of the maximum amplitude
  Bool_t SeparationFlag=0;//SeparationFlag=1 => check for well separated signals, SeparationFlag=0 => don't check

  Bool_t separationcheck;//1 => well separated, 0 =>not well separated
  Int_t notseparatedno = 0;  
  Double_t singlesignal[NPOINTS];//a single pulse for a particle with arrival time 0 and amplitude 1
  Double_t signoise[NPOINTS];//holds the amplitude of noise for each sample in a signal  
  Double_t sigfinal[NPOINTS];//holds the sum of all individual signals in a single event
  Double_t currentsig[NPOINTS];//temporary container for individual hits, rewritten each hit
  
  Int_t TotalToyNo=0;//total number of hits generated
  Int_t lambda=3;//mean number of hits per event
  

  TRandom3 *myRNG=new TRandom3();
  if(SeparationFlag==0)  myRNG->SetSeed(0);//Set seed as clock time
  else myRNG->SetSeed(2);//to create a reproduceable sample for the separable signals, always use the same seed
 
  int zerohiteventnumber=0;

  std::vector<Double_t> arrivaltimens;//time of arrival of signal in ns, rewritten every event
  std::vector<Int_t>    arrivaltimebins;//time of arrival of signal in bins, rewritten every event
  std::vector<Double_t> sortedarrivaltimens;//time of arrival of signal in ns, rewritten every event
  std::vector<Int_t>    sortedarrivaltimebins;//time of arrival of signal in bins, rewritten every event
    
  std::vector<Double_t> hitamp;//Amplitude of signals
  std::vector<Double_t> sortedhitamp;//Amplitude of signals
  
  Int_t eventnumber=0;
  Int_t nhits;
  Double_t CalibrationConst[VetoChannelTot];
  Int_t EventChannel;
  
  //Set up ROOT Tree
  tree->Branch("HitTimeBins",&sortedarrivaltimebins);//Hit arrival time in bins, sorted by order of increasing time
  tree->Branch("HitTimeNs",&sortedarrivaltimens);//Hit arrival time in bins, sorted by order of increasing time
  tree->Branch("HitAmplitude",&sortedhitamp);//Hit amplitude, sorted by order of increasing time
  tree->Branch("SigFinal",&sigfinal,"SigFinal[1024]/D");//Final signal - what the DAQ would register, so the sum of all individual hits
  tree->Branch("EventNo",&eventnumber,"EventNo/I");//Index of event
  tree->Branch("HitsPerEvent",&nhits,"HitsPerEvent/I");//Number of hits per event
  tree->Branch("TBins",&tbins,"TBins[1024]/I");//Time in bins
  tree->Branch("EventChannel",&EventChannel,"EventChannel/I");//Finger of PVeto in which event occurs
  tree->Branch("CalibrationConst",&CalibrationConst,"CalibrationConst[90]/D");//Calibration constant for PVeto channel

  if(UseCalibrationFlag==1)  for(int channel=0;channel<VetoChannelTot;channel++) CalibrationConst[channel]=(myRNG->Uniform(0.7,1.3));
  
  for (eventnumber=0;eventnumber<eventtot;eventnumber++){
    nhits=myRNG->Poisson(lambda);//number of hits follows a Poisson distribution with lambda hits per event window on average
    //Uniform(1,3);//
    //    std::cout<<"nhits "<<nhits<<std::endl;
    if(nhits!=0)   u_expo(twindow,singlesignal);
    //if(nhits!=0)   u_gauss(twindow,singlesignal);
    //if(nhits!=0)   u_scint(twindow,singlesignal);
    
    arrivaltimens.clear();//time of arrival of signal in ns
    arrivaltimebins.clear();
    hitamp.clear();
    sortedarrivaltimens.clear();
    sortedarrivaltimebins.clear();
    sortedhitamp.clear();
    if(SeparationFlag==1&&nhits>1) separationcheck=1;//assume signals are well separated until proven otherwise
    else if(SeparationFlag==1&&nhits<2) {
      separationcheck=0;
      notseparatedno++;
      continue;//if there are fewer than 2 hits in an event, continue to the next event
    }
    
    std::fill(sigfinal, sigfinal+NPOINTS, 0);

    EventChannel=myRNG->Uniform(0,90); //assigns a random number between 0 and 90 to the event to represent the finger of the veto that the event occurs in
    
    if(eventnumber%100==0)
      std::cout<<"eventnumber "<<eventnumber<<" nhits = "<<nhits<<std::endl;

    for(int i=0; i<nhits; i++) {
      TotalToyNo++;
      std::fill(currentsig, currentsig+NPOINTS, 0);
      hitamp.push_back(myRNG->Landau(50,12));//very approximate estimation of average amplitude in mv
      if(UseCalibrationFlag==1)	hitamp[i]=hitamp[i]*CalibrationConst[EventChannel];
     
      arrivaltimens.push_back(myRNG->Uniform(80,280));
      arrivaltimebins.push_back(int(arrivaltimens[i]*BinsPerNs));
      /*      if(SeparationFlag==1&&nhits>1&&i>0&&TMath::Abs(arrivaltimens[i]-arrivaltimens[i-1])<20){
	notseparatedno++;
	separationcheck=0;//if there are fewer than two hits or if the hits are not well separated, make separation check = 0 and continue to the next event
	break;
      }
      if(eventnumber==2) std::cout<<"Event 2, hit "<<i<<" arrivaltimens "<<arrivaltimens[i]<<" arrivaltimebins "<<arrivaltimebins[i]<<" multiplication "<<arrivaltimens[i]*BinsPerNs<<std::endl;
      if(eventnumber==2&&i==5) std::cout<<"Time diff = "<<TMath::Abs(arrivaltimens[2]-arrivaltimens[0])<<std::endl;*/
      for(int j=0;j<NPOINTS;j++){
	tbins[j]=j;
	if(i==0&&j<arrivaltimebins[i]) {//before the arrival of the first hit, the signal is 0
	  sigfinal[j]=0;
	  currentsig[j]=0;
	}
	else if(i==0&&j>=arrivaltimebins[0]){//once the first signal has arrived, the signal is a singlesignal shifted by arrivaltimebins[0]
	  sigfinal[j]=(hitamp[i]*singlesignal[j-arrivaltimebins[0]]);
	  currentsig[j]=(hitamp[i]*singlesignal[j-arrivaltimebins[0]]);
	}
	else if(i>0&&j>=arrivaltimebins[i]){//after the arrival of subsequent signals, the signal is the sum of the new signal and all previous signals
	  sigfinal[j]+=hitamp[i]*singlesignal[j-arrivaltimebins[i]];
	  currentsig[j]=(hitamp[i]*singlesignal[j-arrivaltimebins[i]]);
	}
	else if(i>0&&j<arrivaltimebins[i]) 	  currentsig[j]=0;
	if(i==nhits-1){ //on the final hit, the final signal array will be complete 
	  if(NoiseFlag==0)  signoise[j]=0;
	  if(NoiseFlag==1)  signoise[j]=myRNG->Uniform(-1.5,1.5);
	  sigfinal[j]+=signoise[j];
	}
      }
    } //end of hit loop
    if(SeparationFlag==0&&NoiseFlag==0&&nhits==0) for(int kk=0;kk<NPOINTS;kk++) sigfile<<0<<std::endl;//for any events with 0 hits, print out an event of 0s.*/
    if(SeparationFlag==0&&NoiseFlag==1&&nhits==0) for(int kk=0;kk<NPOINTS;kk++){
	sigfinal[kk]=myRNG->Uniform(-1.5,1.5);//for any events with 0 hits, the event contains only noise.
	sigfile<<sigfinal[kk]<<std::endl;
      }

    std::vector<int> index(arrivaltimebins.size(), 0);
  
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }

    sort(index.begin(), index.end(),
	 [&](const int& a, const int& b) {
	   return (arrivaltimebins[a] < arrivaltimebins[b]);
	 }
	 );

    for (int ii = 0 ; ii != index.size() ; ++ii) {
      sortedhitamp.push_back(hitamp[index[ii]]);
      sortedarrivaltimens.push_back(arrivaltimens[index[ii]]);
      sortedarrivaltimebins.push_back(arrivaltimebins[index[ii]]);
      if(SeparationFlag==1&&(ii>0&&sortedarrivaltimebins[ii]-sortedarrivaltimebins[ii-1]<20*BinsPerNs)){
	std::cout<<"Event number "<<eventnumber<<" hits "<<ii-1<<" and "<<ii<<" times "<<sortedarrivaltimebins[ii-1]<<" and "<<sortedarrivaltimebins[ii]<<" separationcheck = "<<separationcheck<<std::endl;
	separationcheck=0;
      }
    }
    if(separationcheck==0) notseparatedno++;
    
    if(SeparationFlag==0||(SeparationFlag==1&&separationcheck==1)){ //if either the SeparationFlag = 0 or the SeparationFlag = 1 and separation check = 1, write the files
      if(SeparationFlag==1) std::cout<<"Flag "<<SeparationFlag<<" check "<<separationcheck<<std::endl;
      hitfile<<nhits<<std::endl;
      tree->Fill();
      if(nhits>0)  for(int sample=0; sample<NPOINTS; sample++) sigfile<<sigfinal[sample]<<std::endl; //if nhits ==0, the sigfile will already have been written above
    }
  
  }// end of event loop
  fileout->Write();
  delete fileout;

  hitfile.close();
  sigfile.close();
  if(SeparationFlag==1)  std::cout<<"Number of not well separated events = "<<notseparatedno<<std::endl;
}
