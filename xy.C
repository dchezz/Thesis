
// Description:
//      
//       Analysis macro
//
//
// Environment:
//      MPDROOT
//
// Author List:
//       Luis Valenzuela-Cazares          (original author)
//   
//-----------------------------------------------------------
#if !defined(__CINT__) && !defined(__CLING__)
#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"

#include "FairRunSim.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairTrajFilter.h"
#include "FairUrqmdGenerator.h"
#include "FairPrimaryGenerator.h"
#include "FairCave.h"
#include "FairPipe.h"
#include "FairMagnet.h"

#include "MpdSts.h"
#include "TpcDetector.h"
#include "MpdEtof.h"
#include "MpdFsa.h"
#include "MpdBbc.h"
#include "MpdCpc.h"
#include "MpdTof.h"
#include "MpdStrawendcap.h"
#include "MpdZdc.h"
#include "MpdFfd.h"
#include "MpdGetNumEvents.h"

#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

#endif

#include "/home/dario/NICA/mpdroot/macro/mpd/mpdloadlibs.C"

/////////////////////////////////////////
  //Function to get the ring.
   Int_t GetRing(Int_t detID){

  if(detID >= 1 && detID <= 12) return 1;
  else if (detID>=13 && detID<=30) return 2;
  else if (detID>=31 && detID<= 54) return 3;
  else if (detID>=55 && detID<=84) return 4;
  else if (detID>=85 && detID<=120) return 5;
  else if (detID>=121 && detID<=162) return 6;
  
  else if (detID>=163 && detID<=174) return 1;
  else if (detID>=175 && detID<=192) return 2;
  else if (detID>=193 && detID<=216) return 3;
  else if (detID>=217 && detID<=246) return 4;
  else if (detID>=247 && detID<=282) return 5;
  else if (detID>=283 && detID<=324) return 6;

   else return -1;

}
        
Int_t xy(TString inputFile="/home/dario/NICA/mpdroot/macro/mpd/bb_urqmd_10000.root",TString outputFile="xy_urqmd_10000.root")
{
    gROOT->Reset();
    gStyle->SetOptDate(0);      //show day and time
    gStyle->SetOptStat(1);      //show statistic
    gStyle->SetPalette(1,0);    
      
    mpdloadlibs(); 
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////Y vs X both BMD///////////Primary charged Particles//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////X vs Y BMD A/////////////////////////////////////////
  TH2D *hBmdPointXYA  = new TH2D("hBmdPointXYA","BMD A particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",10000,-50,50,1000,-50,50);
  hBmdPointXYA->SetXTitle("Position X (cm)");
  hBmdPointXYA->SetYTitle("Position Y (cm)");

     //Histograms for particle distribution in the plane X-Y in each ring of both BMD detectors.

 ///////////////////////////////////////X vs Y BMD A/////////////////////////////////////////
  TH2D *hBmdPointRing1XYA  = new TH2D("hBmdPointRing1XYA","BMDA ring 1 particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",1000,-50,50,1000,-50,50);
  hBmdPointRing1XYA->SetXTitle("Position X (cm)");
  hBmdPointRing1XYA->SetYTitle("Position Y (cm)");

  TH2D *hBmdPointRing2XYA  = new TH2D("hBmdPointRing2XYA","BMDA ring 2 particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",1000,-50 ,50,1000,-50,50);
  hBmdPointRing2XYA->SetXTitle("Position X (cm)");
  hBmdPointRing2XYA->SetYTitle("Position Y (cm)");  
  
  TH2D *hBmdPointRing3XYA  = new TH2D("hBmdPointRing3XYA","BMDA ring 3 particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",1000,-50,50,1000,-50,50);
  hBmdPointRing3XYA->SetXTitle("Position X (cm)");
  hBmdPointRing3XYA->SetYTitle("Position Y (cm)");  
 
  TH2D *hBmdPointRing4XYA  = new TH2D("hBmdPointRing4XYA","BMDA ring 4 particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",1000,-50,50,1000,-50,50);
  hBmdPointRing4XYA->SetXTitle("Position X (cm)");
  hBmdPointRing4XYA->SetYTitle("Position Y (cm)"); 

  TH2D *hBmdPointRing5XYA  = new TH2D("hBmdPointRing5XYA","BMDA ring 5 particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",1000,-50,50,1000,-50,50);
  hBmdPointRing5XYA->SetXTitle("Position X (cm)");
  hBmdPointRing5XYA->SetYTitle("Position Y (cm)"); 
   
  TH2D *hBmdPointRing6XYA  = new TH2D("hBmdPointRing6XYA","BMDA particle distribution in plane X-Y. 10000 Au+Au @11GeV URQMD. MPDROOT.",1000,-50,50,1000,-50,50);
  hBmdPointRing6XYA->SetXTitle("Position X (cm)");
  hBmdPointRing6XYA->SetYTitle("Position Y (cm)");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////Primary charged particles (PCP) multiplicity/////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////Multiplicity per ring BMD A/////////////////////////////////////////////////////// 
  TH1D *hBmdR1Mul = new TH1D("hBmdR1Mul","Primary charged particles multiplicity BMD A ring 1. 10000 Au+Au @11GeV UrQMD. MPDROOT.",50,0,50);
  hBmdR1Mul->SetXTitle("Ring");
  hBmdR1Mul->SetYTitle("Multiplicity");

  TH1D *hBmdR2Mul = new TH1D("hBmdR2Mul","Primary charged particles multiplicity BMD A ring 2. 10000 Au+Au @11GeV UrQMD. MPDROOT.",50,0,50);
  hBmdR2Mul->SetXTitle("Ring");
  hBmdR2Mul->SetYTitle("Multiplicity");

  TH1D *hBmdR3Mul = new TH1D("hBmdR3Mul","Primary charged particles multiplicity BMD A ring 3. 10000 Au+Au @11GeV UrQMD. MPDROOT.",50,0,50);
  hBmdR3Mul->SetXTitle("Ring");
  hBmdR3Mul->SetYTitle("Multiplicity");

  TH1D *hBmdR4Mul = new TH1D("hBmdR4Mul","Primary charged particles multiplicity BMD A ring 4. 10000 Au+Au @11GeV UrQMD. MPDROOT.",50,0,50);
  hBmdR4Mul->SetXTitle("Ring");
  hBmdR4Mul->SetYTitle("Multiplicity");

  TH1D *hBmdR5Mul = new TH1D("hBmdR5Mul","Primary charged particles multiplicity BMD A ring 5. 10000 Au+Au @11GeV UrQMD. MPDROOT.",50,0,50);
  hBmdR5Mul->SetXTitle("Ring");
  hBmdR5Mul->SetYTitle("Multiplicity");

  TH1D *hBmdR6Mul = new TH1D("hBmdR6Mul","Primary charged particles multiplicity BMD A ring 6. 10000 Au+Au @11GeV UrQMD. MPDROOT.",50,0,50);
  hBmdR6Mul->SetXTitle("Ring");
  hBmdR6Mul->SetYTitle("Multiplicity");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Defining pointers to data
    
    TFile fileInput(inputFile.Data());
    
    TTree *simEvent = (TTree*) fileInput.Get("mpdsim");
    
    //Arrays/////nombre de la variable///////////nombre en la memoria interna
    TClonesArray *bmdPoints = (TClonesArray*) fileInput.FindObjectAny("BmdPoint");
    simEvent->SetBranchAddress("BmdPoint",&bmdPoints);
    
 //   TClonesArray *bbcPoints = (TClonesArray*) fileInput.FindObjectAny("BBCPoint");
 //    simEvent->SetBranchAddress("BBCPoint",&bbcPoints);

    TClonesArray* mcTracks = (TClonesArray*) fileInput.FindObjectAny("MCTrack");
    simEvent->SetBranchAddress("MCTrack",&mcTracks);
  
 //  TClonesArray* tpcPoints = (TClonesArray*) fileInput.FindObjectAny("TpcPoint");
 //   simEvent->SetBranchAddress("TpcPoint",&tpcPoints);
  
    FairMCEventHeader* mcHeader = (FairMCEventHeader*) fileInput.FindObjectAny("MCEventHeader.");
    simEvent->SetBranchAddress("MCEventHeader.",&mcHeader);
        
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


//MCTrack kinematics
//MCEventHeader kind of collision  

    //Declaring variable for the impact parameter.
    Double_t impactParameter; //! impact parameter 

     TFile* fileOutput = new TFile(Form("%s",outputFile.Data()),"RECREATE");

     //Defining the tree and a branch for the impact parameter.
     TTree *bmdTree = new TTree("bmdTree","BMD");
     bmdTree->Branch("impactParameter",&impactParameter);
    
     Int_t events  = simEvent->GetEntries();
  
    //Begins loop for the events
    for (Int_t i = 0000; i < events; ++i) {
    simEvent->GetEntry(i);

     MpdMCEventHeader *extraEventHeader = dynamic_cast<MpdMCEventHeader*> (mcHeader);

      bmdTree->Fill();
  
    Int_t nMCTracks  =  mcTracks->GetEntriesFast();
    Int_t nprimary   =  mcHeader->GetNPrim(); 
  //  cout<<"nMCTracks "<<nMCTracks<<endl;

    Double_t B  = mcHeader->GetB();

    if (bmdPoints != 0 ) {

    Int_t nbmdPoints =  bmdPoints->GetEntriesFast();
    if( nbmdPoints < 1 ) continue;
    
    //cout<<"nbmdPoints: "<<nbmdPoints<<endl;
       
    //Loop for particles in both BMD.
   
      for (Int_t j = 0; j < nbmdPoints; j++) {
      BmdPoint *p1    = (BmdPoint*) bmdPoints->At(j);
      TVector3 recMom;
      Int_t    trkID = p1->GetTrackID();
      Double_t time  = p1->GetTime();
      Double_t eLoss = p1->GetEnergyLoss();
      p1->Momentum(recMom);
      Double_t currPt  = recMom.Pt();
      Double_t currEta = recMom.Eta();
      
      Int_t detID =  p1->GetDetectorID();
     // Int_t statusCode = p1->GetStatusCode();
 
    //    cout<<"trkID: "<<trkID<<endl;  
    //    cout<<"detID: "<<detID<<endl;            

     Int_t x = p1->GetX();
     Int_t y = p1->GetY();
     Int_t z = p1->GetZ();

   //  cout<<"x: "<<x<<"y: "<<y<<"z: "<<z<<endl;   

 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////Analysis per ring////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t ring = -1;
  
ring = GetRing(detID);
 //cout<<"ring: "<<ring<<endl;

        MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(trkID);

          Double_t Xv = mcTr->GetStartX();
          Double_t Yv = mcTr->GetStartY();
          Double_t Zv = mcTr->GetStartZ();
      
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////Primary Particles Condition////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////
  if( mcTr->GetMotherId() < 0 ) {

   //Charged primary particles
         if ( TMath::Abs(mcTr->GetPdgCode() )  != 211  &&
                  TMath::Abs(mcTr->GetPdgCode() )  != 11   &&
                  TMath::Abs(mcTr->GetPdgCode() )  != 2212 && //Protons
                  TMath::Abs(mcTr->GetPdgCode() )  != 321  &&
                  TMath::Abs(mcTr->GetPdgCode() )  != 13   &&
                  (TMath::Abs(Xv) > 0.001 || TMath::Abs(Yv) > 0.001 || TMath::Abs(Zv) > 0.01 )) continue;

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Primary Charged Particles BMD A/////////////////////////
////////////////////////////////////////////////////////////////////////////////////////   

/////////////////////////////////Rings////////////////////////////////////////////        
if(z >=199){ //BMD A    

   //if( TMath::Abs(recMom.Eta()) >  3.5 && TMath::Abs(recMom.Eta()) <= 4.2 ) {
    if( ring==1){  //Ring 1     
             // cout<<"z: "<<z<<endl;                    
                  hBmdPointRing1XYA->Fill(p1->GetX(),p1->GetY());
                  hBmdPointXYA->Fill(p1->GetX(),p1->GetY());
     }    

  // if( TMath::Abs(recMom.Eta()) >  3.2 && TMath::Abs(recMom.Eta()) <= 3.5 ) {                                 
    if( ring==2){ //Ring 2
              //  cout<<"z: "<<z<<endl;       
                 hBmdPointRing2XYA->Fill(p1->GetX(),p1->GetY());  
                 hBmdPointXYA->Fill(p1->GetX(),p1->GetY());

     }    

 //if( TMath::Abs(recMom.Eta()) >  3.0 && TMath::Abs(recMom.Eta()) <= 3.2 ) {      
    if( ring==3){  //Ring 3
                 hBmdPointRing3XYA->Fill(p1->GetX(),p1->GetY());
                 hBmdPointXYA->Fill(p1->GetX(),p1->GetY());

     }    
 
 //if( TMath::Abs(recMom.Eta()) >  2.8 && TMath::Abs(recMom.Eta()) <= 3.0 ) {                      
    if( ring==4){  //Ring 4 
           hBmdPointRing4XYA->Fill(p1->GetX(),p1->GetY());
           hBmdPointXYA->Fill(p1->GetX(),p1->GetY());
     }    
  
//if( TMath::Abs(recMom.Eta()) >  2.6 && TMath::Abs(recMom.Eta()) <= 2.8 ) {                        
    if( ring==5){   //Ring 5    
              hBmdPointRing5XYA->Fill(p1->GetX(),p1->GetY());
              hBmdPointXYA->Fill(p1->GetX(),p1->GetY());
    }    

 //if( TMath::Abs(recMom.Eta()) >=  2.2 && TMath::Abs(recMom.Eta()) <= 2.6 ) {        
    if( ring==6){ //Ring 6
               hBmdPointRing6XYA->Fill(p1->GetX(),p1->GetY());
               hBmdPointXYA->Fill(p1->GetX(),p1->GetY());
     }    


   } //BMD A

///////////////////////////////////////////////////////////////////////////////////////////
     
}


}  //Primary Particles Condition 

}
 

  } //End of the loop's event. 
  
 ///////////////////////////////////////// 

  bmdTree->Write();
  
///////////////////////////////////////////////  
   cout<<"Saving histograms"<<endl;
  
   
   fileOutput->mkdir("Mc");
   fileOutput->cd("Mc");
   
  //gStyle->SetOptTitle(0); //No title for histograms

///////////////////////////////////////////BMD A/////////////////////////////////////////////////////////////////////////

   TCanvas *c20 = new TCanvas("c20","Rings primary charged particles distribution BMD side A");
   //c20->cd(1); 
   hBmdPointXYA->SetStats(kFALSE);
   hBmdPointXYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointXYA->GetYaxis()->SetTitle("Y(cm)");
  // hBmdPointXYA->Draw();
   hBmdPointXYA->Draw("colz");

   c20->SaveAs("xy.pdf");

//////////////////////////////////////////BMD A Rings////////////////////////////////////////////////////////////////
   TCanvas *c21 = new TCanvas("c21","Rings primary charged particles distribution BMD side A");
   c21->Divide(3,2);
   c21->cd(1); 
   hBmdPointRing1XYA->SetStats(kFALSE);
   hBmdPointRing1XYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointRing1XYA->GetYaxis()->SetTitle("Y(cm)");
   hBmdPointRing1XYA->Draw("colz");

   c21->cd(2);
   hBmdPointRing2XYA->SetStats(kFALSE);
   hBmdPointRing2XYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointRing2XYA->GetYaxis()->SetTitle("Y(cm)");
   hBmdPointRing2XYA->Draw("colz");

   c21->cd(3);
   hBmdPointRing3XYA->SetStats(kFALSE);
   hBmdPointRing3XYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointRing3XYA->GetYaxis()->SetTitle("Y(cm)");
   hBmdPointRing3XYA->Draw("colz");

   c21->cd(4); 
   hBmdPointRing4XYA->SetStats(kFALSE);
   hBmdPointRing4XYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointRing4XYA->GetYaxis()->SetTitle("Y(cm)");
   hBmdPointRing4XYA->Draw("colz");

   c21->cd(5);
   hBmdPointRing5XYA->SetStats(kFALSE);
   hBmdPointRing5XYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointRing5XYA->GetYaxis()->SetTitle("Y(cm)");
   hBmdPointRing5XYA->Draw("colz");

   c21->cd(6);
   hBmdPointRing6XYA->SetStats(kFALSE);
   hBmdPointRing6XYA->GetXaxis()->SetTitle("X(cm)");
   hBmdPointRing6XYA->GetYaxis()->SetTitle("Y(cm)");
   hBmdPointRing6XYA->Draw("colz");


   //c21->SaveAs("xyrrings.pdf");


  cout<<"End histograms"<<endl;
  return 0;
   
}

