
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
#include "FairMCTrack.h"

#include "MpdSts.h"
#include "MpdFsa.h"
#include "MpdCpc.h"
#include "MpdStrawendcap.h"
#include "MpdGetNumEvents.h"
#include "MbbDetector.h"

#include <iostream>

#include <TPad.h>
#include <TH1F.h>
#include <THD.h>

#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

#endif


R__ADD_INCLUDE_PATH($VMCWORKDIR)
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
        

Int_t readOutputMBB(TString inputFile="/home/dario/NICA/mpdroot/macro/mpd/mbb_phsd_10000.root",TString outputFile="/home/dario/NICA/bebes/MacrosMBB/mbbresult_mbb_urqmd_10000.root")
{
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptDate(0);      //show day and time
    gStyle->SetOptStat(1);      //show statistic
    gStyle->SetPalette(1,0);
    
    //Double_t epsilon = 0.00001;
      
    mpdloadlibs();
  
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////MPDROOT Analysis//////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////Primary charged particles mbb ///////////////////////////////
    TH1D *hEta = new TH1D("hEta","",240,-6,6);
    hEta->SetXTitle("#eta");
    hEta->SetYTitle("entries");
    
    
    TH1D *hPt = new TH1D("hPt","",500,0,5);
    hPt->SetXTitle("Pt [GeV/C]");
    hPt->SetYTitle("entries");

    TH1D *hEnergy = new TH1D("hEnergy","",500,0,5);
    hEnergy->SetXTitle("E [GeV]");
    hEnergy->SetYTitle("entries");                 
     
    TH2D *hXY = new TH2D("hXY","",600,-30,30,600,-30,30);
    hXY->SetXTitle("X (cm)");
    hXY->SetYTitle("Y (cm)");  
     
     
    TH2D *hPhiVsZ = new TH2D("hPhiVsZ","",600,-3.1416,3.1416,500,-50,50);
    hPhiVsZ->SetXTitle("#phi (rad)");
    hPhiVsZ->SetYTitle("Z (cm)");
    
    
    TH2D *hYVsZ = new TH2D("hYVsZ","",500,-50,50,500,-50,50);
    hYVsZ->SetXTitle("Y (cm)");
    hYVsZ->SetYTitle("Z (cm)");
    
    
    TH3D *hXVsYVsZ = new TH3D("hXVsYVsZ","",140,-35,35,140,-35,35,140,-35,35);
    hXVsYVsZ->SetXTitle("X (cm)");
    hXVsYVsZ->SetYTitle("Y (cm)");
    hXVsYVsZ->SetZTitle("Z (cm)");
    
    
    TH2D *hPhiVsThEta = new TH2D("hPhiVsThEta","",600,-3.1416,3.1416,360,0,180);
    hPhiVsThEta->SetXTitle("#phi");
    hPhiVsThEta->SetYTitle("#theta");
   

    TH1D *hTime = new TH1D("hTime"," ",1000,0,10);
    hTime->SetXTitle("Time-of-fligth (ns)");
    hTime->SetYTitle("entries");
        
    
    TH1D *hMultiplicityVsCell = new TH1D("hMultiplicityVsCell"," ",320,-1,321);
    hMultiplicityVsCell->SetXTitle("Cell ID");
    hMultiplicityVsCell->SetYTitle("Multiplicity");
   
    
    TFile* file = TFile::Open("/home/dario/NICA/bebes/MacrosMBB/urqmd_10000/outtestf140.root");
    file->GetName();


    //Defining pointers to data
    
    TFile fileInput(inputFile.Data());
    
    TTree *simEvent = (TTree*) fileInput.Get("mpdsim");
    
    //Arrays/////nombre de la variable///////////nombre en la memoria interna
     TClonesArray *mbbPoints = (TClonesArray*) fileInput.FindObjectAny("MbbPoint");
     simEvent->SetBranchAddress("MbbPoint",&mbbPoints);
   
     TClonesArray* mcTracks = (TClonesArray*) fileInput.FindObjectAny("MCTrack");
     simEvent->SetBranchAddress("MCTrack",&mcTracks);
  
     FairMCEventHeader* mcHeader = (FairMCEventHeader*) fileInput.FindObjectAny("MCEventHeader.");
     simEvent->SetBranchAddress("MCEventHeader.",&mcHeader);
    
   
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//MCTrack kinematics
//MCEventHeader  kind of collision

    //Declaring variable for the impact parameter.
    

    TFile* fileOutput = new TFile(Form("%s",outputFile.Data()),"RECREATE");       
    
    Int_t events  = simEvent->GetEntries();    

    //Begins loop for the events
    for (Int_t i = 0; i < events; ++i) {
       
         simEvent->GetEntry(i);
         MpdMCEventHeader *extraEventHeader = dynamic_cast<MpdMCEventHeader*> (mcHeader);
    
        Int_t nMCTracks  =  mcTracks->GetEntriesFast();
        Int_t nprimary   =  mcHeader->GetNPrim(); 


    if (mbbPoints != 0 ) {

    Int_t nmbbPoints =  mbbPoints->GetEntriesFast();
    if( nmbbPoints < 1 ) continue;
          
    //cout<<"nmbbPoints: "<<nmbbPoints<<endl;
    //cout<<"nprimary:   "<<nprimary<<endl;
       
    //Loop for particles in MBB
    Int_t nRecTracks = 0;
    Int_t nRecPrimaryTracks = 0;

    for (Int_t j = 0; j < nmbbPoints; j++) {
      MbbPoint *p1    = (MbbPoint*) mbbPoints->At(j);

      TVector3 recMom;
      Int_t    trkID = p1->GetTrackID();
      Double_t time  = p1->GetTime();
      Double_t eLoss = p1->GetEnergyLoss();
      p1->Momentum(recMom);
      Double_t currPt  = recMom.Pt();
      Double_t currEta = recMom.Eta();
      
      Int_t detID =  p1->GetDetectorID();
     // Int_t statusCode = p1->GetStatusCode();
 
     //   cout<<"trkID: "<<trkID<<endl;  
      //  cout<<"detID: "<<detID<<endl;            
           
                  
     Float_t x2 = p1->GetX();
     Float_t y2 = p1->GetY();
     Float_t z2 = p1->GetZ();
     
     TVector3 pos(x2,y2,z2);

      
//MCTrack
  
         TVector3 mcMom;
        
         MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(j);

  if( mcTr->GetMotherId() < 0 ) {
         
        if (  TMath::Abs(mcTr->GetPdgCode() )  != 2212  &&  //Protons
                 TMath::Abs(mcTr->GetPdgCode() )  != 3222 &&  //Sigma+
                 TMath::Abs(mcTr->GetPdgCode() )  != 3112 &&  //Sigma-
                 TMath::Abs(mcTr->GetPdgCode() )  != 3312 &&  //Xi-
                    TMath::Abs(mcTr->GetPdgCode() )  != 211  && //Pions 
                    TMath::Abs(mcTr->GetPdgCode() )  != 11   && //Electrons
                    TMath::Abs(mcTr->GetPdgCode() )  != 321  && //Kaons
                    TMath::Abs(mcTr->GetPdgCode() )  != 13   && //Muons
                 TMath::Abs(mcTr->GetPdgCode() )  != 3334 ) continue; //Omega-
            
            //////////////////Kind of particle////////////////////////
           // if ( TMath::Abs(mcTr->GetPdgCode()) == 3334) {

                cout<<"pdgcode: "<< mcTr->GetPdgCode()<<endl; 

         // if (detID < 10){

                 hPt->Fill(recMom.Pt());
                 hEnergy->Fill(p1->GetEnergyLoss());
                 hEta->Fill(recMom.Eta());
                 hXY->Fill(x2,y2);
                 hTime->Fill(time);
                 hPhiVsThEta->Fill(pos.Phi(),pos.Theta()*TMath::RadToDeg());
                 hPhiVsZ->Fill(pos.Phi(),z2);
                 hYVsZ->Fill(y2,z2);
                 hXVsYVsZ->Fill(x2,y2,z2);
                 hMultiplicityVsCell->Fill(detID);
         // }

          // } //Particle selection

  } //Primary particles condition

  
 } //MBB Points

    
/////////////////////Loop for particles MCTrack//////////////////////////////////////////////////////////
  /*
      for(Int_t k = 0;  k< nMCTracks; k++){
        

         TVector3 mcMom;
        
         FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(k);
         
         mcTr->GetMomentum(mcMom);

           if( mcTr->GetMotherId() < 0 ) {  //Primary particles condition

              if (  TMath::Abs(mcTr->GetPdgCode() )  != 2212  &&  //Protons
                 TMath::Abs(mcTr->GetPdgCode() )  != 3222 &&  //Sigma+
                 TMath::Abs(mcTr->GetPdgCode() )  != 3112 &&  //Sigma-
                 TMath::Abs(mcTr->GetPdgCode() )  != 3312 &&  //Xi-
                    TMath::Abs(mcTr->GetPdgCode() )  != 211  && //Pions
                    TMath::Abs(mcTr->GetPdgCode() )  != 11   && //Electrons
                    TMath::Abs(mcTr->GetPdgCode() )  != 321  && //Kaons
                    TMath::Abs(mcTr->GetPdgCode() )  != 13   && //Muons
                 TMath::Abs(mcTr->GetPdgCode() )  != 3334 ) continue; //Omega-

                     if ( TMath::Abs(mcTr->GetPdgCode()) == 11) { //Particles condition
                       
                      if( TMath::Abs(mcMom.Eta()) < 13) {  
                         hEta->Fill(mcMom.Eta());  
                         hPt->Fill(mcMom.Pt());     
                         hEnergy->Fill(mcTr->GetEnergy()); 

                       cout<<"pdgcode: "<< mcTr->GetPdgCode()<<endl; 

                     }

                     } //Particles condition

                   } //Primary particles condition

               }  //End loop for particles mctrack                
 */  
////////////////////////////////////////////////////////////////////////////////////////////////////////
  

        
    }
} //End loop of events

 ///////////////////////////////////////// 


 
///////////////////////////////////////////////  
   cout<<"Saving histograms"<<endl;
     
   fileOutput->mkdir("Datos");
   fileOutput->cd("Datos");
   hPt->Write();
   hEnergy->Write();
   hEta->Write();
   hXY->Write();
   hTime->Write();
   hPhiVsThEta->Write();
   hMultiplicityVsCell->Write();
   hPhiVsZ->Write();
   hYVsZ->Write();
   hXVsYVsZ->Write();
 

  cout<<"End histograms"<<endl;
  return 0;
 
}

