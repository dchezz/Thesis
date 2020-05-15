
// Description:
//      
//       Analysis macro. This macro is to study several distributions such as pT, energy, and multiplicity.
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
#include "MpdFsa.h"
#include "MpdBbc.h"
#include "MpdCpc.h"
#include "MpdStrawendcap.h"
#include "MpdGetNumEvents.h"

#include <iostream>

#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>
 #include <TPad.h>
  #include <TH1F.h>
  #include <THD.h>


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
        
Int_t readOutput(TString inputFile="~/NICA/mpdroot/macro/mpd/bb_phsd_10000.root",TString outputFile="results_bb_phsd_10000.root")
//Int_t readOutput(TString inputFile="/home/luis/Analysis-BEBE/10000AuAu11GeV16fm.root", TString outputFile="/home/luis/Analysis-BEBE/salida-evetest.root")
{
    gROOT->Reset();
    gStyle->SetOptDate(0);      //show day and time
    gStyle->SetOptStat(1);      //show statistic
    gStyle->SetPalette(1,0);
    
    
    Double_t epsilon = 0.00001;
    
    mpdloadlibs();
    

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////MPDROOT Analysis//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////Primary charged particles bmd ///////////////////////////////
    TH1D *hEtaPrimaryBMD = new TH1D("hEtaPrimaryBMD","@Eta primary charged particles BMD . 10000 Au+Au @11GeV PHSD. MPDROOT.",240,-8,8);
    hEtaPrimaryBMD->SetXTitle("#eta");
    hEtaPrimaryBMD->SetYTitle("#frac{1}{#eta}#frac{dN}{d#eta}");

    TH1D *hPtPrimaryBMD = new TH1D("hPtPrimaryBMD","Pt primary charged particles. BMD  . 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryBMD->SetXTitle("p_T (GeV/#it{c})");
    hPtPrimaryBMD->SetYTitle("entries");

     ///////////////////////////////////////Pt vs Pseudorapidity////////////////////////////////////////////////////////////////
  TH2D *hBmdPtVsEta  = new TH2D("hBmdPtVsEta","BMD pT vs Pseudorapidity. 10000 Au+Au @11GeV PHSD. MPDROOT.",240,-8,8,60,-0,1.5);
  hBmdPtVsEta->SetXTitle("#Eta");
  hBmdPtVsEta->SetYTitle("p_T (GeV/#it{c})");
   
   TH1D *hEnergyPrimaryBMD = new TH1D("hEnergyPrimaryBMD","Energy",1000,0,10);
    hEnergyPrimaryBMD->SetXTitle("Energy (GeV)");

    TH1D *htimeBMD = new TH1D("htimeBMD","Time of flight. BMD  . 10000 Au+Au @11GeV PHSD. MPDROOT.",100,0,10);
    htimeBMD->SetXTitle("t(ns)");
    htimeBMD->SetYTitle("entries");
//////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////Primary charged particles BMD side A per ring/////////////////////////
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////Eta /////////////////////////////////////////////////////////
    TH1D *hEtaPrimaryChargedBMDARing1 = new TH1D("hEtaPrimaryChargedBMDARing1","Charged particles @Eta BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDARing1->SetXTitle("#eta");
    
    TH1D *hEtaPrimaryChargedBMDARing2 = new TH1D("hEtaPrimaryChargedBMDARing2","Charged particles @Eta BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDARing2->SetXTitle("#eta");
   
    TH1D *hEtaPrimaryChargedBMDARing3 = new TH1D("hEtaPrimaryChargedBMDARing3","Charged particles @Eta BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDARing3->SetXTitle("#eta");
  
    TH1D *hEtaPrimaryChargedBMDARing4 = new TH1D("hEtaPrimaryChargedBMDARing4","Charged particles @Eta BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDARing4->SetXTitle("#eta");
  
    TH1D *hEtaPrimaryChargedBMDARing5 = new TH1D("hEtaPrimaryChargedBMDARing5","Charged particles @Eta BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDARing5->SetXTitle("#eta");
   
    TH1D *hEtaPrimaryChargedBMDARing6 = new TH1D("hEtaPrimaryChargedBMDARing6","Charged particles @Eta BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDARing6->SetXTitle("#eta");

///////////////////////////////////////////////Pt/////////////////////////////////////////////////////////////////
    TH1D *hPtPrimaryChargedBMDARing1 = new TH1D("hPtPrimaryChargedBMDARing1","Charged particles Pt BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDARing1 ->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDARing2 = new TH1D("hPtPrimaryChargedBMDARing2","Charged particles Pt BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDARing2->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDARing3 = new TH1D("hPtPrimaryChargedBMDARing3","Charged particles Pt BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDARing3->SetXTitle("p_T (GeV/#it{c})");
  
    TH1D *hPtPrimaryChargedBMDARing4 = new TH1D("hPtPrimaryChargedBMDARing4","Charged particles Pt BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDARing4->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDARing5 = new TH1D("hPtPrimaryChargedBMDARing5","Charged particles Pt BMDA 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDARing5->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDARing6 = new TH1D("hPtPrimaryChargedBMDARing6","Charged particles Pt BMDA 5000 Au+Au @9GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDARing6->SetXTitle("p_T (GeV/#it{c})");
    
    TH1D *hMultVsCellID = new TH1D("hMultVsCellID", " ",162,-1,163);
    hMultVsCellID->SetXTitle("Cell ID");
    hMultVsCellID->SetYTitle("Multiplicity");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////Primary charged particles BMD side C per ring/////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    ////////////////////////////////////////eta /////////////////////////////////////////////////////////
    TH1D *hEtaPrimaryChargedBMDCRing1 = new TH1D("hEtaPrimaryChargedBMDCRing1","@Eta BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDCRing1->SetXTitle("#eta");
    
    TH1D *hEtaPrimaryChargedBMDCRing2 = new TH1D("hEtaPrimaryChargedBMDCRing2","@Eta BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDCRing2->SetXTitle("#eta");
   
    TH1D *hEtaPrimaryChargedBMDCRing3 = new TH1D("hEtaPrimaryChargedBMDCRing3","@Eta BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDCRing3->SetXTitle("#eta");
  
    TH1D *hEtaPrimaryChargedBMDCRing4 = new TH1D("hEtaPrimaryChargedBMDCRing4","@Eta BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDCRing4->SetXTitle("#eta");
  
    TH1D *hEtaPrimaryChargedBMDCRing5 = new TH1D("hEtaPrimaryChargedBMDCRing5","@Eta BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDCRing5->SetXTitle("#eta");
   
    TH1D *hEtaPrimaryChargedBMDCRing6 = new TH1D("hEtaPrimaryChargedBMDCRing6","@Eta BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",200,-5,5);
    hEtaPrimaryChargedBMDCRing6->SetXTitle("#eta");

///////////////////////////////////////////////Pt/////////////////////////////////////////////////////////////////
    TH1D *hPtPrimaryChargedBMDCRing1 = new TH1D("hPtPrimaryChargedBMDCRing1","Pt BMDC 10000 Au+Au @11GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDCRing1 ->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDCRing2 = new TH1D("hPtPrimaryChargedBMDCRing2","Pt BMDC 5000 Au+Au @9GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDCRing2->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDCRing3 = new TH1D("hPtPrimaryChargedBMDCRing3","Pt BMDC 5000 Au+Au @9GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDCRing3->SetXTitle("p_T (GeV/#it{c})");
  
    TH1D *hPtPrimaryChargedBMDCRing4 = new TH1D("hPtPrimaryChargedBMDCRing4","Pt BMDC 5000 Au+Au @9GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDCRing4->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDCRing5 = new TH1D("hPtPrimaryChargedBMDCRing5","Pt BMDC 5000 Au+Au @9GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDCRing5->SetXTitle("p_T (GeV/#it{c})");

    TH1D *hPtPrimaryChargedBMDCRing6 = new TH1D("hPtPrimaryChargedBMDCRing6","Pt BMDC 5000 Au+Au @9GeV PHSD. MPDROOT.",60,0,1.5);
    hPtPrimaryChargedBMDCRing6->SetXTitle("p_T (GeV/#it{c})");
 
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////Monte Carlo Histograms/////////////////////////////////////////////////////////
///////////////////////////////////////////////////MC Track////////////////////////////////////////////////////////////////

  ///////////////////////////////////////Primary charged particles in bmd region//////////////////////////////////
    TH1D *hMCTrackEtaPrimaryBMD = new TH1D("hMCTrackEtaPrimaryBMD","@Eta primary charged particles BMD A & C. 5000 Au+Au @9GeV PHSD. MCTrack.",160,-8,8);
    hMCTrackEtaPrimaryBMD->SetXTitle("#eta");
    hMCTrackEtaPrimaryBMD->SetYTitle("entries");

    TH1D *hMCTrackPtPrimaryBMD = new TH1D("hMCTrackPtPrimaryBMD","Pt primary charged particles. BMD A & C. 5000 Au+Au @9GeV PHSD. MCTrack.",60,0,1.5);
    hMCTrackPtPrimaryBMD->SetXTitle("p_T (GeV/#it{c})");
    hMCTrackPtPrimaryBMD->SetXTitle("entries");

  ///////////////////////////////////////MC all space//////////////////////////////////
    TH1D *hMCTrackPtPrimary = new TH1D("hMCTrackPtPrimary","Pt primary charged particles. BMD A & C. 5000 Au+Au @9GeV PHSD. MCTrack.",60,0,1.5);
    hMCTrackPtPrimary->SetXTitle("p_T (GeV/#it{c})");
    hMCTrackPtPrimary->SetXTitle("entries");

 ///////////////////////////////////////Pt vs Pseudorapidity////////////////////////////////////////////////////////////////
  TH2D *hBmdPtVsEtaMCTrack  = new TH2D("hBmdPtVsEtaMCTrack","BMD pT vs Pseudorapidity. 5000 Au+Au @9GeV PHSD. MCTrack.",240,-8,8,60,-0,1.5);
  hBmdPtVsEtaMCTrack->SetXTitle("#Eta");
  hBmdPtVsEtaMCTrack->SetYTitle("p_T (GeV/#it{c})");

    TH1D *hMCTrackEnergyPrimary = new TH1D("hMCTrackEnergyPrimary","Energy",100,0,10);
    hMCTrackEnergyPrimary->SetXTitle("Energy (GeV)");

    TH1D *hMCTrackEnergyCell = new TH1D("hMCTrackEnergyCell","Energy in cell",100,0,10);
    hMCTrackEnergyCell->SetXTitle("Energy (GeV)");

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////Multiplicity charged particles////////////////////////////////////////////////////////////// 
    TH1D *hBmdChargedMultiplicity = new TH1D("hBmdChargedMultiplicity","Charged particles multiplicity BEBE. 10000 Au+Au @11GeV PHSD. MPDROOT.",1500,0,150);
   hBmdChargedMultiplicity->SetXTitle("Multiplicity");  

    //Defining pointers to data
    
    TFile fileInput(inputFile.Data());
    
    TTree *simEvent = (TTree*) fileInput.Get("mpdsim");
    
    //Arrays/////nombre de la variable///////////nombre en la memoria interna
    TClonesArray *bmdPoints = (TClonesArray*) fileInput.FindObjectAny("BmdPoint");
    simEvent->SetBranchAddress("BmdPoint",&bmdPoints);

   // TClonesArray *bbcPoints = (TClonesArray*) fileInput.FindObjectAny("BBCPoint");
   // simEvent->SetBranchAddress("BBCPoint",&bbcPoints);
    
    TClonesArray* mcTracks = (TClonesArray*) fileInput.FindObjectAny("MCTrack");
    simEvent->SetBranchAddress("MCTrack",&mcTracks);
  
  //  TClonesArray* tpcPoints = (TClonesArray*) fileInput.FindObjectAny("TpcPoint");
  //  simEvent->SetBranchAddress("TpcPoint",&tpcPoints);
  
    FairMCEventHeader* mcHeader = (FairMCEventHeader*) fileInput.FindObjectAny("MCEventHeader.");
    simEvent->SetBranchAddress("MCEventHeader.",&mcHeader);
       
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


 Int_t nChargedMultiplicityBMD=0;

//MCTrack kinematics
//MCEventHeader  kind of collision

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

     MpdMCEventHeader *extraEventHeader = dynamic_cast<MpdMCEventHeader*> (mcHeader); /*the hell is this*/

    // impactParameter = extraEventHeader->GetB();
    // bmdTree->Fill() 
    // cout<<"Impact parameter: "<<extraEventHeader->GetB()<<endl;
       
    Int_t nMCTracks  =  mcTracks->GetEntriesFast();
  //  cout<<"nMCTracks "<<nMCTracks<<endl;
   
    if (bmdPoints != 0 ) {

    Int_t nbmdPoints =  bmdPoints->GetEntriesFast();
    if( nbmdPoints < 1 ) continue;

    /*   
       if (bbcPoints != 0 ) {
       Int_t nbbcPoints =  bbcPoints->GetEntriesFast();
       if( nbbcPoints < 1 ) continue;
  */      
    //cout<<"nbmdPoints: "<<nbmdPoints<<endl;
       
    //Loop for particles in both BMD.
    Int_t nRecTracks = 0;
    Int_t nRecPrimaryTracks = 0; 

       for (Int_t j = 0; j < nbmdPoints; j++) {     
       BmdPoint *p1    = (BmdPoint*) bmdPoints->At(j);

//   for (Int_t j = 0; j < nbbcPoints; j++) {
//    MpdBbcPoint *p1    = (MpdBbcPoint*) bbcPoints->At(j);

      TVector3 recMom;
      Int_t    trkID = p1->GetTrackID();
      Double_t time  = p1->GetTime();
      Double_t eLoss = p1->GetEnergyLoss();
      p1->Momentum(recMom);
      Double_t currPt  = recMom.Pt();
      Double_t currEta = recMom.Eta();
      
      Int_t detID =  p1->GetDetectorID();

      //  cout<<"trkID: "<<trkID<<endl;  
      //  cout<<"detID: "<<detID<<endl;            

     Int_t x = p1->GetX();
     Int_t y = p1->GetY();
     Int_t z = p1->GetZ(); 
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////Analysis per ring////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Int_t ring = -1;
  
//ring = GetRing(detID);
 //cout<<"ring: "<<ring<<endl;

 MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(trkID);
   
    Double_t Xv = mcTr->GetStartX();
    Double_t Yv = mcTr->GetStartY();
    Double_t Zv = mcTr->GetStartZ();

   //  cout<<"x: "<<x<<"y: "<<y<<"z: "<<z<<endl; 

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////Primary Particles Condition////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////////
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

                 cout<<"pdgcode: "<< mcTr->GetPdgCode()<<endl; 



               hMultVsCellID->Fill(detID);
             
             
   //Charged primary particles
         if ( TMath::Abs(mcTr->GetPdgCode() )  != 211  && //Pions
                  TMath::Abs(mcTr->GetPdgCode() )  != 11  && //Electrons
	                TMath::Abs(mcTr->GetPdgCode() )  != 2212 && //Protons
                 TMath::Abs(mcTr->GetPdgCode() )  != 321  && //Kaons
                 TMath::Abs(mcTr->GetPdgCode() )  != 13  && //Muons
	      (TMath::Abs(Xv) > 0.001 || TMath::Abs(Yv) > 0.001 || TMath::Abs(Zv) > 0.01 )) continue;
                            // cout<<"Xv: "<<Xv<<"Yv: "<<Yv<<"Zv: "<<Zv<<endl; 
                            // cout<<"PdgCode: "<<mcTr->GetPdgCode()<<endl;       
                            // cout<<"Detector ID: "<<p1->GetDetectorID()<<endl;       

       htimeBMD->Fill(p1->GetTime()); //Time of flight
      // cout<<"time: "<<p1->GetTime()<<endl;   

//////////////////////////////////BE-BE cells////////////////////////////////////////////
if (p1->GetDetectorID() ==1) {
    hMCTrackEnergyCell->Fill(mcTr->GetEnergy());
}
//////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Primary Charged Particles BMD///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
if(z >=199){ //BMD A


// Rings are turn off
//Ring 1
//  if( ring==1 && z >=199 ){       
// cout<<"z: "<<z<<endl;   
/////////////////////////////////Eta region ring 1///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  3.5 && TMath::Abs(recMom.Eta()) < 4.2 ) {
              // cout<<"PdgCode: "<<mcTr->GetPdgCode()<<endl;       
                                
              //Charged primary particles
              hEtaPrimaryChargedBMDARing1->Fill(recMom.Eta());
              hPtPrimaryChargedBMDARing1->Fill(recMom.Pt());  
         
         //     hEtaPrimaryBMD->Fill(recMom.Eta());     
        //      hPtPrimaryBMD->Fill(recMom.Pt());    
  }    

  //Ring 2
 //  if( ring==2 && z >=199){
 //  cout<<"z: "<<z<<endl;    
/////////////////////////////////Eta region ring 2///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  3.2 && TMath::Abs(recMom.Eta()) < 3.5 ) {                 


         //Charged primary particles
         hEtaPrimaryChargedBMDARing2->Fill(recMom.Eta());  
         hPtPrimaryChargedBMDARing2->Fill(recMom.Pt()); 

        //  hEtaPrimaryBMD->Fill(recMom.Eta());     
        //  hPtPrimaryBMD->Fill(recMom.Pt());        
   } 

//Ring 3
 //  if( ring==3 && z >=199){  
/////////////////////////////////Eta region ring 3///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  3.0 && TMath::Abs(recMom.Eta()) < 3.2 ) {      

     
         //Charged primary particles
         hEtaPrimaryChargedBMDARing3->Fill(recMom.Eta());
         hPtPrimaryChargedBMDARing3->Fill(recMom.Pt());  
         
        //  hEtaPrimaryBMD->Fill(recMom.Eta());     
       //   hPtPrimaryBMD->Fill(recMom.Pt());           
   }    


 //Ring 4
 //  if( ring==4 && z >=199){
/////////////////////////////////Eta region ring 4///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  2.8 && TMath::Abs(recMom.Eta()) < 3.0 ) {      
    
  
         //Charged primary particles
         hEtaPrimaryChargedBMDARing4->Fill(recMom.Eta());
         hPtPrimaryChargedBMDARing4->Fill(recMom.Pt());  

        //  hEtaPrimaryBMD->Fill(recMom.Eta());     
        //  hPtPrimaryBMD->Fill(recMom.Pt());  
     }    


//Ring 5
//  if( ring==5 && z >=199){       
/////////////////////////////////Eta region ring 5///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  2.6 && TMath::Abs(recMom.Eta()) < 2.8 ) {        

         //Charged primary particles
         hEtaPrimaryChargedBMDARing5->Fill(recMom.Eta());
         hPtPrimaryChargedBMDARing5->Fill(recMom.Pt());  

        //  hEtaPrimaryBMD->Fill(recMom.Eta());     
      //    hPtPrimaryBMD->Fill(recMom.Pt());  
        
     }    

//Ring 6
//  if( ring==6 && z >=199){
/////////////////////////////////Eta region ring 6///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  2.2 && TMath::Abs(recMom.Eta()) < 2.6 ) {        
        
         //Charged primary particles
         hEtaPrimaryChargedBMDARing6->Fill(recMom.Eta());
         hPtPrimaryChargedBMDARing6->Fill(recMom.Pt());  

       //   hEtaPrimaryBMD->Fill(recMom.Eta());     
       //   hPtPrimaryBMD->Fill(recMom.Pt());            
     }         

} //BMD A Z
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////Primary Charged Particles BMD C/////////////////////////
////////////////////////////////////////////////////////////////////////////////////////    

/////////////////////////////////Eta region ring 1/////////////////////////////////////////// 
 if( TMath::Abs(recMom.Eta()) >  3.5 && TMath::Abs(recMom.Eta()) < 4.2 ) {            
     //Ring 1
 //   if( ring==1 && z <=-199){       
        if( z <=-199){
  
         //Charged primary particles
         hEtaPrimaryChargedBMDCRing1->Fill(recMom.Eta());
         hPtPrimaryChargedBMDCRing1->Fill(recMom.Pt()); 

         // hEtaPrimaryBMD->Fill(recMom.Eta());     
          // hPtPrimaryBMD->Fill(recMom.Pt());  
        }    
}        
/////////////////////////////////Eta region ring 2///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  3.2 && TMath::Abs(recMom.Eta()) < 3.5 ) {                 
                 //Ring 2
  //  if( ring==2 && z <=-199){
               if( z <=-199){

 
         //Charged primary particles
         hEtaPrimaryChargedBMDCRing2->Fill(recMom.Eta());
         hPtPrimaryChargedBMDCRing2->Fill(recMom.Pt());  

         // hEtaPrimaryBMD->Fill(recMom.Eta());     
        //  hPtPrimaryBMD->Fill(recMom.Pt());  
       }    
}
/////////////////////////////////Eta region ring 3///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  3.0 && TMath::Abs(recMom.Eta()) < 3.2 ) {      
                      //Ring 3
   // if( ring==3 && z <=-199){       
        if( z <=-199){
  
         //Charged primary particles
         hEtaPrimaryChargedBMDCRing3->Fill(recMom.Eta());
         hPtPrimaryChargedBMDCRing3->Fill(recMom.Pt());  

         // hEtaPrimaryBMD->Fill(recMom.Eta());     
        //  hPtPrimaryBMD->Fill(recMom.Pt());    
     }   
 }     

/////////////////////////////////Eta region ring 4///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  2.8 && TMath::Abs(recMom.Eta()) < 3.0 ) {      
                      //Ring 4
    //if( ring==4 && z <=-199){
            if( z <=-199){

         //Charged primary particles
         hEtaPrimaryChargedBMDCRing4->Fill(recMom.Eta());
         hPtPrimaryChargedBMDCRing4->Fill(recMom.Pt());  

         // hEtaPrimaryBMD->Fill(recMom.Eta());     
        //  hPtPrimaryBMD->Fill(recMom.Pt());   
     }    
}

/////////////////////////////////Eta region ring 5///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  2.6 && TMath::Abs(recMom.Eta()) < 2.8 ) {        
                      //Ring 5
   // if( ring==5 && z <=-199){
           if( z <=-199){       

         //Charged primary particles
         hEtaPrimaryChargedBMDCRing5->Fill(recMom.Eta());
         hPtPrimaryChargedBMDCRing5->Fill(recMom.Pt());  

       //  hEtaPrimaryBMD->Fill(recMom.Eta());     
        //  hPtPrimaryBMD->Fill(recMom.Pt());    
     }    
}

/////////////////////////////////Eta region ring 6 ///////////////////////////////////////////
 if( TMath::Abs(recMom.Eta()) >  2.2 && TMath::Abs(recMom.Eta()) < 2.6 ) {        
                      //Ring 6
   // if( ring==6 && z <=-199){
        if( z <=-199){
  
   
         //Charged primary particles
         hEtaPrimaryChargedBMDCRing6->Fill(recMom.Eta());
         hPtPrimaryChargedBMDCRing6->Fill(recMom.Pt());  

     //       hEtaPrimaryBMD->Fill(recMom.Eta());     
    //        hPtPrimaryBMD->Fill(recMom.Pt());    
    }    
 }   


////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////BMD A and C region////////////////////////////////////////
     if( TMath::Abs(recMom.Eta()) >  2.2 && TMath::Abs(recMom.Eta()) < 4.2 ) {  

   //        if( TMath::Abs(recMom.Eta()) >=  3.0 && TMath::Abs(recMom.Eta()) <= 4.2 ) {   //Rings 4, 5 and 6.
                      
                      hEnergyPrimaryBMD->Fill(p1->GetEnergyLoss());
                      hEtaPrimaryBMD->Fill(recMom.Eta());     
                      hPtPrimaryBMD->Fill(recMom.Pt());   
                      hBmdPtVsEta->Fill(recMom.Eta(),recMom.Pt());
                     

                      nChargedMultiplicityBMD++; //BMD multiplicity
                      
    // cout<<"EnergyLoss:" <<p1->GetEnergyLoss()<<endl;
                 } 
 ////////////////////////////////////////////////////////////////////////////////////////////////////

     }  //Primary Particles Condition 

   } 

  }  //bmdpoints 


///////////////////////////////////////////////////////////////////////  

////////////////////////////////Monte Carlo////////////////////////////////////////////////
//////////////////////////////////MCTrack/////////////////////////////////////////////////


    //Loop for particles
    for(Int_t k = 0;  k< nMCTracks; k++){
        
         TVector3 mcMom;
        
         MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(k);
         
         mcTr->GetMomentum(mcMom);

           Double_t Xv = mcTr->GetStartX();
           Double_t Yv = mcTr->GetStartY();
           Double_t Zv = mcTr->GetStartZ();

  if( mcTr->GetMotherId() < 0 ) {
             

         //Charged primary particles

        if ( TMath::Abs(mcTr->GetPdgCode() )  != 211  &&
                  TMath::Abs(mcTr->GetPdgCode() )  != 11   &&
	                TMath::Abs(mcTr->GetPdgCode() )  != 2212 && //Protons
                  TMath::Abs(mcTr->GetPdgCode() )  != 321  &&
                  TMath::Abs(mcTr->GetPdgCode() )  != 13   &&
	     (TMath::Abs(Xv) > 0.001 || TMath::Abs(Yv) > 0.001 || TMath::Abs(Zv) > 0.01 )) continue;
  
         //Charged primary particles
         //hMCTrackEtaPrimary->Fill(mcMom.Eta());
         // hMCTrackPtPrimary->Fill(mcMom.Pt());   //All space

            if( TMath::Abs(mcMom.Eta()) >  2.2 && TMath::Abs(mcMom.Eta()) < 4.2) {          //Pseudorapidity region BE-BE   

	// if( TMath::Abs(mcMom.Eta()) > 3.0 && TMath::Abs(mcMom.Eta()) < 4.2 ) {   //Rings 4, 5 and 6.
         
                     hMCTrackEtaPrimaryBMD->Fill(mcMom.Eta());
                     hMCTrackPtPrimaryBMD->Fill(mcMom.Pt());  

                     hBmdPtVsEtaMCTrack->Fill(mcMom.Eta(),mcMom.Pt());
                     hMCTrackEnergyPrimary->Fill(mcTr->GetEnergy());

                } 
                   
       }

     }
///////////////////////////////////////////////////////////////////////  

///////////////////////Multiplicity ////////////////////////////////////  
    
   //  cout<<"multiplicity:" <<nChargedMultiplicityBMD<<endl;
/////////////////////////////////////////////////////////////////////// 

 } //End of the loop's event. 
  
/*
///////////////////////Multiplicity ////////////////////////////////////  
     hBmdChargedMultiplicity->Fill(nChargedMultiplicityBMD);
     cout<<"multiplicity:" <<nChargedMultiplicityBMD<<endl;
/////////////////////////////////////////////////////////////////////// 
*/

 // bmdTree->Write();

 
///////////////////////////////////////////////  
   cout<<"Saving histograms"<<endl;
     
   fileOutput->mkdir("Mc");
   fileOutput->cd("Mc");
   hMultVsCellID->Write();

  gStyle->SetOptTitle(0); //No title for histograms

  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////MPDROOT Histograms////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   TCanvas *c14 = new TCanvas("c14","Time of flight primary charged particles BMD. MPDROOT.");
   c14->cd(1);
   htimeBMD->SetStats(kFALSE);
   htimeBMD->GetXaxis()->SetTitle("t(ns)");
   htimeBMD->GetYaxis()->SetTitle("1/N dN/dt");
   htimeBMD->SetStats(kFALSE);
   htimeBMD->Scale(1./htimeBMD->Integral());
   htimeBMD->Draw();

//////////////////////////Eta and pT primary charged particles BMD A and C/////////////////////////////////////////
   TCanvas *c15 = new TCanvas("c15","Pt primary charged particles BMD. MPDROOT.");
   c15->cd(1);
   hPtPrimaryBMD->SetStats(kFALSE);
   hPtPrimaryBMD->GetXaxis()->SetTitle("p_T (GeV/#it{c})");
   hPtPrimaryBMD->GetYaxis()->SetTitle("1/N dN/dpT");
   hPtPrimaryBMD->SetStats(kFALSE);
   TH1D *hPtPrimaryBMDEfficiency = (TH1D*) hPtPrimaryBMD->Clone();
   hPtPrimaryBMD->Scale(1./hMCTrackPtPrimaryBMD->Integral());
   hPtPrimaryBMD->Draw();


   TCanvas *c16 = new TCanvas("c16","@Eta primary charged particles BMD. MPDROOT.");
   c16->cd(1);
   hEtaPrimaryBMD->SetStats(kFALSE);
   hEtaPrimaryBMD->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryBMD->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryBMD->Scale(1./hEtaPrimaryBMD->Integral());
   hEtaPrimaryBMD->Draw();

/*
   TCanvas *c18 = new TCanvas("c18","Pt vs Eta primary charged particles BMD. MPDROOT.");
   c18->cd(1);
   hBmdPtVsEta->GetXaxis()->SetTitle("#eta");
   hBmdPtVsEta->GetYaxis()->SetTitle("p_T (GeV/#it{c})");
   hBmdPtVsEta->Scale(1./hBmdPtVsEta->Integral());
//   hBmdPtVsEta->Scale(1./hBmdPtVsEta->Integral());
   hBmdPtVsEta->Draw();
*/

     TCanvas *c19 = new TCanvas("c19","@Eta primary charged particles BMD region.  MPDROOT.");
   hEnergyPrimaryBMD->SetStats(kFALSE);
   hEnergyPrimaryBMD->GetXaxis()->SetTitle("Energy (GeV)");
   hEnergyPrimaryBMD->GetYaxis()->SetTitle("1/N dN/dE");
   hEnergyPrimaryBMD->Scale(1./hEnergyPrimaryBMD->Integral());
   hEnergyPrimaryBMD->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////////////Analysis per ring side A/////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
/////////////////////Eta primary charged particles///////////////////////////////////////////////////////

   TCanvas *c23 = new TCanvas("c23","Pseudorapidity per ring BMD side A. MPDROOT");
   c23->Divide(3,2);
   c23->cd(1); 
   hEtaPrimaryChargedBMDARing1->SetStats(kFALSE);
   hEtaPrimaryChargedBMDARing1->GetXaxis()->SetTitle("#eta");
   hEtaPrimaryChargedBMDARing1->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryChargedBMDARing1->Scale(1./(hEtaPrimaryBMD->GetEntries()));
   hEtaPrimaryChargedBMDARing1->Draw();

   c23->cd(2);
   hEtaPrimaryChargedBMDARing2->SetStats(kFALSE);
   hEtaPrimaryChargedBMDARing2->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing2->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryChargedBMDARing2->Scale(1./(hEtaPrimaryBMD->GetEntries()));
   hEtaPrimaryChargedBMDARing2->Draw();

   c23->cd(3);
   hEtaPrimaryChargedBMDARing3->SetStats(kFALSE);
   hEtaPrimaryChargedBMDARing3->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing3->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryChargedBMDARing3->Scale(1./(hEtaPrimaryBMD->GetEntries()));
   hEtaPrimaryChargedBMDARing3->Draw();

   c23->cd(4);
   hEtaPrimaryChargedBMDARing4->SetStats(kFALSE);
   hEtaPrimaryChargedBMDARing4->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing4->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryChargedBMDARing4->Scale(1./(hEtaPrimaryBMD->GetEntries()));
   hEtaPrimaryChargedBMDARing4->Draw();

   c23->cd(5);
   hEtaPrimaryChargedBMDARing5->SetStats(kFALSE);
   hEtaPrimaryChargedBMDARing5->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing5->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryChargedBMDARing5->Scale(1./(hEtaPrimaryBMD->GetEntries()));
   hEtaPrimaryChargedBMDARing5->Draw();
   
   c23->cd(6);
   hEtaPrimaryChargedBMDARing6->SetStats(kFALSE);
   hEtaPrimaryChargedBMDARing6->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing6->GetYaxis()->SetTitle("1/N dN/d#eta");
   hEtaPrimaryChargedBMDARing6->Scale(1./(hEtaPrimaryBMD->GetEntries()));
   hEtaPrimaryChargedBMDARing6->Draw();

   /////////////////////////////////////////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////
 
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////Pt charged primary particles BMD A/////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////   

 TCanvas *c25 = new TCanvas("c25","Pt primary charged particles");
    c25->Divide(3,2);
    c25->cd(1);
    hPtPrimaryChargedBMDARing1->SetStats(kFALSE);
    hPtPrimaryChargedBMDARing1->GetXaxis()->SetTitle("Pt");
    hPtPrimaryChargedBMDARing1->GetYaxis()->SetTitle("1/N dN/dpT");
    hPtPrimaryChargedBMDARing1->Scale(1./(hPtPrimaryBMD->GetEntries()));
    hPtPrimaryChargedBMDARing1->Draw();

    c25->cd(2);
    hPtPrimaryChargedBMDARing2->SetStats(kFALSE);
    hPtPrimaryChargedBMDARing2->GetXaxis()->SetTitle("Pt");
    hPtPrimaryChargedBMDARing2->GetYaxis()->SetTitle("1/N dN/dpT");
    hPtPrimaryChargedBMDARing2->Scale(1./(hPtPrimaryBMD->GetEntries()));
    hPtPrimaryChargedBMDARing2->Draw();

    c25->cd(3);
    hPtPrimaryChargedBMDARing3->SetStats(kFALSE);
    hPtPrimaryChargedBMDARing3->GetXaxis()->SetTitle("Pt");
    hPtPrimaryChargedBMDARing3->GetYaxis()->SetTitle("1/N dN/dpT");
    hPtPrimaryChargedBMDARing3->Scale(1./(hPtPrimaryBMD->GetEntries()));
    hPtPrimaryChargedBMDARing3->Draw();
  
    c25->cd(4);
    hPtPrimaryChargedBMDARing4->SetStats(kFALSE);
    hPtPrimaryChargedBMDARing4->GetXaxis()->SetTitle("Pt");
    hPtPrimaryChargedBMDARing4->GetYaxis()->SetTitle("1/N dN/dpT");
    hPtPrimaryChargedBMDARing4->Scale(1./(hPtPrimaryBMD->GetEntries()));
    hPtPrimaryChargedBMDARing4->Draw();

    c25->cd(5);
    hPtPrimaryChargedBMDARing5->SetStats(kFALSE);
    hPtPrimaryChargedBMDARing5->GetXaxis()->SetTitle("Pt");
    hPtPrimaryChargedBMDARing5->GetYaxis()->SetTitle("1/N dN/dpT");
    hPtPrimaryChargedBMDARing5->Scale(1./(hPtPrimaryBMD->GetEntries()));
    hPtPrimaryChargedBMDARing5->Draw();

    c25->cd(6);
    hPtPrimaryChargedBMDARing6->SetStats(kFALSE);
    hPtPrimaryChargedBMDARing6->GetXaxis()->SetTitle("Pt");
    hPtPrimaryChargedBMDARing6->GetYaxis()->SetTitle("1/N dN/dpT");
    hPtPrimaryChargedBMDARing6->Scale(1./(hPtPrimaryBMD->GetEntries()));
    hPtPrimaryChargedBMDARing6->Draw();

//////////////////////////////////////////////////////////////////////////////

///////////////////////////////Monte Carlo histograms////////////////////////////
///////////////////////Primary charged particles in Be-Be region/////////////////
/////////////////////////////////////////////////////////////////////////////////
   TCanvas *c12 = new TCanvas("c12","@Eta primary charged particles BMD.  MCTrack ");
   hMCTrackEtaPrimaryBMD->SetStats(kFALSE);
   hMCTrackEtaPrimaryBMD->GetXaxis()->SetTitle("Pseudorapidity");
   hMCTrackEtaPrimaryBMD->GetYaxis()->SetTitle("1/N dN/d#eta");
   hMCTrackEtaPrimaryBMD->Scale(1./hMCTrackEtaPrimaryBMD->Integral());
   hMCTrackEtaPrimaryBMD->Draw();

   TCanvas *c13 = new TCanvas("c13","Pt primary charged particles BMD. MCTrack");
   c13->cd(1);
   hMCTrackPtPrimaryBMD->SetStats(kFALSE);
   hMCTrackPtPrimaryBMD->GetXaxis()->SetTitle("p_T (GeV/#it{c})");
   hMCTrackPtPrimaryBMD->GetYaxis()->SetTitle("1/N dN/dpT");
   TH1D *hMCTrackPtPrimaryBMDEfficiency = (TH1D*) hMCTrackPtPrimaryBMD->Clone();
  // TH1F *h2 = (TH1F*) hist->Clone();
   hMCTrackPtPrimaryBMD->Scale(1./hMCTrackPtPrimaryBMD->Integral());
   hMCTrackPtPrimaryBMD->Draw();

   /*
      TCanvas *c14 = new TCanvas("c14","Pt vs Eta primary charged particles BMD.  MCTrack");
      c14->cd(1);
      hBmdPtVsEtaMCTrack->SetStats(kFALSE);
      hBmdPtVsEtaMCTrack->GetXaxis()->SetTitle("#Eta");
      hBmdPtVsEtaMCTrack->GetYaxis()->SetTitle("#pT");
      hBmdPtVsEtaMCTrack->Scale(1./hBmdPtVsEtaMCTrack->Integral());
      hBmdPtVsEtaMCTrack->Draw();
  */           

 ///////////////////////////Multiplicity-Efficiency////////////////////////////////////  
      TCanvas *c22 = new TCanvas("c22","Efficiency");
      c22->cd(1);
      hPtPrimaryBMDEfficiency->GetXaxis()->SetTitle("p_T (GeV/#it{c})");
      hPtPrimaryBMDEfficiency->GetYaxis()->SetTitle("Efficiency");
      hPtPrimaryBMDEfficiency->SetName("Efficiency");
      hPtPrimaryBMDEfficiency->Sumw2();
      hMCTrackPtPrimaryBMDEfficiency->Sumw2();
      hPtPrimaryBMDEfficiency->Rebin();
      hMCTrackPtPrimaryBMDEfficiency->Rebin();
      hPtPrimaryBMDEfficiency->Divide(hMCTrackPtPrimaryBMDEfficiency);
      hPtPrimaryBMDEfficiency->Draw();


   TCanvas *c54 = new TCanvas("c54","Charged particles multiplicity BEBE");
   hBmdChargedMultiplicity->SetStats(kFALSE);
   hBmdChargedMultiplicity->GetXaxis()->SetTitle("Multiplicity");
   hBmdChargedMultiplicity->GetYaxis()->SetTitle("Events number");
   hBmdChargedMultiplicity->Scale(1./hBmdChargedMultiplicity->GetEntries());
   hBmdChargedMultiplicity->Draw();
  

   TCanvas *c17 = new TCanvas("c17","Energy primary charged particles BMD region.  MCTrack.");
   hMCTrackEnergyPrimary->SetStats(kFALSE);
   hMCTrackEnergyPrimary->GetXaxis()->SetTitle("Energy (GeV)");
   hMCTrackEnergyPrimary->GetYaxis()->SetTitle("1/N dN/dE");
   hMCTrackEnergyPrimary->Scale(1./hMCTrackEnergyPrimary->Integral());
   hMCTrackEnergyPrimary->Draw();

      TCanvas *c20 = new TCanvas("c20","Energy primary charged particles BMD cell.  MCTrack.");
   hMCTrackEnergyCell->SetStats(kFALSE);
   hMCTrackEnergyCell->GetXaxis()->SetTitle("Energy (GeV)");
   hMCTrackEnergyCell->GetYaxis()->SetTitle("1/N dN/dE");
   hMCTrackEnergyCell->Scale(1./hMCTrackEnergyCell->Integral());
   hMCTrackEnergyCell->Draw();


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////    

  cout<<"End histograms"<<endl;
  return 0;


}

