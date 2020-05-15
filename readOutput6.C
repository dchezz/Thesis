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
#include "TpcDetector.h"
#include "MpdFsa.h"
#include "MpdBbc.h"
#include "MpdCpc.h"
#include "MpdStrawendcap.h"
#include "MpdGetNumEvents.h"
#include "MpdGetNumEvents.h"
#include <iostream>

#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <THD.h>
#include "Rtypes.h"

using namespace std;

#endif

#include "/home/dario/NICA/mpdroot/macro/mpd/mpdloadlibs.C"

/////////////////////////////////////////
//Función que elige los 20 anillos con 16 detectores cada uno
/*
Int_t GetRing(Int_t detID){
  
  if(detID == 1 || detID == 21 || detID == 41 || detID == 61 || detID == 81 || detID == 101 || detID == 121 || detID == 141 || detID == 161 || detID == 181 || detID == 201 || detID == 221 || detID == 241 || detID == 261 || detID == 281 || detID == 301) return 1;
  else if (detID == 2 || detID == 22 || detID == 42 || detID == 62 || detID == 82 || detID == 102 || detID == 122 || detID == 142 || detID == 162 || detID == 182 || detID == 202 || detID == 222 || detID == 242 || detID == 262 || detID == 282 || detID == 302) return 2;
  else if (detID == 3 || detID == 23 || detID == 43 || detID == 63 || detID == 83 || detID == 103 || detID == 123 || detID == 143 || detID == 163 || detID == 183 || detID == 203 || detID == 223 || detID == 243 || detID == 263 || detID == 283 || detID == 303) return 3;
  else if (detID == 4 || detID == 24 || detID == 44 || detID == 64 || detID == 84 || detID == 104 || detID == 124 || detID == 144 || detID == 164 || detID == 184 || detID == 204 || detID == 224 || detID == 244 || detID == 264 || detID == 284 || detID == 304) return 4;
  else if (detID == 5 || detID == 25 || detID == 45 || detID == 65 || detID == 85 || detID == 105 || detID == 125 || detID == 145 || detID == 165 || detID == 185 || detID == 205 || detID == 225 || detID == 245 || detID == 265 || detID == 285 || detID == 305) return 5;
  else if (detID == 6 || detID == 26 || detID == 46 || detID == 66 || detID == 86 || detID == 106 || detID == 126 || detID == 146 || detID == 166 || detID == 186 || detID == 206 || detID == 226 || detID == 246 || detID == 266 || detID == 286 || detID == 306) return 6;
  else if (detID == 7 || detID == 27 || detID == 47 || detID == 67 || detID == 87 || detID == 107 || detID == 127 || detID == 147 || detID == 167 || detID == 187 || detID == 207 || detID == 227 || detID == 247 || detID == 267 || detID == 287 || detID == 307) return 7;
  else if (detID == 8 || detID == 28 || detID == 48 || detID == 68 || detID == 88 || detID == 108 || detID == 128 || detID == 148 || detID == 168 || detID == 188 || detID == 208 || detID == 228 || detID == 248 || detID == 268 || detID == 288 || detID == 308) return 8;
  else if (detID == 9 || detID == 29 || detID == 49 || detID == 69 || detID == 89 || detID == 109 || detID == 129 || detID == 149 || detID == 169 || detID == 189 || detID == 209 || detID == 229 || detID == 249 || detID == 269 || detID == 289 || detID == 309) return 9;
  else if (detID == 10 || detID == 30 || detID == 50 || detID == 70 || detID == 90 || detID == 110 || detID == 130 || detID == 150 || detID == 170 || detID == 190 || detID == 210 || detID == 230 || detID == 250 || detID == 270 || detID == 290 || detID == 310) return 10;
  else if (detID == 11 || detID == 31 || detID == 51 || detID == 71 || detID == 91 || detID == 111 || detID == 131 || detID == 151 || detID == 171 || detID == 191 || detID == 211 || detID == 231 || detID == 251 || detID == 271 || detID == 291 || detID == 311) return 11;
  else if (detID == 12 || detID == 32 || detID == 52 || detID == 72 || detID == 92 || detID == 112 || detID == 132 || detID == 152 || detID == 172 || detID == 192 || detID == 212 || detID == 232 || detID == 252 || detID == 272 || detID == 292 || detID == 312) return 12;
  else if (detID == 13 || detID == 33 || detID == 53 || detID == 73 || detID == 93 || detID == 113 || detID == 133 || detID == 153 || detID == 173 || detID == 193 || detID == 213 || detID == 233 || detID == 253 || detID == 273 || detID == 293 || detID == 313) return 13;
  else if (detID == 14 || detID == 34 || detID == 54 || detID == 74 || detID == 94 || detID == 114 || detID == 134 || detID == 154 || detID == 174 || detID == 194 || detID == 214 || detID == 234 || detID == 254 || detID == 274 || detID == 294 || detID == 314) return 14;
  else if (detID == 15 || detID == 35 || detID == 55 || detID == 75 || detID == 95 || detID == 115 || detID == 135 || detID == 155 || detID == 175 || detID == 195 || detID == 215 || detID == 235 || detID == 255 || detID == 275 || detID == 295 || detID == 315) return 15;
  else if (detID == 16 || detID == 36 || detID == 56 || detID == 76 || detID == 96 || detID == 116 || detID == 136 || detID == 156 || detID == 176 || detID == 196 || detID == 216 || detID == 236 || detID == 256 || detID == 276 || detID == 296 || detID == 316) return 16;
  else if (detID == 17 || detID == 37 || detID == 57 || detID == 77 || detID == 97 || detID == 117 || detID == 137 || detID == 157 || detID == 177 || detID == 197 || detID == 217 || detID == 237 || detID == 257 || detID == 277 || detID == 297 || detID == 317) return 17;
  else if (detID == 18 || detID == 38 || detID == 58 || detID == 78 || detID == 98 || detID == 118 || detID == 138 || detID == 158 || detID == 178 || detID == 198 || detID == 218 || detID == 238 || detID == 258 || detID == 278 || detID == 298 || detID == 318) return 18;
  else if (detID == 19 || detID == 39 || detID == 59 || detID == 79 || detID == 99 || detID == 119 || detID == 139 || detID == 159 || detID == 179 || detID == 199 || detID == 219 || detID == 239 || detID == 259 || detID == 279 || detID == 299 || detID == 319) return 19;
  else if (detID == 20 || detID == 40 || detID == 60 || detID == 80 || detID == 100 || detID == 120 || detID == 140 || detID == 160 || detID == 180 || detID == 200 || detID == 220 || detID == 240 || detID == 260 || detID == 280 || detID == 300 || detID == 320) return 20;
  
  
  else return 0;
  
}
*/
/////////////////////////////////////eLoss para los 20 anillos del mbb, cada uno con 16 detectores//////////////////////////////////////////////////////////
Int_t nv=100;

TH1D *heLossA1 = new TH1D("heLossA1","Energy loss in ring 1",nv,0,0.010);
TH1D *heLossA2 = new TH1D("heLossA2","Energy loss in ring 2",nv,0,0.010);
TH1D *heLossA3 = new TH1D("heLossA3","Energy loss in ring 3",nv,0,0.010);
TH1D *heLossA4 = new TH1D("heLossA4","Energy loss in ring 4",nv,0,0.010);
TH1D *heLossA5 = new TH1D("heLossA5","Energy loss in ring 5",nv,0,0.010);
TH1D *heLossA6 = new TH1D("heLossA6","Energy loss in ring 6",nv,0,0.010);
TH1D *heLossA7 = new TH1D("heLossA7","Energy loss in ring 7",nv,0,0.010);
TH1D *heLossA8 = new TH1D("heLossA8","Energy loss in ring 7",nv,0,0.010);
TH1D *heLossA9 = new TH1D("heLossA9","Energy loss in ring 8",nv,0,0.010);
TH1D *heLossA10 = new TH1D("heLossA10","Energy loss in ring 10",nv,0,0.010);
TH1D *heLossA11 = new TH1D("heLossA11","Energy loss in ring 11",nv,0,0.010);
TH1D *heLossA12 = new TH1D("heLossA12","Energy loss in ring 12",nv,0,0.010);
TH1D *heLossA13 = new TH1D("heLossA13","Energy loss in ring 13",nv,0,0.010);
TH1D *heLossA14 = new TH1D("heLossA14","Energy loss in ring 14",nv,0,0.010);
TH1D *heLossA15 = new TH1D("heLossA15","Energy loss in ring 15",nv,0,0.010);
TH1D *heLossA16 = new TH1D("heLossA16","Energy loss in ring 16",nv,0,0.010);
TH1D *heLossA17 = new TH1D("heLossA17","Energy loss in ring 17",nv,0,0.010);
TH1D *heLossA18 = new TH1D("heLossA18","Energy loss in ring 18",nv,0,0.010);
TH1D *heLossA19 = new TH1D("heLossA19","Energy loss in ring 19",nv,0,0.010);
TH1D *heLossA20 = new TH1D("heLossA20","Energy loss in ring 20",nv,0,0.010);


char name[50];
/*for (Int_t nf = 0; nf < 1000 ; nf++){
  sprintf(name,"outtestf14%d",nf);
  cout <<" Openning file" << name << endl;*/

Int_t readOutput6(){
  Int_t contador=0;
  
  //char name[50];
  for (Int_t nf = 0; nf < 10 ; nf++){      // los archivos se llaman outtestf14#.root, dónde lee los nf número de archivos 
    sprintf(name,"~/NICA/bebes/MacrosMBB/phsd_10000/outtestf14%d.root",nf);
    cout <<" Openning file " << name << endl;
    TString inputFile=name;
    TString outputFile="AuAu11_phsd_10000.root";
    
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetOptDate(0);      //show day and time
    gStyle->SetOptStat(1);      //show statistic
    gStyle->SetPalette(1,0);
    
    
    Double_t epsilon = 0.00001;
    
    //    mpdloadlibs();
    TFile fileInput(inputFile.Data());
    if (fileInput.IsZombie()) continue; // Se salta los archivos que no existen
    ++contador;  //cuenta los archivos que si existen
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////MPDROOT Analysis//////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //Defining pointers to data
    
    TTree *simEvent = (TTree*) fileInput.Get("mpdsim");
    
    //Arrays/////nombre de la variable///////////nombre en la memoria interna
  //  TClonesArray *mbbPoints = (TClonesArray*) fileInput.FindObjectAny("MbbPoint");
   // simEvent->SetBranchAddress("MbbPoint",&mbbPoints);
    
TClonesArray *mbbPoints = (TClonesArray*) fileInput.FindObjectAny("MbbPoint");
     simEvent->SetBranchAddress("MbbPoint",&mbbPoints);

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
    TTree *mbbTree = new TTree("mbbTree","MBB");
    mbbTree->Branch("impactParameter",&impactParameter);
    
    Int_t events  = simEvent->GetEntries(); 
    
    //Begins loop for the events
    for (Int_t i = 0000; i < events; ++i) {
      simEvent->GetEntry(i);
      
      MpdMCEventHeader *extraEventHeader = dynamic_cast<MpdMCEventHeader*> (mcHeader);
      
      // impactParameter = extraEventHeader->GetB();
      // bmdTree->Fill() 
      // cout<<"Impact parameter: "<<extraEventHeader->GetB()<<endl;
      
      Int_t nMCTracks  =  mcTracks->GetEntriesFast();
      //  cout<<"nMCTracks "<<nMCTracks<<endl;
      
      if (mbbPoints != 0 ) { //if mbbpoints!=0
	
	Int_t nmbbPoints =  mbbPoints->GetEntriesFast();
	if( nmbbPoints < 1 ) continue;
	
	//cout<<"nmbbPoints: "<<nbmdPoints<<endl;
	
	//Loop for particles in both MBB.
	Int_t nRecTracks = 0;
	Int_t nRecPrimaryTracks = 0; 
	
	for (Int_t j = 0; j < nmbbPoints; j++) { // inicia loop mbbpoints    
	  MbbPoint *p1    = (MbbPoint*) mbbPoints->At(j);
	  
	  	  
	  TVector3 recMom;
	  Int_t    trkID = p1->GetTrackID();
	  Double_t time  = p1->GetTime();
	  Double_t eLoss = p1->GetEnergyLoss();
	  p1->Momentum(recMom);
	  Double_t currPt  = recMom.Pt();
	  Double_t currEta = recMom.Eta();
	  Int_t x = p1->GetX();
	  Int_t y = p1->GetY();
	  Int_t z = p1->GetZ();
	  
	  Double_t detID =  p1->GetDetectorID();
    //cout<<detID<<endl;
	  	  
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////Analysis per ring////////////////////////////////////////////////////////
	  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  
	  	  
	  MpdMCTrack* mcTr = (MpdMCTrack*) mcTracks->UncheckedAt(j);
//	  Double_t Xv = mcTr->GetStartX();
//	  Double_t Yv = mcTr->GetStartY();
//	  Double_t Zv = mcTr->GetStartZ();
	  Int_t ring=-1;

	  if((TMath::Abs(mcTr->GetPdgCode())== 211 || TMath::Abs(mcTr->GetPdgCode())== 11 || TMath::Abs(mcTr->GetPdgCode())== 2212 || TMath::Abs(mcTr->GetPdgCode())== 321 || TMath::Abs(mcTr->GetPdgCode())== 13) && mcTr->GetMotherId() < 0){
	    
	  //ring=GetRing(detID);
	  /////////////eLoss  de los 20 anillos /////////////////////////////////////////         
	  if(detID == 1 || detID == 21 || detID == 41 || detID == 61 || detID == 81 || detID == 101 || detID == 121 || detID == 141 || detID == 161 || detID == 181 || detID == 201 || detID == 221 || detID == 241 || detID == 261 || detID == 281 || detID == 301){
	    heLossA1->Fill(eLoss);
	  }
	  if(detID == 2 || detID == 22 || detID == 42 || detID == 62 || detID == 82 || detID == 102 || detID == 122 || detID == 142 || detID == 162 || detID == 182 || detID == 202 || detID == 222 || detID == 242 || detID == 262 || detID == 282 || detID == 302){
	    heLossA2->Fill(eLoss);
	  }
	  if(detID == 3 || detID == 23 || detID == 43 || detID == 63 || detID == 83 || detID == 103 || detID == 123 || detID == 143 || detID == 163 || detID == 183 || detID == 203 || detID == 223 || detID == 243 || detID == 263 || detID == 283 || detID == 303){
	    heLossA3->Fill(eLoss);
	  }
	  if(detID == 4 || detID == 24 || detID == 44 || detID == 64 || detID == 84 || detID == 104 || detID == 124 || detID == 144 || detID == 164 || detID == 184 || detID == 204 || detID == 224 || detID == 244 || detID == 264 || detID == 284 || detID == 304){
	    heLossA4->Fill(eLoss);
	  }
	  if(detID == 5 || detID == 25 || detID == 45 || detID == 65 || detID == 85 || detID == 105 || detID == 125 || detID == 145 || detID == 165 || detID == 185 || detID == 205 || detID == 225 || detID == 245 || detID == 265 || detID == 285 || detID == 305){
	    heLossA5->Fill(eLoss);
	  }
	  if(detID == 6 || detID == 26 || detID == 46 || detID == 66 || detID == 86 || detID == 106 || detID == 126 || detID == 146 || detID == 166 || detID == 186 || detID == 206 || detID == 226 || detID == 246 || detID == 266 || detID == 286 || detID == 306){
	    heLossA6->Fill(eLoss);
	  }
	  if(detID == 7 || detID == 27 || detID == 47 || detID == 67 || detID == 87 || detID == 107 || detID == 127 || detID == 147 || detID == 167 || detID == 187 || detID == 207 || detID == 227 || detID == 247 || detID == 267 || detID == 287 || detID == 307){
            heLossA7->Fill(eLoss);
          }
	  if(detID == 8 || detID == 28 || detID == 48 || detID == 68 || detID == 88 || detID == 108 || detID == 128 || detID == 148 || detID == 168 || detID == 188 || detID == 208 || detID == 228 || detID == 248 || detID == 268 || detID == 288 || detID == 308){
            heLossA8->Fill(eLoss);
          }
	  if(detID == 9 || detID == 29 || detID == 49 || detID == 69 || detID == 89 || detID == 109 || detID == 129 || detID == 149 || detID == 169 || detID == 189 || detID == 209 || detID == 229 || detID == 249 || detID == 269 || detID == 289 || detID == 309){
            heLossA9->Fill(eLoss);
          }
	  if(detID == 10 || detID == 30 || detID == 50 || detID == 70 || detID == 90 || detID == 110 || detID == 130 || detID == 150 || detID == 170 || detID == 190 || detID == 210 || detID == 230 || detID == 250 || detID == 270 || detID == 290 || detID == 310){
            heLossA10->Fill(eLoss);
          }
	  if(detID == 11 || detID == 31 || detID == 51 || detID == 71 || detID == 91 || detID == 111 || detID == 131 || detID == 151 || detID == 171 || detID == 191 || detID == 211 || detID == 231 || detID == 251 || detID == 271 || detID == 291 || detID == 311){
            heLossA11->Fill(eLoss);
          }
	  if(detID == 12 || detID == 32 || detID == 52 || detID == 72 || detID == 92 || detID == 112 || detID == 132 || detID == 152 || detID == 172 || detID == 192 || detID == 212 || detID == 232 || detID == 252 || detID == 272 || detID == 292 || detID == 312){
            heLossA12->Fill(eLoss);
          }
	  if(detID == 13 || detID == 33 || detID == 53 || detID == 73 || detID == 93 || detID == 113 || detID == 133 || detID == 153 || detID == 173 || detID == 193 || detID == 213 || detID == 233 || detID == 253 || detID == 273 || detID == 293 || detID == 313){
            heLossA13->Fill(eLoss);
          }
	  if(detID == 14 || detID == 34 || detID == 54 || detID == 74 || detID == 94 || detID == 114 || detID == 134 || detID == 154 || detID == 174 || detID == 194 || detID == 214 || detID == 234 || detID == 254 || detID == 274 || detID == 294 || detID == 314){
            heLossA14->Fill(eLoss);
          }
	  if(detID == 15 || detID == 35 || detID == 55 || detID == 75 || detID == 95 || detID == 115 || detID == 135 || detID == 155 || detID == 175 || detID == 195 || detID == 215 || detID == 235 || detID == 255 || detID == 275 || detID == 295 || detID == 315){
            heLossA15->Fill(eLoss);
          }
	  if(detID == 16 || detID == 36 || detID == 56 || detID == 76 || detID == 96 || detID == 116 || detID == 136 || detID == 156 || detID == 176 || detID == 196 || detID == 216 || detID == 236 || detID == 256 || detID == 276 || detID == 296 || detID == 316){
            heLossA16->Fill(eLoss);
          }
	  if(detID == 17 || detID == 37 || detID == 57 || detID == 77 || detID == 97 || detID == 117 || detID == 137 || detID == 157 || detID == 177 || detID == 197 || detID == 217 || detID == 237 || detID == 257 || detID == 277 || detID == 297 || detID == 317){
            heLossA17->Fill(eLoss);
          }
	  if(detID == 18 || detID == 38 || detID == 58 || detID == 78 || detID == 98 || detID == 118 || detID == 138 || detID == 158 || detID == 178 || detID == 198 || detID == 218 || detID == 238 || detID == 258 || detID == 278 || detID == 298 || detID == 318){
            heLossA18->Fill(eLoss);
          }
	  if(detID == 19 || detID == 39 || detID == 59 || detID == 79 || detID == 99 || detID == 119 || detID == 139 || detID == 159 || detID == 179 || detID == 199 || detID == 219 || detID == 239 || detID == 259 || detID == 279 || detID == 299 || detID == 319){
            heLossA19->Fill(eLoss);
          }
	  if(detID == 20 || detID == 40 || detID == 60 || detID == 80 || detID == 100 || detID == 120 || detID == 140 || detID == 160 || detID == 180 || detID == 200 || detID == 220 || detID == 240 || detID == 260 || detID == 280 || detID == 300 || detID == 320){
            heLossA20->Fill(eLoss);
          }
	  }//end if getpdgcode

	  /////////////////////////////////////////////////////////////////////7

                    cout<<"pdgcode: "<< mcTr->GetPdgCode()<<endl; 

	  
	}  //termina loop mbbpoints 
	
      } // termina if loop mbbpoints!=0
    }  //termina el loop de los events
    
    cout<<"Saving histograms"<<endl;
    
    fileOutput->mkdir("Mc");
    fileOutput->cd("Mc");
    
    
    gStyle->SetOptTitle(0); //No title for histograms
    
    
    
    ////////////////////////////// detector  eLoss //////////////////
    
  } //Para el loop de agregar los outtestf14#.root
  Double_t ndetect [20]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}; // 20 anillos
  Double_t norm = contador*10000*16;  //Normalización del número de archivos, 100 eventos por cada archivo, 16 detectores por anillo
  Double_t numdet=16; //número de detectores por anillo
  
  Double_t elossA[20]={(heLossA1->GetMean())/numdet,(heLossA2->GetMean())/numdet,(heLossA3->GetMean())/numdet,(heLossA4->GetMean())/numdet,(heLossA5->GetMean())/numdet,(heLossA6->GetMean())/numdet,(heLossA7->GetMean())/numdet,(heLossA8->GetMean())/numdet,(heLossA9->GetMean())/numdet,(heLossA10->GetMean())/numdet,(heLossA11->GetMean())/numdet,(heLossA12->GetMean())/numdet,(heLossA13->GetMean())/numdet,(heLossA14->GetMean())/numdet,(heLossA15->GetMean())/numdet,(heLossA16->GetMean())/numdet,(heLossA17->GetMean())/numdet,(heLossA18->GetMean())/numdet,(heLossA19->GetMean())/numdet,(heLossA20->GetMean())/numdet};  //Valor medio de eLoss por anillo / 16 detectores en cada anillo
  
  Double_t elosserrA[20]={(heLossA1->GetMeanError())/numdet,(heLossA2->GetMeanError())/numdet,(heLossA3->GetMeanError())/numdet,(heLossA4->GetMeanError())/numdet,(heLossA5->GetMeanError())/numdet,(heLossA6->GetMeanError())/numdet,(heLossA7->GetMeanError())/numdet,(heLossA8->GetMeanError())/numdet,(heLossA9->GetMeanError())/numdet,(heLossA10->GetMeanError())/numdet,(heLossA11->GetMeanError())/numdet,(heLossA12->GetMeanError())/numdet,(heLossA13->GetMeanError())/numdet,(heLossA14->GetMeanError())/numdet,(heLossA15->GetMeanError())/numdet,(heLossA16->GetMeanError())/numdet,(heLossA17->GetMeanError())/numdet,(heLossA18->GetMeanError())/numdet,(heLossA19->GetMeanError())/numdet,(heLossA20->GetMeanError())/numdet};  //error de cada anillo / 16 detectores de cada anillo
  
  Double_t entriesA[20]={(heLossA1->GetEntries())/norm,(heLossA2->GetEntries())/norm,(heLossA3->GetEntries())/norm,(heLossA4->GetEntries())/norm,(heLossA5->GetEntries())/norm,(heLossA6->GetEntries())/norm,(heLossA7->GetEntries())/norm,(heLossA8->GetEntries())/norm,(heLossA9->GetEntries())/norm,(heLossA10->GetEntries())/norm,(heLossA11->GetEntries())/norm,(heLossA12->GetEntries())/norm,(heLossA13->GetEntries())/norm,(heLossA14->GetEntries())/norm,(heLossA15->GetEntries())/norm,(heLossA16->GetEntries())/norm,(heLossA17->GetEntries())/norm,(heLossA18->GetEntries())/norm,(heLossA19->GetEntries())/norm,(heLossA20->GetEntries())/norm};  // Hits de cada anillo / normalización 
  
  Double_t errorhA[20]={(heLossA1->GetStdDev())/(Sqrt(heLossA1->GetEntries())*numdet),(heLossA2->GetStdDev())/(Sqrt(heLossA2->GetEntries())*numdet),(heLossA3->GetStdDev())/(Sqrt(heLossA3->GetEntries())*numdet),(heLossA4->GetStdDev())/(Sqrt(heLossA4->GetEntries())*numdet),(heLossA5->GetStdDev())/(Sqrt(heLossA5->GetEntries())*numdet),(heLossA6->GetStdDev())/(Sqrt(heLossA6->GetEntries())*numdet),(heLossA7->GetStdDev())/(Sqrt(heLossA7->GetEntries())*numdet),(heLossA8->GetStdDev())/(Sqrt(heLossA8->GetEntries())*numdet),(heLossA9->GetStdDev())/(Sqrt(heLossA9->GetEntries())*numdet),(heLossA10->GetStdDev())/(Sqrt(heLossA10->GetEntries())*numdet),(heLossA11->GetStdDev())/(Sqrt(heLossA11->GetEntries())*numdet),(heLossA12->GetStdDev())/(Sqrt(heLossA12->GetEntries())*numdet),(heLossA13->GetStdDev())/(Sqrt(heLossA13->GetEntries())*numdet),(heLossA14->GetStdDev())/(Sqrt(heLossA14->GetEntries())*numdet),(heLossA15->GetStdDev())/(Sqrt(heLossA15->GetEntries())*numdet),(heLossA16->GetStdDev())/(Sqrt(heLossA16->GetEntries())*numdet),(heLossA17->GetStdDev())/(Sqrt(heLossA17->GetEntries())*numdet),(heLossA18->GetStdDev())/(Sqrt(heLossA18->GetEntries())*numdet),(heLossA19->GetStdDev())/(Sqrt(heLossA19->GetEntries())*numdet),(heLossA20->GetStdDev())/(Sqrt(heLossA20->GetEntries())*numdet)}; // error desviación estandar / (sqrt(N) * 16 detectores de cada anillo)

  
   
  TCanvas *c5 = new TCanvas("c5","Energy deposited");
  //c5->SetWindowSize(1000,1000);
  TGraph* gr1A = new TGraphErrors(20,ndetect,elossA,0,elosserrA); // datos con error
  TMultiGraph *mg1 = new TMultiGraph();
  gr1A->SetMarkerStyle(20);
  //gr1A->SetMarkerColor(2);
  //gr1A->SetLineColor(2);
  //gr1C->SetMarkerColor(4);
  //gr1C->SetLineColor(4);
  mg1->Add(gr1A); // agrega los datos con el error
  mg1->Draw("AP");
  mg1->GetYaxis()->SetNdivisions(6);
  TGaxis::SetMaxDigits(3); 
  mg1->GetXaxis()->SetTitle("Ring");
  mg1->GetYaxis()->SetTitle("< E > [GeV]");
  mg1->GetYaxis()->SetTitleOffset(1.25);
  mg1->GetXaxis()->SetTitleOffset(1.25);
  mg1->GetXaxis()->CenterTitle(true);
  mg1->GetYaxis()->CenterTitle(true);
  mg1->Write("mg1");
  TLegend *leg1 = new TLegend(0.45,0.7,0.80,0.85);
  leg1->SetTextFont(62);
  //leg2->SetTextSize(0.04);                                    
  leg1->SetLineColor(0);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(1001);
  leg1->AddEntry("","PHSD Au-Au #sqrt{s_{NN}} = 11 GeV","");
  leg1->Draw();
  TCanvas *c6 = new TCanvas("c6","Number Hits");
  TGraph* gr2A = new TGraphErrors(20,ndetect,entriesA,0,errorhA); // datos con error
  TMultiGraph *mg2 = new TMultiGraph();
  gr2A->SetMarkerStyle(20);
  //gr2A->SetMarkerColor(2);
  //gr2A->SetLineColor(2);
  //gr2C->SetMarkerColor(4);
  //gr2C->SetLineColor(4);
  mg2->Add(gr2A); // agrega los datos con el error
  //void TAxis::SetRangeUser(0,6);
  //mg2->GetYaxis()->
  mg2->Draw("AP");
  //TGaxis::SetMaxDigits(1);
  //mg2->GetYaxis()->SetMaxDigits(6);
  //mg2->GetYaxis()->SetMaxDigits(3);
  mg2->GetYaxis()->SetTitleOffset(1.25);
  mg2->GetXaxis()->SetTitleOffset(1.25);
  mg2->GetXaxis()->SetTitle("Ring");
  mg2->GetYaxis()->SetTitle("Hits");
  mg2->GetXaxis()->CenterTitle(true);
  mg2->GetYaxis()->CenterTitle(true);
  mg2->Write("mg2");
  
  TLegend *leg2 = new TLegend(0.45,0.7,0.80,0.85);
  leg2->SetTextFont(62);
  //leg2->SetTextSize(0.04);                                                            
  leg2->SetLineColor(0);
  leg2->SetLineStyle(0);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(1001);
  leg2->AddEntry("","PHSD Au-Au #sqrt{s_{NN}} = 11 GeV","");
  leg2->Draw();
  
  
  TCanvas *c3 = new TCanvas("c3","Energy loss of detector A"); // Para probar eLoss de los primeros 6 anillos
  c3->Divide(3,2);
  c3->cd(1);
  heLossA1->GetXaxis()->SetTitle("E_{Loss} (GeV)");
  heLossA1->GetYaxis()->SetTitle("1/N dN/dE");
  heLossA1->Scale(1./heLossA1->Integral());
  heLossA1->Draw();
  
  c3->cd(2);
  heLossA2->GetXaxis()->SetTitle("E_{Loss} (GeV)");
  heLossA2->GetYaxis()->SetTitle("1/N dN/dE");
  heLossA2->Scale(1./heLossA2->Integral());
  heLossA2->Draw();
  
  c3->cd(3);
  heLossA3->GetXaxis()->SetTitle("E_{Loss} (GeV)");
  heLossA3->GetYaxis()->SetTitle("1/N dN/dE");
  heLossA3->Scale(1./heLossA3->Integral());
  heLossA3->Draw();
  
  c3->cd(4);
  heLossA4->GetXaxis()->SetTitle("E_{Loss} (GeV)");
  heLossA4->GetYaxis()->SetTitle("1/N dN/dE");
  heLossA4->Scale(1./heLossA4->Integral());
  heLossA4->Draw();
  
  c3->cd(5);
  heLossA5->GetXaxis()->SetTitle("E_{Loss} (GeV)");
  heLossA5->GetYaxis()->SetTitle("1/N dN/dE");
  heLossA5->Scale(1./heLossA5->Integral());
  heLossA5->Draw();
  
  c3->cd(6);
  heLossA6->GetXaxis()->SetTitle("E_{Loss} (GeV)");
  heLossA6->GetYaxis()->SetTitle("1/N dN/dE");
  heLossA6->Scale(1./heLossA6->Integral());
  heLossA6->Draw();
  cout<<"Contador="<<contador<<endl; 
  cout<<"End histograms"<<endl;


  return 0;
  
  
} //readoutput4

