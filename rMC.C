
// Description:
//      
//       Analysis macro. This macro is to study several distributions such as pT, energy, and multiplicity generated on URQMD.
//
//
// Environment:
//      ROOT
//
// Author List:
//       Luis Valenzuela-Cazares          (original author)
//   
//-----------------------------------------------------------

#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>


Int_t rMC(){
  
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////All space distributions/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////Pseudorapity distributions/////////////////////////////////////////////////////////////////////
//////////////////////////////////Charged particles in all space///////////////////////////////////////////////////////////////////////////////////////////
   //TH1D *hMCEtaCharged = new TH1D("hMCEtaCharged","Pseudorapity charged particles. 10000 Au+Au @11GeV UrQMD.",240,-8,8);
         //hMCEtaCharged->SetYTitle("1/N dN/d#eta");
         //hMCEtaCharged->SetXTitle("#eta");
         //hMCEtaCharged->SetMarkerColor(kBlack);
         //hMCEtaCharged->SetMarkerStyle(kDot);
         //hMCEtaCharged->GetXaxis()->CenterTitle(true);
         //hMCEtaCharged->GetXaxis()->SetTitleSize(0.04);
         //hMCEtaCharged->GetXaxis()->SetLabelSize(0.03);
         //hMCEtaCharged->GetXaxis()->SetTitleOffset(1.4);

TH1F *heta1 = new TH1F("h1b","UrQMD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#eta;# entries",277,-8,8);
TH1F *heta2 = new TH1F("h2b","eta",277,-8,8);
TH1F *heta3 = new TH1F("h3b","eta",277,-8,8);
TH1F *heta4 = new TH1F("h4b","eta",277,-8,8);


  //////////////////////////////////Charged pions in all space///////////////////////////////////////////////////////////////////////////////////////////
    //TH1D *etapionsH = new TH1D("etapionsH","@Eta pions . MC", 240,-8,8);
          //etapionsH->SetMarkerColor(kRed);
          //etapionsH->SetMarkerStyle(kPlus);

TH1F *pheta1 = new TH1F("h1c","UrQMD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#eta;# entries",277,-8,8);
TH1F *pheta2 = new TH1F("h2c","eta",277,-8,8);

  //////////////////////////////////Protons in all space///////////////////////////////////////////////////////////////////////////////////////////
    //TH1D *etaprotonsH = new TH1D("etaprotonsH","Protons . MC", 240,-8,8);
          //etaprotonsH->SetMarkerColor(kGreen);
          //etaprotonsH->SetMarkerStyle(kStar);

    //TH1D *etaSpectatorprotonsH = new TH1D("etaSpectatorprotonsH","Spectator protons . MC", 240,-8,8);
          //etaSpectatorprotonsH->SetMarkerColor(kMagenta);
          //etaSpectatorprotonsH->SetMarkerStyle(kStar);

   //TH1D  *etaParticipantprotonsH = new TH1D("etaParticipantprotonsH","Participant protons . MC", 240,-8,8);
          //etaParticipantprotonsH->SetMarkerColor(kYellow);
          //etaParticipantprotonsH->SetMarkerStyle(kPlus);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////BEBE region distributions//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////Some charged particles in bmd region//////////////////////////////////
    //TH1D *hMCEtaBMDCharged = new TH1D("hMCEtaBMDCharged","@Eta charged particles BEBE region. 10000 Au+Au @11GeV UrQMD. MC.",240,-8,8);
         //hMCEtaBMDCharged->SetYTitle("1/N dN/d#eta");
         //hMCEtaBMDCharged->SetXTitle("#eta");
         //hMCEtaBMDCharged->SetMarkerColor(kBlack);
         //hMCEtaBMDCharged->SetMarkerStyle(kDot);
         //hMCEtaBMDCharged->GetXaxis()->CenterTitle(true);
         //hMCEtaBMDCharged->GetXaxis()->SetTitleSize(0.04);
         //hMCEtaBMDCharged->GetXaxis()->SetLabelSize(0.03);
         //hMCEtaBMDCharged->GetXaxis()->SetTitleOffset(1.4);
    
    //TH1D *etapionsHBMD = new TH1D("etapionsHBMD","@Eta pions BEBE region. MC", 240,-8,8);
          //etapionsHBMD->SetMarkerColor(kRed);
          //etapionsHBMD->SetMarkerStyle(kPlus);

   //////////////////////////////////Protons ///////////////////////////////////////////////////////////////////////////////////////////
    //TH1D *etaprotonsHBMD = new TH1D("etaprotonsHBMD","Protons . MC", 240,-8,8);
          //etaprotonsHBMD->SetMarkerColor(kGreen);
          //etaprotonsHBMD->SetMarkerStyle(kStar);

    //TH1D *etaSpectatorprotonsHBMD = new TH1D("etaSpectatorprotonsHBMD","Spectator protons . MC", 240,-8,8);
          //etaSpectatorprotonsHBMD->SetMarkerColor(kMagenta);
          //etaSpectatorprotonsHBMD->SetMarkerStyle(kStar);

   //TH1D  *etaParticipantprotonsHBMD = new TH1D("etaParticipantprotonsHBMD","Participant protons . MC", 240,-8,8);
          //etaParticipantprotonsHBMD->SetMarkerColor(kYellow);
          //etaParticipantprotonsHBMD->SetMarkerStyle(kPlus);        

 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         

/////////////////////////////////////////////////////////////////////////////Individual plots//////////////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////Primary charged particles in bmd region//////////////////////////////////
    //TH1D *hMCEtaPrimaryBMD = new TH1D("hMCEtaPrimaryBMD","@Eta primary charged particles BEBE region. 5000 Au+Au @11GeV UrQMD. MC.",240,-8,8);
    //hMCEtaPrimaryBMD->SetXTitle("#eta");

    //TH1D *hMCPtPrimaryBMDA = new TH1D("hMCPtPrimaryBMDA","Pt primary charged particles BEBE region. 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hMCPtPrimaryBMDA->SetXTitle("p_T (GeV/#it{c})");

 ///////////////////////////////////////Pt vs Pseudorapidity////////////////////////////////////////////////////////////////
  //TH2D *hBmdPtVsEta  = new TH2D("hBmdPtVsEta","BMD pT vs Pseudorapidity. 5000 Au+Au @9GeV UrQMD.",240,-8,8,60,-0,1.5);
  //hBmdPtVsEta->SetXTitle("#Eta");
  //hBmdPtVsEta->SetYTitle("p_T (GeV/#it{c})");



   
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////Primary charged particles BMD side A per ring/////////////////////////
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////Eta /////////////////////////////////////////////////////////
    //TH1D *hEtaPrimaryChargedBMDARing1 = new TH1D("hEtaPrimaryChargedBMDARing1","Charged particles @Eta BMDA 5000 Au+Au @9GeV UrQMD. MC.",200,-5,5);
    //hEtaPrimaryChargedBMDARing1->SetXTitle("#eta");
    
    //TH1D *hEtaPrimaryChargedBMDARing2 = new TH1D("hEtaPrimaryChargedBMDARing2","Charged particles @Eta BMDA 5000 Au+Au @9GeV UrQMD. MC.",200,-5,5);
    //hEtaPrimaryChargedBMDARing2->SetXTitle("#eta");
   
    //TH1D *hEtaPrimaryChargedBMDARing3 = new TH1D("hEtaPrimaryChargedBMDARing3","Charged particles @Eta BMDA 5000 Au+Au @9GeV UrQMD. MC.",200,-5,5);
    //hEtaPrimaryChargedBMDARing3->SetXTitle("#eta");
  
    //TH1D *hEtaPrimaryChargedBMDARing4 = new TH1D("hEtaPrimaryChargedBMDARing4","Charged particles @Eta BMDA 5000 Au+Au @9GeV UrQMD. MC.",200,-5,5);
    //hEtaPrimaryChargedBMDARing4->SetXTitle("#eta");
  
    //TH1D *hEtaPrimaryChargedBMDARing5 = new TH1D("hEtaPrimaryChargedBMDARing5","Charged particles @Eta BMDA 5000 Au+Au @9GeV UrQMD. MC.",200,-5,5);
    //hEtaPrimaryChargedBMDARing5->SetXTitle("#eta");
   
    //TH1D *hEtaPrimaryChargedBMDARing6 = new TH1D("hEtaPrimaryChargedBMDARing6","Charged particles @Eta BMDA 5000 Au+Au @9GeV UrQMD. MC.",200,-5,5);
    //hEtaPrimaryChargedBMDARing6->SetXTitle("#eta");

    ///////////////////////////////////////////////Pt/////////////////////////////////////////////////////////////////
    //TH1D *hPtPrimaryChargedBMDARing1 = new TH1D("hPtPrimaryChargedBMDARing1","Charged particles Pt BMDA 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hPtPrimaryChargedBMDARing1 ->SetXTitle("p_T (GeV/#it{c})");

    //TH1D *hPtPrimaryChargedBMDARing2 = new TH1D("hPtPrimaryChargedBMDARing2","Charged particles Pt BMDA 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hPtPrimaryChargedBMDARing2->SetXTitle("p_T (GeV/#it{c})");

    //TH1D *hPtPrimaryChargedBMDARing3 = new TH1D("hPtPrimaryChargedBMDARing3","Charged particles Pt BMDA 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hPtPrimaryChargedBMDARing3->SetXTitle("p_T (GeV/#it{c})");
  
    //TH1D *hPtPrimaryChargedBMDARing4 = new TH1D("hPtPrimaryChargedBMDARing4","Charged particles Pt BMDA 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hPtPrimaryChargedBMDARing4->SetXTitle("p_T (GeV/#it{c})");

    //TH1D *hPtPrimaryChargedBMDARing5 = new TH1D("hPtPrimaryChargedBMDARing5","Charged particles Pt BMDA 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hPtPrimaryChargedBMDARing5->SetXTitle("p_T (GeV/#it{c})");

    //TH1D *hPtPrimaryChargedBMDARing6 = new TH1D("hPtPrimaryChargedBMDARing6","Charged particles Pt BMDA 5000 Au+Au @9GeV UrQMD. MC.",60,0,1.5);
    //hPtPrimaryChargedBMDARing6->SetXTitle("p_T (GeV/#it{c})");

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  gROOT->Reset();
  TChain mychain("T");
  mychain.Add("urqmd_10000.root");
  Int_t nlines = 0;
  struct particula_t 
  {
    Float_t time,X,Y,Z,E,Px,Py,Pz,Pt,P,Eta,m,id,isoespin,charge,lastcoll,numbercoll,history,frezetime,frezeX,frezeY,frezeZ,frezeE,frezePx;
    Float_t frezePy,frezePz,frezePt,frezeP,frezeEta;
  } PARTICLE;
  
  particula_t  particle;
  mychain.SetBranchAddress("particle",&particle);
  Int_t nevent = mychain.GetEntries();
  for (Int_t i=0;i<nevent;i++) 
    {
      mychain.GetEvent(i);
      //cout << "pid=" << particle.id << " charge=" << particle.charge << " mass=" << particle.m << endl;

//        cout<<"pT: "<<particle.Pt<<endl;  
                 
                       heta1->Fill(particle.Eta);

                          // cout<<"N coll : "<<particle.numbercoll<<endl;  

////////////////////////////////////////All space/////////////////////////////////////////////////

       if( particle.charge!=0){
                      heta4->Fill(particle.Eta);  
                     //hMCEtaCharged->Fill(particle.Eta);
        //cout << "pid=" << particle.id << " charge=" << particle.charge << " mass=" << particle.m << endl;
       }

       if(particle.id==101 && particle.charge==-1){//charged pions
                     pheta1->Fill(particle.Eta);
    //   cout << "eta pions:" << particle.Eta <<endl;
        }      
 
       if(particle.id==101 && particle.charge==1){
                     pheta2->Fill(particle.Eta);
        }

       if(particle.id==1 && particle.charge!=0){//protons
                      heta2->Fill(particle.Eta);
               
               //if(particle.numbercoll<1) {  //spectators
                      //etaSpectatorprotonsH->Fill(particle.Eta);
               //}

               //if(particle.numbercoll>=1) {  //Participants
                      //etaParticipantprotonsH->Fill(particle.Eta);

               }

       if(particle.id==1 && particle.charge==0){//neutrons
                    heta3->Fill(particle.Eta);
        }

     //   cout << "pid=" << particle.id << " charge=" << particle.charge << " mass=" << particle.m << endl;

       }      
                                           
///////////////////////////////////////////////////////////////////////////////////////////////////
                     //if(particle.charge!=0){

                     //hBmdPtVsEta->Fill(particle.Eta,particle.Pt);                 
                     //}
///////////////////////////////////////////Be-Be pseudorapidity region/////////////////////////////////////////////////////////////////////
                //if( TMath::Abs(particle.Eta) >  2.2 && TMath::Abs(particle.Eta) < 4.2 ) {
                  //        cout << "pid=" << particle.id << " charge=" << particle.charge << " mass=" << particle.m << endl;
                     

                     //if(particle.charge!=0){
                      //    if(particle.charge!=0 && particle.id!=1){ //without protons
                          //hMCEtaPrimaryBMD->Fill(particle.Eta);
                          
                        //   hBmdPtVsEta->Fill(particle.Eta,particle.Pt);


                          //}                     

                                          //Charged particles 
                                         //if(particle.charge!=0){
                                         //Charged primary particles  in pseudorapidity region of BEBE                                
                                          //hMCEtaBMDCharged->Fill(particle.Eta);
                                          //}        

                                           //if(particle.id==101 && particle.charge!=0){//charged pions
                                           //etapionsHBMD->Fill(particle.Eta);
                                            //}

                                           //if(particle.id==1 && particle.charge!=0){//protons
                                           //etaprotonsHBMD->Fill(particle.Eta);

                                                      //if(particle.numbercoll<1) {  //spectators
                                                      //etaSpectatorprotonsHBMD->Fill(particle.Eta);
                                                      //}

                                                      //if(particle.numbercoll>=1) {  //Participants
                                                      //etaParticipantprotonsHBMD->Fill(particle.Eta);
                                                      //}

                                            //}    
                                          
                   //} //Pseudorapidity region BEBE

                  

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////Pt BMD A pseudorapidity region///////////////////////////////////////////////////////////////////////////////////////////////

 //if( particle.Eta >  2.2 && particle.Eta < 4.2 ) {
                  //        cout << "pid=" << particle.id << " charge=" << particle.charge << " mass=" << particle.m << endl;
                         
                  //       if(particle.charge!=0 && particle.id!=1){ //without protons
                         //Charged particles 
                         //if( particle.charge!=0){
    
                            //hMCPtPrimaryBMDA->Fill(particle.Pt);
                            // cout<<"pT: "<<particle.Pt<<endl;  
                          //}
 //}                
               
//////////////////////////////////////////////////Rings analysis//////////////////////////////////////////////////////
  //Charged particles


  /*((if(particle.charge!=0 && particle.id!=1){ //Charged particles except protons

  //if( particle.charge!=0){ //Charged particles*

   
    ((if(particle.Eta <  4.2 && particle.Eta > 3.5 ) {
                        //Charged primary particles  in pseudorapidity region of BEBE                                
                         hEtaPrimaryChargedBMDARing1->Fill(particle.Eta);
                         hPtPrimaryChargedBMDARing1->Fill(particle.Pt);

     }

    if( particle.Eta <  3.5 && particle.Eta > 3.2 ) {
                        //Charged primary particles  in pseudorapidity region of BEBE                                
                         hEtaPrimaryChargedBMDARing2->Fill(particle.Eta);
                         hPtPrimaryChargedBMDARing2->Fill(particle.Pt);

     }

    if( particle.Eta <  3.2 && particle.Eta > 3.0 ) {
                        //Charged primary particles  in pseudorapidity region of BEBE                                
                         hEtaPrimaryChargedBMDARing3->Fill(particle.Eta);
                         hPtPrimaryChargedBMDARing3->Fill(particle.Pt);
                       
    }


    if( particle.Eta <  3.0 && particle.Eta > 2.8 ) {
                        //Charged primary particles  in pseudorapidity region of BEBE                                
                         hEtaPrimaryChargedBMDARing4->Fill(particle.Eta);
                         hPtPrimaryChargedBMDARing4->Fill(particle.Pt);

    }


    if( particle.Eta <  2.8 && particle.Eta > 2.6 ) {
                        //Charged primary particles  in pseudorapidity region of BEBE                                
                         hEtaPrimaryChargedBMDARing5->Fill(particle.Eta);
                         hPtPrimaryChargedBMDARing5->Fill(particle.Pt);

    }


    if( particle.Eta <  2.6 && particle.Eta > 2.2 ) {
                        //Charged primary particles  in pseudorapidity region of BEBE                                
                         hEtaPrimaryChargedBMDARing6->Fill(particle.Eta);
                         hPtPrimaryChargedBMDARing6->Fill(particle.Pt);

        }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////7


             
      }

  ///////////////////////////////////Eta all space MC charged particles//////////////////////////////////////////////////////////////////////
    TCanvas *c1 = new TCanvas("c1","UrQMD test example",800,800);
   gStyle->SetOptStat(false);                                   
   //gStyle->SetPalette(1);                                       
   c1->SetRightMargin(0.0465116);
   c1->SetTopMargin(0.1);
   c1->SetFillColor(0);*/



//////////////////////////////////////////////////////////////////////////
   

    /*hMCEtaCharged->Scale(1.0/hMCEtaCharged->GetEntries());
    hMCEtaCharged->Draw();
    
    etapionsH->Scale(1.0/hMCEtaCharged->GetEntries());
    etapionsH->Draw("sames");

    etaprotonsH->Scale(1.0/hMCEtaCharged->GetEntries());
    etaprotonsH->Draw("sames");   

    etaSpectatorprotonsH->Scale(1.0/hMCEtaCharged->GetEntries());
    etaSpectatorprotonsH->Draw("sames");  

    etaParticipantprotonsH->Scale(1.0/hMCEtaCharged->GetEntries());
    etaParticipantprotonsH->Draw("sames");*/
//////////////////////////////////////////////////////////////////////////////
  
   //Int_t IntegralCharged = hMCEtaCharged->Integral();
        //cout<<"Integral: "<<IntegralCharged<<endl;  
    /*      Int_t IntegralPions = etapionsH->Integral();
         cout<<"IntegralPions: "<<IntegralPions<<endl;
  */

   /*TLegend *leg = new TLegend(0.65,0.8,0.92,0.89);
   leg->SetTextFont(62);
   //leg2->SetTextSize(0.04);                                    
   leg->SetLineColor(0);
   leg->SetLineStyle(0);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry("","UrQMD Au-Au #sqrt{s_{NN}} = 11 GeV","");
   leg->AddEntry("hMCEtaCharged","all","p");
   leg->AddEntry("etapionsH","charged pions","p");
   leg->AddEntry("etaprotonsH","protons","p");
   leg->AddEntry("etaSpectatorprotonsH","spectator protons","p"); 
   leg->AddEntry("etaParticipantprotonsH","participant protons","p");
   leg->Draw();*/

TCanvas *c2=new TCanvas("c2","c2",1080,720);
c2->cd();

heta1->SetMarkerColor(kBlue);
heta1->SetMarkerStyle(kFullCross);
heta1-> Scale(1./(heta1->GetEntries()));
heta1->Draw();
heta2->SetMarkerColor(kGreen);
heta2->SetMarkerStyle(kFullStar);
heta2-> Scale(1./(heta1->GetEntries()));
heta2->Draw("Same");
heta3->SetMarkerColor(kMagenta);
heta3->SetMarkerStyle(kCircle);
heta3-> Scale(1./(heta1->GetEntries()));
heta3->Draw("Same");
heta4->SetMarkerColor(kBlack);
heta4->SetMarkerStyle(kFullDotLarge);
heta4-> Scale(1./(heta1->GetEntries()));
heta4->Draw("Same");

TLegend *legend2=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend2->SetHeader("");
legend2->AddEntry(heta1,"all particles","p");
legend2->AddEntry(heta2,"protons","p");
legend2->AddEntry(heta3,"neutrons","p");
legend2->AddEntry(heta4,"charged particles","p");
legend2->Draw();

TCanvas *c3=new TCanvas("c3","c3",1080,720);
c3->cd();

pheta1->SetMarkerColor(kBlue);
pheta1->SetMarkerStyle(kFullDotLarge);
pheta1->Scale(1./(heta1->GetEntries()));
pheta1->Draw();
pheta2->SetMarkerColor(kBlack);
pheta2->SetMarkerStyle(kFullStar);
pheta2->Scale(1./(heta1->GetEntries()));
pheta2->Draw("Same");

TLegend *legend3=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend3->SetHeader("");
legend3->AddEntry(pheta1,"#pi^{-}","p");
legend3->AddEntry(pheta2,"#pi^{+}","p");
legend3->Draw();

   //c1->SaveAs("MC-eta-9GeV.pdf");
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 ///////////////////////////////////Eta BE-BE region. MC charged particles//////////////////////////////////////////////////////////////////////
    /*TCanvas *c2 = new TCanvas("c2","UrQMD test example",800,800);
   gStyle->SetOptStat(false);                                   
   //gStyle->SetPalette(1);                                       
   c2->SetRightMargin(0.0465116);
   c2->SetTopMargin(0.1);
   c2->SetFillColor(0);

  double norm1 = 1.;

    hMCEtaBMDCharged->Scale(1.0/hMCEtaBMDCharged->GetEntries());
    hMCEtaBMDCharged->Draw();
    
     etapionsHBMD->Scale(1.0/hMCEtaBMDCharged->GetEntries());
     etapionsHBMD->Draw("sames");

     etaprotonsHBMD->Scale(1.0/hMCEtaBMDCharged->GetEntries());
     etaprotonsHBMD->Draw("sames");   

    etaSpectatorprotonsHBMD->Scale(1.0/hMCEtaBMDCharged->GetEntries());
    etaSpectatorprotonsHBMD->Draw("sames");  

    etaParticipantprotonsHBMD->Scale(1.0/hMCEtaBMDCharged->GetEntries());
    etaParticipantprotonsHBMD->Draw("sames"); 
//hMCEtaBMDCharged
//hMCEtaCharged    

   TLegend *leg2 = new TLegend(0.65,0.8,0.92,0.89);
   leg2->SetTextFont(62);
   //leg2->SetTextSize(0.04);                                    
   leg2->SetLineColor(0);
   leg2->SetLineStyle(0);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(1001);
   leg2->AddEntry("","UrQMD Au-Au #sqrt{s_{NN}} = 11 GeV","");
   leg2->AddEntry("hMCEtaBMDCharged","all","p");
   leg2->AddEntry("etapionsHBMD","charged pions","p");
   leg2->AddEntry("etaprotonsHBMD","protons","p");
   leg2->AddEntry("etaSpectatorprotonsHBMD","spectator protons","p"); 
   leg2->AddEntry("etaParticipantprotonsHBMD","participant protons","p");
   leg2->Draw();

   c2->SaveAs("MC-eta-BEBE-9GeV.pdf");
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

///////////////////////Primary charged particles in Be-Be region/////////////////

//////////////////////////////////////Pt///////////////////////////////////////////
   TCanvas *c13 = new TCanvas("c13","Pt primary charged particles BMD A");
   c13->cd(1);
   hMCPtPrimaryBMDA->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
   hMCPtPrimaryBMDA->GetYaxis()->SetTitle("Entries");
   hMCPtPrimaryBMDA->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
   hMCPtPrimaryBMDA->Draw();
   
   c13->SaveAs("MC-pT-BEBE-A-9GeV.pdf");

     TCanvas *c18 = new TCanvas("c18","Pt primary charged particles BMD");
     c18->cd(1);
      hBmdPtVsEta->GetXaxis()->SetTitle("#Eta");
      hBmdPtVsEta->GetYaxis()->SetTitle("#pT");
      hBmdPtVsEta->Scale(1./(hBmdPtVsEta->Integral()));
      hBmdPtVsEta->Draw();


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   //////////////////////Analysis per ring side A/////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    
/////////////////////Eta primary charged particles///////////////////////////////////////////////////////

   TCanvas *c23 = new TCanvas("c23","Pseudorapidity per ring BMD side A");
   c23->Divide(3,2);
   c23->cd(1); 
   hEtaPrimaryChargedBMDARing1->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing1->GetYaxis()->SetTitle("Entries");
   hEtaPrimaryChargedBMDARing1->Scale(1.0/hMCEtaCharged->GetEntries());
   hEtaPrimaryChargedBMDARing1->Draw();

   c23->cd(2);
   hEtaPrimaryChargedBMDARing2->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing2->GetYaxis()->SetTitle("Entries");
   hEtaPrimaryChargedBMDARing2->Scale(1.0/hMCEtaCharged->GetEntries());
   hEtaPrimaryChargedBMDARing2->Draw();

   c23->cd(3);
   hEtaPrimaryChargedBMDARing3->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing3->GetYaxis()->SetTitle("Entries");
   hEtaPrimaryChargedBMDARing3->Scale(1.0/hMCEtaCharged->GetEntries());
   hEtaPrimaryChargedBMDARing3->Draw();

   c23->cd(4); 
   hEtaPrimaryChargedBMDARing4->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing4->GetYaxis()->SetTitle("Entries");
   hEtaPrimaryChargedBMDARing4->Scale(1.0/hMCEtaCharged->GetEntries());
   hEtaPrimaryChargedBMDARing4->Draw();

   c23->cd(5);
   hEtaPrimaryChargedBMDARing5->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing5->GetYaxis()->SetTitle("Entries");
   hEtaPrimaryChargedBMDARing5->Scale(1.0/hMCEtaCharged->GetEntries());
   hEtaPrimaryChargedBMDARing5->Draw();

   c23->cd(6);
   hEtaPrimaryChargedBMDARing6->GetXaxis()->SetTitle("Pseudorapidity");
   hEtaPrimaryChargedBMDARing6->GetYaxis()->SetTitle("Entries");
   hEtaPrimaryChargedBMDARing6->Scale(1.0/hMCEtaCharged->GetEntries());
   hEtaPrimaryChargedBMDARing6->Draw();

   c23->SaveAs("MC-eta-BEBE-A-Rings-9GeV.pdf");
   ///////////////////////////////////////////////////////////////////////////////////
 

 
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////Pt charged primary particles BMD A/////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////   

 TCanvas *c25 = new TCanvas("c25","Pt primary charged particles");
    c25->Divide(3,2);
    c25->cd(1);
    hPtPrimaryChargedBMDARing1->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
    hPtPrimaryChargedBMDARing1->GetYaxis()->SetTitle("Entries");
    hPtPrimaryChargedBMDARing1->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
    hPtPrimaryChargedBMDARing1->Draw();

    c25->cd(2);
    hPtPrimaryChargedBMDARing2->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
    hPtPrimaryChargedBMDARing2->GetYaxis()->SetTitle("Entries");
    hPtPrimaryChargedBMDARing2->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
    hPtPrimaryChargedBMDARing2->Draw();

    c25->cd(3);
    hPtPrimaryChargedBMDARing3->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
    hPtPrimaryChargedBMDARing3->GetYaxis()->SetTitle("Entries");
    hPtPrimaryChargedBMDARing3->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
    hPtPrimaryChargedBMDARing3->Draw();

    c25->cd(4);
    hPtPrimaryChargedBMDARing4->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
    hPtPrimaryChargedBMDARing4->GetYaxis()->SetTitle("Entries");
    hPtPrimaryChargedBMDARing4->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
    hPtPrimaryChargedBMDARing4->Draw();

    c25->cd(5);
    hPtPrimaryChargedBMDARing5->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
    hPtPrimaryChargedBMDARing5->GetYaxis()->SetTitle("Entries");
    hPtPrimaryChargedBMDARing5->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
    hPtPrimaryChargedBMDARing5->Draw();

    c25->cd(6);
    hPtPrimaryChargedBMDARing6->GetXaxis()->SetTitle("Pt (GeV/#it{c})");
    hPtPrimaryChargedBMDARing6->GetYaxis()->SetTitle("Entries");
    hPtPrimaryChargedBMDARing6->Scale(1.0/hMCPtPrimaryBMDA->GetEntries());
    hPtPrimaryChargedBMDARing6->Draw();

    c25->SaveAs("MC-pT-BEBE-A-Rings-9GeV.pdf");*/

//////////////////////////////////////////////////////////////////////////////

 

return 0;

}

