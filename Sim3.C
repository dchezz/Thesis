#define Sim3_cxx
#include "Sim3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>

void Sim3::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Sim3.C
//      root> Sim3 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
      //Long64_t ientry = LoadTree(jentry);
      //if (ientry < 0) break;
      //nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

//All particles, neutrons, protons and all charged particles
TH1F *heta1 = new TH1F("h1b","PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#eta;# entries",277,-8,8);
TH1F *heta2 = new TH1F("h2b","eta",277,-8,8);
TH1F *heta3 = new TH1F("h3b","eta",277,-8,8);
TH1F *heta4 = new TH1F("h4b","eta",277,-8,8);

TH1F *htheta1 = new TH1F("h1a","PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#theta;# entries",277,-0.1,3.15);
TH1F *htheta2= new TH1F("h2a","theta",277,-0.1,3.15);
TH1F *htheta3= new TH1F("h3a","theta",277,-0.1,3.15);
TH1F *htheta4= new TH1F("h4a","theta",277,-0.1,3.15);

//Pions
TH1F *pheta1 = new TH1F("h1c","PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#eta;# entries",277,-8,8);
TH1F *pheta2 = new TH1F("h2c","eta",277,-8,8);

TH1F *phtheta1 = new TH1F("h1d","PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#theta;# entries",277,-0.1,3.15);
TH1F *phtheta2= new TH1F("h2d","theta",277,-0.1,3.15);

TH1F *hImpact = new TH1F("hb_10000", "PHSD Au+Au, 11 GeV, 10000 events, bmin=0 fm, bmax=14 fm;#b;# entries",100,-0.1,14.2);




for(Long64_t jentry=0;jentry<nentries;jentry++) {
GetEntry(jentry); //Loop over events

hImpact->Fill(b);

for(Int_t j=0;j<n;j++) { //Loop over particles

Float_t pt[2164];
Float_t theta[2164];
Float_t eta[2164];
pt[j] = sqrt(px[j]*px[j]+py[j]*py[j]);
theta[j] = atan2(pt[j],pz[j]);
eta[j] = -log(tan(theta[j]/2));

  heta1->Fill(eta[j]);
  htheta1->Fill(theta[j]);
if(id[j]==2212) {
  heta2->Fill(eta[j]);
  htheta2->Fill(theta[j]);
     }
if(id[j]==2112) {
  heta3->Fill(eta[j]);
  htheta3->Fill(theta[j]);
     }
if(q[j]!=0) {
  heta4->Fill(eta[j]);
  htheta4->Fill(theta[j]);
     }
if(id[j]==-211)  {
  pheta1->Fill(eta[j]);
  phtheta1->Fill(theta[j]);
     }
if(id[j]==211)  {
  pheta2->Fill(eta[j]);
  phtheta2->Fill(theta[j]);
     }
   
   }
 }

TCanvas *c1=new TCanvas("c1","c1",1080,720);
c1->cd();

htheta1->SetMarkerColor(kBlue);
htheta1->SetMarkerStyle(kFullCross);
htheta1-> Scale(1./(htheta1->GetEntries()));
htheta1->Draw();
htheta2->SetMarkerColor(kGreen);
htheta2->SetMarkerStyle(kFullStar);
htheta2-> Scale(1./(htheta1->GetEntries()));
htheta2->Draw("Same");
htheta3->SetMarkerColor(kMagenta);
htheta3->SetMarkerStyle(kCircle);
htheta3-> Scale(1./(htheta1->GetEntries()));
htheta3->Draw("Same");
htheta4->SetMarkerColor(kBlack);
htheta4->SetMarkerStyle(kFullDotLarge);
htheta4-> Scale(1./(htheta1->GetEntries()));
htheta4->Draw("Same");

TLegend *legend=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend->SetHeader("");
legend->AddEntry(htheta1,"all particles","p");
legend->AddEntry(htheta2,"protons","p");
legend->AddEntry(htheta3,"neutrons","p");
legend->AddEntry(htheta4,"charged particles","p");
legend->Draw();


TCanvas *c2=new TCanvas("c2","c2",1080,720);
c2->cd();

heta1->SetMarkerColor(kBlue);
heta1->SetMarkerStyle(kFullCross);
heta1-> Scale(1./(htheta1->GetEntries()));
heta1->Draw();
heta2->SetMarkerColor(kGreen);
heta2->SetMarkerStyle(kFullStar);
heta2-> Scale(1./(htheta1->GetEntries()));
heta2->Draw("Same");
heta3->SetMarkerColor(kMagenta);
heta3->SetMarkerStyle(kCircle);
heta3-> Scale(1./(htheta1->GetEntries()));
heta3->Draw("Same");
heta4->SetMarkerColor(kBlack);
heta4->SetMarkerStyle(kFullDotLarge);
heta4-> Scale(1./(htheta1->GetEntries()));
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
pheta1->Scale(1./(htheta1->GetEntries()));
pheta1->Draw();
pheta2->SetMarkerColor(kBlack);
pheta2->SetMarkerStyle(kFullStar);
pheta2->Scale(1./(htheta1->GetEntries()));
pheta2->Draw("Same");

TLegend *legend3=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend3->SetHeader("");
legend3->AddEntry(pheta1,"#pi^{-}","p");
legend3->AddEntry(pheta2,"#pi^{+}","p");
legend3->Draw();


TCanvas *c4=new TCanvas("c4","c4",1080,720);
c4->cd();

phtheta1->SetMarkerColor(kBlue);
phtheta1->SetMarkerStyle(kFullDotLarge);
phtheta1->Scale(1./(htheta1->GetEntries()));
phtheta1->Draw();
phtheta2->SetMarkerColor(kBlack);
phtheta2->SetMarkerStyle(kFullStar);
phtheta2->Scale(1./(htheta1->GetEntries()));
phtheta2->Draw("Same");

TLegend *legend4=new TLegend(0.779055,0.572874,0.979566,0.77327);
legend4->SetHeader("");
legend4->AddEntry(phtheta1,"#pi^{-}","p");
legend4->AddEntry(phtheta2,"#pi^{+}","p");
legend4->Draw();

TCanvas *c5=new TCanvas("c5","c5",1080,720);
c5->cd();

hImpact->SetLineColor(kBlack);
hImpact->SetLineWidth(1);
hImpact->Scale(1./1000);
hImpact->Draw();


}
