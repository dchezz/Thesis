//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 29 15:27:20 2020 by ROOT version 6.16/00
// from TTree phsd/
// found on file: phsd.root
//////////////////////////////////////////////////////////

#ifndef Sim3_h
#define Sim3_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Sim3 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           isub;
   Int_t           irun;
   Float_t         b;
   Float_t         ibw;
   Float_t         psi[4];
   Float_t         epsilon[4];
   Int_t           np;
   Int_t           n;
   Int_t           id[2164];   //[n]
   Short_t         q[2164];   //[n]
   Float_t         e[2164];   //[n]
   Float_t         px[2164];   //[n]
   Float_t         py[2164];   //[n]
   Float_t         pz[2164];   //[n]
   Int_t           code1[2164];   //[n]
   Int_t           code2[2164];   //[n]

   // List of branches
   TBranch        *b_isub;   //!
   TBranch        *b_irun;   //!
   TBranch        *b_b;   //!
   TBranch        *b_ibw;   //!
   TBranch        *b_psi;   //!
   TBranch        *b_epsilon;   //!
   TBranch        *b_np;   //!
   TBranch        *b_n;   //!
   TBranch        *b_id;   //!
   TBranch        *b_q;   //!
   TBranch        *b_e;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_code1;   //!
   TBranch        *b_code2;   //!

   Sim3(TTree *tree=0);
   virtual ~Sim3();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Sim3_cxx
Sim3::Sim3(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("phsd.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("phsd.root");
      }
      f->GetObject("phsd",tree);

   }
   Init(tree);
}

Sim3::~Sim3()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Sim3::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Sim3::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Sim3::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isub", &isub, &b_isub);
   fChain->SetBranchAddress("irun", &irun, &b_irun);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("ibw", &ibw, &b_ibw);
   fChain->SetBranchAddress("psi", psi, &b_psi);
   fChain->SetBranchAddress("epsilon", epsilon, &b_epsilon);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("n", &n, &b_n);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("code1", code1, &b_code1);
   fChain->SetBranchAddress("code2", code2, &b_code2);
   Notify();
}

Bool_t Sim3::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Sim3::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Sim3::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Sim3_cxx
