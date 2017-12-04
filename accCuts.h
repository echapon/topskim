//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  4 14:31:31 2017 by ROOT version 6.04/02
// from TTree skimTree/skimTree
// found on file: emuSkim_ttbarMC_0.root
//////////////////////////////////////////////////////////

#ifndef accCuts_h
#define accCuts_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class accCuts {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       evt;
   UInt_t          lumi;
   Int_t           hiBin;
   Float_t         vz;
   Int_t           passed;
   Int_t           nLep;
   Int_t           lepID[4];
   Float_t         lepPt[4];
   Float_t         lepPhi[4];
   Float_t         lepEta[4];
   Int_t           lepChg[4];
   Float_t         lepIso[4];
   Float_t         lepIsoTk[4];
   Float_t         lepDRjet[4];
   Float_t         lepDRgen[4];
   Float_t         lepPtrel[4];
   Float_t         lepPtrel2[4];
   Float_t         lepJetPt[4];
   Float_t         emuPt;
   Float_t         emuRap;
   Float_t         emuPhi;
   Float_t         emuMass;
   Int_t           emuSgn;
   Float_t         emuPhistar;
   Float_t         emuDpt;
   Float_t         emuDeta;
   Float_t         emuDphi;
   Float_t         emuDR;
   Float_t         emuMinLepPt;
   Float_t         emuSumLepPt;
   Float_t         emuMinDRjet;
   Float_t         emuMinPtrel;
   Float_t         emuMinPtrel2;
   Float_t         emuMinJetPt;
   Float_t         MHT;
   Float_t         genMuPt;
   Float_t         genMuEta;
   Float_t         genMuPhi;
   Float_t         genElePt;
   Float_t         genEleEta;
   Float_t         genElePhi;
   Int_t           nJt;
   Float_t         jtPt[60];   //[nJt]
   Float_t         jtPhi[60];   //[nJt]
   Float_t         jtEta[60];   //[nJt]
   Float_t         jtM[60];   //[nJt]
   Float_t         discr_csvV1[60];   //[nJt]
   Int_t           refparton_flavorForB[60];   //[nJt]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_passed;   //!
   TBranch        *b_nLep;   //!
   TBranch        *b_lepID;   //!
   TBranch        *b_lepPt;   //!
   TBranch        *b_lepPhi;   //!
   TBranch        *b_lepEta;   //!
   TBranch        *b_lepChg;   //!
   TBranch        *b_lepIso;   //!
   TBranch        *b_lepIsoTk;   //!
   TBranch        *b_lepDRjet;   //!
   TBranch        *b_lepDRgen;   //!
   TBranch        *b_lepPtrel;   //!
   TBranch        *b_lepPtrel2;   //!
   TBranch        *b_lepJetPt;   //!
   TBranch        *b_emuPt;   //!
   TBranch        *b_emuRap;   //!
   TBranch        *b_emuPhi;   //!
   TBranch        *b_emuMass;   //!
   TBranch        *b_emuSgn;   //!
   TBranch        *b_emuPhistar;   //!
   TBranch        *b_emuDpt;   //!
   TBranch        *b_emuDeta;   //!
   TBranch        *b_emuDphi;   //!
   TBranch        *b_emuDR;   //!
   TBranch        *b_emuMinLepPt;   //!
   TBranch        *b_emuSumLepPt;   //!
   TBranch        *b_emuMinDRjet;   //!
   TBranch        *b_emuMinPtrel;   //!
   TBranch        *b_emuMinPtrel2;   //!
   TBranch        *b_emuMinJetPt;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_genMuPt;   //!
   TBranch        *b_genMuEta;   //!
   TBranch        *b_genMuPhi;   //!
   TBranch        *b_genElePt;   //!
   TBranch        *b_genEleEta;   //!
   TBranch        *b_genElePhi;   //!
   TBranch        *b_nJt;   //!
   TBranch        *b_jtPt;   //!
   TBranch        *b_jtPhi;   //!
   TBranch        *b_jtEta;   //!
   TBranch        *b_jtM;   //!
   TBranch        *b_discr_csvV1;   //!
   TBranch        *b_refparton_flavorForB;   //!

   accCuts(TTree *tree=0);
   virtual ~accCuts();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef accCuts_cxx
accCuts::accCuts(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("emuSkim_ttbarMC_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("emuSkim_ttbarMC_0.root");
      }
      f->GetObject("skimTree",tree);

   }
   Init(tree);
}

accCuts::~accCuts()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t accCuts::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t accCuts::LoadTree(Long64_t entry)
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

void accCuts::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("passed", &passed, &b_passed);
   fChain->SetBranchAddress("nLep", &nLep, &b_nLep);
   fChain->SetBranchAddress("lepID", lepID, &b_lepID);
   fChain->SetBranchAddress("lepPt", lepPt, &b_lepPt);
   fChain->SetBranchAddress("lepPhi", lepPhi, &b_lepPhi);
   fChain->SetBranchAddress("lepEta", lepEta, &b_lepEta);
   fChain->SetBranchAddress("lepChg", lepChg, &b_lepChg);
   fChain->SetBranchAddress("lepIso", lepIso, &b_lepIso);
   fChain->SetBranchAddress("lepIsoTk", lepIsoTk, &b_lepIsoTk);
   fChain->SetBranchAddress("lepDRjet", lepDRjet, &b_lepDRjet);
   fChain->SetBranchAddress("lepDRgen", lepDRgen, &b_lepDRgen);
   fChain->SetBranchAddress("lepPtrel", lepPtrel, &b_lepPtrel);
   fChain->SetBranchAddress("lepPtrel2", lepPtrel2, &b_lepPtrel2);
   fChain->SetBranchAddress("lepJetPt", lepJetPt, &b_lepJetPt);
   fChain->SetBranchAddress("emuPt", &emuPt, &b_emuPt);
   fChain->SetBranchAddress("emuRap", &emuRap, &b_emuRap);
   fChain->SetBranchAddress("emuPhi", &emuPhi, &b_emuPhi);
   fChain->SetBranchAddress("emuMass", &emuMass, &b_emuMass);
   fChain->SetBranchAddress("emuSgn", &emuSgn, &b_emuSgn);
   fChain->SetBranchAddress("emuPhistar", &emuPhistar, &b_emuPhistar);
   fChain->SetBranchAddress("emuDpt", &emuDpt, &b_emuDpt);
   fChain->SetBranchAddress("emuDeta", &emuDeta, &b_emuDeta);
   fChain->SetBranchAddress("emuDphi", &emuDphi, &b_emuDphi);
   fChain->SetBranchAddress("emuDR", &emuDR, &b_emuDR);
   fChain->SetBranchAddress("emuMinLepPt", &emuMinLepPt, &b_emuMinLepPt);
   fChain->SetBranchAddress("emuSumLepPt", &emuSumLepPt, &b_emuSumLepPt);
   fChain->SetBranchAddress("emuMinDRjet", &emuMinDRjet, &b_emuMinDRjet);
   fChain->SetBranchAddress("emuMinPtrel", &emuMinPtrel, &b_emuMinPtrel);
   fChain->SetBranchAddress("emuMinPtrel2", &emuMinPtrel2, &b_emuMinPtrel2);
   fChain->SetBranchAddress("emuMinJetPt", &emuMinJetPt, &b_emuMinJetPt);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("genMuPt", &genMuPt, &b_genMuPt);
   fChain->SetBranchAddress("genMuEta", &genMuEta, &b_genMuEta);
   fChain->SetBranchAddress("genMuPhi", &genMuPhi, &b_genMuPhi);
   fChain->SetBranchAddress("genElePt", &genElePt, &b_genElePt);
   fChain->SetBranchAddress("genEleEta", &genEleEta, &b_genEleEta);
   fChain->SetBranchAddress("genElePhi", &genElePhi, &b_genElePhi);
   fChain->SetBranchAddress("nJt", &nJt, &b_nJt);
   fChain->SetBranchAddress("jtPt", jtPt, &b_jtPt);
   fChain->SetBranchAddress("jtPhi", jtPhi, &b_jtPhi);
   fChain->SetBranchAddress("jtEta", jtEta, &b_jtEta);
   fChain->SetBranchAddress("jtM", jtM, &b_jtM);
   fChain->SetBranchAddress("discr_csvV1", discr_csvV1, &b_discr_csvV1);
   fChain->SetBranchAddress("refparton_flavorForB", refparton_flavorForB, &b_refparton_flavorForB);
   Notify();
}

Bool_t accCuts::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void accCuts::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t accCuts::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef accCuts_cxx
