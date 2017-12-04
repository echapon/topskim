#define accCuts_cxx
#include "accCuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void accCuts::Loop()
{
//   In a ROOT session, you can do:
//      root> .L accCuts.C
//      root> accCuts t
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
   fChain->SetBranchStatus("*",0);  // disable all branches
   fChain->SetBranchStatus("genMuPt",1);  // activate branchname
   fChain->SetBranchStatus("genElePt",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH2F *h1 = new TH2F("h1",";p_{T}(e);p_{T}(#mu)",20,0,100,20,0,100);
   TH2F *h2 = new TH2F("h2",";p_{T}(l);p_{T}(#mu)+p_{T}(e)",20,0,100,20,0,100);

   int ntot=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if (genMuPt<0 || genElePt<0) continue;
      ntot++;

      // fill h1
      for (int i=0; i<20; i++)
         for (int j=0; j<20; j++) 
            if (genElePt>i*5&&genMuPt>j*5)
               h1->Fill(i*5+1,j*5+1);

      // fill h2
      for (int i=0; i<20; i++)
         for (int j=0; j<20; j++) 
            if (genElePt>i*5&&genMuPt>i*5&&genElePt+genMuPt>j*5)
               h2->Fill(i*5+1,j*5+1);
   }

   h1->Scale(1./ntot);
   h2->Scale(1./ntot);

   TCanvas *c1 = new TCanvas(); h1->Draw("COLZ");
   TCanvas *c2 = new TCanvas(); h2->Draw("COLZ");
}
