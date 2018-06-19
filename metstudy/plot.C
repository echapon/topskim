#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"

const double mumass = 105.6583745e-3;

void plot(const char* file) {
   TFile *f = TFile::Open(file);
   TTree *tr = (TTree*) f->Get("particleTree");
   Int_t           nMu = 0;
   vector<float>   *muPt = 0;
   vector<float>   *muEta = 0;
   vector<float>   *muPhi = 0;
   vector<int>     *muCharge = 0;
   tr->SetBranchAddress("nMu",&nMu);
   tr->SetBranchAddress("muPt",&muPt);
   tr->SetBranchAddress("muPhi",&muPhi);
   tr->SetBranchAddress("muEta",&muEta);
   tr->SetBranchAddress("muCharge",&muCharge);

   TH1F *h = new TH1F("hmass","dimuon mass;mass[GeV];Entries",100,0,200);

   int nentries = tr->GetEntries();
   for (int i=0; i<nentries; i++) {
      tr->GetEntry(i);
      if (nMu<2) continue;
      if (muCharge->at(0) == muCharge->at(1)) continue;
      TLorentzVector mu0; mu0.SetPtEtaPhiM(muPt->at(0), muEta->at(0), muPhi->at(0), mumass);
      TLorentzVector mu1; mu1.SetPtEtaPhiM(muPt->at(1), muEta->at(1), muPhi->at(1), mumass);
      TLorentzVector dimu = mu0+mu1;
      h->Fill(dimu.M());
   }

   h->Draw();
}
