#ifndef LepJetsSkimTree_h
#define LepJetsSkimTree_h

#include "TTree.h"

#include <iostream>

TTree* skimTree_p = 0;

const Int_t nLep = 4;
const Int_t eleID = 11;
const Int_t muID = 13;
const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;

UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t vz_;
Int_t passed_;

Int_t nLep_;
Int_t lepID_[nLep];
Float_t lepPt_[nLep];
Float_t lepPhi_[nLep];
Float_t lepEta_[nLep];
Int_t lepChg_[nLep];
Float_t lepIso_[nLep];
Float_t lepIsoTk_[nLep];
Float_t lepDRjet_[nLep];
Float_t lepDRgen_[nLep];
Float_t lepPtrel_[nLep];
Float_t lepPtrel2_[nLep];
Float_t lepJetPt_[nLep];

Int_t   emuFlav_;
Float_t emuPt_;
Float_t emuRap_;
Float_t emuPhi_;
Float_t emuMass_;
Int_t emuSgn_;
Float_t emuPhistar_;
Float_t emuDpt_;
Float_t emuDeta_;
Float_t emuDphi_;
Float_t emuDR_;
Float_t emuMinLepPt_;
Float_t emuSumLepPt_;
Float_t emuMinDRjet_;
Float_t emuMinPtrel_;
Float_t emuMinPtrel2_;
Float_t emuMinJetPt_;
Float_t MHT_;

Float_t genMuPt_;
Float_t genMuEta_;
Float_t genMuPhi_;
Float_t genElePt_;
Float_t genEleEta_;
Float_t genElePhi_;

const Int_t nMaxJets = 500;
Int_t nJt_;
Float_t jtPt_[nMaxJets];
Float_t jtPhi_[nMaxJets];
Float_t jtEta_[nMaxJets];
Float_t jtM_[nMaxJets];
Float_t jtDRlep_[nMaxJets];
Float_t discr_csvV1_[nMaxJets];
int     refparton_flavorForB_[nMaxJets];

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
TBranch        *b_lepIso; //!
TBranch        *b_lepIsoTk; //!
TBranch        *b_lepDRjet; //!
TBranch        *b_lepDRgen; //!
TBranch        *b_lepPtrel; //!
TBranch        *b_lepPtrel2; //!
TBranch        *b_lepJetPt; //!
TBranch        *b_emuFlav; //!
TBranch        *b_emuPt; //!
TBranch        *b_emuRap; //!
TBranch        *b_emuPhi; //!
TBranch        *b_emuMass; //!
TBranch        *b_emuSgn; //!
TBranch        *b_emuPhistar; //!
TBranch        *b_emuDpt; //!
TBranch        *b_emuDeta; //!
TBranch        *b_emuDphi; //!
TBranch        *b_emuDR; //!
TBranch        *b_emuMinLepPt; //!
TBranch        *b_emuSumLepPt; //!
TBranch        *b_emuMinDRjet; //!
TBranch        *b_emuMinPtrel; //!
TBranch        *b_emuMinPtrel2; //!
TBranch        *b_emuMinJetPt; //!
TBranch        *b_genMuPt; //!
TBranch        *b_genMuEta; //!
TBranch        *b_genMuPhi; //!
TBranch        *b_genElePt; //!
TBranch        *b_genEleEta; //!
TBranch        *b_genElePhi; //!
TBranch        *b_MHT; //!
TBranch        *b_nJt;   //!
TBranch        *b_jtPt;   //!
TBranch        *b_jtPhi;   //!
TBranch        *b_jtEta;   //!
TBranch        *b_jtM;   //!
TBranch        *b_jtDRlep;   //!
TBranch        *b_discr_csvV1;//!
TBranch        *b_refparton_flavorForB;//!

void BookTree()
{
  if(skimTree_p == NULL){
    std::cout << "BOOKTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }

  skimTree_p->Branch("run", &run_, "run/i");
  skimTree_p->Branch("evt", &evt_, "evt/l");
  skimTree_p->Branch("lumi", &lumi_, "lumi/i");
  skimTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  skimTree_p->Branch("vz", &vz_, "vz/F");
  skimTree_p->Branch("passed", &passed_, "passed/I");

  skimTree_p->Branch("nLep", &nLep_, "nLep/I");
  skimTree_p->Branch("lepID", lepID_, Form("lepID[%d]/I",nLep));
  skimTree_p->Branch("lepPt", lepPt_, Form("lepPt[%d]/F", nLep));
  skimTree_p->Branch("lepPhi", lepPhi_, Form("lepPhi[%d]/F", nLep));
  skimTree_p->Branch("lepEta", lepEta_, Form("lepEta[%d]/F", nLep));
  skimTree_p->Branch("lepChg", lepChg_, Form("lepChg[%d]/I", nLep));
  skimTree_p->Branch("lepIso", lepIso_, Form("lepIso[%d]/F", nLep)); 
  skimTree_p->Branch("lepIsoTk", lepIsoTk_, Form("lepIsoTk[%d]/F", nLep)); 
  skimTree_p->Branch("lepDRjet", lepDRjet_, Form("lepDRjet[%d]/F", nLep)); 
  skimTree_p->Branch("lepDRgen", lepDRgen_, Form("lepDRgen[%d]/F", nLep)); 
  skimTree_p->Branch("lepPtrel", lepPtrel_, Form("lepPtrel[%d]/F", nLep)); 
  skimTree_p->Branch("lepPtrel2", lepPtrel2_, Form("lepPtrel2[%d]/F", nLep)); 
  skimTree_p->Branch("lepJetPt", lepJetPt_, Form("lepJetPt[%d]/F", nLep)); 

  skimTree_p->Branch("emuFlav", &emuFlav_, "emuFlav/I"); 
  skimTree_p->Branch("emuPt", &emuPt_, "emuPt/F"); 
  skimTree_p->Branch("emuRap", &emuRap_, "emuRap/F"); 
  skimTree_p->Branch("emuPhi", &emuPhi_, "emuPhi/F"); 
  skimTree_p->Branch("emuMass", &emuMass_, "emuMass/F"); 
  skimTree_p->Branch("emuSgn", &emuSgn_, "emuSgn/I"); 
  skimTree_p->Branch("emuPhistar", &emuPhistar_, "emuPhistar/F"); 
  skimTree_p->Branch("emuDpt", &emuDpt_, "emuDpt/F"); 
  skimTree_p->Branch("emuDeta", &emuDeta_, "emuDeta/F"); 
  skimTree_p->Branch("emuDphi", &emuDphi_, "emuDphi/F"); 
  skimTree_p->Branch("emuDR", &emuDR_, "emuDR/F"); 
  skimTree_p->Branch("emuMinLepPt", &emuMinLepPt_, "emuMinLepPt/F"); 
  skimTree_p->Branch("emuSumLepPt", &emuSumLepPt_, "emuSumLepPt/F"); 
  skimTree_p->Branch("emuMinDRjet", &emuMinDRjet_, "emuMinDRjet/F"); 
  skimTree_p->Branch("emuMinPtrel", &emuMinPtrel_, "emuMinPtrel/F"); 
  skimTree_p->Branch("emuMinPtrel2", &emuMinPtrel2_, "emuMinPtrel2/F"); 
  skimTree_p->Branch("emuMinJetPt", &emuMinJetPt_, "emuMinJetPt/F"); 
  skimTree_p->Branch("MHT", &MHT_, "MHT/F"); 

  skimTree_p->Branch("genMuPt", &genMuPt_, "genMuPt/F"); 
  skimTree_p->Branch("genMuEta", &genMuEta_, "genMuEta/F"); 
  skimTree_p->Branch("genMuPhi", &genMuPhi_, "genMuPhi/F"); 
  skimTree_p->Branch("genElePt", &genElePt_, "genElePt/F"); 
  skimTree_p->Branch("genEleEta", &genEleEta_, "genEleEta/F"); 
  skimTree_p->Branch("genElePhi", &genElePhi_, "genElePhi/F"); 
 
  skimTree_p->Branch("nJt", &nJt_, "nJt/I");
  skimTree_p->Branch("jtPt", jtPt_, "jtPt[nJt]/F");
  skimTree_p->Branch("jtPhi", jtPhi_, "jtPhi[nJt]/F");
  skimTree_p->Branch("jtEta", jtEta_, "jtEta[nJt]/F");
  skimTree_p->Branch("jtM", jtM_, "jtM[nJt]/F");
  skimTree_p->Branch("jtDRlep", jtDRlep_, "jtDRlep[nJt]/F");
  skimTree_p->Branch("discr_csvV1", discr_csvV1_, "discr_csvV1[nJt]/F");
  skimTree_p->Branch("refparton_flavorForB", refparton_flavorForB_, "refparton_flavorForB[nJt]/I");

  return;
}


void ReadTree()
{
  if(skimTree_p == NULL){
    std::cout << "READTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }

  /* skimTree_p->SetBranchAddress("run", &run, &b_run); */
  /* skimTree_p->SetBranchAddress("evt", &evt, &b_evt); */
  /* skimTree_p->SetBranchAddress("lumi", &lumi, &b_lumi); */
  /* skimTree_p->SetBranchAddress("hiBin", &hiBin, &b_hiBin); */
  /* skimTree_p->SetBranchAddress("vz", &vz, &b_vz); */
  /* skimTree_p->SetBranchAddress("nLep", &nLep, &b_nLep); */
  /* skimTree_p->SetBranchAddress("lepID", lepID, &b_lepID); */
  /* skimTree_p->SetBranchAddress("lepPt", lepPt, &b_lepPt); */
  /* skimTree_p->SetBranchAddress("lepPhi", lepPhi, &b_lepPhi); */
  /* skimTree_p->SetBranchAddress("lepEta", lepEta, &b_lepEta); */
  /* skimTree_p->SetBranchAddress("lepChg", lepChg, &b_lepChg); */
  /* skimTree_p->SetBranchAddress("nJt", &nJt, &b_nJt); */
  /* skimTree_p->SetBranchAddress("jtPt", jtPt, &b_jtPt); */
  /* skimTree_p->SetBranchAddress("jtPhi", jtPhi, &b_jtPhi); */
  /* skimTree_p->SetBranchAddress("jtEta", jtEta, &b_jtEta); */
  /* skimTree_p->SetBranchAddress("jtM", jtM, &b_jtM); */


  return;
}

#endif
