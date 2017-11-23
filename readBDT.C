#include "TMVA/Reader.h"
#include "TTree.h"
#include "TFile.h"

void readBDT() {
   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used

   Float_t emuPt_;
   Float_t emuRap_;
   Float_t absemuRap_;
   Float_t emuPhi_;
   Float_t emuMass_;
   Int_t emuSgn_;
   Float_t emuPhistar_;
   Float_t emuDpt_;
   Float_t emuDptrel_;
   Float_t emuDeta_;
   Float_t emuDphi_;
   Float_t emuDR_;
   Int_t nJt_;
   Float_t discr_csvV1_[500];
   reader->AddVariable( "emuPt", &emuPt_ );
   reader->AddVariable( "abs(emuRap)", &absemuRap_ );
   reader->AddVariable( "emuMass", &emuMass_ );
   // reader->AddVariable( "emuPhistar", &emuPhistar_ );
   reader->AddVariable( "emuDpt/emuPt", &emuDptrel_ );
   reader->AddVariable( "emuDeta", &emuDeta_ );
   reader->AddVariable( "emuDphi", &emuDphi_ );
   // reader->AddVariable( "emuDR", &emuDR_ );
   reader->AddSpectator( "emuSgn", &emuSgn_ );
   reader->AddSpectator( "nJt", &nJt_ );
   reader->AddVariable( "discr_csvV1[0]", discr_csvV1_ );

   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   TString methodNameShort[5] {"BDTG","MLP","LikelihoodPCA","LikelihoodMIX", "LikelihoodKDE"};
   for (int i=0; i<5; i++) {
      TString methodName = methodNameShort[i] + TString("_bjet") + TString(" method");
      TString weightfile =  TString("weights/TMVAClassification_") + methodNameShort[i] + TString("_bjet") + TString(".weights.xml");
      reader->BookMVA( methodName, weightfile ); 
   }

   TFile *f = new TFile("emuSkim_data.root");
   TTree *tr = (TTree*) f->Get("skimTree");
   tr->SetBranchAddress( "emuPt", &emuPt_ );
   tr->SetBranchAddress( "emuRap", &emuRap_ );
   tr->SetBranchAddress( "emuMass", &emuMass_ );
   tr->SetBranchAddress( "emuPhistar", &emuPhistar_ );
   tr->SetBranchAddress( "emuDpt", &emuDpt_ );
   tr->SetBranchAddress( "emuDeta", &emuDeta_ );
   tr->SetBranchAddress( "emuDphi", &emuDphi_ );
   tr->SetBranchAddress( "emuDR", &emuDR_ );
   tr->SetBranchAddress( "emuSgn", &emuSgn_ );
   tr->SetBranchAddress( "nJt", &nJt_ );
   tr->SetBranchAddress( "discr_csvV1", discr_csvV1_ );

   TFile *f2 = new TFile("tree_dataAndBDT.root","RECREATE");
   TH1F *hbdtOS = new TH1F("hbdtOS","",40,-1,1);
   TH1F *hbdtSS1 = new TH1F("hbdtSS1","",40,-1,1);
   TH1F *hbdtSS2 = new TH1F("hbdtSS2","",40,-1,1);

   TTree *tr2 = tr->CloneTree(0);
   float bdtg, mlp, lhpca, lhmix, lhkde;
   tr2->Branch("bdtg",&bdtg,"bdtg/F");
   tr2->Branch("mlp",&mlp,"mlp/F");
   tr2->Branch("lhpca",&lhpca,"lhpca/F");
   tr2->Branch("lhmix",&lhmix,"lhmix/F");
   tr2->Branch("lhkde",&lhkde,"lhkde/F");
   for (int i=0; i<tr->GetEntries(); i++) {
      tr->GetEntry(i);
      if (emuPt_<=0) continue;
      if (nJt_<=0) continue;
      if (discr_csvV1_[0]<=-1) continue;
      
      absemuRap_ = fabs(emuRap_);
      emuDptrel_ = emuDpt_/emuPt_;
      bdtg = reader->EvaluateMVA("BDTG_bjet method");
      mlp = reader->EvaluateMVA("MLP_bjet method");
      lhpca = reader->EvaluateMVA("LikelihoodPCA_bjet method");
      lhmix = reader->EvaluateMVA("LikelihoodMIX_bjet method");
      lhkde = reader->EvaluateMVA("LikelihoodKDE_bjet method");
      
      float mvaval = lhpca;
      if (emuSgn_==0) hbdtOS->Fill(mvaval);
      else if(i<tr->GetEntries()/2) hbdtSS1->Fill(mvaval);
      else hbdtSS2->Fill(mvaval);

      tr2->Fill();
   }
   hbdtSS1->SetLineColor(kRed);;
   hbdtSS1->SetMarkerColor(kRed);;
   hbdtSS1->Draw("hist");
   hbdtSS2->SetLineColor(kBlue);;
   hbdtSS2->SetMarkerColor(kBlue);;
   hbdtSS2->Draw("hist same");
   hbdtOS->Draw("same");

   f2->Write();
   f2->Close();
}
