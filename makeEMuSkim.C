#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include "LepJetsSkimTree.h"
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "ForestTreeHeaders/ForestElectrons.h"
#include "ForestTreeHeaders/ForestMuons.h"

#include "Helpers/EnergyRegression.h"

#define eleMass 0.5109989461e-3
#define muMass  105.6583745e-3

using namespace std;

const bool isDebug = false;

const float jetPtCut = 25.;
const float jetEtaCut = 2.;

const float lepPtrelCut = -1;//5;
const float lepDRjetCut = 0.2;

// FIXME: Need to check lepton selections for PbPb
const float muEtaCut = 2.4;
const float muPtCut = 15;

// Run-II tight cuts (https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon)
const float muChi2NDFCut   = 10;
const float muInnerD0Cut   = 0.3;// same as HIN-16-004, was 0.2
const float muInnerDzCut   = 20;// same as HIN-16-004, was 0.5
const int   muMuonHitsCut  = 0;
const int   muStationsCut  = 1;
const int   muTrkLayersCut = 5;
const int   muPixelHitsCut = 0;

//PbPb: https://twiki.cern.ch/twiki/bin/view/CMS/ElectronPbPb5TeV#3_Selection_for_different_centra
//pp:   https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
//WARNING: this code has PbPb selections activated
const float eleEtaCut = 2.4;
const float elePtCut = 15;

const float barrelEndcapEta = 1.479;
const int nBarrelEndcap = 2; // barrel = 0, endcap = 1
//For PbPb electron ID is centrality dependent
//centrality bins 0-20%, 20-40%, 40-70%, 70-100%
const int nCentEleId = 4;
double centMinEleId[4] = {0.,20.,40.,70.};
double centMaxEleId[4] = {20.,40.,70.,100.};

float eleSigmaIEtaIEta_VetoCut[nBarrelEndcap][nCentEleId];
float eleDEtaIn_VetoCut[nBarrelEndcap][nCentEleId];
float eleDPhiIn_VetoCut[nBarrelEndcap][nCentEleId];
float eleHOverE_VetoCut[nBarrelEndcap][nCentEleId];
float eleRelIsoWithEA_VetoCut[nBarrelEndcap][nCentEleId];
float eleOOEmooP_VetoCut[nBarrelEndcap][nCentEleId];
float eleD0_VetoCut[nBarrelEndcap][nCentEleId];
float eleDZ_VetoCut[nBarrelEndcap][nCentEleId];
float eleMissingInnerHits_VetoCut[nBarrelEndcap][nCentEleId];
float eleEoverPInv_VetoCut[nBarrelEndcap][nCentEleId];

// const float eleSigmaIEtaIEta_VetoCut[nBarrelEndcap] = {0.0114, 0.0352};
// const float eleDEtaIn_VetoCut[nBarrelEndcap] = {0.0152, 0.0113};
// const float eleDPhiIn_VetoCut[nBarrelEndcap] = {0.216, 0.237};
// const float eleHOverE_VetoCut[nBarrelEndcap] = {0.181, 0.116};
// const float eleRelIsoWithEA_VetoCut[nBarrelEndcap] = {0.126, 0.144};
// const float eleOOEmooP_VetoCut[nBarrelEndcap] = {0.207, 0.174};
// const float eleD0_VetoCut[nBarrelEndcap] = {0.0564, 0.222};
// const float eleDZ_VetoCut[nBarrelEndcap] = {0.472, 0.921};
// const float eleMissingInnerHits_VetoCut[nBarrelEndcap] = {2, 3};

int getCentBinEleId(double cent) {

  int centBin = -1;
  for(int i = 0; i<nCentEleId; ++i) {
    if(cent>=centMinEleId[i] && cent<centMaxEleId[i]) centBin = i;
  }
  return centBin;
}

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi, std::vector<int> *pfId, bool tkonly=false);

double phistar(TLorentzVector v1, TLorentzVector v2);

float calcLeptonMinDRjet(float lepPt, float lepEta, float lepPhi, int njt, float *jtpt, float *jteta, float *jtphi, float &lepPtrel, float &lepPtrel2, float &jetpt);


void makeEMuSkim(const std::string outFileName = "", const std::string inFileName = "", bool isMC = false)
{
   if(!strcmp(inFileName.c_str(), "")){
      std::cout << "No inputs specified. return" << std::endl;
      return;
   }

   if(isDebug) std::cout << __LINE__ << std::endl;

   TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
   skimTree_p = new TTree("skimTree", "skimTree");
   BookTree();

   std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
   inFileNames_p->push_back(inFileName);
   // if(strcmp(inFileName.c_str(), "") != 0) inFileNames_p->push_back(inFileName);

   //  Printf("Create EnergyRegression object");
   EnergyRegression energyRegression;
   //Printf("Create EnergyRegression object done. calling initreader");
   if(!isMC) energyRegression.initreader();

   //initialize electron Id variables [barrel/endcap][centBin]
   eleSigmaIEtaIEta_VetoCut[0][0] = 0.01325;
   eleSigmaIEtaIEta_VetoCut[1][0] = 0.04272;
   eleSigmaIEtaIEta_VetoCut[0][1] = 0.01098;
   eleSigmaIEtaIEta_VetoCut[1][1] = 0.03569;
   eleSigmaIEtaIEta_VetoCut[0][2] = 0.01078;
   eleSigmaIEtaIEta_VetoCut[1][2] = 0.03155;
   eleSigmaIEtaIEta_VetoCut[0][3] = 0.01038;
   eleSigmaIEtaIEta_VetoCut[1][3] = 0.02955;

   eleDEtaIn_VetoCut[0][0] = 0.04524;
   eleDEtaIn_VetoCut[1][0] = 0.42269;
   eleDEtaIn_VetoCut[0][1] = 0.02799;
   eleDEtaIn_VetoCut[1][1] = 0.01592;
   eleDEtaIn_VetoCut[0][2] = 0.01275;
   eleDEtaIn_VetoCut[1][2] = 0.01074;
   eleDEtaIn_VetoCut[0][3] = 0.00595;
   eleDEtaIn_VetoCut[1][3] = 0.00927;

   eleDPhiIn_VetoCut[0][0] = 0.15133;
   eleDPhiIn_VetoCut[1][0] = 0.33656;
   eleDPhiIn_VetoCut[0][1] = 0.10697;
   eleDPhiIn_VetoCut[1][1] = 0.18786;
   eleDPhiIn_VetoCut[0][2] = 0.11228;
   eleDPhiIn_VetoCut[1][2] = 0.12940;
   eleDPhiIn_VetoCut[0][3] = 0.22180;
   eleDPhiIn_VetoCut[1][3] = 0.20464;

   eleHOverE_VetoCut[0][0] = 0.12879;
   eleHOverE_VetoCut[1][0] = 0.14855;
   eleHOverE_VetoCut[0][1] = 0.09844;
   eleHOverE_VetoCut[1][1] = 0.11125;
   eleHOverE_VetoCut[0][2] = 0.02355;
   eleHOverE_VetoCut[1][2] = 0.05202;
   eleHOverE_VetoCut[0][3] = 0.02997;
   eleHOverE_VetoCut[1][3] = 0.01670;

   eleD0_VetoCut[0][0] = 0.18115;
   eleD0_VetoCut[1][0] = 0.12069;
   eleD0_VetoCut[0][1] = 0.06520;
   eleD0_VetoCut[1][1] = 0.16610;
   eleD0_VetoCut[0][2] = 0.05574;
   eleD0_VetoCut[1][2] = 0.14651;
   eleD0_VetoCut[0][3] = 0.05396;
   eleD0_VetoCut[1][3] = 0.12990;

   eleDZ_VetoCut[0][0] = 0.11313;
   eleDZ_VetoCut[1][0] = 0.28022;
   eleDZ_VetoCut[0][1] = 0.06983;
   eleDZ_VetoCut[1][1] = 0.24015;
   eleDZ_VetoCut[0][2] = 0.02011;
   eleDZ_VetoCut[1][2] = 0.18170;
   eleDZ_VetoCut[0][3] = 0.06513;
   eleDZ_VetoCut[1][3] = 0.24262;

   eleMissingInnerHits_VetoCut[0][0] = 1.00005;
   eleMissingInnerHits_VetoCut[1][0] = 1.00005;
   eleMissingInnerHits_VetoCut[0][1] = 1.00005;
   eleMissingInnerHits_VetoCut[1][1] = 1.00005;
   eleMissingInnerHits_VetoCut[0][2] = 1.00005;
   eleMissingInnerHits_VetoCut[1][2] = 1.00005;
   eleMissingInnerHits_VetoCut[0][3] = 1.00005;
   eleMissingInnerHits_VetoCut[1][3] = 1.00005;

   eleEoverPInv_VetoCut[0][0] = 0.20965;
   eleEoverPInv_VetoCut[1][0] = 0.20556;
   eleEoverPInv_VetoCut[0][1] = 0.12291;
   eleEoverPInv_VetoCut[1][1] = 0.29670;
   eleEoverPInv_VetoCut[0][2] = 0.30592;
   eleEoverPInv_VetoCut[1][2] = 0.20473;
   eleEoverPInv_VetoCut[0][3] = 0.23732;
   eleEoverPInv_VetoCut[1][3] = 0.11148;

   TChain *lepTree_p = new TChain("ggHiNtuplizer/EventTree");
   TChain *jetTree_p = new TChain("akCs2PFJetAnalyzer/t");
   TChain *hiTree_p = new TChain("hiEvtAnalyzer/HiTree");
   TChain *hltTree_p = new TChain("hltanalysis/HltTree");
   TChain *pfTree_p = new TChain("pfcandAnalyzerCS/pfTree");
   TChain *skimAnaTree_p = new TChain("skimanalysis/HltTree");

   const int nFiles = (int)inFileNames_p->size();

   for(int fileIter = 0; fileIter < nFiles; fileIter++){
      std::cout << "On file: " << fileIter << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
      lepTree_p->Add(inFileNames_p->at(fileIter).c_str());
      jetTree_p->Add(inFileNames_p->at(fileIter).c_str());
      hiTree_p->Add(inFileNames_p->at(fileIter).c_str());
      hltTree_p->Add(inFileNames_p->at(fileIter).c_str());
      pfTree_p->Add(inFileNames_p->at(fileIter).c_str());
      skimAnaTree_p->Add(inFileNames_p->at(fileIter).c_str());
   }
   // TFile *inFile_p = TFile::Open(inFileNames_p->at(fileIter).c_str(), "READ");
   // TTree *lepTree_p = dynamic_cast<TTree*>(inFile_p->Get("ggHiNtuplizer/EventTree"));
   // TTree *jetTree_p = dynamic_cast<TTree*>(inFile_p->Get("akCs2PFJetAnalyzer/t"));
   // TTree *hiTree_p = dynamic_cast<TTree*>(inFile_p->Get("hiEvtAnalyzer/HiTree"));
   // TTree *hltTree_p = dynamic_cast<TTree*>(inFile_p->Get("hltanalysis/HltTree"));
   // TTree *pfTree_p = dynamic_cast<TTree*>(inFile_p->Get("pfcandAnalyzerCS/pfTree"));
   // TTree *skimAnaTree_p = dynamic_cast<TTree*>(inFile_p->Get("skimanalysis/HltTree"));

   ForestMuons fForestMu;
   ForestElectrons fForestEle;

   const int maxJets = 5000;
   int           nref;
   float         jtpt[maxJets];   //[nref]
   float         jteta[maxJets];   //[nref]
   float         jtphi[maxJets];   //[nref]
   float         jtm[maxJets];   //[nref]
   float         discr_csvV1[maxJets]; //[nref]
   int           refparton_flavorForB[maxJets]; //[nref]

   //pf particles pfId, pfPt, pfEta, pfPhi
   std::vector<int>           *pfId = 0;
   std::vector<float>         *pfPt = 0;
   std::vector<float>         *pfEta = 0;
   std::vector<float>         *pfPhi = 0;

   // MC
   Int_t           nMC;
   vector<int>     *mcPID = 0;
   vector<float>   *mcPt = 0;
   vector<float>   *mcEta = 0;
   vector<float>   *mcPhi = 0;
   vector<float>   *mcMass = 0;
   vector<int>     *mcMomPID = 0;

   int trig = 1;

   //event selections
   int phfCoincFilter = 1;
   int HBHENoiseFilterResult = 1;
   int pprimaryVertexFilter = 1;
   int pcollisionEventSelection = 1;

   lepTree_p->SetBranchStatus("mu*", 1);
   lepTree_p->SetBranchStatus("mc*", 1);

   lepTree_p->SetBranchAddress("nMC", &nMC);
   lepTree_p->SetBranchAddress("mcPID", &mcPID);
   lepTree_p->SetBranchAddress("mcPt", &mcPt);
   lepTree_p->SetBranchAddress("mcEta", &mcEta);
   lepTree_p->SetBranchAddress("mcPhi", &mcPhi);
   lepTree_p->SetBranchAddress("mcMass", &mcMass);
   lepTree_p->SetBranchAddress("mcMomPID", &mcMomPID);

   lepTree_p->SetBranchAddress("muPt", &fForestMu.muPt);
   lepTree_p->SetBranchAddress("muPhi", &fForestMu.muPhi);
   lepTree_p->SetBranchAddress("muEta", &fForestMu.muEta);
   lepTree_p->SetBranchAddress("muCharge", &fForestMu.muCharge);
   lepTree_p->SetBranchAddress("muChi2NDF", &fForestMu.muChi2NDF);
   lepTree_p->SetBranchAddress("muInnerD0", &fForestMu.muInnerD0);
   lepTree_p->SetBranchAddress("muInnerDz", &fForestMu.muInnerDz);
   lepTree_p->SetBranchAddress("muMuonHits", &fForestMu.muMuonHits);
   lepTree_p->SetBranchAddress("muStations", &fForestMu.muStations);
   lepTree_p->SetBranchAddress("muTrkLayers", &fForestMu.muTrkLayers);
   lepTree_p->SetBranchAddress("muPixelHits", &fForestMu.muPixelHits);    


   lepTree_p->SetBranchStatus("ele*", 1);

   lepTree_p->SetBranchAddress("elePt", &fForestEle.elePt);
   lepTree_p->SetBranchAddress("elePhi", &fForestEle.elePhi);
   lepTree_p->SetBranchAddress("eleEta", &fForestEle.eleEta);
   lepTree_p->SetBranchAddress("eleCharge", &fForestEle.eleCharge);
   lepTree_p->SetBranchAddress("eleSigmaIEtaIEta", &fForestEle.eleSigmaIEtaIEta);
   lepTree_p->SetBranchAddress("eledEtaAtVtx", &fForestEle.eledEtaAtVtx);
   lepTree_p->SetBranchAddress("eledPhiAtVtx", &fForestEle.eledPhiAtVtx);
   lepTree_p->SetBranchAddress("eleHoverE", &fForestEle.eleHoverE);
   lepTree_p->SetBranchAddress("eleD0", &fForestEle.eleD0);
   lepTree_p->SetBranchAddress("eleDz", &fForestEle.eleDz);
   lepTree_p->SetBranchAddress("eleEoverPInv", &fForestEle.eleEoverPInv);

   lepTree_p->SetBranchAddress("eleSCPhi", &fForestEle.eleSCPhi);
   lepTree_p->SetBranchAddress("eleSCEta", &fForestEle.eleSCEta);
   lepTree_p->SetBranchAddress("eleSigmaIEtaIEta_2012", &fForestEle.eleSigmaIEtaIEta_2012);
   lepTree_p->SetBranchAddress("eleSigmaIPhiIPhi", &fForestEle.eleSigmaIPhiIPhi);
   lepTree_p->SetBranchAddress("eleSCEtaWidth", &fForestEle.eleSCEtaWidth);
   lepTree_p->SetBranchAddress("eleSCPhiWidth", &fForestEle.eleSCPhiWidth);

   jetTree_p->SetBranchStatus("*", 0);
   jetTree_p->SetBranchStatus("nref", 1);
   jetTree_p->SetBranchStatus("jtpt", 1);
   jetTree_p->SetBranchStatus("jtphi", 1);
   jetTree_p->SetBranchStatus("jteta", 1);
   jetTree_p->SetBranchStatus("jtm", 1);
   jetTree_p->SetBranchStatus("discr_csvV1", 1);
   if (isMC) jetTree_p->SetBranchStatus("refparton_flavorForB", 1);

   jetTree_p->SetBranchAddress("nref", &nref);
   jetTree_p->SetBranchAddress("jtpt", jtpt);
   jetTree_p->SetBranchAddress("jtphi", jtphi);
   jetTree_p->SetBranchAddress("jteta", jteta);
   jetTree_p->SetBranchAddress("jtm", jtm);
   jetTree_p->SetBranchAddress("discr_csvV1", discr_csvV1);
   if (isMC) jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);

   hiTree_p->SetBranchStatus("*", 0);
   hiTree_p->SetBranchStatus("run", 1);
   hiTree_p->SetBranchStatus("evt", 1);
   hiTree_p->SetBranchStatus("lumi", 1);
   hiTree_p->SetBranchStatus("hiBin", 1);
   hiTree_p->SetBranchStatus("vz", 1);

   hiTree_p->SetBranchAddress("run", &run_);
   hiTree_p->SetBranchAddress("evt", &evt_);
   hiTree_p->SetBranchAddress("lumi", &lumi_);
   hiTree_p->SetBranchAddress("hiBin", &hiBin_);
   hiTree_p->SetBranchAddress("vz", &vz_);

   pfTree_p->SetBranchAddress("pfId", &pfId);
   pfTree_p->SetBranchAddress("pfPt", &pfPt);
   pfTree_p->SetBranchAddress("pfEta", &pfEta);
   pfTree_p->SetBranchAddress("pfPhi", &pfPhi);

   hltTree_p->SetBranchStatus("HLT_HIL2Mu15_v2",1);
   hltTree_p->SetBranchAddress("HLT_HIL2Mu15_v2",&trig);

   skimAnaTree_p->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
   skimAnaTree_p->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
   skimAnaTree_p->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
   skimAnaTree_p->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);

   if(isDebug) std::cout << __LINE__ << std::endl;

   if(isDebug) std::cout << __LINE__ << std::endl;

   int nEntries = (int)lepTree_p->GetEntries();
   //nEntries = 100;
   int entryDiv = ((int)(nEntries/20));

   if(isDebug) std::cout << __LINE__ << std::endl;

   map<string,int> cnt;
   cnt["tot"]=0;
   cnt["filt1"]=0;
   cnt["filt2"]=0;
   cnt["filt3"]=0;
   cnt["filt4"]=0;
   cnt["filt5"]=0;
   cnt["filt6"]=0;
   cnt["gene"]=0;
   cnt["genm"]=0;
   cnt["genem"]=0;
   cnt["rece"]=0;
   cnt["recm"]=0;
   cnt["recem"]=0;
   cnt["recee"]=0;
   cnt["recmm"]=0;
   cnt["recj"]=0;
   cnt["recbj"]=0;
   cnt["recbq"]=0;
   cnt["mu0"]=0;
   cnt["mu1"]=0;
   cnt["mu2"]=0;
   cnt["mu3"]=0;
   cnt["mu4"]=0;
   cnt["mu5"]=0;
   cnt["mu6"]=0;
   cnt["mu7"]=0;
   cnt["mu8"]=0;
   cnt["mu9"]=0;
   cnt["mu10"]=0;
   cnt["ele0"]=0;
   cnt["ele1"]=0;
   cnt["ele2"]=0;
   cnt["ele3"]=0;
   cnt["ele4"]=0;
   cnt["ele5"]=0;
   cnt["ele6"]=0;
   cnt["ele7"]=0;
   cnt["ele8"]=0;
   cnt["ele9"]=0;
   cnt["ele10"]=0;
   cnt["ele11"]=0;
   cnt["matchm0"]=0;
   cnt["matchm1"]=0;
   cnt["matchm2"]=0;
   cnt["matchm3"]=0;
   cnt["matchm4"]=0;
   cnt["matchm5"]=0;
   cnt["matchm6"]=0;
   cnt["matchm7"]=0;
   cnt["matchm8"]=0;
   cnt["matchm9"]=0;
   cnt["matchm10"]=0;
   cnt["matche0"]=0;
   cnt["matche1"]=0;
   cnt["matche2"]=0;
   cnt["matche3"]=0;
   cnt["matche4"]=0;
   cnt["matche5"]=0;
   cnt["matche6"]=0;
   cnt["matche7"]=0;
   cnt["matche8"]=0;
   cnt["matche9"]=0;
   cnt["matche10"]=0;
   cnt["matche11"]=0;

   for(int entry = 0; entry < nEntries; entry++){
      if(isDebug) std::cout << __LINE__ << std::endl;

      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      if(isDebug) std::cout << __LINE__ << std::endl;

      lepTree_p->GetEntry(entry);
      jetTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      hltTree_p->GetEntry(entry);
      pfTree_p->GetEntry(entry);
      skimAnaTree_p->GetEntry(entry);
      cnt["tot"]++;
      passed_+=true;

      // if(!isMC && !trig) continue;
      if(!trig) passed_=false;
      if(passed_) cnt["filt1"]++;
      if(!phfCoincFilter) passed_=false;
      if(passed_) cnt["filt2"]++;
      if(!HBHENoiseFilterResult) passed_=false;
      if(passed_) cnt["filt3"]++;
      if(!pcollisionEventSelection) passed_=false;
      if(passed_) cnt["filt4"]++;
      if(!pprimaryVertexFilter) passed_=false;
      if(passed_) cnt["filt5"]++;

      if(TMath::Abs(vz_) > 15) passed_=false;
      if(passed_) cnt["filt6"]++;
      
      if (!isMC && !passed_) continue;

      if(isDebug) std::cout << __LINE__ << std::endl;

      if (isMC) {
         genMuPt_=-999.;
         genMuEta_=-999.;
         genMuPhi_=-999.;
         genElePt_=-999.;
         genEleEta_=-999.;
         genElePhi_=-999.;

         for (int i=0; i<nMC; i++) {
            if (abs(mcMomPID->at(i))==24) {
               if (abs(mcPID->at(i))==11 && mcPt->at(i)>genElePt_) {
                  genElePt_ = mcPt->at(i);
                  genEleEta_ = mcEta->at(i);
                  genElePhi_ = mcPhi->at(i);
               } else if (abs(mcPID->at(i))==13 && mcPt->at(i)>genMuPt_) {
                  genMuPt_ = mcPt->at(i);
                  genMuEta_ = mcEta->at(i);
                  genMuPhi_ = mcPhi->at(i);
               }
            } // if mcMom is W
         } // genParticles loop
         if (genMuPt_>muPtCut&&fabs(genMuEta_)<muEtaCut) cnt["genm"]++;
         if (genElePt_>elePtCut&&fabs(genEleEta_)<eleEtaCut) cnt["gene"]++;
         if (genMuPt_>muPtCut&&fabs(genMuEta_)<muEtaCut&&genElePt_>elePtCut&&fabs(genEleEta_)<eleEtaCut) cnt["genem"]++;
      } // if (isMC)

      int centEleId = getCentBinEleId((double)hiBin_/2.);

      float tempMuPt_[nLep];
      float tempMuPhi_[nLep];
      float tempMuEta_[nLep];
      int tempMuChg_[nLep];
      float tempMuIso_[nLep];     
      float tempMuIsoTk_[nLep];     
      float tempMuDRjet_[nLep];     
      float tempMuDRgen_[nLep];     
      float tempMuPtrel_[nLep];     
      float tempMuPtrel2_[nLep];     
      float tempMuJetPt_[nLep];     

      float tempElePt_[nLep];
      float tempElePhi_[nLep];
      float tempEleEta_[nLep];
      int tempEleChg_[nLep];
      float tempEleIso_[nLep];            
      float tempEleIsoTk_[nLep];            
      float tempEleDRjet_[nLep];            
      float tempEleDRgen_[nLep];            
      float tempElePtrel_[nLep];            
      float tempElePtrel2_[nLep];            
      float tempEleJetPt_[nLep];            

      for(int lepIter = 0; lepIter < nLep; lepIter++){
         lepPt_[lepIter] = -999;
         lepPhi_[lepIter] = -999;
         lepEta_[lepIter] = -999;
         lepChg_[lepIter] = -999;
         lepID_[lepIter] = -999;
         lepIso_[lepIter] = -999;
         lepIsoTk_[lepIter] = -999;
         lepDRjet_[lepIter] = -999;
         lepDRgen_[lepIter] = -999;
         lepPtrel_[lepIter] = -999;
         lepPtrel2_[lepIter] = -999;
         lepJetPt_[lepIter] = -999;
      }

      for(int lepIter = 0; lepIter < 2; lepIter++){
         tempMuPt_[lepIter] = -999;
         tempMuPhi_[lepIter] = -999;
         tempMuEta_[lepIter] = -999;
         tempMuChg_[lepIter] = -999;
         tempMuIso_[lepIter] = -999;
         tempMuIsoTk_[lepIter] = -999;
         tempMuDRjet_[lepIter] = -999;
         tempMuDRgen_[lepIter] = -999;
         tempMuPtrel_[lepIter] = -999;
         tempMuPtrel2_[lepIter] = -999;
         tempMuJetPt_[lepIter] = -999;

         tempElePt_[lepIter] = -999;
         tempElePhi_[lepIter] = -999;
         tempEleEta_[lepIter] = -999;
         tempEleChg_[lepIter] = -999;
         tempEleIso_[lepIter] = -999;
         tempEleDRjet_[lepIter] = -999;
         tempEleDRgen_[lepIter] = -999;
         tempElePtrel_[lepIter] = -999;
         tempElePtrel2_[lepIter] = -999;
         tempEleJetPt_[lepIter] = -999;
      }

      for(int ij = 0; ij<nMaxJets; ++ij) {
         jtPt_[ij] = -999.;
         jtEta_[ij] = -999.;
         jtPhi_[ij] = -999.;
         jtM_[ij] = -999.;
         discr_csvV1_[ij] = -999.;
         refparton_flavorForB_[ij] = -999.;
      }

      if(isDebug) std::cout << __LINE__ << std::endl;

      int njets = 0;
      TLorentzVector mht;
      for(int jetIter = 0; jetIter < nref; jetIter++){
         if(jtpt[jetIter]<jetPtCut) continue;
         if(fabs(jteta[jetIter])>jetEtaCut) continue;
         jtPt_[njets]  = jtpt[jetIter];
         jtEta_[njets] = jteta[jetIter];
         jtPhi_[njets] = jtphi[jetIter];
         jtM_[njets]   = jtm[jetIter]; 
         TLorentzVector tjet; tjet.SetPtEtaPhiM(jtpt[jetIter],jteta[jetIter],jtphi[jetIter],jtm[jetIter]);
         mht = mht+tjet;
         discr_csvV1_[njets] = discr_csvV1[jetIter];
         if (isMC) refparton_flavorForB_[njets] = refparton_flavorForB[jetIter];
         ++njets;
      }
      nJt_ = njets;

      // not interested in events with no jet passing cuts
      if (nJt_==0) passed_=false;
      if(passed_) cnt["recj"]++;
      if (!isMC && !passed_) continue;
      if (discr_csvV1_[0]>0.7) cnt["recbj"]++;
      if (isMC && refparton_flavorForB_[0]) cnt["recbq"]++;

      //Find two leading muons
      for(unsigned int muIter = 0; muIter < fForestMu.muPt->size(); ++muIter) {
         cnt["mu0"]++;
         bool matchm = (sqrt(pow(fForestMu.muEta->at(muIter)-genMuEta_,2)+pow(TVector2::Phi_mpi_pi(fForestMu.muPhi->at(muIter)-genMuPhi_),2))<0.1);
         if (matchm) cnt["matchm0"]++;
         if(TMath::Abs(fForestMu.muEta->at(muIter)) > muEtaCut) continue;
         cnt["mu1"]++;
         if (matchm) cnt["matchm1"]++;
         if(fForestMu.muPt->at(muIter) < muPtCut) continue;
         cnt["mu2"]++;
         if (matchm) cnt["matchm2"]++;

         if(fForestMu.muChi2NDF->at(muIter) >= muChi2NDFCut) continue;
         cnt["mu3"]++;
         if (matchm) cnt["matchm3"]++;
         if(TMath::Abs(fForestMu.muInnerD0->at(muIter)) > muInnerD0Cut) continue;
         cnt["mu4"]++;
         if (matchm) cnt["matchm4"]++;
         if(TMath::Abs(fForestMu.muInnerDz->at(muIter)) > muInnerDzCut) continue;
         cnt["mu5"]++;
         if (matchm) cnt["matchm5"]++;
         if(fForestMu.muMuonHits->at(muIter) <= muMuonHitsCut) continue;
         cnt["mu6"]++;
         if (matchm) cnt["matchm6"]++;
         if(fForestMu.muStations->at(muIter) <= muStationsCut) continue;
         cnt["mu7"]++;
         if (matchm) cnt["matchm7"]++;
         if(fForestMu.muTrkLayers->at(muIter) <= muTrkLayersCut) continue;
         cnt["mu8"]++;
         if (matchm) cnt["matchm8"]++;
         if(fForestMu.muPixelHits->at(muIter) <= muPixelHitsCut) continue;
         cnt["mu9"]++;
         if (matchm) cnt["matchm9"]++;
         float muDRjet, muPtrel, muPtrel2, muJetPt;
         muDRjet = calcLeptonMinDRjet(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),nJt_,jtPt_,jtEta_,jtPhi_,muPtrel,muPtrel2,muJetPt);
         // if (muDRjet<lepDRjetCut) continue;
         cnt["mu10"]++;
         if (matchm) cnt["matchm10"]++;
         // if (muPtrel<lepPtrelCut) continue;

         float muDRgen = isMC ? sqrt(pow(fForestMu.muEta->at(muIter)-genMuEta_,2)+pow(TVector2::Phi_mpi_pi(fForestMu.muPhi->at(muIter)-genMuPhi_),2)) : -999;

         if(fForestMu.muPt->at(muIter) > lepPt_[0]){
            tempMuPt_[1]  = tempMuPt_[0];
            tempMuPhi_[1] = tempMuPhi_[0];
            tempMuEta_[1] = tempMuEta_[0];
            tempMuChg_[1] = tempMuChg_[0];
            tempMuIso_[1] = tempMuIso_[0]; 
            tempMuIsoTk_[1] = tempMuIsoTk_[0]; 
            tempMuDRjet_[1] = tempMuDRjet_[0]; 
            tempMuDRgen_[1] = tempMuDRgen_[0]; 
            tempMuPtrel_[1] = tempMuPtrel_[0]; 
            tempMuPtrel2_[1] = tempMuPtrel2_[0]; 
            tempMuJetPt_[1] = tempMuJetPt_[0]; 

            tempMuPt_[0]  = fForestMu.muPt->at(muIter);
            tempMuPhi_[0] = fForestMu.muPhi->at(muIter);
            tempMuEta_[0] = fForestMu.muEta->at(muIter);
            tempMuChg_[0] = fForestMu.muCharge->at(muIter);
            tempMuIso_[0] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi,pfId); //iso;
            tempMuIsoTk_[0] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi,pfId,true); //tk iso;
            tempMuDRjet_[0] = muDRjet;
            tempMuDRgen_[0] = muDRgen;
            tempMuPtrel_[0] = muPtrel;
            tempMuPtrel2_[0] = muPtrel2;
            tempMuJetPt_[0] = muJetPt;
         }
         else if(fForestMu.muPt->at(muIter) > tempMuPt_[1]){
            tempMuPt_[1]  = fForestMu.muPt->at(muIter);
            tempMuPhi_[1] = fForestMu.muPhi->at(muIter);
            tempMuEta_[1] = fForestMu.muEta->at(muIter);
            tempMuChg_[1] = fForestMu.muCharge->at(muIter);
            tempMuIso_[1] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi,pfId); //iso;
            tempMuIsoTk_[1] = calcLeptonIsolation(fForestMu.muPt->at(muIter),fForestMu.muEta->at(muIter),fForestMu.muPhi->at(muIter),pfPt,pfEta,pfPhi,pfId,true);//tk iso;
            tempMuDRjet_[1] = muDRjet;
            tempMuDRgen_[1] = muDRgen;
            tempMuPtrel_[1] = muPtrel;
            tempMuPtrel_[1] = muPtrel;
            tempMuPtrel2_[1] = muPtrel2;
            tempMuJetPt_[1] = muJetPt;
         }
      }
      if (tempMuPt_[0]>0) cnt["recm"]++;

      //Find two leading electrons
      for(unsigned int eleIter = 0; eleIter < fForestEle.elePt->size(); ++eleIter) {
         cnt["ele0"]++;
         bool matche = (sqrt(pow(fForestEle.eleEta->at(eleIter)-genEleEta_,2)+pow(TVector2::Phi_mpi_pi(fForestEle.elePhi->at(eleIter)-genElePhi_),2))<0.1);
         if (matche) cnt["matche0"]++;
         if(TMath::Abs(fForestEle.eleEta->at(eleIter)) > eleEtaCut) continue;
         cnt["ele1"]++;
         if (matche) cnt["matche1"]++;
         if(fForestEle.elePt->at(eleIter) < elePtCut) continue;	
         cnt["ele2"]++;
         if (matche) cnt["matche2"]++;

         int eleEtaCutPos = 0;
         if(TMath::Abs(fForestEle.eleEta->at(eleIter)) > barrelEndcapEta) eleEtaCutPos = 1;
         cnt["ele3"]++;
         if (matche) cnt["matche3"]++;

         if(fForestEle.eleSigmaIEtaIEta->at(eleIter) > eleSigmaIEtaIEta_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele4"]++;
         if (matche) cnt["matche4"]++;
         if(TMath::Abs(fForestEle.eledEtaAtVtx->at(eleIter)) > eleDEtaIn_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele5"]++;
         if (matche) cnt["matche5"]++;
         if(TMath::Abs(fForestEle.eledPhiAtVtx->at(eleIter)) > eleDPhiIn_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele6"]++;
         if (matche) cnt["matche6"]++;
         if(fForestEle.eleHoverE->at(eleIter) > eleHOverE_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele7"]++;
         if (matche) cnt["matche7"]++;
         if(TMath::Abs(fForestEle.eleD0->at(eleIter)) > eleD0_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele8"]++;
         if (matche) cnt["matche8"]++;
         if(TMath::Abs(fForestEle.eleDz->at(eleIter)) > eleDZ_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele9"]++;
         if (matche) cnt["matche9"]++;
         if(TMath::Abs(fForestEle.eleEoverPInv->at(eleIter)) > eleEoverPInv_VetoCut[eleEtaCutPos][centEleId]) continue;
         cnt["ele10"]++;
         if (matche) cnt["matche10"]++;
         float eleDRjet, elePtrel, elePtrel2, eleJetPt;
         eleDRjet = calcLeptonMinDRjet(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),nJt_,jtPt_,jtEta_,jtPhi_,elePtrel,elePtrel2,eleJetPt);
         // if (eleDRjet<lepDRjetCut) continue;
         cnt["ele11"]++;
         if (matche) cnt["matche11"]++;
         // if (elePtrel<lepPtrelCut) continue;

         if(isDebug) std::cout << __LINE__ << std::endl;
         float ptEle = fForestEle.elePt->at(eleIter);
         //Printf("ptEle: %f",ptEle);
         float eleCorr = 1.;
         if(!isMC) eleCorr = energyRegression.ElectronRegressionTMVA(fForestEle.eleSCPhi->at(eleIter),
               fForestEle.eleSCEta->at(eleIter),
               fForestEle.eleSigmaIEtaIEta->at(eleIter),
               fForestEle.eleSigmaIPhiIPhi->at(eleIter),
               fForestEle.eleSCEtaWidth->at(eleIter),
               fForestEle.eleSCPhiWidth->at(eleIter),
               fForestEle.eleHoverE->at(eleIter));


         fForestEle.elePt->at(eleIter) = ptEle * eleCorr;
         //Printf("ptEle: %f  eleCorr: %f elePt: %f",ptEle,eleCorr,fForestEle.elePt->at(eleIter));

         float eleDRgen = isMC ? sqrt(pow(fForestEle.eleEta->at(eleIter)-genEleEta_,2)+pow(TVector2::Phi_mpi_pi(fForestEle.elePhi->at(eleIter)-genElePhi_),2)) : -999;

         if(fForestEle.elePt->at(eleIter) > tempElePt_[0]){
            tempElePt_[1] = tempElePt_[0];
            tempElePhi_[1] = tempElePhi_[0];
            tempEleEta_[1] = tempEleEta_[0];
            tempEleChg_[1] = tempEleChg_[0];
            tempEleIso_[1] = tempEleIso_[0];	 
            tempEleIsoTk_[1] = tempEleIsoTk_[0];	 
            tempEleDRjet_[1] = tempEleDRjet_[0];	 
            tempEleDRgen_[1] = tempEleDRgen_[0];	 
            tempElePtrel_[1] = tempElePtrel_[0];	 
            tempElePtrel2_[1] = tempElePtrel2_[0];	 
            tempEleJetPt_[1] = tempEleJetPt_[0];	 

            tempElePt_[0] = fForestEle.elePt->at(eleIter);
            tempElePhi_[0] = fForestEle.elePhi->at(eleIter);
            tempEleEta_[0] = fForestEle.eleEta->at(eleIter);
            tempEleChg_[0] = fForestEle.eleCharge->at(eleIter);
            tempEleIso_[0] = calcLeptonIsolation(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),pfPt,pfEta,pfPhi,pfId); //iso;
            tempEleIsoTk_[0] = calcLeptonIsolation(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),pfPt,pfEta,pfPhi,pfId,true); //tk iso;
            tempEleDRjet_[0] = eleDRjet;
            tempEleDRgen_[0] = eleDRgen;
            tempElePtrel_[0] = elePtrel;
            tempElePtrel2_[0] = elePtrel2;
            tempEleJetPt_[0] = eleJetPt;
         }
         else if(fForestEle.elePt->at(eleIter) > tempElePt_[1]){
            tempElePt_[1]  = fForestEle.elePt->at(eleIter);
            tempElePhi_[1] = fForestEle.elePhi->at(eleIter);
            tempEleEta_[1] = fForestEle.eleEta->at(eleIter);
            tempEleChg_[1] = fForestEle.eleCharge->at(eleIter);
            tempEleIso_[1] = calcLeptonIsolation(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),pfPt,pfEta,pfPhi,pfId);//iso;
            tempEleIsoTk_[1] = calcLeptonIsolation(fForestEle.elePt->at(eleIter),fForestEle.eleEta->at(eleIter),fForestEle.elePhi->at(eleIter),pfPt,pfEta,pfPhi,pfId,true);//tk iso;
            tempEleDRjet_[1] = eleDRjet;
            tempEleDRgen_[1] = eleDRgen;
            tempElePtrel_[1] = elePtrel;
            tempElePtrel2_[1] = elePtrel2;
            tempEleJetPt_[1] = eleJetPt;
         }
      }
      if (tempElePt_[0]>0) cnt["rece"]++;

      //store electrons and muons in out tree
      int lepIter = 0;
      for(int muIter = 0; muIter < 2; muIter++){
         if(tempMuPt_[muIter]<0.) continue;
         lepPt_[lepIter] = tempMuPt_[muIter];
         lepPhi_[lepIter] = tempMuPhi_[muIter];
         lepEta_[lepIter] = tempMuEta_[muIter];
         lepChg_[lepIter] = tempMuChg_[muIter];
         lepID_[lepIter] = muID;
         lepIso_[lepIter] = tempMuIso_[muIter];
         lepIsoTk_[lepIter] = tempMuIsoTk_[muIter];
         lepDRjet_[lepIter] = tempMuDRjet_[muIter];
         lepDRgen_[lepIter] = tempMuDRgen_[muIter];
         lepPtrel_[lepIter] = tempMuPtrel_[muIter];
         lepPtrel2_[lepIter] = tempMuPtrel2_[muIter];
         lepJetPt_[lepIter] = tempMuJetPt_[muIter];
         ++lepIter;
      }
      for(int eleIter = 0; eleIter < 2; eleIter++){
         if(tempElePt_[eleIter]<0.) continue;
         lepPt_[lepIter] = tempElePt_[eleIter];
         lepPhi_[lepIter] = tempElePhi_[eleIter];
         lepEta_[lepIter] = tempEleEta_[eleIter];
         lepChg_[lepIter] = tempEleChg_[eleIter];
         lepID_[lepIter] = eleID;
         lepIso_[lepIter] = tempEleIso_[eleIter];
         lepIsoTk_[lepIter] = tempEleIsoTk_[eleIter];
         lepDRjet_[lepIter] = tempEleDRjet_[eleIter];
         lepDRgen_[lepIter] = tempEleDRgen_[eleIter];
         lepPtrel_[lepIter] = tempElePtrel_[eleIter];
         lepPtrel2_[lepIter] = tempElePtrel2_[eleIter];
         lepJetPt_[lepIter] = tempEleJetPt_[eleIter];
         ++lepIter;
      }
      if(lepIter<2) passed_=false;
      if (!isMC && lepIter<1) continue;
      nLep_ = lepIter;

      // emu variables
      TLorentzVector tele, tmu, temu;
      int chg1=-999,chg2=-999;
      float drj1,drj2,pr1,pr2,pr21,pr22,jp1,jp2;
      if (tempMuPt_[0]>0 && tempElePt_[0]>0) {
         tele.SetPtEtaPhiM(tempElePt_[0],tempEleEta_[0],tempElePhi_[0],eleMass);
         tmu.SetPtEtaPhiM(tempMuPt_[0],tempMuEta_[0],tempMuPhi_[0],muMass);
         chg1=tempEleChg_[0];
         chg2=tempMuChg_[0];
         drj1=tempEleDRjet_[0];drj2=tempMuDRjet_[0];
         pr1=tempElePtrel_[0];pr2=tempMuPtrel_[0];
         pr21=tempElePtrel2_[0];pr22=tempMuPtrel2_[0];
         jp1=tempEleJetPt_[0];jp2=tempMuJetPt_[0];
         emuFlav_=0;
         cnt["recem"]++;
      } else if (tempElePt_[0]>0 && tempElePt_[1]>0) {
         tele.SetPtEtaPhiM(tempElePt_[0],tempEleEta_[0],tempElePhi_[0],eleMass);
         tmu.SetPtEtaPhiM(tempElePt_[1],tempEleEta_[1],tempElePhi_[1],eleMass);
         chg1=tempEleChg_[0];
         chg2=tempEleChg_[1];
         drj1=tempEleDRjet_[0];drj2=tempEleDRjet_[1];
         pr1=tempElePtrel_[0];pr2=tempElePtrel_[1];
         pr21=tempElePtrel2_[0];pr22=tempElePtrel2_[1];
         jp1=tempEleJetPt_[0];jp2=tempEleJetPt_[1];
         emuFlav_=1;
         cnt["recee"]++;
      } else if (tempMuPt_[0]>0 && tempMuPt_[1]>0) {
         tele.SetPtEtaPhiM(tempMuPt_[1],tempMuEta_[1],tempMuPhi_[1],muMass);
         tmu.SetPtEtaPhiM(tempMuPt_[0],tempMuEta_[0],tempMuPhi_[0],muMass);
         chg1=tempMuChg_[1];
         chg2=tempMuChg_[0];
         drj1=tempMuDRjet_[1];drj2=tempMuDRjet_[0];
         pr1=tempMuPtrel_[1];pr2=tempMuPtrel_[0];
         pr21=tempMuPtrel2_[1];pr22=tempMuPtrel2_[0];
         jp1=tempMuJetPt_[1];jp2=tempMuJetPt_[0];
         emuFlav_=2;
         cnt["recmm"]++;
      }
      if (chg1!=-999) {
         temu = tele+tmu;
         emuPt_ = temu.Pt();
         emuRap_ = temu.Rapidity();
         emuPhi_ = temu.Phi();
         emuMass_ = temu.M();
         emuSgn_ = chg1+chg2;
         emuPhistar_ = phistar(tele,tmu);
         emuDpt_ = fabs(tmu.Pt()-tele.Pt());
         emuDeta_ = fabs(tmu.Eta()-tele.Eta());
         emuDphi_ = fabs(tmu.DeltaPhi(tele));
         emuDR_ = tmu.DeltaR(tele);
         emuMinLepPt_=min(tmu.Pt(),tele.Pt());
         emuSumLepPt_=tmu.Pt()+tele.Pt();
         emuMinDRjet_=min(drj1,drj2);
         emuMinPtrel_=min(pr1,pr2);
         emuMinPtrel2_=min(pr21,pr22);
         emuMinJetPt_=min(jp1,jp2);
      } else {
         emuFlav_=-999;
         emuPt_=-999;
         emuRap_=-999;
         emuPhi_=-999;
         emuMass_=-999;
         emuSgn_=-999;
         emuPhistar_=-999;
         emuDpt_=-999;
         emuDeta_=-999;
         emuDphi_=-999;
         emuDR_=-999;
         emuMinLepPt_=-999;
         emuSumLepPt_=-999;
         emuMinDRjet_=-999;
         emuMinPtrel_=-999;
         emuMinPtrel2_=-999;
         emuMinJetPt_=-999;
      }
      MHT_ = mht.Pt();

      skimTree_p->Fill();

   }//entries loop

   //    inFile_p->Close();
   // delete inFile_p;    

   // print counters
   for (map<string,int>::const_iterator it=cnt.begin(); it!=cnt.end(); it++) {
      cout << it->first << ": " << it->second << endl;
   }

   outFile_p->cd();
   TNamed pathStr1("pathStr1", outFileName.c_str());
   pathStr1.Write("", TObject::kOverwrite);
   TNamed pathStr2("pathStr2", inFileName.c_str());
   pathStr2.Write("", TObject::kOverwrite);
   skimTree_p->Write("", TObject::kOverwrite);

   outFile_p->Close();
   delete outFile_p;

   return;
}

double calcLeptonIsolation(float lepPt, float lepEta, float lepPhi, std::vector<float> *pfPt, std::vector<float> *pfEta, std::vector<float> *pfPhi, std::vector<int> *pfId, bool tkonly) {
  //calculate lepton isolation from pf candidates.
  //Isolation cone R=0.3
  
  double conePt = 0.;
  for(unsigned int i = 0; i<pfPt->size(); ++i) {
     // if (!tkonly && pfId->at(i)!=1 && pfId->at(i)!=4 && pfId->at(i)!=5) continue; // keep only tracks, photons and neutral hadrons
     if (tkonly && pfId->at(i)!=1) continue; // keep only tracks for tkonly

    double deltaR = sqrt(pow(acos(cos(lepPhi-pfPhi->at(i))),2)+pow(lepEta-pfEta->at(i),2));

    if(deltaR<0.03 || deltaR>0.3) continue;
    // if(deltaR>0.3) continue;

    conePt+=pfPt->at(i);
  }
  double relIso = conePt;
  if(lepPt>0.) relIso = conePt/lepPt;
  
  return relIso;
}

double phistar(TLorentzVector v1, TLorentzVector v2) {
   double eta1 = v1.Eta();
   double eta2 = v2.Eta();
   double phi1 = v1.Phi();
   double phi2 = v2.Phi();
   return tan((TMath::Pi()-fabs(TVector2::Phi_mpi_pi(phi2-phi1)))/2.)*sin(acos(tanh(fabs(eta2-eta1)/2.)));
}

float calcLeptonMinDRjet(float lepPt, float lepEta, float lepPhi, int njt, float *jtpt, float *jteta, float *jtphi, float &lepPtrel, float &lepPtrel2, float &jetpt) {
   TVector3 tlep; tlep.SetPtEtaPhi(lepPt,lepEta,lepPhi);

   double mindr=999;
   for (int i=0; i<njt; i++) {
      TVector3 tjet; tjet.SetPtEtaPhi(jtpt[i],jteta[i],jtphi[i]);
      double dr = tlep.DeltaR(tjet);
      if (dr==0) {
         lepPtrel=0;
         lepPtrel2=0;
         jetpt=jtpt[i];
         return 0;
      } else if (dr<mindr) {
         mindr=dr;
         jetpt=jtpt[i];
         // from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/src/MuonMvaEstimator.cc#L39
         TVector3 tjet2 = tjet-tlep;
         float dot = tlep.Dot( tjet2 );
         lepPtrel = tlep.Perp2() - dot*dot/tjet2.Perp2();
         lepPtrel = lepPtrel>0 ? sqrt(lepPtrel) : 0.0;

         // my estimate from https://cds.cern.ch/record/1356198/files/ATLAS-CONF-2011-089.pdf:
         // the variable pTrel [...] is defined as the momentum of the muon transverse to the combined muon plus jet axis.
         TVector3 tlj = tjet+tlep;
         lepPtrel2 = fabs((tlj.Orthogonal()).Dot(tlep));
      }
   }
   return mindr;
}
