#!/bin/bash

# MC signal
for i in `seq 0 3`; do 
   root -l -b -q makeEMuSkim.C+'("emuSkim_ttbarMC_'${i}'.root","/eos/cms/store/group/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT2l_172v5_PowhegV2_hvq/Forest/v1/merge/HiForest_'${i}'.root",true)'; 
done

# MC bkg
for i in `seq 0 2`; do 
   root -l -b -q makeEMuSkim.C+'("emuSkim_ZmmMC_'${i}'.root","/eos/cms/store/group/cmst3/group/hintt/mverweij/CS/MC/PbPb/Pythia8_Zmu10mu10_Hydjet_MB/crab_HiForestZMu10Mu10V1/160625_200912/merge/HiForest_'${i}'.root",true)'; 
done

# data
for i in `seq 0 10`; do 
   root -l -b -q makeEMuSkim.C+'("emuSkim_data_'${i}'.root","/eos/cms/store/cmst3/group/hintt/mverweij/PbPb5TeV/data/HIEWQExo/crab_FilteredSingleMuHighPt_v2/160421_135925/mergePartial/HiForest_'${i}'.root",false)' 
done
