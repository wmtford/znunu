#!/bin/bash

# source batch/exportProd.sh

# these arrays must all have the same length
SAMPLES=(
Run2017B-31Mar2018-v1.HTMHT_0_RA2AnalysisTree.root \
Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root \
tree_SingleMuon_2017B.root \
tree_SingleElectron_2017B.root \
tree_DYJetsToLL_M-50_HT-100to200_MC2017.root \
)
DIRS=(
root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV17
root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV17
/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV17/tree_DYm_CleanVars
/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV17/tree_DYe_CleanVars
/mnt/hadoop/store/user/aperloff/SusyRA2Analysis2015/Skims/Run2ProductionV17/tree_DYm_CleanVars
)
TREES=(
TreeMaker2/PreSelection \
TreeMaker2/PreSelection \
tree \
tree \
tree \
)

HEADER="
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <vector>
#include <string>
"
