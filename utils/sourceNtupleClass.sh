#!/bin/bash

# source batch/exportProd.sh

# these arrays must all have the same length
SAMPLES=(
Run2017B-31Mar2018-v1.MET_0_RA2AnalysisTree.root \
Autumn18.GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_0_RA2AnalysisTree.root \
tree_EGamma_2018A.root \
tree_GJets_DR-0p4_HT-100to200_MC2018.root \
)
DIRS=(
root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV18
root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV18
root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV18/tree_GJet_CleanVars
root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Skims/Run2ProductionV18/tree_GJet_CleanVars
)
TREES=(
TreeMaker2/PreSelection \
TreeMaker2/PreSelection \
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
