//
//  Make histograms for Zinv background prediction in the RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#include "RA2bZinvAnalysis.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TCut.h>
#include <TTreeCache.h>

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::stringstream;

#include <stdio.h>
// using std::sprintf;

#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

ClassImp(RA2bZinvAnalysis)

// ======================================================================================

RA2bZinvAnalysis::RA2bZinvAnalysis() {
  Init();
}

RA2bZinvAnalysis::RA2bZinvAnalysis(const std::string& cfg_filename, const std::string& runBlock) {
  Init(cfg_filename);
  if (!runBlock.empty()) runBlock_ = runBlock;
  cout << "The runBlock is " << runBlock_ << endl;
}

void
RA2bZinvAnalysis::Init(const std::string& cfg_filename) {

  // Set up configuration, using boost/program_options.
  po::options_description desc("Config");
  desc.add_options()
    ("verbosity", po::value<int>(&verbosity_))
    ("era", po::value<std::string>(&era_))
    ("runBlock", po::value<std::string>(&runBlock_))  // May be overridden by constructor
    ("tree path", po::value<std::string>(&treeLoc_))
    ("root file index", po::value<std::string>(&fileListsFile_))
    ("delta phi sample", po::value<std::string>(&deltaPhi_), "nominal, hdp, ldp")
    ("integrated luminosity", po::value<double>(&intLumi_))
    ("apply Z mass cut", po::value<bool>(&applyMassCut_))
    ("apply Z/gamma Pt cut", po::value<bool>(&applyPtCut_))
    ("use analysis bin from tree", po::value<bool>(&useTreeCCbin_))
    ("use DeepCSV", po::value<bool>(&useDeepCSV_))
    ("apply b-tag SF", po::value<bool>(&applyBTagSF_))
    ("apply pileup weight", po::value<bool>(&applyPuWeight_))
    ("use custom pileup weight", po::value<bool>(&customPuWeight_))
    ("apply HEM jet veto", po::value<bool>(&applyHEMjetVeto_))
    ("apply Z Pt weight", po::value<bool>(&applyZptWt_))
    ("apply double-ratio fit weight", po::value<bool>(&applyDRfitWt_))
    ("apply scale factors to MC", po::value<bool>(&applySFwtToMC_))

    ;
  po::variables_map vm;
  std::ifstream cfi_file("RA2bZinvAnalysis.cfi");
  po::store(po::parse_config_file(cfi_file , desc), vm);
  po::notify(vm);
  if (!cfg_filename.empty()) {
    std::ifstream cfg_file(cfg_filename);
    vm = po::variables_map();  // Clear the map.
    po::store(po::parse_config_file(cfg_file , desc), vm);
    po::notify(vm);
  }

  // useTreeCCbin_  only in skims
  // applyBTagSF_  overridden false if !isMC_
  // applyPuWeight_  overridden false if !isMC_
  // customPuWeight_  Substitute Kevin P recipe for the PuWeight in the tree
  puWeight = 1;  // overridden from tree if isMC_
  Weight = 1;  // overridden from tree if isMC_
  TrueNumInteractions = 20;  // overridden from tree if isMC_
  NonPrefiringProb = 1;  // overridden from tree if isMC_
  RA2bin = 0;  // overridden from tree if isSkim_
  NElectrons = 0;  // overridden from tree if >=V15
  NMuons = 0;  // overridden from tree if >=V15
  ecalBadCalibFilter = 0;  // overridden from tree if >=V15

  if (!isMC_) {
    applyBTagSF_ = false;
    applyPuWeight_ = false;
  }
  if (!isSkim_) {
    useTreeCCbin_ = false;
    treeName_ = "TreeMaker2/PreSelection";  // For ntuple
  } else {
    treeName_ = "tree";  // For skims
  }
  if (applyBTagSF_
      || (ntupleVersion_ == "V15" && runBlock_.find("2016")==std::string::npos)
      || useDeepCSV_)
    useTreeCCbin_ = false;
  if (ntupleVersion_ == "V12") useDeepCSV_ = false;
  if (ntupleVersion_ == "V15" && runBlock_.find("2016")==std::string::npos) csvMthreshold_ = 0.8838;

  // Needed branches
  activeBranches_.push_back("NJets");
  activeBranches_.push_back("BTags");
  activeBranches_.push_back("HT");
  activeBranches_.push_back("HT5");
  activeBranches_.push_back("MHT");
  activeBranches_.push_back("JetID");
  activeBranches_.push_back("Jets");
  activeBranches_.push_back("Jets_hadronFlavor");
  activeBranches_.push_back("Jets_HTMask");
  activeBranches_.push_back("isoElectronTracks");
  activeBranches_.push_back("isoMuonTracks");
  activeBranches_.push_back("isoPionTracks");
  activeBranches_.push_back("DeltaPhi1");
  activeBranches_.push_back("DeltaPhi2");
  activeBranches_.push_back("DeltaPhi3");
  activeBranches_.push_back("DeltaPhi4");
  if (isSkim_) {
    activeBranches_.push_back("RA2bin");
  } else {
    activeBranches_.push_back("NJetsclean");
    activeBranches_.push_back("BTagsclean");
    activeBranches_.push_back("BTagsDeepCSVclean");
    activeBranches_.push_back("HTclean");
    activeBranches_.push_back("HT5clean");
    activeBranches_.push_back("MHTclean");
    activeBranches_.push_back("JetIDclean");
    activeBranches_.push_back("Jetsclean");
    activeBranches_.push_back("Jetsclean_hadronFlavor");
    activeBranches_.push_back("Jetsclean_HTMask");
    activeBranches_.push_back("isoElectronTracksclean");
    activeBranches_.push_back("isoMuonTracksclean");
    activeBranches_.push_back("isoPionTracksclean");
    activeBranches_.push_back("DeltaPhi1clean");
    activeBranches_.push_back("DeltaPhi2clean");
    activeBranches_.push_back("DeltaPhi3clean");
    activeBranches_.push_back("DeltaPhi4clean");
  }
  if (ntupleVersion_ != "V12") {
    activeBranches_.push_back("NMuons");
    activeBranches_.push_back("NElectrons");
    activeBranches_.push_back("BTagsDeepCSV");
    activeBranches_.push_back("ecalBadCalibFilter");
  }
  activeBranches_.push_back("RunNum");
  activeBranches_.push_back("LumiBlockNum");
  activeBranches_.push_back("EvtNum");
  activeBranches_.push_back("Jets_bDiscriminatorCSV");
  activeBranches_.push_back("Muons");
  activeBranches_.push_back("Electrons");
  activeBranches_.push_back("ZCandidates");
  activeBranches_.push_back("Photons");
  activeBranches_.push_back("Photons_nonPrompt");
  activeBranches_.push_back("Photons_fullID");
  activeBranches_.push_back("Photons_hasPixelSeed");
  activeBranches_.push_back("Photons_isEB");
  activeBranches_.push_back("NVtx");
  activeBranches_.push_back("TriggerNames");
  activeBranches_.push_back("TriggerPass");
  activeBranches_.push_back("TriggerPrescales");
  activeBranches_.push_back("HBHENoiseFilter");
  activeBranches_.push_back("HBHEIsoNoiseFilter");
  activeBranches_.push_back("eeBadScFilter");
  activeBranches_.push_back("EcalDeadCellTriggerPrimitiveFilter");
  activeBranches_.push_back("globalTightHalo2016Filter");
  activeBranches_.push_back("globalSuperTightHalo2016Filter");
  activeBranches_.push_back("BadChargedCandidateFilter");
  activeBranches_.push_back("BadPFMuonFilter");
  activeBranches_.push_back("PFCaloMETRatio");
  activeBranches_.push_back("METRatioFilter");
  activeBranches_.push_back("EcalNoiseJetFilter");
  activeBranches_.push_back("HTRatioDPhiFilter");
  activeBranches_.push_back("HTRatioFilter");

  activeBranches_.push_back("nAllVertices");
  if (isMC_) {
    activeBranches_.push_back("puWeight");
    activeBranches_.push_back("Weight");
    activeBranches_.push_back("TrueNumInteractions");
    activeBranches_.push_back("madMinPhotonDeltaR");
    activeBranches_.push_back("GenMuons");
    activeBranches_.push_back("GenElectrons");
    activeBranches_.push_back("GenHT");
    activeBranches_.push_back("GenTaus");
    activeBranches_.push_back("GenParticles");
    activeBranches_.push_back("GenParticles_PdgId");
    activeBranches_.push_back("GenParticles_Status");
    activeBranches_.push_back("NonPrefiringProb");
    activeBranches_.push_back("NonPrefiringProbUp");
    activeBranches_.push_back("NonPrefiringProbDn");
  }

  // For now we have only 2016 pileup and BTag correction files
  if (runBlock_.empty() || runBlock_.find("2016")!=std::string::npos || runBlock_.find("2017")!=std::string::npos || runBlock_.find("2018")!=std::string::npos) {
    if (applyPuWeight_ && customPuWeight_) {
      TFile* pufile = TFile::Open("../../Analysis/corrections/PileupHistograms_0121_69p2mb_pm4p6.root", "READ");
      puHist_ = (TH1*) pufile->Get("pu_weights_down");
    }
    BTagSFfile_ = "../../Analysis/btag/CSVv2_Moriond17_B_H_mod.csv";

  } // runBlock 2016

  CCbins_ = new CCbinning(era_, deltaPhi_);

  cout << "After initialization," << endl;
  cout << "The verbosity level is " << verbosity_ << endl;
  cout << "The ntuple version is " << ntupleVersion_ << endl;
  cout << "The MC flag is " << isMC_ << endl;
  cout << "The input-files-are-skims flag is " << isSkim_ << endl;
  cout << "The era is " << era_ << endl;
  cout << "The integrated luminosity = " << intLumi_ << endl;
  cout << "The path to input files is " << treeLoc_ << endl;
  cout << "The minDeltaPhi cuts are " << deltaPhi_ << endl;
  cout << "Apply Z/gamma Pt cut is " << applyPtCut_ << endl;
  cout << "Use analysis bin from tree is " << useTreeCCbin_ << endl;
  cout << "Use DeepCSV is " << useDeepCSV_ << endl;
  cout << "Apply b-tag scale factors is " << applyBTagSF_ << endl;
  cout << "Apply pileup weight is " << applyPuWeight_ << endl;
  cout << "Use custom if pileup weight is " << customPuWeight_ << endl;
  cout << "Apply HEM jet veto for 2018HEM is " << applyHEMjetVeto_ << endl;
  cout << "Apply Z Pt weight for 2017 Z MC is " << applyZptWt_ << endl;
  cout << "Apply double-ratio fit weight is " << applyDRfitWt_ << endl;
  cout << "Apply scale factors to MC for non-DR histograms is " << applySFwtToMC_ << endl;

  fillCutMaps();
  effPurCorr_.openFiles();

}  // ======================================================================================

TChain*
RA2bZinvAnalysis::getChain(const char* sample, Int_t* fCurrent, bool makeClass) {

  bool activateAllBranches = false;  // Can be set true for debugging
  cout << "activateAllBranches = " << activateAllBranches << endl;

  TString theSample(sample);
  TString key;
  if      (theSample.Contains("zinv")) key = TString("zinv");
  else if (theSample.Contains("ttzvv")) key = TString("ttzvv");
  else if (theSample.Contains("dymm")) key = TString("dymm");
  else if (theSample.Contains("dyee")) key = TString("dyee");
  else if (theSample.Contains("ttzmm")) key = TString("ttzmm");
  else if (theSample.Contains("ttzee")) key = TString("ttzee");
  else if (theSample.Contains("VVmm")) key = TString("VVmm");
  else if (theSample.Contains("VVee")) key = TString("VVee");
  else if (theSample.Contains("ttmm")) key = TString("ttmm");
  else if (theSample.Contains("ttee")) key = TString("ttee");
  else if (theSample.Contains("zmm")) key = TString("zmm");
  else if (theSample.Contains("zee")) key = TString("zee");
  else if (theSample.Contains("photon")) key = TString("photon");
  else if (theSample.Contains("gjetsqcd")) key = TString("gjetsqcd");
  else if (theSample.Contains("gjets")) key = TString("gjets");
  if (deltaPhi_ == "ldp" && isSkim_) key += "ldp";
  if (!runBlock_.empty()) key += runBlock_;  key("HEM") = "";

  TChain* chain = new TChain(treeName_.data());
  std::vector<TString> files = fileList(key);
  for (auto file : files) {
    if (verbosity_ >= 2) cout << file << endl;
    chain->Add(file);
  }
  if (fCurrent != nullptr) *fCurrent = -1;
  if (makeClass) return chain;

  cout << "Initial size of cache for chain = " << chain->GetCacheSize() << endl;
  TTreeCache::SetLearnEntries(1);
  chain->SetCacheSize(200*1024*1024);
  chain->SetCacheEntryRange(0, chain->GetEntries());
  if (activateAllBranches) {
    chain->AddBranchToCache("*", true);
  } else {
    for (auto theBranch : activeBranches_) chain->AddBranchToCache(theBranch, true);
  }
  chain->StopCacheLearningPhase();
  cout << "Reset size of cache for chain = " << chain->GetCacheSize() << endl;

  setBranchAddress(chain);

  if (!activateAllBranches) {
    chain->SetBranchStatus("*", 0);  // disable all branches
    for (auto theBranch : activeBranches_) chain->SetBranchStatus(theBranch, 1);
  }

  return chain;

}  // ======================================================================================

std::vector<TString>
RA2bZinvAnalysis::fileList(TString sampleKey) {

  std::vector<TString> files;
  TString dir(treeLoc_);
  TString key, toReserve;
  const char* fileName = fileListsFile_.data();
  ifstream dataStream;
  cout << "Getting root file list from " << fileName << endl;
  dataStream.open(fileName); // open the data file
  if (!dataStream.good()) {
    cout << "Open failed for file " << fileName << endl;
    return files; // exit if file not found
  }

  TString buf;
  Ssiz_t from;
  do {  // Look for matching sample key
    buf.ReadLine(dataStream);
    if (dataStream.eof()) {
      cout << "EOF found prematurely while searching for key " << sampleKey << endl;
      return files;
    }
    Ssiz_t pos = buf.First('#');
    if (pos <= 0) continue;  // Ignore comment lines
    buf = buf.Remove(pos, buf.Length()-pos);  // Remove comment following the tokens
    from = 0;
    buf.Tokenize(key, from);  // First token is the key
  }
  while (key != sampleKey);

  int nReserve = 0;  // Second token is the optional size to reserve in the file list vector
  bool next = buf.Tokenize(toReserve, from);
  if (next && toReserve.IsDec()) sscanf(toReserve.Data(), "%d", &nReserve);

  do {  // Read the file names for this sample
    buf.ReadLine(dataStream);
    if (dataStream.eof()) {
      cout << "EOF found prematurely while searching for end of file list" << sampleKey << endl;
      return files;
    }
    if (buf.Contains("#*")) {  // Ignore #* ... *# comment block
      do {
	buf.ReadLine(dataStream);
	if (dataStream.eof()) {
	  cout << "EOF found prematurely while searching for comment block terminator" << sampleKey << endl;
	  return files;
	}
      }
      while (!buf.Contains("*#"));
    }
    if (buf.Contains("#")) continue;
    buf = buf.Strip();  // Remove any trailing whitespace
    files.push_back(dir+buf);
  }
  while (!buf.Contains("##"));  // Dataset terminator is ##

  return files;

}  // ======================================================================================

TCut
RA2bZinvAnalysis::getCuts(const TString sample) {

  TCut cuts;
  TString sampleKey;
  try {sampleKey = sampleKeyMap_.at(sample);}
  catch (const std::out_of_range& oor) {
    std::cerr << oor.what() << " getCuts called with invalid sample = " << sample << endl;
    return cuts;
  }

  // // commonCuts_ = "(JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0 && BadPFMuonFilter && PFCaloMETRatio < 5)";  // Troy revision+
  // if (trigger.empty()) {
  //   commonCuts_ = "JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";  // Troy revision-;  moved JetID to !isSkim_
  // } else {
  //   commonCuts_ = "JetID==1&& globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";  // Troy revision-;  moved JetID to !isSkim_
  // }

  if (ntupleVersion_ == "V12" || ntupleVersion_ == "V15")
    commonCuts_ =    "globalTightHalo2016Filter==1";
  else
    commonCuts_ =    "globalSuperTightHalo2016Filter==1";
  commonCuts_ += " && HBHENoiseFilter==1";
  commonCuts_ += " && HBHEIsoNoiseFilter==1";
  commonCuts_ += " && EcalDeadCellTriggerPrimitiveFilter==1";
  commonCuts_ += " && BadChargedCandidateFilter";
  commonCuts_ += " && BadPFMuonFilter";
  commonCuts_ += " && NVtx > 0";
  if (ntupleVersion_ != "V12")
    // commonCuts_ += " && ecalBadCalibFilter==1";  // Added for 94X
  if (!isMC_)
    commonCuts_ += " && eeBadScFilter==1";
  // Kevin Pedro, re V16:  Please note that I have included only the JetID "event cleaning" cut in these skims.
  // The PFCaloMETRatio, HT5/HT, and noMuonJet cuts are left out,
  // so the impact of these cuts can be tested and refined if necessary.
  if (isSkim_ && ntupleVersion_ == "V16") {
    commonCuts_ += " && METRatioFilter";
    if (!isMC_ && runBlock_.find("2017") != std::string::npos)
      commonCuts_ += " && EcalNoiseJetFilter";
    if (era_ == "2016")
      commonCuts_ += " && HTRatioFilter";
    else
      commonCuts_ += " && HTRatioDPhiFilter";
      // Commoncuts_ += " && HT5/HT <= (DeltaPhi1 - (-0.5875))/1.025";
      // above implements "DeltaPhi1 >= 1.025*HT5/HT - 0.5875", avoiding "+" and "*" for regexp;
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming (xV16), single lepton
  }
  if (!isSkim_) {
    commonCuts_ += " && JetIDclean";
    commonCuts_ += " && PFCaloMETRatio < 5";
    if (era_ == "2016")
      commonCuts_ += " && HT5clean/HTclean <= 2";
    else
      commonCuts_ += " && HT5clean/HTclean <= (DeltaPhi1clean - (-0.5875))/1.025";
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming (xV16), single lepton
    // commonCuts_ += " && noFakeJet";  // Defined in loop, applied in skimming FastSim
  }

  massCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyMassCut_)
    massCut_ = "@ZCandidates.size()==1 && ZCandidates.M()>=76.188 && ZCandidates.M()<=106.188";

  ptCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyPtCut_)
    ptCut_ = "@ZCandidates.size()==1 && ZCandidates.Pt()>=200.";

  photonDeltaRcut_ = "1";  // For cut histograms only; incorporated into objCutMap
  if (sampleKey == "photon") {
    if (applyPtCut_) ptCut_ = "@Photons.size()==1 && Photons.Pt()>=200.";
    if (isMC_) photonDeltaRcut_ = "madMinPhotonDeltaR>=0.4";
  }

  if (isSkim_) {
    HTcut_ = std::string("HT>=") + std::to_string(CCbins_->htThreshold(0, 0));
    NJetscut_ = std::string("NJets>=") + std::to_string(CCbins_->nJetThreshold(0));
  } else {
    HTcut_ = std::string("HTclean>=") + std::to_string(CCbins_->htThreshold(0, 0));
    NJetscut_ = std::string("NJetsclean>=") + std::to_string(CCbins_->nJetThreshold(0));
  }
  MHTcut_ = MHTCutMap_.at(deltaPhi_);
  objcut_ = objCutMap_.at(sampleKey);
  minDphicut_ = minDphiCutMap_.at(deltaPhi_);
  // HEMgeomcut_ = "(-3.0 < Electrons[0].Eta() && Electrons[0].Eta() < -1.4 && -1.57 < Electrons[0].Phi() && Electrons[0].Phi() < -0.87)
  //             || (-3.0 < Electrons[1].Eta() && Electrons[1].Eta() < -1.4 && -1.57 < Electrons[1].Phi() && Electrons[1].Phi() < -0.87)";
  // Extended HEM cut η[-3.2, -1.2] φ[-1.77, -0.67] 


  // For test of same-flavor isoLeptonTracks veto
  isoSFlepTksVeto_ = "1";
  if (isSkim_) {
    if (sampleKey == "zmm") {isoSFlepTksVeto_ = "isoMuonTracks!=0";  isoSFlepTksCut_ = "isoMuonTracks==0";}
    if (sampleKey == "zee") {isoSFlepTksVeto_ = "isoElectronTracks!=0";  isoSFlepTksCut_ = "isoElectronTracks==0";}
  } else {
    if (sampleKey == "zmm") {isoSFlepTksVeto_ = "isoMuonTracksclean!=0";  isoSFlepTksCut_ = "isoMuonTracksclean==0";}
    if (sampleKey == "zee") {isoSFlepTksVeto_ = "isoElectronTracksclean!=0";  isoSFlepTksCut_ = "isoElectronTracksclean==0";}
  }
  photonVeto_ = "@Photons.size()==0";
  photonCut_ = "@Photons.size()!=0";

  cuts += objcut_;
  cuts += HTcut_;
  cuts += NJetscut_;
  cuts += MHTcut_;
  if (!isSkim_) cuts += minDphicut_;
  cuts += ptCut_;
  cuts += massCut_;
  cuts += commonCuts_;

  return cuts;

}  // ======================================================================================

void
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample, std::vector<histConfig*>& histograms, TCut baselineCuts) {
  //
  // Define N - 1 (or N - multiple) cuts, book histograms.  Traverse the chain and fill.
  //
  TString sampleKey = sampleKeyMap_.count(sample) > 0 ? sampleKeyMap_.at(sample) : TString("none");
  int currentYear = -1;
  Int_t fCurrent;  // current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);
  TObjArray* forNotify = new TObjArray;
  forNotify->SetOwner();  // so that TreeFormulas will be deleted

  TTreeFormula* baselineFormula = new TTreeFormula("baselineCuts", baselineCuts, chain);
  forNotify->Add(baselineFormula);

  // For B-tagging corrections
  if (applyBTagSF_) {
    btagcorr_ = new BTagCorrector;
    btagcorr_->SetCalib(BTagSFfile_);
  } else {
    btagcorr_ = nullptr;
  }

  // Z Pt weights
  Double_t ptBins[297] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.,64.,65.,66.,67.,68.,69.,70.,71.,72.,73.,74.,75.,76.,77.,78.,79.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,200.,202.,204.,206.,208.,210.,212.,214.,216.,218.,220.,222.,224.,226.,228.,230.,232.,234.,236.,238.,240.,242.,244.,246.,248.,250.,252.,254.,256.,258.,260.,262.,264.,266.,268.,270.,272.,274.,276.,278.,280.,282.,284.,286.,288.,290.,292.,294.,296.,298.,300.,304.,308.,312.,316.,320.,324.,328.,332.,336.,340.,344.,348.,352.,356.,360.,364.,368.,372.,376.,380.,384.,388.,392.,396.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,520.,540.,560.,580.,600.,650.,700.,750.,800.,900.,1000.};
  Double_t ptWgts[297] = {0.618, 0.629, 0.653, 0.693, 0.742, 0.801, 0.856, 0.904, 0.941, 0.975, 1.000, 1.021, 1.040, 1.055, 1.069, 1.085, 1.099, 1.109, 1.125, 1.134, 1.147, 1.155, 1.164, 1.167, 1.172, 1.174, 1.175, 1.173, 1.175, 1.173, 1.172, 1.171, 1.168, 1.164, 1.163, 1.157, 1.156, 1.147, 1.143, 1.143, 1.138, 1.137, 1.130, 1.122, 1.119, 1.116, 1.122, 1.112, 1.104, 1.103, 1.104, 1.096, 1.095, 1.092, 1.084, 1.082, 1.075, 1.074, 1.069, 1.074, 1.069, 1.061, 1.060, 1.062, 1.057, 1.053, 1.055, 1.055, 1.042, 1.043, 1.039, 1.038, 1.034, 1.036, 1.034, 1.040, 1.017, 1.022, 1.029, 1.019, 1.020, 1.004, 1.014, 1.021, 0.994, 1.008, 1.000, 0.997, 1.000, 0.992, 0.992, 0.994, 0.979, 0.976, 0.983, 0.977, 0.980, 0.970, 0.991, 0.980, 0.971, 0.974, 0.975, 0.975, 0.967, 0.973, 0.950, 0.958, 0.952, 0.968, 0.936, 0.964, 0.942, 0.950, 0.958, 0.951, 0.945, 0.939, 0.954, 0.947, 0.951, 0.930, 0.919, 0.920, 0.929, 0.913, 0.926, 0.945, 0.944, 0.928, 0.937, 0.908, 0.938, 0.926, 0.912, 0.919, 0.906, 0.919, 0.909, 0.899, 0.912, 0.917, 0.929, 0.887, 0.887, 0.884, 0.901, 0.905, 0.901, 0.915, 0.881, 0.877, 0.912, 0.872, 0.861, 0.869, 0.875, 0.867, 0.903, 0.884, 0.880, 0.854, 0.864, 0.881, 0.848, 0.898, 0.856, 0.887, 0.893, 0.853, 0.880, 0.844, 0.885, 0.822, 0.844, 0.834, 0.864, 0.869, 0.852, 0.819, 0.866, 0.857, 0.893, 0.812, 0.797, 0.877, 0.814, 0.827, 0.900, 0.853, 0.843, 0.820, 0.850, 0.859, 0.887, 0.846, 0.812, 0.855, 0.820, 0.785, 0.812, 0.845, 0.830, 0.829, 0.808, 0.791, 0.878, 0.798, 0.770, 0.794, 0.815, 0.833, 0.826, 0.795, 0.818, 0.782, 0.818, 0.779, 0.797, 0.785, 0.828, 0.755, 0.791, 0.740, 0.833, 0.804, 0.816, 0.842, 0.720, 0.789, 0.735, 0.805, 0.735, 0.748, 0.760, 0.842, 0.782, 0.756, 0.708, 0.760, 0.802, 0.772, 0.803, 0.777, 0.646, 0.710, 0.771, 0.736, 0.711, 0.756, 0.701, 0.743, 0.772, 0.721, 0.726, 0.703, 0.730, 0.755, 0.723, 0.692, 0.696, 0.717, 0.661, 0.736, 0.756, 0.701, 0.758, 0.739, 0.717, 0.643, 0.588, 0.679, 0.662, 0.710, 0.759, 0.702, 0.791, 0.726, 0.591, 0.681, 0.752, 0.651, 0.670, 0.667, 0.639, 0.595, 0.581, 0.634, 0.592, 0.659, 0.628, 0.536, 0.484, 0.568, 0.530, 0.676, 0.536};

  cutHistos cutHistFiller(chain, forNotify);  // for cutFlow histograms

  // Book histograms
  if (verbosity_ >= 1) cout << endl << "baseline = " << endl << baselineCuts << endl << endl;
  for (auto & hg : histograms) {
    if (hg->NbinsY > 0) {
      hg->hist = new TH2F(hg->name, hg->title, hg->NbinsX, hg->rangeX.first, hg->rangeX.second,
			  hg->NbinsY, hg->rangeY.first, hg->rangeY.second);
      hg->hist->SetOption("colz");
    } else {
      if (hg->binsX == nullptr)
	hg->hist = new TH1D(hg->name, hg->title, hg->NbinsX, hg->rangeX.first, hg->rangeX.second);
      else
	hg->hist = new TH1D(hg->name, hg->title, hg->NbinsX, hg->binsX);
      hg->hist->SetOption("hist");
      hg->hist->SetMarkerSize(0);
    }
    hg->hist->Sumw2();
    hg->hist->GetXaxis()->SetTitle(hg->axisTitles.first);
    hg->hist->GetYaxis()->SetTitle(hg->axisTitles.second);
    if (hg->name.Contains(TString("Cut"))) cutHistFiller.setAxisLabels((TH1D*) hg->hist);
    if (hg->name.Contains(TString("hCut")) || hg->name.Contains(TString("hgen"))) {
      hg->NminusOneCuts = "1";
    } else {
      hg->NminusOneCuts = baselineCuts;
      for (auto cutToOmit : hg->omitCuts) hg->NminusOneCuts(*cutToOmit) = "1";
      if (strlen(hg->addCuts) != 0) hg->NminusOneCuts += TString(" && ") + hg->addCuts;
    }
    if (verbosity_ >= 1) {
      cout << "\n For sample " << sample << ", histo " << hg->name << ", hg->omitCuts = ";
      for (auto cutToOmit : hg->omitCuts) cout << *cutToOmit << " ";
      cout << "; hg->addCuts = " << hg->addCuts;
      cout << ", cuts = " << endl << hg->NminusOneCuts << endl;
    }
    hg->NminusOneFormula = new TTreeFormula(hg->name, hg->NminusOneCuts, chain);
    forNotify->Add(hg->NminusOneFormula);
  }
  chain->SetNotify(forNotify);

  // Traverse the tree and fill histograms
  double MCwtCorr = 1.;
  Long64_t Nentries = chain->GetEntries();
  if (verbosity_ >= 1) cout << "Nentries in tree = " << Nentries << endl;
  int count = 0, countInFile = 0, countInSel = 0, countNegWt = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (verbosity_ >= 1 && count % 100000 == 0) cout << "Entry number " << count << endl;

    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      // New input root file encountered
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) {
    	if (btagcorr_) btagcorr_->SetEffs(thisFile);
	TString path = thisFile->GetName();
	if (verbosity_ >= 1) cout << "Current file in chain: " << path << endl;
	// Set MCwtCorr for this file
	if (path.Contains("V16") && path.Contains("MC2017") &&
	    (path.Contains("DYJetsToLL") || path.Contains("ZJetsToNuNu"))) {
	  if (applyZptWt_) {
	    if      (path.Contains("HT-100to200")) MCwtCorr = 1.05713;
	    else if (path.Contains("HT-200to400")) MCwtCorr = 1.20695;
	    else if (path.Contains("HT-400to600")) MCwtCorr = 1.30533;
	    else if (path.Contains("HT-600to800")) MCwtCorr = 1.38453;
	    else if (path.Contains("HT-800to1200")) MCwtCorr = 1.40301;
	    else if (path.Contains("HT-1200to2500")) MCwtCorr = 1.42145;
	    else if (path.Contains("HT-2500toInf")) MCwtCorr = 1.11697;
	  } else {
	    if      (path.Contains("HT-100to200")) MCwtCorr = 1.09226;
	    else if (path.Contains("HT-200to400")) MCwtCorr = 1.18517;
	    else if (path.Contains("HT-400to600")) MCwtCorr = 1.22966;
	    else if (path.Contains("HT-600to800")) MCwtCorr = 1.27798;
	    else if (path.Contains("HT-800to1200")) MCwtCorr = 1.27728;
	    else if (path.Contains("HT-1200to2500")) MCwtCorr = 1.27279;
	    else if (path.Contains("HT-2500toInf")) MCwtCorr = 0.975599;
	  }
	}
	// Following is the original HT reweighting to match 2016 MC
	// if (applyZptWt_ && path.Contains("V16") && path.Contains("MC2017") &&
	//     (path.Contains("DYJetsToLL") || path.Contains("ZJetsToNuNu"))) {
	//   if      (path.Contains("HT-100to200")) MCwtCorr = 0.968;
	//   else if (path.Contains("HT-200to400")) MCwtCorr = 1.018;
	//   else if (path.Contains("HT-400to600")) MCwtCorr = 1.062;
	//   else if (path.Contains("HT-600to800")) MCwtCorr = 1.083;
	//   else if (path.Contains("HT-800to1200")) MCwtCorr = 1.098;
	//   else if (path.Contains("HT-1200to2500")) MCwtCorr = 1.117;
	//   else if (path.Contains("HT-2500toInf")) MCwtCorr = 1.145;
	// }
      }
      chain->GetEntry(entry);  // Pull in tree variables for reinitialization
      setTriggerIndexList(sample);
      if (isMC_ && verbosity_ >= 1) cout << "MC weight for this file is " << Weight << " times correction " << MCwtCorr << endl;
      // Load the year-dependent correction factors.
      int theYear = -1;
      if (isMC_) {
	if (     runBlock_.find("2016") != std::string::npos) theYear = Year2016;
	else if (runBlock_.find("2017") != std::string::npos) theYear = Year2017;
	else if (runBlock_.find("2018") != std::string::npos) theYear = Year2018;
      } else {
	if (RunNum < Start2017) theYear = Year2016;
	else if (RunNum < Start2018) theYear = Year2017;
	else theYear = Year2018;
      }
      if (theYear != currentYear) {
	currentYear = theYear;
	cout << "currentYear = " << currentYear << endl;
	effPurCorr_.getHistos(sample, currentYear);  // For purity, Fdir, trigger eff, reco eff
      }
    }
    chain->GetEntry(entry);
    countInFile++;

    cleanVars();  // If unskimmed input, copy <var>clean to <var>
    Int_t BTagsOrig = setBTags();
    // if (countInFile <= 100) cout << "BTagsOrig, BTags, BTagsDeepCSV = " << BTagsOrig << ", " << BTags << ", " << BTagsDeepCSV << endl;

    if (ZCandidates->size() > 1 && verbosity_ >= 2) cout << ZCandidates->size() << " Z candidates found" << endl;
    // baselineFormula->GetNdata();
    // double baselineWt = baselineFormula->EvalInstance(0);

    // Compute event weight factors
    Double_t eventWt = 1, MCwt = 1, PUweight = 1, NoPrefireWt = 1, ZPtWt = 1;
    if (isMC_) {
      MCwt = 1000*intLumi_*Weight*MCwtCorr;  // MCwtCorr for 2017 MC
      if (applyPuWeight_) {
	// Pileup weight for 2016
	if (customPuWeight_ && puHist_ != nullptr) {
	  // This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
	  PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
	} else
	  PUweight = puWeight;  // Take puWeight directly from the tree
	MCwt *= PUweight;
      }

      if (currentYear == Year2016 || currentYear == Year2017) {
	// Apply L1 prefire weight for 2016 and 2017
	if (sampleKey.Contains("ee")) {
	  for (unsigned j = 0; j < Jets->size(); ++j) {
	    NoPrefireWt *= prefiring_weight_jet(j);
	  }
	  double eeNoPFwt = 1;
	  for (unsigned e = 0; e < Electrons->size(); ++e) {
	    double w = prefiring_weight_electron(e);
	    if (w < eeNoPFwt) eeNoPFwt = w;
	  }
	  NoPrefireWt *= eeNoPFwt;
	} else {
	NoPrefireWt = NonPrefiringProb;
	}
	MCwt *= NoPrefireWt;
      }  // 2016 or 2017

      if (applyZptWt_ && (TString(sample).Contains("dy") || TString(sample).Contains("zinv"))
	  && (currentYear == Year2017 || currentYear == Year2018)) {
	// Apply Z Pt weight for 2017 MC
	double ptZ = getPtZ();
	ZPtWt = 1.0;
	if (ptZ > 0.0) {
	  int iptbin;
	  for (iptbin = 1; iptbin<297; iptbin++) {
	    if (ptZ < ptBins[iptbin]) break;
	  }
	  ZPtWt *= ptWgts[iptbin-1];
	}
	MCwt *= ZPtWt;
      }  // DY or Zinv in 2017 or 2018

      eventWt *= MCwt;
    }  // isMC_

    pair<double, double> efficiency = effPurCorr_.weight(CCbins_, NJets, BTags, MHT, HT, *ZCandidates, *Photons,
							 *Electrons, *Muons, *Photons_isEB, applyDRfitWt_, currentYear);
    effWt_ = efficiency.first;  effSys_ = efficiency.second;

    // Trigger requirements
    bool passTrg = true;
    if (!isMC_) {
      passTrg = false;
      for (auto trgIndex : triggerIndexList_)
	if (TriggerPass->at(trgIndex)) passTrg = true;
    }
    // HEM veto for 2018HEM
    bool passHEM = true;
    if (applyHEMjetVeto_ && !passHEMjetVeto()) passHEM = false;

    for (auto & hg : histograms) {
      if (hg->name.Contains(TString("hCut"))) {
	double cutHistWt = 1;
	if (hg->name.Contains(TString("Wt"))) {
	  cutHistWt = eventWt;
	  if (applySFwtToMC_ && isMC_) cutHistWt *= effWt_;
	}
	cutHistFiller.fill((TH1D*) hg->hist, cutHistWt, passTrg, passHEM);
	continue;
      }
      if (!passTrg) break;
      if (!passHEM) break;
      // bool keep = false; for (auto & theE : *Electrons) {if (!passHEMobjVeto(theE)) keep = true;}  if (!keep) break;

      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt == 0) continue;

      if (hg->name.Contains("hCC_")) {
	countInSel++;
	if (eventWt < 0) countNegWt++;
      }

      double eventWt0 = eventWt;
      if ((applySFwtToMC_ && isMC_) || hg->name.Contains(TString("_DR"))) {
	// For MC, or double ratio, apply weights for trigger eff, reco eff.
	if (hg->name.Contains(TString("_DR"))) {
	  if (BTags > 0) continue;
	  int CCbin = CCbins_->jbk(CCbins_->jbin(NJets), CCbins_->bbin(NJets, BTags), CCbins_->kinBin(HT, MHT));
	  if (CCbin <= 0) continue;
	  if ((UInt_t) CCbin != RA2bin) cout << "CCbin = " << CCbin << ", != " << RA2bin << endl;
	}
	eventWt *= effWt_;
      }

      if (hg->dvalue != nullptr) {
	hg->hist->Fill(*(hg->dvalue), selWt*eventWt);
      }
      else if (hg->ivalue != nullptr) {
	hg->hist->Fill(Double_t(*(hg->ivalue)), selWt*eventWt);
      }
      else if (hg->filler1D != nullptr) (this->*(hg->filler1D))((TH1D*) hg->hist, selWt*eventWt);
      else if (hg->filler2D != nullptr) (this->*(hg->filler2D))((TH2F*) hg->hist, selWt*eventWt);
      else cerr << "No method to fill histogram provided for " << hg->name << endl;
      eventWt = eventWt0;  // Restore event weight after processing efficiency, purity weighted histograms

    }  // loop over histograms
  }  // loop over entries
  cout << "At end, count = " << countInSel << ", with negative weights = " << countNegWt << endl;

  delete forNotify;
  delete chain->GetCurrentFile();
  if (btagcorr_) delete btagcorr_;

}  // ======================================================================================

std::vector<TH1*>
RA2bZinvAnalysis::makeHistograms(const char* sample) {
  //
  // Define histograms, variables to fill, and cuts to be modified.
  // Return a vector of histograms.
  // For an event variable of type double (int), set member dvalue (ivalue)
  // to point to the tree variable.  For other types of tree branches, supply
  // a function to fill the histogram, and set member filler to point to it.
  //
  // When drawing, control number of digits in axis labels by
  // TGaxis myTGaxis;
  // myTGaxis.SetMaxDigits(4);

  std::vector<histConfig*> histograms;
  TCut baselineCuts = getCuts(sample);
  TString sampleKey = sampleKeyMap_.count(sample) > 0 ? sampleKeyMap_.at(sample) : TString("none");
  bool isZll = (sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll");
  bool isPhoton = (sampleKey == "photon" || sampleKey =="photonqcd");

  histConfig hCutFlow;
  hCutFlow.name = TString("hCutFlow_") + sample;  hCutFlow.title = "Cut flow unweighted";
  hCutFlow.NbinsX = 12;  hCutFlow.rangeX.first = 0;  hCutFlow.rangeX.second = 12;
  hCutFlow.axisTitles.first = "";  hCutFlow.axisTitles.second = "Events surviving";
  histograms.push_back(&hCutFlow);

  histConfig hCutFlowWt;
  hCutFlowWt.name = TString("hCutFlowWt_") + sample;  hCutFlowWt.title = "Cut flow with weights";
  hCutFlowWt.NbinsX = 12;  hCutFlowWt.rangeX.first = 0;  hCutFlowWt.rangeX.second = 12;
  hCutFlowWt.axisTitles.first = "";  hCutFlowWt.axisTitles.second = "Weighted events surviving";
  histograms.push_back(&hCutFlowWt);

  histConfig hCuts;
  hCuts.name = TString("hCuts_") + sample;  hCuts.title = "Cuts passed";
  hCuts.NbinsX = 12;  hCuts.rangeX.first = 0;  hCuts.rangeX.second = 12;
  hCuts.axisTitles.first = "";  hCuts.axisTitles.second = "Events passing";
  histograms.push_back(&hCuts);
  // These hCut* histograms are filled before the trigger cut.

  histConfig hFilterCuts;
  hFilterCuts.name = TString("hFilterCuts_") + sample;  hFilterCuts.title = "Filter cuts";
  hFilterCuts.NbinsX = 16;  hFilterCuts.rangeX.first = 0;  hFilterCuts.rangeX.second = 16;
  hFilterCuts.axisTitles.first = "";  hFilterCuts.axisTitles.second = "Events failing";
  hFilterCuts.filler1D = &RA2bZinvAnalysis::fillFilterCuts;
  hFilterCuts.omitCuts.push_back(&commonCuts_);
  histograms.push_back(&hFilterCuts);

  histConfig hCC;
  hCC.name = TString("hCC_") + sample;  hCC.title = "Cut & count analysis";
  Int_t MaxBins = CCbins_->bins();
  cout << "MaxBins = " << MaxBins << endl;
  hCC.NbinsX = MaxBins;  hCC.rangeX.first = 0.5;  hCC.rangeX.second = MaxBins+0.5;
  hCC.axisTitles.first = "Njets, Nb, (HT, MHT)";  hCC.axisTitles.second = "Events / bin";
  hCC.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCC);

  histConfig hCCjb;
  hCCjb.name = TString("hCCjb_") + sample;  hCCjb.title = "Cut & count, kinematics integrated";
  Int_t MaxBinsjb = CCbins_->binsjb();
  cout << "MaxBinsjb = " << MaxBinsjb << endl;
  hCCjb.NbinsX = MaxBinsjb;  hCCjb.rangeX.first = 0.5;  hCCjb.rangeX.second = MaxBinsjb+0.5;
  hCCjb.axisTitles.first = "Njets, Nb";  hCCjb.axisTitles.second = "Events / bin";
  hCCjb.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCjb);

  histConfig hCCspl;
  hCCspl.name = TString("hCCspl_") + sample;  hCCspl.title = "Cut & count, NJet split";
  Int_t MaxBinsSpl = CCbins_->binsSpl();
  cout << "MaxBinsSpl = " << MaxBinsSpl << endl;
  hCCspl.NbinsX = MaxBinsSpl;  hCCspl.rangeX.first = 0.5;  hCCspl.rangeX.second = MaxBinsSpl+0.5;
  hCCspl.axisTitles.first = "Njets, Nb, (HT, MHT)";  hCCspl.axisTitles.second = "Events / bin";
  hCCspl.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCspl);

  histConfig hCCJb;
  hCCJb.name = TString("hCCJb_") + sample;  hCCJb.title = "Cut & count, kinematics integrated, NJet split";
  Int_t MaxBinsJb = CCbins_->binsJb();
  cout << "MaxBinsJb = " << MaxBinsJb << endl;
  hCCJb.NbinsX = MaxBinsJb;  hCCJb.rangeX.first = 0.5;  hCCJb.rangeX.second = MaxBinsJb+0.5;
  hCCJb.axisTitles.first = "Njets, Nb";  hCCJb.axisTitles.second = "Events / bin";
  hCCJb.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCJb);

  histConfig hCCjk;
  hCCjk.name = TString("hCCjk_") + sample;  hCCjk.title = "Cut & count, Nb = 0";
  Int_t MaxBinsjk = CCbins_->binsjk();
  cout << "MaxBinsjk = " << MaxBinsjk << endl;
  hCCjk.NbinsX = MaxBinsjk;  hCCjk.rangeX.first = 0.5;  hCCjk.rangeX.second = MaxBinsjk+0.5;
  hCCjk.axisTitles.first = "Njets, (HT, MHT)";  hCCjk.axisTitles.second = "Events / bin";
  hCCjk.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCjk);

  histConfig hgenHT;
  hgenHT.name = TString("hgenHT_") + sample;  hgenHT.title = "Generated HT";
  hgenHT.NbinsX = 60;  hgenHT.rangeX.first = 0;  hgenHT.rangeX.second = 3000;
  hgenHT.axisTitles.first = "HT [GeV]";  hgenHT.axisTitles.second = "Events / 50 GeV";
  hgenHT.dvalue = &HT;
  if (isMC_) histograms.push_back(&hgenHT);

  histConfig hPrefire;
  hPrefire.name = TString("hPrefire_") + sample;  hPrefire.title = "Probability to survive trigger prefire";
  hPrefire.NbinsX = 120;  hPrefire.rangeX.first = 0;  hPrefire.rangeX.second = 1.2;
  hPrefire.axisTitles.first = "NonPrefireProb";  hPrefire.axisTitles.second = "Events / 0.01";
  hPrefire.dvalue = &NonPrefiringProb;
  if (isMC_) histograms.push_back(&hPrefire);

  histConfig hHT;
  hHT.name = TString("hHT_") + sample;  hHT.title = "HT";
  hHT.NbinsX = 60;  hHT.rangeX.first = 0;  hHT.rangeX.second = 3000;
  hHT.axisTitles.first = "HT [GeV]";  hHT.axisTitles.second = "Events / 50 GeV";
  hHT.dvalue = &HT;  hHT.omitCuts.push_back(&HTcut_);
  histograms.push_back(&hHT);

  histConfig hMHT;
  hMHT.name = TString("hMHT_") + sample;  hMHT.title = "MHT";
  hMHT.NbinsX = 60;  hMHT.rangeX.first = 0;  hMHT.rangeX.second = 3000;
  hMHT.axisTitles.first = "MHT [GeV]";  hMHT.axisTitles.second = "Events / 50 GeV";
  hMHT.dvalue = &MHT;  hMHT.omitCuts.push_back(&MHTcut_);  hMHT.omitCuts.push_back(&ptCut_);
  histograms.push_back(&hMHT);

  histConfig hNJets;
  hNJets.name = TString("hNJets_") + sample;  hNJets.title = "NJets";
  hNJets.NbinsX = 20;  hNJets.rangeX.first = 0;  hNJets.rangeX.second = 20;
  hNJets.axisTitles.first = "N (jets)";  hNJets.axisTitles.second = "Events / bin";
  hNJets.ivalue = &NJets;  hNJets.omitCuts.push_back(&NJetscut_);
  histograms.push_back(&hNJets);

  histConfig hBTags;
  hBTags.name = TString("hBTags_") + sample;  hBTags.title = "BTags";
  hBTags.NbinsX = 20;  hBTags.rangeX.first = 0;  hBTags.rangeX.second = 20;
  hBTags.axisTitles.first = "N (b jets)";  hBTags.axisTitles.second = "Events / bin";
  hBTags.ivalue = &BTags;
  histograms.push_back(&hBTags);

  histConfig hVertices;
  hVertices.name = TString("hVertices_") + sample;  hVertices.title = "Number of reco vertices";
  hVertices.NbinsX = 100;  hVertices.rangeX.first = 0;  hVertices.rangeX.second = 100;
  hVertices.axisTitles.first = "No. of vertices";  hVertices.axisTitles.second = "Events / bin";
  hVertices.ivalue = &nAllVertices;
  histograms.push_back(&hVertices);

  histConfig hTrueNumInt;
  hTrueNumInt.name = TString("hTrueNumInt_") + sample;  hTrueNumInt.title = "Number of generated interactions";
  hTrueNumInt.NbinsX = 100;  hTrueNumInt.rangeX.first = 0;  hTrueNumInt.rangeX.second = 100;
  hTrueNumInt.axisTitles.first = "No. of interactions";  hTrueNumInt.axisTitles.second = "Events / bin";
  hTrueNumInt.dvalue = &TrueNumInteractions;
  if (isMC_) histograms.push_back(&hTrueNumInt);

  histConfig hPUwtvsNint;
  hPUwtvsNint.name = TString("hPUwtvsNint_") + sample;  hPUwtvsNint.title = "PU wt vs Number of generated interactions";
  hPUwtvsNint.NbinsX = 100;  hPUwtvsNint.rangeX.first = 0;  hPUwtvsNint.rangeX.second = 100;
  hPUwtvsNint.NbinsY = 120;  hPUwtvsNint.rangeY.first = 0;  hPUwtvsNint.rangeY.second = 6;
  hPUwtvsNint.axisTitles.first = "No. of interactions";  hPUwtvsNint.axisTitles.second = "PU weight";
  hPUwtvsNint.filler2D = &RA2bZinvAnalysis::fillPUwtvsNint;
  if (isMC_) histograms.push_back(&hPUwtvsNint);

  // For double ratio

  Double_t hHT_DR_bins[9] = {300, 400, 450, 550, 650, 750, 1000, 1200, 1600};  // Check if CCbinning changes
  histConfig hHT_DR;
  hHT_DR.name = TString("hHT_DR_") + sample;  hHT_DR.title = "HT";
  hHT_DR.NbinsX = 8;
  hHT_DR.binsX = hHT_DR_bins;
  hHT_DR.axisTitles.first = "HT [GeV]";  hHT_DR.axisTitles.second = "Events / bin";
  hHT_DR.dvalue = &HT;
  histograms.push_back(&hHT_DR);

  histConfig hHT_DR_xWt;  // for weighted centers of bins
  hHT_DR_xWt.name = TString("hHT_DR_xWt_") + sample;  hHT_DR_xWt.title = "HT bin values";
  hHT_DR_xWt.NbinsX = hHT_DR.NbinsX;
  hHT_DR_xWt.binsX = hHT_DR.binsX;
  hHT_DR_xWt.axisTitles.first = "HT [GeV]";  hHT_DR_xWt.axisTitles.second = "Bin value / bin";
  hHT_DR_xWt.filler1D = &RA2bZinvAnalysis::fillHT_DR_xWt;
  if (isPhoton) histograms.push_back(&hHT_DR_xWt);

  std::vector<double> HTthresh1;
  // Check if CCbinning changes [
  for (Size_t i = 0; i < CCbins_->binsht(0); ++i) HTthresh1.push_back(CCbins_->htThreshold(0, i));
  HTthresh1.push_back(2500);
  // cout << "HTthresh1 = " ; for (Size_t i = 0; i < HTthresh1.size(); ++i) cout << HTthresh1[i] << ", "; cout << endl;
  //  ]
  Double_t* hHT1_DRCC_bins = &HTthresh1[0];
  histConfig hHT1_DRCC;
  hHT1_DRCC.name = TString("hHT1_DRCC_") + sample;  hHT1_DRCC.title = "HT for DR fit 1";
  hHT1_DRCC.NbinsX = HTthresh1.size() - 1;
  hHT1_DRCC.binsX = hHT1_DRCC_bins;
  hHT1_DRCC.axisTitles.first = "HT [GeV]";  hHT1_DRCC.axisTitles.second = "Events / bin";
  hHT1_DRCC.dvalue = &HT;  hHT1_DRCC.omitCuts.push_back(&HTcut_);
  if (isPhoton) histograms.push_back(&hHT1_DRCC);
  histConfig hHT1_DRCC_xWt;  // for weighted centers of bins
  hHT1_DRCC_xWt.name = TString("hHT1_DRCC_xWt_") + sample;  hHT1_DRCC_xWt.title = "HT CC bin values";
  hHT1_DRCC_xWt.NbinsX = hHT1_DRCC.NbinsX;
  hHT1_DRCC_xWt.binsX = hHT1_DRCC.binsX;
  hHT1_DRCC_xWt.axisTitles.first = "HT [GeV]";  hHT1_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hHT1_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillHT_DR_xWt;
  hHT1_DRCC_xWt.omitCuts.push_back(&HTcut_);
  if (isPhoton) histograms.push_back(&hHT1_DRCC_xWt);

  std::vector<double> HTthresh2;
  // Check if CCbinning changes [
  HTthresh2.push_back(CCbins_->htThreshold(0, 0));
  HTthresh2.push_back(CCbins_->htThreshold(1, 0));
  HTthresh2.push_back(CCbins_->htThreshold(2, 0));
  HTthresh2.push_back(CCbins_->htThreshold(3, 0));
  HTthresh2.push_back(CCbins_->htThreshold(3, 1));
  HTthresh2.push_back(2500);
  // cout << "HTthresh2 = " ; for (Size_t i = 0; i < HTthresh2.size(); ++i) cout << HTthresh2[i] << ", "; cout << endl;
  //  ]
  Double_t* hHT2_DRCC_bins = &HTthresh2[0];
  histConfig hHT2_DRCC;
  hHT2_DRCC.name = TString("hHT2_DRCC_") + sample;  hHT2_DRCC.title = "HT for DR fit 2";
  hHT2_DRCC.NbinsX = HTthresh2.size() - 1;
  hHT2_DRCC.binsX = hHT2_DRCC_bins;
  hHT2_DRCC.axisTitles.first = "HT [GeV]";  hHT2_DRCC.axisTitles.second = "Events / bin";
  hHT2_DRCC.dvalue = &HT;  hHT2_DRCC.omitCuts.push_back(&HTcut_);
  if (isPhoton) histograms.push_back(&hHT2_DRCC);
  histConfig hHT2_DRCC_xWt;  // for weighted centers of bins
  hHT2_DRCC_xWt.name = TString("hHT2_DRCC_xWt_") + sample;  hHT2_DRCC_xWt.title = "HT CC bin values";
  hHT2_DRCC_xWt.NbinsX = hHT2_DRCC.NbinsX;
  hHT2_DRCC_xWt.binsX = hHT2_DRCC.binsX;
  hHT2_DRCC_xWt.axisTitles.first = "HT [GeV]";  hHT2_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hHT2_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillHT_DR_xWt;
  hHT2_DRCC_xWt.omitCuts.push_back(&HTcut_);
  if (isPhoton) histograms.push_back(&hHT2_DRCC_xWt);

  Double_t hMHT_DR_bins[7] = {300., 350., 400., 450., 600., 750, 900.};  // Check if CCbinning changes
  histConfig hMHT_DR;
  hMHT_DR.name = TString("hMHT_DR_") + sample;  hMHT_DR.title = "MHT for DR fit";
  hMHT_DR.NbinsX = 6;
  hMHT_DR.binsX = hMHT_DR_bins;
  hMHT_DR.axisTitles.first = "MHT [GeV]";  hMHT_DR.axisTitles.second = "Events / bin";
  hMHT_DR.dvalue = &MHT;
  histograms.push_back(&hMHT_DR);

  histConfig hMHT_DR_xWt;  // for weighted centers of bins
  hMHT_DR_xWt.name = TString("hMHT_DR_xWt_") + sample;  hMHT_DR_xWt.title = "MHT bin values";
  hMHT_DR_xWt.NbinsX = hMHT_DR.NbinsX;
  hMHT_DR_xWt.binsX = hMHT_DR.binsX;
  hMHT_DR_xWt.axisTitles.first = "MHT [GeV]";  hMHT_DR_xWt.axisTitles.second = "Bin value / bin";
  hMHT_DR_xWt.filler1D = &RA2bZinvAnalysis::fillMHT_DR_xWt;
  if (isPhoton) histograms.push_back(&hMHT_DR_xWt);

  std::vector<double> MHTthresh;
  for (Size_t m = 0; m < CCbins_->binsmht(); ++m) MHTthresh.push_back(CCbins_->mhtThreshold(m));
  MHTthresh.insert(MHTthresh.begin(), MHTthresh.back());
  MHTthresh.pop_back();  MHTthresh.push_back(900);  // Check if CCbinning changes
  Double_t* hMHT_DRCC_bins = &MHTthresh[0];
  histConfig hMHT_DRCC;
  hMHT_DRCC.name = TString("hMHT_DRCC_") + sample;  hMHT_DRCC.title = "MHT for DR";
  hMHT_DRCC.NbinsX = MHTthresh.size() - 1;
  hMHT_DRCC.binsX = hMHT_DRCC_bins;
  hMHT_DRCC.axisTitles.first = "MHT [GeV]";  hMHT_DRCC.axisTitles.second = "Events / bin";
  hMHT_DRCC.dvalue = &MHT;  hMHT_DRCC.omitCuts.push_back(&MHTcut_);  hMHT_DRCC.omitCuts.push_back(&ptCut_);
  if (isPhoton) histograms.push_back(&hMHT_DRCC);
  histConfig hMHT_DRCC_xWt;  // for weighted centers of bins
  hMHT_DRCC_xWt.name = TString("hMHT_DRCC_xWt_") + sample;  hMHT_DRCC_xWt.title = "MHT CC bin values";
  hMHT_DRCC_xWt.NbinsX = hMHT_DRCC.NbinsX;
  hMHT_DRCC_xWt.binsX = hMHT_DRCC.binsX;
  hMHT_DRCC_xWt.axisTitles.first = "MHT [GeV]";  hMHT_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hMHT_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillMHT_DR_xWt;
  hMHT_DRCC_xWt.omitCuts.push_back(&MHTcut_);  hMHT_DRCC_xWt.omitCuts.push_back(&ptCut_);
  if (isPhoton) histograms.push_back(&hMHT_DRCC_xWt);

  histConfig hNJets_DR;
  hNJets_DR.name = TString("hNJets_DR_") + sample;  hNJets_DR.title = "NJets";
  hNJets_DR.NbinsX = 8;  hNJets_DR.rangeX.first = 1.5;  hNJets_DR.rangeX.second = 9.5;
  hNJets_DR.axisTitles.first = "N (jets)";  hNJets_DR.axisTitles.second = "Events / bin";
  hNJets_DR.ivalue = &NJets;
  histograms.push_back(&hNJets_DR);

  std::vector<double> NJetsthresh;
  for (auto thresh : CCbins_->nJetThresholds()) NJetsthresh.push_back(thresh);
  NJetsthresh.push_back(12);  // Check if CCbinning changes
  Double_t* hNJets_DRCC_bins = &NJetsthresh[0];
  histConfig hNJets_DRCC;
  hNJets_DRCC.name = TString("hNJets_DRCC_") + sample;  hNJets_DRCC.title = "NJets for DR";
  hNJets_DRCC.NbinsX = NJetsthresh.size() - 1;
  hNJets_DRCC.binsX = hNJets_DRCC_bins;
  hNJets_DRCC.axisTitles.first = "N (jets)";  hNJets_DRCC.axisTitles.second = "Events / bin";
  hNJets_DRCC.ivalue = &NJets;
  histograms.push_back(&hNJets_DRCC);
  histConfig hNJets_DRCC_xWt;  // for weighted centers of bins
  hNJets_DRCC_xWt.name = TString("hNJets_DRCC_xWt_") + sample;  hNJets_DRCC_xWt.title = "NJets CC bin values";
  hNJets_DRCC_xWt.NbinsX = hNJets_DRCC.NbinsX;
  hNJets_DRCC_xWt.binsX = hNJets_DRCC.binsX;
  hNJets_DRCC_xWt.axisTitles.first = "N (jets)";  hNJets_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hNJets_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillNJets_DR_xWt;
  if (isPhoton) histograms.push_back(&hNJets_DRCC_xWt);

  histConfig hSFwt_DR;
  hSFwt_DR.name = TString("hSFwt_DR_") + sample;  hSFwt_DR.title = "Weights for efficiency, purity";
  hSFwt_DR.NbinsX = 80;  hSFwt_DR.rangeX.first = 0.8;  hSFwt_DR.rangeX.second = 1.2;
  hSFwt_DR.axisTitles.first = "SF weight";  hSFwt_DR.axisTitles.second = "Events / bin";
  hSFwt_DR.filler1D = &RA2bZinvAnalysis::fillSFwt_DR;
  histograms.push_back(&hSFwt_DR);

  histConfig hSFsys_DR;
  hSFsys_DR.name = TString("hSFsys_DR_") + sample;  hSFsys_DR.title = "Syst errors for efficiency, purity";
  hSFsys_DR.NbinsX = 100;  hSFsys_DR.rangeX.first = 0;  hSFsys_DR.rangeX.second = 0.2;
  hSFsys_DR.axisTitles.first = "SF error";  hSFsys_DR.axisTitles.second = "Events / bin";
  hSFsys_DR.filler1D = &RA2bZinvAnalysis::fillSFsys_DR;
  histograms.push_back(&hSFsys_DR);

  histConfig hnZcand;
  hnZcand.name = TString("hnZcand_") + sample;  hnZcand.title = "Number of Z candidates";
  hnZcand.NbinsX = 10;  hnZcand.rangeX.first = 0;  hnZcand.rangeX.second = 10;
  hnZcand.axisTitles.first = "N(Z candidates)";  hnZcand.axisTitles.second = "Events / bin";
  hnZcand.filler1D = &RA2bZinvAnalysis::fillnZcand;  hnZcand.omitCuts.push_back(&massCut_);
  histograms.push_back(&hnZcand);

  histConfig hgenZpt;
  hgenZpt.name = TString("hgenZpt_") + sample;  hgenZpt.title = "Generated Zpt";
  hgenZpt.NbinsX = 60;  hgenZpt.rangeX.first = 0;  hgenZpt.rangeX.second = 3000;
  hgenZpt.axisTitles.first = "Gen P_t(Z) [GeV]";  hgenZpt.axisTitles.second = "Events / 50 GeV";
  hgenZpt.filler1D = &RA2bZinvAnalysis::fillgenZpt;
  if (isMC_) histograms.push_back(&hgenZpt);

  histConfig hZpt;
  hZpt.name = TString("hZpt_") + sample;  hZpt.title = "Z Pt";
  hZpt.NbinsX = 60;  hZpt.rangeX.first = 0;  hZpt.rangeX.second = 3000;
  hZpt.axisTitles.first = "Pt(Z) [GeV]";  hZpt.axisTitles.second = "Events / 50 GeV";
  hZpt.filler1D = &RA2bZinvAnalysis::fillZpt;  hZpt.omitCuts.push_back(&ptCut_);  hZpt.omitCuts.push_back(&MHTcut_);
  if (isZll) histograms.push_back(&hZpt);

  histConfig hPhotonPt;
  hPhotonPt.name = TString("hPhotonPt_") + sample;  hPhotonPt.title = "Photon Pt";
  hPhotonPt.NbinsX = 60;  hPhotonPt.rangeX.first = 0;  hPhotonPt.rangeX.second = 3000;
  hPhotonPt.axisTitles.first = "Pt(Photon) [GeV]";  hPhotonPt.axisTitles.second = "Events / 50 GeV";
  hPhotonPt.filler1D = &RA2bZinvAnalysis::fillPhotonPt;  hPhotonPt.omitCuts.push_back(&ptCut_);  hPhotonPt.omitCuts.push_back(&MHTcut_);
  if (isPhoton) histograms.push_back(&hPhotonPt);

  histConfig hPhotonEta;
  hPhotonEta.name = TString("hPhotonEta_") + sample;  hPhotonEta.title = "Photon Eta";
  hPhotonEta.NbinsX = 60;  hPhotonEta.rangeX.first = -3;  hPhotonEta.rangeX.second = 3;
  hPhotonEta.axisTitles.first = "Eta(photon) [GeV]";  hPhotonEta.axisTitles.second = "Events / 0.1";
  hPhotonEta.filler1D = &RA2bZinvAnalysis::fillPhotonEta;
  if (isPhoton) histograms.push_back(&hPhotonEta);

  histConfig hMuonEta;
  hMuonEta.name = TString("hMuonEta_") + sample;  hMuonEta.title = "Muon Eta";
  hMuonEta.NbinsX = 60;  hMuonEta.rangeX.first = -3;  hMuonEta.rangeX.second = 3;
  hMuonEta.axisTitles.first = "Eta(muon) [GeV]";  hMuonEta.axisTitles.second = "Events / 0.1";
  hMuonEta.filler1D = &RA2bZinvAnalysis::fillMuonEta;
  if (sampleKey == "zmm") histograms.push_back(&hMuonEta);

  histConfig hElectronEta;
  hElectronEta.name = TString("hElectronEta_") + sample;  hElectronEta.title = "Electron Eta";
  hElectronEta.NbinsX = 60;  hElectronEta.rangeX.first = -3;  hElectronEta.rangeX.second = 3;
  hElectronEta.axisTitles.first = "Eta(electron) [GeV]";  hElectronEta.axisTitles.second = "Events / 0.1";
  hElectronEta.filler1D = &RA2bZinvAnalysis::fillElectronEta;
  if (sampleKey == "zee") histograms.push_back(&hElectronEta);

  histConfig hZmass;
  hZmass.name = TString("hZmass_") + sample;  hZmass.title = "Z mass";
  hZmass.NbinsX = 30;  hZmass.rangeX.first = 60;  hZmass.rangeX.second = 120;
  hZmass.axisTitles.first = "M(Z) [GeV]";  hZmass.axisTitles.second = "Events / 2 GeV";
  hZmass.filler1D = &RA2bZinvAnalysis::fillZmass;  hZmass.omitCuts.push_back(&massCut_);
  if (isZll) histograms.push_back(&hZmass);

  // Z mass in Njet, Nb bins
  histConfig hZmass_2j0b(hZmass);
  hZmass_2j0b.name = TString("hZmass_2j0b_") + sample;  hZmass_2j0b.title = "Z mass, 2 jets & 0 b jets";
  hZmass_2j0b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_2j0b);

  histConfig hZmass_2j1b(hZmass);
  hZmass_2j1b.name = TString("hZmass_2j1b_") + sample;  hZmass_2j1b.title = "Z mass, 2 jets & 1 b jet";
  hZmass_2j1b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_2j1b);

  histConfig hZmass_2j2b(hZmass);
  hZmass_2j2b.name = TString("hZmass_2j2b_") + sample;  hZmass_2j2b.title = "Z mass, 2 jets & >=2 b jets";
  hZmass_2j2b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_2j2b);
  //
  histConfig hZmass_3j0b(hZmass);
  hZmass_3j0b.name = TString("hZmass_3j0b_") + sample;  hZmass_3j0b.title = "Z mass, 3-4 jets & 0 b jets";
  hZmass_3j0b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_3j0b);

  histConfig hZmass_3j1b(hZmass);
  hZmass_3j1b.name = TString("hZmass_3j1b_") + sample;  hZmass_3j1b.title = "Z mass, 3-4 jets & 1 b jet";
  hZmass_3j1b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_3j1b);

  histConfig hZmass_3j2b(hZmass);
  hZmass_3j2b.name = TString("hZmass_3j2b_") + sample;  hZmass_3j2b.title = "Z mass, 3-4 jets & >=2 b jets";
  hZmass_3j2b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_3j2b);
  //
  histConfig hZmass_5j0b(hZmass);
  hZmass_5j0b.name = TString("hZmass_5j0b_") + sample;  hZmass_5j0b.title = "Z mass, >=5 jets & 0 B jets";
  hZmass_5j0b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_5j0b);

  histConfig hZmass_5j1b(hZmass);
  hZmass_5j1b.name = TString("hZmass_5j1b_") + sample;  hZmass_5j1b.title = "Z mass, >=5 jets & 1 B jet";
  hZmass_5j1b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_5j1b);

  histConfig hZmass_5j2b(hZmass);
  hZmass_5j2b.name = TString("hZmass_5j2b_") + sample;  hZmass_5j2b.title = "Z mass, >=5 jets & >=2 B jets";
  hZmass_5j2b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  if (isZll) histograms.push_back(&hZmass_5j2b);

  // Special studies

  // histConfig hZmass_sfLepTksVeto(hZmass);
  // hZmass_sfLepTksVeto.name = TString("hZmass_sfLepTksVeto_") + sample;  hZmass_sfLepTksVeto.title = "Z mass, SF lepton vetoed";
  // hZmass_sfLepTksVeto.omitCuts.push_back(&isoSFlepTksCut_);  hZmass_sfLepTksVeto.addCuts = isoSFlepTksVeto_.Data();
  // histograms.push_back(&hZmass_sfLepTksVeto);

  // histConfig hZmass_photonVeto(hZmass);
  // hZmass_photonVeto.name = TString("hZmass_photonVeto_") + sample;  hZmass_photonVeto.title = "Z mass, photon vetoed";
  // hZmass_photonVeto.omitCuts.push_back(&photonVeto_);  hZmass_photonVeto.addCuts = photonCut_.Data();
  // histograms.push_back(&hZmass_photonVeto);

  // histConfig hGpt;
  // hGpt.name = TString("hGpt_") + sample;  hGpt.title = "Photon Pt";
  // hGpt.NbinsX = 60;  hGpt.rangeX.first = 0;  hGpt.rangeX.second = 3000;
  // hGpt.axisTitles.first = "Pt(gamma) [GeV]";  hGpt.axisTitles.second = "Events / 50 GeV";
  // hGpt.filler1D = &RA2bZinvAnalysis::fillGpt;
  // hGpt.omitCuts.push_back(&photonVeto_);  hGpt.addCuts = photonCut_.Data();
  // histograms.push_back(&hGpt);

  // histConfig hZGmass;
  // hZGmass.name = TString("hZGmass_") + sample;  hZGmass.title = "Z-gamma mass";
  // hZGmass.NbinsX = 100;  hZGmass.rangeX.first = 0;  hZGmass.rangeX.second = 2000;
  // hZGmass.axisTitles.first = "M(Z gamma) [GeV]";  hZGmass.axisTitles.second = "Events / 20 GeV";
  // hZGmass.filler1D = &RA2bZinvAnalysis::fillZGmass;
  // hZGmass.omitCuts.push_back(&photonVeto_);  hZGmass.addCuts = photonCut_.Data();
  // histograms.push_back(&hZGmass);

  // histConfig hZGdRvsM;
  // hZGdRvsM.name = TString("hZGdRvsM_")+sample;
  // hZGdRvsM.title = "Min Delta R vs M(Zgamma) for Z(ll) leptons";
  // hZGdRvsM.NbinsX = 100;  hZGdRvsM.rangeX.first = 0;  hZGdRvsM.rangeX.second = 200;
  // hZGdRvsM.NbinsY = 40;  hZGdRvsM.rangeY.first = 0;  hZGdRvsM.rangeY.second = 0.02;
  // hZGdRvsM.axisTitles.first = "M(Z gamma) [GeV]";  hZGdRvsM.axisTitles.second = "DR(l gamma)";
  // hZGdRvsM.filler2D = &RA2bZinvAnalysis::fillZGdRvsM;
  // hZGdRvsM.omitCuts.push_back(&photonVeto_);  hZGdRvsM.addCuts = photonCut_.Data();
  // histograms.push_back(&hZGdRvsM);

  // histConfig hGJdR;
  // hGJdR.name = TString("hGJdR_") + sample;  hGJdR.title = "min DR(photon-jet)";
  // hGJdR.NbinsX = 200;  hGJdR.rangeX.first = 0;  hGJdR.rangeX.second = 4;
  // hGJdR.axisTitles.first = "Delta R";  hGJdR.axisTitles.second = "Events / 0.02";
  // hGJdR.filler1D = &RA2bZinvAnalysis::fillGJdR;
  // hGJdR.omitCuts.push_back(&photonVeto_);  hGJdR.addCuts = photonCut_.Data();
  // histograms.push_back(&hGJdR);

  // histConfig hGLdRnoPixelSeed;
  // hGLdRnoPixelSeed.name = TString("hGLdRnoPixelSeed_") + sample;  hGLdRnoPixelSeed.title = "min DR(photon-lepton), noPixelSeed";
  // hGLdRnoPixelSeed.NbinsX = 200;  hGLdRnoPixelSeed.rangeX.first = 0;  hGLdRnoPixelSeed.rangeX.second = 4;
  // hGLdRnoPixelSeed.axisTitles.first = "Delta R";  hGLdRnoPixelSeed.axisTitles.second = "Events / 0.02";
  // hGLdRnoPixelSeed.filler1D = &RA2bZinvAnalysis::fillGLdRnoPixelSeed;
  // hGLdRnoPixelSeed.omitCuts.push_back(&photonVeto_);  hGLdRnoPixelSeed.addCuts = photonCut_.Data();
  // histograms.push_back(&hGLdRnoPixelSeed);

  // histConfig hGLdRpixelSeed;
  // hGLdRpixelSeed.name = TString("hGLdRpixelSeed_") + sample;  hGLdRpixelSeed.title = "min DR(photon-lepton), pixelSeed";
  // hGLdRpixelSeed.NbinsX = 200;  hGLdRpixelSeed.rangeX.first = 0;  hGLdRpixelSeed.rangeX.second = 4;
  // hGLdRpixelSeed.axisTitles.first = "Delta R";  hGLdRpixelSeed.axisTitles.second = "Events / 0.02";
  // hGLdRpixelSeed.filler1D = &RA2bZinvAnalysis::fillGLdRpixelSeed;
  // hGLdRpixelSeed.omitCuts.push_back(&photonVeto_);  hGLdRpixelSeed.addCuts = photonCut_.Data();
  // histograms.push_back(&hGLdRpixelSeed);

  bookAndFillHistograms(sample, histograms, baselineCuts);

  std::vector<TH1*> theHists;
  for (auto & thisHist : histograms) theHists.push_back(thisHist->hist);

  return theHists;

}  // ======================================================================================

void
RA2bZinvAnalysis::fillCC(TH1D* h, double wt) {

  // Filler for the analysis-bin histograms and variants
  TString hName(h->GetName());
  Int_t binCC = 0, binCCjb = 0;

  int binKin = CCbins_->kinBin(HT, MHT);
  if (binKin < 0) return;
  int binNjets = (hName.Contains("spl") || hName.Contains("Jb")) ? CCbins_->Jbin(NJets) : CCbins_->jbin(NJets);
  if (binNjets < 0) return;
  if (!applyBTagSF_) {
    int binNb = CCbins_->bbin(NJets, BTags);
    binCC = (hName.Contains("spl") || hName.Contains("Jb")) ?
      CCbins_->Jbk(binNjets, binNb, binKin) :
      useTreeCCbin_ ? RA2bin : CCbins_->jbk(binNjets, binNb, binKin);
    if (binCC <= 0) return;
    if (verbosity_ >= 4) cout << "j = " << binNjets << ", b = " << binNb << ", k = " << binKin << ", binCC = "
			      << binCC << ", RA2bin = " << RA2bin << endl;
    if (hName.Contains("jb") || hName.Contains("Jb")) {
      // Above test on j, b, k needed even here, to exclude j = 3,4, k = 0,3
      binCCjb = hName.Contains("jb") ? CCbins_->jb(binNjets, binNb) : CCbins_->Jb(binNjets, binNb);
      if (binCCjb > 0) h->Fill(Double_t(binCCjb), wt);
    } else if (hName.Contains("jk")) {
      if (BTags == 0) {
	int binCCjk = CCbins_->jk(binNjets, binKin);
	if (binCCjk > 0) h->Fill(Double_t(binCCjk), wt);
      }
    } else {
      h->Fill(Double_t(binCC), wt);
    }
  } else {
    // apply BTagSF to all Nb bins
    if (verbosity_ >= 4) cout << "Size of input Jets = " << Jets->size()
			      << ", Jets_hadronFlavor = " << Jets_hadronFlavor->size()
			      << " Jets_HTMask = " << Jets_HTMask->size() << endl;
    vector<double> probNb = btagcorr_->GetCorrections(Jets, Jets_hadronFlavor, Jets_HTMask);
    for (int binNb = 0; binNb < CCbins_->binsb(CCbins_->jbin(NJets)); ++binNb) {
      binCC = hName.Contains("spl") ? CCbins_->Jbk(binNjets, binNb, binKin) : CCbins_->jbk(binNjets, binNb, binKin);
      if (binCC <= 0) return;
      if (verbosity_ >= 4) cout << "j = " << binNjets << ", NbTags = " << BTags << ", b = " << binNb
				<< ", b wt = " << probNb[binNb] << ", k = " << binKin << ", binCC = " << binCC << endl;
      if (hName.Contains("jb") || hName.Contains("Jb")) {
	binCCjb = hName.Contains("jb") ? CCbins_->jb(binNjets, binNb) : CCbins_->Jb(binNjets, binNb);
	if (binCCjb <= 0) return;
	h->Fill(Double_t(binCCjb), wt*probNb[binNb]);
      } else if (hName.Contains("jk")) {
	int binCCjk = CCbins_->jk(binNjets, binKin);
	if (binCCjk > 0) h->Fill(Double_t(binCCjk), wt*probNb[0]);
      } else {
	h->Fill(Double_t(binCC), wt*probNb[binNb]);
      }
    }
  }  // if apply BTagSF

}  // ======================================================================================

void
RA2bZinvAnalysis::fillZmassjb(TH1D* h, double wt) {
  if (ZCandidates->size() == 0) return;
  Int_t j, b;
  TString hName(h->GetName());
  j = hName(hName.First('j')-1) - '0';
  b = hName(hName.First('b')-1) - '0';
  // j = 2, 3, 5 in the histogram name really means
  // NJets in the first, second, or third and higher NJets bin, respectively.
  if ((j == 2 && NJets >= CCbins_->nJetThreshold(0) && NJets < CCbins_->nJetThreshold(1)) ||
      (j == 3 && (NJets >= CCbins_->nJetThreshold(1) && NJets < CCbins_->nJetThreshold(2))) ||
      (j == 5 && NJets >= CCbins_->nJetThreshold(2))) {
    if ((b < 2 && BTags == b) || (b == 2 && BTags >= 2)) {
      h->Fill(ZCandidates->at(0).M(), wt);
    }
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::fillFilterCuts(TH1D* h, double wt) {
  h->Fill(0.5, wt);
  if (!(globalTightHalo2016Filter==1)) h->Fill(1.5, wt);
  if (ntupleVersion_ != "V12" && ntupleVersion_ != "V15" && !(globalSuperTightHalo2016Filter==1)) h->Fill(2.5, wt);
  if (!(HBHENoiseFilter==1)) h->Fill(3.5, wt);
  if (!(HBHEIsoNoiseFilter==1)) h->Fill(4.5, wt);
  if (!(EcalDeadCellTriggerPrimitiveFilter==1)) h->Fill(5.5, wt);
  if (!BadChargedCandidateFilter) h->Fill(6.5, wt);
  if (!BadPFMuonFilter) h->Fill(7.5, wt);
  if (!(NVtx > 0)) h->Fill(8.5, wt);
  if (!(eeBadScFilter==1)) h->Fill(9.5, wt);
  if (!isMC_ && !(ecalBadCalibFilter==1)) h->Fill(10.5, wt);
  if (!(JetID)) h->Fill(12.5, wt);
  if (!(PFCaloMETRatio < 5)) h->Fill(13.5, wt);
  if (era_ == "2016" && !(HT5/HT <= 2))
    h->Fill(14.5, wt);
  // else if (!(DeltaPhi1 >= 1.025*HT5/HT - 0.5875))
  else if (!HTRatioDPhiFilter)
    h->Fill(14.5, wt);
  if (isSkim_ && ntupleVersion_ == "V16" && runBlock_.find("2017") != std::string::npos && !EcalNoiseJetFilter)
    h->Fill(15.5, wt);

}  // ======================================================================================

void
RA2bZinvAnalysis::fillZGmass(TH1D* h, double wt) {
  for (auto & theZ : *ZCandidates) {
    for (auto & aPhoton : *Photons) {
      TLorentzVector Zg(theZ);  Zg += aPhoton;
      h->Fill(Zg.M(), wt);
    }
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGJdR(TH1D* h, double wt) {
  TLorentzVector thePhoton = Photons->at(0);
  if (Jets->size() > 0) {
    Double_t dR = 999.;
    for (auto & thisJet : *Jets)
	dR = thePhoton.DeltaR(thisJet) < dR ? thePhoton.DeltaR(thisJet) : dR;
    h->Fill(dR, wt);
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillZGdRvsM(TH2F* h, double wt) {
  if (ZCandidates->size() == 0) return;
  TLorentzVector theZ = ZCandidates->at(0);
  TLorentzVector thePhoton = Photons->at(0);
  TLorentzVector Zg(theZ);  Zg += thePhoton;
  if (Muons->size() > 0) {
    Double_t dR = 999.;
    for (auto & thisMuon : *Muons)
	dR = thePhoton.DeltaR(thisMuon) < dR ? thePhoton.DeltaR(thisMuon) : dR;
    h->Fill(Zg.M(), dR, wt);
  }
  if (Electrons->size() > 0) {
    Double_t dR = 999.;
    for (auto & thisElectron : *Electrons)
	dR = thePhoton.DeltaR(thisElectron) < dR ? thePhoton.DeltaR(thisElectron) : dR;
    h->Fill(Zg.M(), dR, wt);
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGLdRnoPixelSeed(TH1D* h, double wt) {
  for (size_t i = 0; i < Photons->size(); ++i) {
    // if (Photons_fullID->size()) cout << "Photons_fullID = " << Photons_fullID->at(i) << endl;
    // if (Photons_passElectronVeto->size()) cout << "Photons_passElectronVeto = " << Photons_passElectronVeto->at(i) << endl;
    if (!Photons_fullID->at(i)) continue;
    if (Photons_hasPixelSeed->at(i)) continue;
    if (Muons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisMuon : *Muons)
	dR = Photons->at(i).DeltaR(thisMuon) < dR ? Photons->at(i).DeltaR(thisMuon) : dR;
      h->Fill(dR, wt);
      // if (Photons_isEB->size() && dR < 0.02) cout << "isEB photon = " << Photons_isEB->at(i) << endl;
      // if (Photons_hasPixelSeed->size() && dR < 0.02) cout << "Photons_hasPixelSeed, fake = " << Photons_hasPixelSeed->at(i) << endl;
      // if (Photons_hasPixelSeed->size() && dR >= 0.02) cout << "Photons_hasPixelSeed, real = " << Photons_hasPixelSeed->at(i) << endl;
    }
    if (Electrons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisElectron : *Electrons)
	dR = Photons->at(i).DeltaR(thisElectron) < dR ? Photons->at(i).DeltaR(thisElectron) : dR;
      h->Fill(dR, wt);
      // if (Photons_isEB->size() && dR < 0.02) cout << "isEB photon = " << Photons_isEB->at(i) << endl;
      // if (Photons_hasPixelSeed->size() && dR < 0.02) cout << "Photons_hasPixelSeed, fake = " << Photons_hasPixelSeed->at(i) << endl;
      // if (Photons_hasPixelSeed->size() && dR >= 0.02) cout << "Photons_hasPixelSeed, real = " << Photons_hasPixelSeed->at(i) << endl;
    }
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGLdRpixelSeed(TH1D* h, double wt) {
  for (size_t i = 0; i < Photons->size(); ++i) {
    if (!Photons_fullID->at(i)) continue;
    if (!Photons_hasPixelSeed->at(i)) continue;
    if (Muons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisMuon : *Muons)
	dR = Photons->at(i).DeltaR(thisMuon) < dR ? Photons->at(i).DeltaR(thisMuon) : dR;
      h->Fill(dR, wt);
    }
    if (Electrons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisElectron : *Electrons)
	dR = Photons->at(i).DeltaR(thisElectron) < dR ? Photons->at(i).DeltaR(thisElectron) : dR;
      h->Fill(dR, wt);
    }
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillCutMaps() {

  sampleKeyMap_["sig"] = "sig";
  sampleKeyMap_["T1bbbb"] = "sig";
  sampleKeyMap_["T1qqqq"] = "sig";
  sampleKeyMap_["T1tttt"] = "sig";
  sampleKeyMap_["zinv"] = "sig";
  sampleKeyMap_["topW"] = "sig";
  sampleKeyMap_["topWsle"] = "sle";
  sampleKeyMap_["topWslm"] = "slm";
  sampleKeyMap_["topsle"] = "sle";
  sampleKeyMap_["topslm"] = "slm";
  sampleKeyMap_["wjetssle"] = "sle";
  sampleKeyMap_["wjetsslm"] = "slm";
  sampleKeyMap_["qcd"] = "sig";
  sampleKeyMap_["qcdIDP"] = "sig";
  sampleKeyMap_["zinvIDP"] = "sig";
  sampleKeyMap_["wjets"] = "sig";
  sampleKeyMap_["other"] = "sig";
  sampleKeyMap_["ttzvv"] = "ttz";
  sampleKeyMap_["zmm"] = "zmm";
  sampleKeyMap_["zmm20"] = "zmm";
  sampleKeyMap_["dymm"] = "zmm";
  sampleKeyMap_["ttmm"] = "zmm";
  sampleKeyMap_["ttzmm"] = "zmm";
  sampleKeyMap_["VVmm"] = "zmm";
  sampleKeyMap_["tribosonmm"] = "zmm";
  sampleKeyMap_["zee"] = "zee";
  sampleKeyMap_["zee20"] = "zee";
  sampleKeyMap_["dyee"] = "zee";
  sampleKeyMap_["ttzee"] = "zee";
  sampleKeyMap_["ttee"] = "zee";
  sampleKeyMap_["VVee"] = "zee";
  sampleKeyMap_["tribosonee"] = "zee";
  sampleKeyMap_["zll"] = "zll";
  sampleKeyMap_["zll20"] = "zll";
  sampleKeyMap_["dyll"] = "zll";
  sampleKeyMap_["ttzll"] = "zll";
  sampleKeyMap_["ttll"] = "zll";
  sampleKeyMap_["VVll"] = "zll";
  sampleKeyMap_["photon"] = "photon";
  sampleKeyMap_["photon20"] = "photon";
  sampleKeyMap_["gjets"] = "photon";
  sampleKeyMap_["gjetsold"] = "photon";
  sampleKeyMap_["ttgjets"] = "photon";
  sampleKeyMap_["gjetsqcd"] = "photonqcd";

  if (isSkim_) {
    if (ntupleVersion_ == "V12") {
      objCutMap_["sig"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
      objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";  // Troy mod+
      // objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
      objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
      objCutMap_["zll"] = "((@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0) || (@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
      objCutMap_["ttz"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "@Muons.size()==1 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";
      objCutMap_["sle"] = "@Muons.size()==0 && @Electrons.size()==1 && isoMuonTracks==0 && isoPionTracks==0";

    } else if (ntupleVersion_ == "V15" || ntupleVersion_ == "V16") {

      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // zmm skim cuts:  NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0
      objCutMap_["zmm"] = "1";
      // zee skim cuts:  NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0
      objCutMap_["zee"] = "1";
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0) || (NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0))";
      // photon skim cuts:  "Sum$(Photons_fullID)==1 && Photons_hasPixelSeed==0 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0"
      objCutMap_["photon"] = "@Photons.size()==1";  // Andrew recommendation, no skim cuts
      if (isMC_) objCutMap_["photon"] += " && !Photons_nonPrompt && madMinPhotonDeltaR>=0.4";  // Andrew recommendation, no skim cuts
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "@Photons.size()==1";
      if (isMC_) objCutMap_["photonqcd"] += " && Photons_nonPrompt";
      objCutMap_["ttz"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "NMuons==1 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0";
      objCutMap_["sle"] = "NMuons==0 && NElectrons==1 && isoMuonTracks==0 && isoPionTracks==0";
    }

    minDphiCutMap_["nominal"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
    minDphiCutMap_["hdp"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
    minDphiCutMap_["ldp"] = "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3 || DeltaPhi4<0.3)";

    MHTCutMap_["nominal"] = "MHT>=300";
    MHTCutMap_["hdp"] = "MHT>=250";
    MHTCutMap_["ldp"] = "MHT>=250";

  } else {  // ntuple

    if (ntupleVersion_ == "V12") {
    } else if (ntupleVersion_ == "V15" || ntupleVersion_ == "V16") {
      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      // objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoMuonTracksclean==0";
      objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0";
      // objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoElectronTracksclean==0";
      objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0) || (NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1 && Photons_hasPixelSeed==0) && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["ttz"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "NMuons==1 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["slm"] += " && TransverseMass(METPt,METPhi,@Muons.at(0).Pt(),@Muons.at(0).Phi()) < 100";  // add'l skim cut
      objCutMap_["sle"] = "NMuons==0 && NElectrons==1 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["sle"] += " && TransverseMass(METPt,METPhi,@Electrons.at(0).Pt(),@Electrons.at(0).Phi()) < 100";  // add'l skim cut
    }

    minDphiCutMap_["nominal"] = "DeltaPhi1clean>0.5 && DeltaPhi2clean>0.5 && DeltaPhi3clean>0.3 && DeltaPhi4clean>0.3";
    minDphiCutMap_["hdp"] = "DeltaPhi1clean>0.5 && DeltaPhi2clean>0.5 && DeltaPhi3clean>0.3 && DeltaPhi4clean>0.3";
    minDphiCutMap_["ldp"] = "(DeltaPhi1clean<0.5 || DeltaPhi2clean<0.5 || DeltaPhi3clean<0.3 || DeltaPhi4clean<0.3)";

    MHTCutMap_["nominal"] = "MHTclean>=300";
    MHTCutMap_["hdp"] = "MHTclean>=250";
    MHTCutMap_["ldp"] = "MHTclean>=250";
  }

  triggerMapByName_["zmm"] = {
    "HLT_IsoMu20_v",  // Sam+
    "HLT_IsoMu22_v",  // Sam+
    "HLT_IsoMu24_v",  // prescaled in late 2017 --Owen
    "HLT_IsoMu27_v",
    "HLT_IsoMu22_eta2p1_v",  // Sam+
    "HLT_IsoMu24_eta2p1_v",  // Sam+
    "HLT_IsoTkMu22_v",  // Sam+
    "HLT_IsoTkMu24_v",
    "HLT_Mu15_IsoVVVL_PFHT350_v",
    "HLT_Mu15_IsoVVVL_PFHT400_v",
    "HLT_Mu15_IsoVVVL_PFHT450_v",  // Sam+
    "HLT_Mu15_IsoVVVL_PFHT600_v",  // Sam+
    "HLT_Mu50_IsoVVVL_PFHT400_v",  // Sam+
    "HLT_Mu50_IsoVVVL_PFHT450_v",  // Sam+
    "HLT_Mu50_v",
    "HLT_Mu55_v"  // Sam+
    // From AN-18-271
    // HLT_IsoMuX_v* (X=20,22,24,27);
    // HLT_IsoMuX_eta2p1_v* (X=22,24) ;
    // HLT_IsoTkMuX_v* (X=22,24) ;
    // HLT_Mu15_IsoVVVL_PFHTX_v* (X=350,400,450,600) ;
    // HLT_Mu50_IsoVVVL_PFHTX_v* (400,450) ;
    // HLT_MuX_v* (X=50,55) .
  };

  triggerMapByName_["zee"] = {
    "HLT_Ele105_CaloIdVT_GsfTrkIdT_v",
    "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
    "HLT_Ele135_CaloIdVT_GsfTrkIdT_v",  // Sam+
    "HLT_Ele145_CaloIdVT_GsfTrkIdT_v",  // Sam+
    "HLT_Ele15_IsoVVVL_PFHT350_v",
    "HLT_Ele15_IsoVVVL_PFHT400_v",
    "HLT_Ele15_IsoVVVL_PFHT450_v",  // Sam+
    "HLT_Ele15_IsoVVVL_PFHT600_v",  // Sam+
    "HLT_Ele27_WPTight_Gsf_v",
    "HLT_Ele35_WPTight_Gsf_v",  // Sam+
    "HLT_Ele20_WPLoose_Gsf_v",  // Sam+
    "HLT_Ele45_WPLoose_Gsf_v",  // Sam+
    "HLT_Ele25_eta2p1_WPTight_Gsf_v",  // Sam+
    "HLT_Ele20_eta2p1_WPLoose_Gsf_v",  // Sam+
    "HLT_Ele27_eta2p1_WPLoose_Gsf_v"
    // From AN-18-271
    // HLT_EleX_CaloIdVT_GsfTrkIdT_v*(X=105, 115, 135, 145);
    // HLT_Ele25_eta2p1_WPTight_Gsf_v*;
    // HLT_EleX_eta2p1_WPLoose_Gsf_v*(X=20, 27);
    // HLT_Ele15_IsoVVVL_PFHTX_v*(X=350, 400, 450, 600);
    // HLT_EleX_WPTight_Gsf_v*(X=27, 35); and
    // HLT_EleX_WPLoose_Gsf_v*(X=20,45).
  };
  triggerMapByName_["photon"] = {
    "HLT_Photon175_v",
    "HLT_Photon200_v"
  };
  triggerMapByName_["sig"] = {
    "HLT_PFMET100_PFMHT100_IDTight_v",
    "HLT_PFMET110_PFMHT110_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_v",
    "HLT_PFMET130_PFMHT130_IDTight_v",  // Sam+
    "HLT_PFMET140_PFMHT140_IDTight_v",  // Sam+
    "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v",
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v",  // Sam+
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",  // Sam+
    "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET110_PFMHT110_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET130_PFMHT130_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMET140_PFMHT140_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",  // Sam+
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_PFHT60_v"  // Sam+
    // From AN-18-271
    // HLT_PFMETX_PFMHTX_IDTight_v* (X=100,110,120,130,140).
    // HLT_PFMETNoMuX_PFMHTNoMuX_IDTight_v* (X=100,110,120,130,140)
    // HLT_PFMETX_PFMHTX_IDTight_PFHT60_v* (X=100,110,120,130,140)
    // HLT_PFMETNoMuX_PFMHTNoMuX_IDTight_PFHT60_v* (X=100,110,120,130)
  };
  triggerMapByName_["zll"].reserve(triggerMapByName_["zmm"].size() + triggerMapByName_["zee"].size());
  triggerMapByName_["zll"] = triggerMapByName_["zmm"];
  triggerMapByName_["zll"].insert(triggerMapByName_["zll"].end(), triggerMapByName_["zee"].begin(), triggerMapByName_["zee"].end());
  triggerMapByName_["sle"] = triggerMapByName_["sig"];
  triggerMapByName_["slm"] = triggerMapByName_["sig"];

}  // ======================================================================================

void
RA2bZinvAnalysis::setTriggerIndexList(const char* sample) {

  triggerIndexList_.clear();
  std::vector<TString> triggers;
  if (triggerMapByName_.count(sample) > 0) {
    triggers = triggerMapByName_.at(sample);
  } else {
    cout << "No matches in triggerMapByName for sample " << sample << endl;
    return;
  }
  for (auto myTrigName : triggers) {
    for (unsigned int ti = 0; ti < TriggerNames->size(); ++ti) {
      if (TString(TriggerNames->at(ti)).Contains(myTrigName)) {
	triggerIndexList_.push_back(ti);
	Int_t prescale = TriggerPrescales->at(ti);
	if (verbosity_ >= 2 || (verbosity_ >= 1 && prescale != 1))
	  cout << "Trigger " << TriggerNames->at(ti) << " (" << ti << ") prescaled by " << prescale << endl;
	continue;
      }
    }
  }

}  // ======================================================================================

Int_t
RA2bZinvAnalysis::setBTags() {
  Int_t BTagsOrig = BTags;
  if (useDeepCSV_) {
    BTags = BTagsDeepCSV;
  } else if (ntupleVersion_ == "V15" && runBlock_.find("2016") == std::string::npos) {
    // Recompute BTags with a different discriminator threshold
    BTags = 0;
    for (size_t j = 0; j < Jets->size(); ++j) {
      if(!Jets_HTMask->at(j)) continue;
      if (Jets_bDiscriminatorCSV->at(j) > csvMthreshold_) BTags++;
    }
  }
  // From Rishi email of 20 Feb 2019, DeepCSV 2018 WP: 0.4184 and the 2017 WP: 0.4941
  return BTagsOrig;
}  // ======================================================================================

RA2bZinvAnalysis::cutHistos::cutHistos(TChain* chain, TObjArray* forNotify) : forNotify_(forNotify) {
  HTcutf_ = new TTreeFormula("HTcut", HTcut_, chain);  forNotify->Add(HTcutf_);
  MHTcutf_ = new TTreeFormula("MHTcut", MHTcut_, chain);  forNotify->Add(MHTcutf_);
  NJetscutf_ = new TTreeFormula("NJetscut", NJetscut_, chain);  forNotify->Add(NJetscutf_);
  minDphicutf_ = new TTreeFormula("minDphicut", minDphicut_, chain);  forNotify->Add(minDphicutf_);
  objcutf_ = new TTreeFormula("objcut", objcut_, chain);  forNotify->Add(objcutf_);
  ptcutf_ = new TTreeFormula("ptcut", ptCut_, chain);  forNotify->Add(ptcutf_);
  masscutf_ = new TTreeFormula("masscut", massCut_, chain);  forNotify->Add(masscutf_);
  photonDeltaRcutf_ = new TTreeFormula("photonDeltaRcut", photonDeltaRcut_, chain);  forNotify->Add(photonDeltaRcutf_);
  commoncutf_ = new TTreeFormula("commoncut", commonCuts_, chain);  forNotify->Add(commoncutf_);
}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::setAxisLabels(TH1D* hcf) {
  if (TString(hcf->GetName()).Contains(TString("hFilterCuts"))) {
    hcf->GetXaxis()->SetBinLabel(1, "1 None");
    hcf->GetXaxis()->SetBinLabel(2, "2 TightHalo");
    hcf->GetXaxis()->SetBinLabel(3, "3 SupTgtHalo");
    hcf->GetXaxis()->SetBinLabel(4, "4 HBENoise");
    hcf->GetXaxis()->SetBinLabel(5, "5 HBEIsoNoise");
    hcf->GetXaxis()->SetBinLabel(6, "6 EcalDeadCell");
    hcf->GetXaxis()->SetBinLabel(7, "7 BadChCand");
    hcf->GetXaxis()->SetBinLabel(8, "8 BadPFMu");
    hcf->GetXaxis()->SetBinLabel(9, "9 NVtx");
    hcf->GetXaxis()->SetBinLabel(10, "10 eeBadSc");
    hcf->GetXaxis()->SetBinLabel(11, "11 ecalBadCal");
    hcf->GetXaxis()->SetBinLabel(12, "12");
    hcf->GetXaxis()->SetBinLabel(13, "13 JetID");
    hcf->GetXaxis()->SetBinLabel(14, "14 METRatio");
    hcf->GetXaxis()->SetBinLabel(15, "15 HTdPhiRatio");
    hcf->GetXaxis()->SetBinLabel(16, "16 EcalNoiseJet");
    hcf->GetXaxis()->LabelsOption("vu");
  }
  if (TString(hcf->GetName()).Contains(TString("hCut"))) {
    hcf->GetXaxis()->SetBinLabel(1, "1 None");
    hcf->GetXaxis()->SetBinLabel(2, "2 HT");
    hcf->GetXaxis()->SetBinLabel(3, "3 MHT");
    hcf->GetXaxis()->SetBinLabel(4, "4 NJets");
    hcf->GetXaxis()->SetBinLabel(5, "5 mnDphi");
    hcf->GetXaxis()->SetBinLabel(6, "6 objects");
    hcf->GetXaxis()->SetBinLabel(7, "7 Z/gamma pt");
    hcf->GetXaxis()->SetBinLabel(8, "8 m(Z)");
    hcf->GetXaxis()->SetBinLabel(9, "9 gamma dR");
    hcf->GetXaxis()->SetBinLabel(10, "10 Trigger");
    hcf->GetXaxis()->SetBinLabel(11, "11 HEM");
    hcf->GetXaxis()->SetBinLabel(12, "12 Filters");
    hcf->GetXaxis()->LabelsOption("vu");
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::fill(TH1D* hcf, Double_t wt, bool passTrg, bool passHEM) {
  HTcutf_->GetNdata();
  MHTcutf_->GetNdata();
  NJetscutf_->GetNdata();
  minDphicutf_->GetNdata();
  objcutf_->GetNdata();
  ptcutf_->GetNdata();
  masscutf_->GetNdata();
  photonDeltaRcutf_->GetNdata();
  commoncutf_->GetNdata();

  hcf->Fill(0.5, wt);
  if (TString(hcf->GetName()).Contains(TString("Flow"))) {
    if (HTcutf_->EvalInstance(0)) {
      hcf->Fill(1.5, wt);
      if (MHTcutf_->EvalInstance(0)) {
	hcf->Fill(2.5, wt);
	if (NJetscutf_->EvalInstance(0)) {
	  hcf->Fill(3.5, wt);
	  if (minDphicutf_->EvalInstance(0)) {
	    hcf->Fill(4.5, wt);
	    if (objcutf_->EvalInstance(0)) {
	      hcf->Fill(5.5, wt);
	      if (ptcutf_->EvalInstance(0)) {
		hcf->Fill(6.5, wt);
		if (masscutf_->EvalInstance(0)) {
		  hcf->Fill(7.5, wt);
		  if (photonDeltaRcutf_->EvalInstance(0)) {
		    hcf->Fill(8.5, wt);
		    if (passTrg) {
		      hcf->Fill(9.5, wt);
		      if (passHEM) {
			hcf->Fill(10.5, wt);
			if (commoncutf_->EvalInstance(0)) {
			  hcf->Fill(11.5, wt);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  } else {
    if (HTcutf_->EvalInstance(0)) hcf->Fill(1.5, wt);
    if (MHTcutf_->EvalInstance(0)) hcf->Fill(2.5, wt);
    if (NJetscutf_->EvalInstance(0)) hcf->Fill(3.5, wt);
    if (minDphicutf_->EvalInstance(0)) hcf->Fill(4.5, wt);
    if (objcutf_->EvalInstance(0)) hcf->Fill(5.5, wt);
    if (ptcutf_->EvalInstance(0)) hcf->Fill(6.5, wt);
    if (masscutf_->EvalInstance(0)) hcf->Fill(7.5, wt);
    if (photonDeltaRcutf_->EvalInstance(0)) hcf->Fill(8.5, wt);
    if (passTrg) hcf->Fill(9.5, wt);
    if (passHEM) hcf->Fill(10.5, wt);
    if (commoncutf_->EvalInstance(0)) hcf->Fill(11.5, wt);
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::efficiencyAndPurity::openFiles() {
  TString plotDir("../plots/histograms/");

  prefiringWeightFile_ = new TFile((plotDir+"L1PrefiringMaps_new.root").Data(), "read");

  purityTrigEffFile_.push_back(new TFile((plotDir+"effHists.root").Data(), "read"));
  purityTrigEffFile_.push_back(new TFile((plotDir+"effHists.root").Data(), "read"));
  purityTrigEffFile_.push_back(new TFile((plotDir+"effHists.root").Data(), "read"));

  photonTrigEffFile_.push_back(new TFile((plotDir+"triggersRa2bRun2_v3_withTEffs.root").Data(), "read"));
  photonTrigEffFile_.push_back(new TFile((plotDir+"triggersRa2bRun2_v3_withTEffs.root").Data(), "read"));
  photonTrigEffFile_.push_back(new TFile((plotDir+"triggersRa2bRun2_v3_withTEffs.root").Data(), "read"));

  photonSFFile_.push_back(new TFile((plotDir+"Photons2016_SF_all.root").Data(), "read"));
  photonSFFile_.push_back(new TFile((plotDir+"Photons2017_SF_all.root").Data(), "read"));
  photonSFFile_.push_back(new TFile((plotDir+"Photons2018_SF_all.root").Data(), "read"));

  // elecSFFile_.push_back(new TFile((plotDir+"ElectronScaleFactors_Run2017.root").Data(), "read"));
  // elecSFFile_.push_back(new TFile((plotDir+"ElectronScaleFactors_Run2017.root").Data(), "read"));
  // elecSFFile_.push_back(new TFile((plotDir+"ElectronScaleFactors_Run2017.root").Data(), "read"));

  elecIDandIsoSFFile_.push_back(new TFile((plotDir+"Electrons2016_SF_ISOandID.root").Data(), "read"));
  elecIDandIsoSFFile_.push_back(new TFile((plotDir+"Electrons2017_SF_ISOandID.root").Data(), "read"));
  elecIDandIsoSFFile_.push_back(new TFile((plotDir+"Electrons2018_SF_ISOandID.root").Data(), "read"));

  elecRecoLowSFFile_.push_back(new TFile((plotDir+"Electrons2016_SF_RECOlow.root").Data(), "read"));
  elecRecoLowSFFile_.push_back(new TFile((plotDir+"Electrons2017_SF_RECOlow.root").Data(), "read"));
  //elecRecoLowSFFile_.push_back(new TFile((plotDir+"Electrons2018_SF_RECOlow.root").Data(), "read"));//not available yet for 2018

  elecRecoHighSFFile_.push_back(new TFile((plotDir+"Electrons2016_SF_RECOhigh.root").Data(), "read"));//pt>20
  elecRecoHighSFFile_.push_back(new TFile((plotDir+"Electrons2017_SF_RECOhigh.root").Data(), "read"));//pt>20
  elecRecoHighSFFile_.push_back(new TFile((plotDir+"Electrons2018_SF_RECOhigh.root").Data(), "read"));//pt>10

  muonIDSFFile_.push_back(new TFile((plotDir+"Muons2016_SF_ID.root").Data(), "read"));
  muonIDSFFile_.push_back(new TFile((plotDir+"Muons2017_SF_ID.root").Data(), "read"));
  muonIDSFFile_.push_back(new TFile((plotDir+"Muons2018_SF_ID.root").Data(), "read"));

  muonIsoSFFile_.push_back(new TFile((plotDir+"Muons2016_SF_ISO.root").Data(), "read"));
  muonIsoSFFile_.push_back(new TFile((plotDir+"Muons2017_SF_ISO.root").Data(), "read"));
  muonIsoSFFile_.push_back(new TFile((plotDir+"Muons2018_SF_ISO.root").Data(), "read"));

  DRfun_ = new TF1("DRfun", "[0] + [1]*min([3], x)");
  // DRfun_ = new TF1("DRfun", "[2] + (1/3)*[1]*(min([3], x) - ([2] - [0]) / [1])");
  DRpars_.push_back({0.8229, 0.0001665, 0.9061, 900});
  DRpars_.push_back({0.8229, 0.0001665, 0.9061, 900});
  DRpars_.push_back({0.8229, 0.0001665, 0.9061, 900});
  // Graph_from_hHT_DR_zmm  Run 2 no PU
  // Function parameter 0:  0.822859680122 +/- 0.0188240089302
  // Function parameter 1:  0.000166574459092 +/- 3.55474346141e-05
  // DRpars_.push_back({0.8378, 0.0001363, 0.9054, 900});
  // Graph_from_hHT_DR_zmm  2016 noPU
  // Function parameter 0:  0.837859796025 +/- 0.0352792948296
  // Function parameter 1:  0.000136284149882 +/- 6.6930863488e-05
  // DRpars_.push_back({0.8430, 0.0001147, 0.9000, 900});
  // Graph_from_hHT_DR_zmm  2016 first iteration
  // Function parameter 0:  0.843020160316 +/- 0.0348336458629
  // Function parameter 1:  0.000114662993277 +/- 6.57638091784e-05
  // Average y0 = 0.900.  x0 = (y0 -p0) / p1
  // DRpars_.push_back({0.8246, 0.0001274, 0.8877, 900});
  // DRpars_.push_back({0.8246, 0.0001274, 0.8877, 900});
  // Graph_from_hHT_DR_zmm  2017 HT17wt ZptWt
  // Function parameter 0:  0.824632395818 +/- 0.0327675716282
  // Function parameter 1:  0.000127461473018 +/- 6.25698022129e-05
  // DRpars_.push_back({0.9611, 0.0002574, 1.0871, 900});
  // DRpars_.push_back({0.9611, 0.0002574, 1.0871, 900});
  // effWt /= (min(HT, 900.0) - 489.2)*(0.00025739) + 1.0858;
  // Graph_from_hHT_DR_zmm    // 2017 HT16Wt ZptWt
  // Function parameter 0:  0.961193586511 +/- 0.0401249241206
  // Function parameter 1:  0.000257387819072 +/- 7.74713707439e-05
  // Average y0 = 1.0871.  x0 = (y0 -p0) / p1

  // effWt /= (min(HT, 900.0) - 497.4)*(0.0002288) + 1.0395;  // Run2 Z Pt weighted
  // Graph_from_hHT_DR_zmm
  // Function parameter 0:  0.92567079283 +/- 0.021512658121
  // Function parameter 1:  0.000228788345936 +/- 4.08070918994e-05
  // Average y0 = 1.0395.  x0 = (y0 -p0) / p1
  // effWt /= (min(MHT, 900.0) - 399.6)*(  -0.00040321 *1/3) + 0.8389;  // New Fdir, 28 Jan, 2019
  // Graph_from_hMHT_DR_zmm
  // Function parameter 0:  0.999926406018 +/- 0.0291554520092
  // Function parameter 1:  -0.000403207439065 +/- 7.19133297508e-05
  // Average y0 = 0.8389.  x0 = (y0 -p0) / p1
  // effWt /= (min(MHT, 900.0) - 400.2)*( -0.00039455 *2/3) + 0.8401;  // gJets_signal.dat 28 Jan, 2019
  // Function parameter 0:  0.997861381979 +/- 0.0292274226244
  // Function parameter 1:  -0.000394553221914 +/- 7.21321517572e-05
  // Average y0 = 0.8401.  x0 = (y0 -p0) / p1
  // effWt /= (min(MHT, 900.0) - 397.7)*( -0.00058276 *2/3) + 0.8401;
  // Function parameter 0:  1.07188729834 +/- 0.// 047180710618
  // Function parameter 1:  -0.000582762306737 +/- 0.000115148917583
  // Average y0 = 0.8401.  x0 = (y0 -p0) / p1
  // effWt /= (min(MHT, 900.0) - 397.4)*( -0.0005098 *2/3) + 0.8134;  // 16 Jan 2019
  // Function parameter 0:  1.01601616809 +/- 0.0459169526921
  // Function parameter 1:  -0.00050983094812 +/- 0.000112511842395
  // Average y0 = 0.8134.  x0 = (y0 -p0) / p1
  // effWt /= min(MHT, 900.0)*( -0.0005284) + 1.056;  // FIXME:  hard-wired correction

}  // ======================================================================================

void
RA2bZinvAnalysis::efficiencyAndPurity::getHistos(const char* sample, int currentYear) {

  std::vector<const char*> year;  year.push_back("2016");  year.push_back("2017");  year.push_back("2018");

  // For L1 prefiring weight
  if (currentYear == Year2016) {
    hPrefiring_photon_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_photonptvseta_2016BtoH");
    if (hPrefiring_photon_ == nullptr) cout << "***** Histogram L1prefiring_photonptvseta_2016BtoH not found *****" << endl;
    hPrefiring_jet_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_jetptvseta_2016BtoH");
    if (hPrefiring_jet_ == nullptr) cout << "***** Histogram L1prefiring_jetptvseta_2016BtoH not found *****" << endl;
  } else if (currentYear == Year2017) {
    hPrefiring_photon_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_photonptvseta_2017BtoF");
    if (hPrefiring_photon_ == nullptr) cout << "***** Histogram L1prefiring_photonptvseta_2017BtoF not found *****" << endl;
    hPrefiring_jet_ = (TH2F*) prefiringWeightFile_->Get("L1prefiring_jetptvseta_2017BtoF");
    if (hPrefiring_jet_ == nullptr) cout << "***** Histogram L1prefiring_jetptvseta_2017BtoF not found *****" << endl;
  }

  // For purity, Fdir, trigger eff, reco eff
  theSample_ = TString(sample);
  // hSFeff_ = nullptr;
  hSFeff_.clear();
  FdirGraph_ = nullptr;
  hPurity_.clear();
  hTrigEff_.clear();
  if (theSample_.Contains("zmm") && !theSample_.Contains("tt")) {
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_pur_m"));
    if (hPurity_.back() == nullptr) cout << "***** Histogram h_pur_m not found *****" << endl;
  } else if (theSample_.Contains("zee") && !theSample_.Contains("tt")) {
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_pur_e"));
    if (hPurity_.back() == nullptr) cout << "***** Histogram h_pur_e not found *****" << endl;
  } else if (theSample_.Contains("photon")) {
    TString hn("h_pur_eb_"); hn += year.at(currentYear);
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get(hn));
    if (hPurity_.back() == nullptr) cout << "***** Histogram " << hn << " not found *****" << endl;
    hn("eb") = "ec";
    hPurity_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get(hn));
    if (hPurity_.back() == nullptr) cout << "***** Histogram " << hn << " not found *****" << endl;
    // FdirHist_ = (TH1D*) purityTrigEffFile_.at(currentYear)->Get("h_bin46_NJets8910");
    // if (FdirHist_ == nullptr) cout << "***** Histogram h_bin46_NJets8910 not found *****" << endl;
    FdirGraph_ = (TGraphErrors*) purityTrigEffFile_.at(currentYear)->Get("h_bin46_NJets8910");
    if (FdirGraph_ == nullptr) cout << "***** Histogram h_bin46_NJets8910 not found *****" << endl;
  } else if (theSample_.Contains("dymm") || theSample_.Contains("ttmm") ||
	     theSample_.Contains("ttzmm") || theSample_.Contains("VVmm")) {
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_m"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_m not found *****" << endl;
    // hSFeff_ = (TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_SFm_MHT");
    // if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;

    //These change based on year
    if (currentYear == Year2016) {
      hSFeff_.push_back((TH2F*) muonIDSFFile_.at(currentYear)->Get("NUM_MediumID_DEN_genTracks_eta_pt"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for muon ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) muonIsoSFFile_.at(currentYear)->Get("NUM_TightRelIso_DEN_MediumID_eta_pt"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for muon iso SFs not found *****" << endl;
    }
    else {
      hSFeff_.push_back((TH2F*) muonIDSFFile_.at(currentYear)->Get("NUM_MediumID_DEN_genTracks_pt_abseta"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for muon ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) muonIsoSFFile_.at(currentYear)->Get("NUM_TightRelIso_DEN_MediumID_pt_abseta"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for muon iso SFs not found *****" << endl;
    }

  } else if (theSample_.Contains("dyee") || theSample_.Contains("ttee") ||
	     theSample_.Contains("ttzee") || theSample_.Contains("VVee")) {
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_e"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_e not found *****" << endl;
    // hSFeff_ = (TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_SFe_MHT");  // Maybe this should be h_NJets
    // if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;

    //These change based on year
    if (currentYear == Year2016) {
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2016_CutBasedVetoNoIso94XV2"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2016_Mini"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron iso SFs not found *****" << endl;
    }
    else if (currentYear == Year2017) {
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2017_CutBasedVetoNoIso94XV2"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2017_MVAVLooseTightIP2DMini"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron iso SFs not found *****" << endl;
    }
    else if (currentYear == Year2018) {
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2018_CutBasedVetoNoIso94XV2"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron ID SFs not found *****" << endl;
      hSFeff_.push_back((TH2F*) elecIDandIsoSFFile_.at(currentYear)->Get("Run2018_Mini"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron iso SFs not found *****" << endl;
    }
    hSFeff_.push_back((TH2F*) elecRecoHighSFFile_.at(currentYear)->Get("EGamma_SF2D"));
    if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron reco (high pT) SFs not found *****" << endl;

    //Electron reco (pt<10 GeV) SFs for 2018 not (yet) available
    if (currentYear != Year2018) {
      hSFeff_.push_back((TH2F*) elecRecoLowSFFile_.at(currentYear)->Get("EGamma_SF2D"));
      if (hSFeff_.back() == nullptr) cout << "***** Histogram for electron reco (low pT) SFs not found *****" << endl;
    }


  } else if (theSample_.Contains("gjets")) {
    // // Troy-style trigger efficiency histograms:
    // hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_eb"));
    // if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_eb not found *****" << endl;
    // hTrigEff_.push_back((TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_trig_ec"));
    // if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_ec not found *****" << endl;
    // // Sam-style trigger efficiency functions:
    // fTrigEff_.push_back((TF1*) purityTrigEffFile_.at(currentYear)->Get("f_trig_eb"));
    // if (fTrigEff_.back() == nullptr) {
    //   cout << "***** Histogram f_trig_eb not found *****" << endl;
    // } else {
    //   cout << "Asymptotic EB trigger efficiency = " << 1 + fTrigEff_[0]->GetParameter(0)
    // 	   << "+/-" << fTrigEff_[0]->GetParError(0) << endl;
    // }
    // fTrigEff_.push_back((TF1*) purityTrigEffFile_.at(currentYear)->Get("f_trig_ec"));
    // if (fTrigEff_.back() == nullptr) {
    //   cout << "***** Histogram f_trig_ec not found *****" << endl;
    // } else {
    //   cout << "Asymptotic EC trigger efficiency = " << 1 + fTrigEff_[1]->GetParameter(0)
    // 	   << "+/-" << fTrigEff_[1]->GetParError(0) << endl;
    // }
    // Sam-style trigger efficiency TEfficiency's
    TString te("teff_SinglePhotonBarrelLooseHdp_hists_Run");  te += year.at(currentYear);  te += "_JetHT";
    eTrigEff_.push_back((TEfficiency*) photonTrigEffFile_.at(currentYear)->Get(te.Data()));
    if (eTrigEff_.back() == nullptr)  cout << "***** Histogram " << te << " not found *****" << endl;
    te("Barrel") = "Endcap";  if (theSample_.Contains("ldp")) te("Hdp") = "Ldp";
    eTrigEff_.push_back((TEfficiency*) photonTrigEffFile_.at(currentYear)->Get(te.Data()));
    if (eTrigEff_.back() == nullptr)  cout << "***** Histogram " << te << " not found *****" << endl;
    // hSFeff_ = (TH1F*) purityTrigEffFile_.at(currentYear)->Get("h_SFg_MHT");  // Maybe this should be h_NJets
    // if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;

    hSFeff_.push_back((TH2F*) photonSFFile_.at(currentYear)->Get("EGamma_SF2D"));
    if (hSFeff_.back() == nullptr) cout << "***** Histogram for photon SFs not found *****" << endl;

    DRfun_->SetParameters(DRpars_.at(currentYear).data());
    cout << "DR intercept = " << DRpars_.at(currentYear)[0] << ", slope = " << DRpars_.at(currentYear)[1]
	 << ", Ave. value = " << DRpars_.at(currentYear)[2] << ", saturation point = " << DRpars_.at(currentYear)[3] << endl;
  }

}  // ======================================================================================

pair<double, double>
RA2bZinvAnalysis::efficiencyAndPurity::weight(CCbinning* CCbins, Int_t NJets, Int_t BTags,
					      Double_t MHT, Double_t HT,
					      vector<TLorentzVector> ZCandidates,
					      vector<TLorentzVector> Photons,
					      vector<TLorentzVector> Electrons,
					      vector<TLorentzVector> Muons,
					      vector<double> EBphoton,
					      bool applyDRfitWt,
                int currentYear) {
  // For double ratio, apply weights for purity, Fdir, trigger eff, reco eff.
  // Systematic error is relative.
  double effWt = 1, effSys = 0;

  if ((theSample_.Contains("zmm") || theSample_.Contains("zee")) && !theSample_.Contains("tt")) {
    int binCCjb = CCbins->jb(CCbins->jbin(NJets), CCbins->bbin(NJets, BTags));
    if (hPurity_[0] != nullptr) {
      double eff = hPurity_[0]->GetBinContent(binCCjb);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[0]->GetBinError(binCCjb) / eff : 0;
    }

  } else if (theSample_.Contains("photon")) {
    if(EBphoton.at(0) == 1 && hPurity_[0] != nullptr) {
      int bin = hPurity_[0]->GetNbinsX();  while (MHT < hPurity_[0]->GetBinLowEdge(bin)) bin--;
      double eff = hPurity_[0]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[0]->GetBinError(bin) / eff : 0;
    }
    if(EBphoton.at(0) == 0 && hPurity_[1] != nullptr) {
      int bin = hPurity_[1]->GetNbinsX();  while (MHT < hPurity_[1]->GetBinLowEdge(bin)) bin--;
      double eff = hPurity_[1]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hPurity_[1]->GetBinError(bin) / eff : 0;
    }
    int CCbinjk = CCbins->jk(CCbins->jbin(NJets), CCbins->kinBin(HT, MHT));
    // if (CCbinjk > 0 && FdirHist_ != nullptr && CCbinjk <= FdirHist_->GetNbinsX()) effWt *= FdirHist_->GetBinContent(CCbinjk);
    if (CCbinjk > 0 && FdirGraph_ != nullptr && CCbinjk < FdirGraph_->GetN()) {
      double eff = FdirGraph_->GetY()[CCbinjk-1];
      effWt *= eff;
      effSys = eff > 0 ? quadSum(effSys, FdirGraph_->GetErrorY(CCbinjk-1)) / eff : 0;
    }
    // FIXME:  For LDP, fill extra bins with value 0.825

  } else if (theSample_.Contains("dymm") || theSample_.Contains("dyee") ||
	     theSample_.Contains("ttmm") || theSample_.Contains("ttee") ||
	     theSample_.Contains("ttzmm") || theSample_.Contains("ttzee") ||
	     theSample_.Contains("VVmm") || theSample_.Contains("VVee")) {
    if (hTrigEff_[0] != nullptr) {
      Double_t zpt = ZCandidates.at(0).Pt();  zpt = max(hTrigEff_[0]->GetBinLowEdge(1), zpt);
      int bin = hTrigEff_[0]->GetNbinsX();  while (zpt < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
      double eff = hTrigEff_[0]->GetBinContent(bin);
      effWt *= eff;
      effSys = eff > 0 ? hTrigEff_[0]->GetBinError(bin) / eff : 0;
    }
    if (hSFeff_[0] != nullptr && hSFeff_[1] != nullptr) {
      float pt = 0; float eta = 0;
      int globalbin_ID = 0; int globalbin_ISO = 0;
      
      double sysCorr = 0;  // Errors are correlated between Z daughter leptons
      if (theSample_.Contains("ee") && hSFeff_[2]!=nullptr) {
	int globalbin_RECO = 0;
	int numElectrons = Electrons.size();
	for (int i=0; i<numElectrons; i++){
	  pt  = Electrons.at(i).Pt(); if (pt>500) pt=499.9;
	  eta = Electrons.at(i).Eta();
	  globalbin_ID  = hSFeff_[0]->FindBin(eta,pt); globalbin_ISO = hSFeff_[1]->FindBin(eta,pt);
	  double effID = hSFeff_[0]->GetBinContent(globalbin_ID);
	  double effISO = hSFeff_[1]->GetBinContent(globalbin_ISO);
	  effWt *= effID * effISO;
	  double sysOne = quadSum(effID > 0 ? hSFeff_[0]->GetBinError(globalbin_ID)/effID : 0,
				  effISO > 0 ? hSFeff_[1]->GetBinError(globalbin_ISO)/effISO : 0);
	  if (currentYear == Year2018) { //2018 only has reco SFs for pt>10
	    if (pt<=10.0) pt = 10.1;
	    globalbin_RECO = hSFeff_[2]->FindBin(eta,pt);
	    double eff = hSFeff_[2]->GetBinContent(globalbin_RECO);
	    effWt *= eff;
	    sysOne = eff > 0 ? quadSum(sysOne, hSFeff_[2]->GetBinError(globalbin_RECO)/eff) : sysOne;
	  }
	  else { //2016 and 2017 have reco SFs for 0-20 GeV and >20 GeV
	    if (pt>20) {
	      globalbin_RECO = hSFeff_[2]->FindBin(eta,pt);
	      double eff = hSFeff_[2]->GetBinContent(globalbin_RECO);
	      effWt *= eff;
	      sysOne = eff > 0 ? quadSum(sysOne, hSFeff_[2]->GetBinError(globalbin_RECO)/eff) : sysOne;
	    }
	    else {
	      globalbin_RECO = hSFeff_[3]->FindBin(eta,pt);
	      double eff = hSFeff_[3]->GetBinContent(globalbin_RECO);
	      effWt *= eff;
	      sysOne = eff > 0 ? quadSum(sysOne, hSFeff_[3]->GetBinError(globalbin_RECO)/eff) : 0;
	    }
	  } //2016 or 2017
	  sysCorr += sysOne;
	} //loop over all electrons (should be two)
      }

      else if (theSample_.Contains("mm")) {
	int numMuons = Muons.size();
	for (int i=0; i<numMuons; i++){
	  pt = Muons.at(i).Pt(); if (pt>120) pt=119.9;
	  if (currentYear == Year2016) {
	    eta = Muons.at(i).Eta();
	    globalbin_ID = hSFeff_[0]->FindBin(eta,pt); globalbin_ISO = hSFeff_[1]->FindBin(eta,pt);
	  }
	  else { //2017 and 2018 have abs(eta), and the x- and y-axis are swapped compared to 2016
	    eta = abs(Muons.at(i).Eta());
	    globalbin_ID = hSFeff_[0]->FindBin(pt,eta); globalbin_ISO = hSFeff_[1]->FindBin(pt,eta);
	  }
	  double effSF = hSFeff_[0]->GetBinContent(globalbin_ID);
	  double effISO = hSFeff_[1]->GetBinContent(globalbin_ISO);
	  effWt *= effSF * effISO;
	  sysCorr += quadSum(effSF > 0 ? hSFeff_[0]->GetBinError(globalbin_ID)/effSF : 0,
			     effISO > 0 ? hSFeff_[1]->GetBinError(globalbin_ISO)/effISO : 0);
	}
      }
      effSys = quadSum(effSys, sysCorr);

    }

  } else if (theSample_.Contains("gjets")) {
    // if(EBphoton.at(0) == 1 && hTrigEff_[0] != nullptr) {
    //   int bin = hTrigEff_[0]->GetNbinsX();  while (MHT < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
    //   effWt *= hTrigEff_[0]->GetBinContent(bin);
    // }
    // if(EBphoton.at(0) == 0 && hTrigEff_[1] != nullptr) {
    //   int bin = hTrigEff_[1]->GetNbinsX();  while (MHT < hTrigEff_[1]->GetBinLowEdge(bin)) bin--;
    //   effWt *= hTrigEff_[1]->GetBinContent(bin);
    // }
    // if(EBphoton.at(0) == 1 && fTrigEff_[0] != nullptr) {
    //   effWt *= fTrigEff_[0]->Eval(max(double(205), Photons.at(0).Pt()));  // hard-wired cutoff
    // }
    // if(EBphoton.at(0) == 0 && fTrigEff_[1] != nullptr) {
    //   effWt *= fTrigEff_[1]->Eval(max(double(205), Photons.at(0).Pt()));  // hard-wired cutoff
    // }
    if(EBphoton.at(0) == 1 && eTrigEff_[0] != nullptr) {
      TH1F* htot = (TH1F*) eTrigEff_[0]->GetTotalHistogram();
      Int_t bin = min(htot->GetNbinsX(), htot->FindBin(Photons.at(0).Pt()));
      double eff = eTrigEff_[0]->GetEfficiency(bin);
      effWt *= eff;
      effSys = eff > 0 ? eTrigEff_[0]->GetEfficiencyErrorLow(bin) / eff : 0;
    } else if(EBphoton.at(0) == 0 && eTrigEff_[1] != nullptr) {
      TH1F* htot = (TH1F*) eTrigEff_[1]->GetTotalHistogram();
      Int_t bin = min(htot->GetNbinsX(), htot->FindBin(Photons.at(0).Pt()));
      double eff = eTrigEff_[1]->GetEfficiency(bin);
      effWt *= eff;
      effSys = eff > 0 ? eTrigEff_[1]->GetEfficiencyErrorLow(bin) / eff : 0;
    }

    if (applyDRfitWt) effWt /= DRfun_->Eval(HT);

    if (hSFeff_[0] != nullptr) {
      float photon_pt = 0; float photon_eta = 0;
      int globalbin_photon = 0;
      int numPhotons = Photons.size();
      for (int i=0; i<numPhotons; i++){
	photon_pt = Photons.at(i).Pt(); photon_eta = Photons.at(i).Eta();
	if (photon_pt>500) photon_pt=499.9;
	globalbin_photon  = hSFeff_[0]->FindBin(photon_eta, photon_pt);
	double eff = hSFeff_[0]->GetBinContent(globalbin_photon);
	effWt *= eff;
	effSys = eff > 0 ? quadSum(effSys, hSFeff_[0]->GetBinError(globalbin_photon)/eff) : effSys;
      }
    }
  }
  return pair<double, double>(effWt, effSys);

}  // ======================================================================================

void
RA2bZinvAnalysis::runMakeClass(const std::string& sample) {
  TChain* chain = getChain(sample.data(), nullptr, true);
  TString templateName("TreeMkrTemplate_");
  if (!isSkim_) templateName += "unskimmed_";
  templateName += isMC_ ? "MC_" : "data_";
  templateName += ntupleVersion_;
  chain->MakeClass(templateName.Data());

}  // ======================================================================================

void
RA2bZinvAnalysis::checkTrigPrescales(const char* sample) {
  Int_t fCurrent;  // current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);
  Long64_t Nentries = chain->GetEntries();
  Int_t countInFile = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) cout << "Current file in chain: " << thisFile->GetName() << endl;
      chain->GetEntry(entry);
      for (unsigned int ti=0; ti<TriggerNames->size(); ++ti) {
	Int_t prescale = TriggerPrescales->at(ti);
	// if (prescale != 1)
	  cout << "Trigger " << TriggerNames->at(ti) << " (" << ti << ") prescaled by " << prescale << endl;
      }
    }
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::dumpSelEvIDs(const char* sample, const char* idFileName) {
  //
  // For debugging, where individual events need to be selected
  //
  FILE* idFile;
  idFile = fopen(idFileName, "w");

  Int_t fCurrent;  // current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);

  TObjArray* forNotify = new TObjArray;
  forNotify->SetOwner();  // so that TreeFormulas will be deleted

  TCut baselineCuts = getCuts(sample);
  if (verbosity_ >= 1) cout << endl << "baseline = " << endl << baselineCuts << endl << endl;
  TTreeFormula* baselineFormula = new TTreeFormula("baselineCuts", baselineCuts, chain);
  forNotify->Add(baselineFormula);
  TTreeFormula* HTcutf_ = new TTreeFormula("HTcut", HTcut_, chain);  forNotify->Add(HTcutf_);
  TTreeFormula* MHTcutf_ = new TTreeFormula("MHTcut", MHTcut_, chain);  forNotify->Add(MHTcutf_);
  TTreeFormula* NJetscutf_ = new TTreeFormula("NJetscut", NJetscut_, chain);  forNotify->Add(NJetscutf_);
  TTreeFormula* minDphicutf_ = new TTreeFormula("minDphicut", minDphicut_, chain);  forNotify->Add(minDphicutf_);
  TTreeFormula* objcutf_ = new TTreeFormula("objcut", objcut_, chain);  forNotify->Add(objcutf_);
  TTreeFormula* ptcutf_ = new TTreeFormula("ptcut", ptCut_, chain);  forNotify->Add(ptcutf_);
  TTreeFormula* masscutf_ = new TTreeFormula("masscut", massCut_, chain);  forNotify->Add(masscutf_);
  TTreeFormula* photonDeltaRcutf_ = new TTreeFormula("photonDeltaRcut", photonDeltaRcut_, chain);  forNotify->Add(photonDeltaRcutf_);
  TTreeFormula* commoncutf_ = new TTreeFormula("commoncut", commonCuts_, chain);  forNotify->Add(commoncutf_);

  chain->SetNotify(forNotify);

  std::vector<ULong64_t> EvtNums = {
    // 721713515,
    // 276225503,
    // 243443474,
    // 210695036,
    // 644860046,
    // 245742780,
    // 265782210,
    // 220373684,
    // 133598095,
    //  12214825,
    // 760192813,
    // 265408355,
    // 743702352
      703085123,
      872764606,
      743702352,
     1106837339,
      528776780,
     1295400923,
     4021888665,
      116551966
  };

  Long64_t Nentries = chain->GetEntries();
  if (verbosity_ >= 1) cout << "Nentries in tree = " << Nentries << endl;
  int count = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (verbosity_ >= 1 && count % 100000 == 0) cout << "Entry number " << count << endl;

    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) {
	if (verbosity_ >= 1) cout << "Current file in chain: " << thisFile->GetName() << endl;
      }
    }
    chain->GetEntry(entry);

    if (std::find(EvtNums.begin(), EvtNums.end(), EvtNum) == EvtNums.end()) continue;
    printf("%15u %15u %15llu\n", RunNum, LumiBlockNum, EvtNum);
    fprintf(idFile, "%15u %15u %15llu\n", RunNum, LumiBlockNum, EvtNum);

    cleanVars();  // If unskimmed input, copy <var>clean to <var>
    Int_t BTagsOrig = setBTags();

    baselineFormula->GetNdata();  cout << "baseline = " << baselineFormula->EvalInstance(0) << endl;
    HTcutf_->GetNdata();  cout << "HTcut = " << HTcutf_->EvalInstance(0) << endl;
    MHTcutf_->GetNdata();  cout << "MHTcut = " << MHTcutf_->EvalInstance(0) << endl;
    NJetscutf_->GetNdata();  cout << "NJetscut = " << NJetscutf_->EvalInstance(0) << endl;
    minDphicutf_->GetNdata();  cout << "minDphicut = " << minDphicutf_->EvalInstance(0) << endl;
    objcutf_->GetNdata();  cout << "objcut = " << objcutf_->EvalInstance(0) << endl;
    ptcutf_->GetNdata();  cout << "ptcut = " << ptcutf_->EvalInstance(0) << endl;
    masscutf_->GetNdata();  cout << "masscut = " << masscutf_->EvalInstance(0) << endl;
    photonDeltaRcutf_->GetNdata();  cout << "photonDeltaRcut = " << photonDeltaRcutf_->EvalInstance(0) << endl;
    commoncutf_->GetNdata();  cout << "commoncut = " << commoncutf_->EvalInstance(0) << endl;
    cout << "globalTightHalo2016Filter = " << globalTightHalo2016Filter << endl;
    cout << "globalSuperTightHalo2016Filter = " << globalTightHalo2016Filter << endl;
    cout << "HBHENoiseFilter = " << HBHENoiseFilter << endl;
    cout << "HBHEIsoNoiseFilter = " << HBHEIsoNoiseFilter << endl;
    cout << "eeBadScFilter = " << eeBadScFilter << endl;
    cout << "EcalDeadCellTriggerPrimitiveFilter = " << EcalDeadCellTriggerPrimitiveFilter << endl;
    cout << "BadChargedCandidateFilter = " << BadChargedCandidateFilter << endl;
    cout << "BadPFMuonFilter = " << BadPFMuonFilter << endl;
    cout << "NVtx = " << NVtx << endl;
    cout << "PFCaloMETRatio = " << PFCaloMETRatio << endl;
    cout << "JetID = " << JetID << endl;
    cout << "HT = " << HT << ", HT5 = " << HT5 << ", HT5/HT = " << HT5/HT << endl;
    // cout << "JetIDclean = " << JetIDclean << endl;
    // cout << "HTclean = " << HTclean << ", HT5clean = " << HT5clean << ", HT5clean/HTclean = " << HT5clean/HTclean << endl;

    // if (baselineFormula->EvalInstance(0)) {
      // if (std::find(EvtNums.begin(), EvtNums.end(), EvtNum) != EvtNums.end()) {
	// cout << "Passes baseline" << endl;
	// printf("%15u %15u %15llu\n", RunNum, LumiBlockNum, EvtNum);
	// fprintf(idFile, "%15u %15u %15llu\n", RunNum, LumiBlockNum, EvtNum);
	// cout << " PFCaloMETRatio = " << PFCaloMETRatio << ", HT5/HT = " << HT5/HT << endl;
      // }
    // }
  }
  fclose(idFile);

}  // ======================================================================================
