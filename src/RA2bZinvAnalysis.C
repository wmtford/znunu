//
//  Zinv background prediction for RA2b analysis
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

#include <TMath.h>
using TMath::Sqrt; using TMath::Power;

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
    ("apply Z Pt cut", po::value<bool>(&applyPtCut_))
    ("apply photon min Delta R cut", po::value<bool>(&applyMinDeltaRCut_))
    ("use analysis bin from tree", po::value<bool>(&useTreeCCbin_))
    ("use DeepCSV", po::value<bool>(&useDeepCSV_))
    ("apply b-tag SF", po::value<bool>(&applyBTagSF_))
    ("apply pileup weight", po::value<bool>(&applyPuWeight_))
    ("use custom pileup weight", po::value<bool>(&customPuWeight_))
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
  activeBranches_.push_back("NVtx");
  activeBranches_.push_back("TriggerPass");
  activeBranches_.push_back("TriggerPrescales");
  activeBranches_.push_back("HBHENoiseFilter");
  activeBranches_.push_back("HBHEIsoNoiseFilter");
  activeBranches_.push_back("eeBadScFilter");
  activeBranches_.push_back("EcalDeadCellTriggerPrimitiveFilter");
  activeBranches_.push_back("globalTightHalo2016Filter");
  activeBranches_.push_back("BadChargedCandidateFilter");
  activeBranches_.push_back("BadPFMuonFilter");
  activeBranches_.push_back("PFCaloMETRatio");
  activeBranches_.push_back("nAllVertices");
  if (isMC_) {
    activeBranches_.push_back("puWeight");
    activeBranches_.push_back("Weight");
    activeBranches_.push_back("TrueNumInteractions");
    activeBranches_.push_back("madMinPhotonDeltaR");
    activeBranches_.push_back("GenMuons");
    activeBranches_.push_back("GenElectrons");
    activeBranches_.push_back("GenTaus");
  }

  if (era_ == TString("2016")) {

    if (applyPuWeight_ && customPuWeight_) {
      TFile* pufile = TFile::Open("../../Analysis/corrections/PileupHistograms_0121_69p2mb_pm4p6.root","READ");
      puHist_ = (TH1*) pufile->Get("pu_weights_down");
    }

    if (isMC_) BTagSFfile_ = "../../Analysis/btag/CSVv2_Moriond17_B_H_mod.csv";

  } // era 2016

  CCbinning CCmaps(era_, deltaPhi_);
  kinThresholds_ = CCmaps.kinThresholds();
  nJet1Thresholds_ = CCmaps.nJet1Thresholds();
  nJetThresholds_ = CCmaps.nJetThresholds();
  nbThresholds_ = CCmaps.nbThresholds();
  toCCbin_ = CCmaps.toCCbin();
  toCCbinjb_ = CCmaps.toCCbinjb();
  toCCbinSpl_ = CCmaps.toCCbinSpl();
  toCCbinJb_ = CCmaps.toCCbinJb();

  fillCutMaps();

  cout << "After initialization," << endl;
  cout << "The verbosity level is " << verbosity_ << endl;
  cout << "The ntuple version is " << ntupleVersion_ << endl;
  cout << "The MC flag is " << isMC_ << endl;
  cout << "The input-files-are-skims flag is " << isSkim_ << endl;
  cout << "The era is " << era_ << endl;
  cout << "The integrated luminosity = " << intLumi_ << endl;
  cout << "The path to input files is " << treeLoc_ << endl;
  cout << "The minDeltaPhi cuts are " << deltaPhi_ << endl;
  cout << "Use analysis bin from tree is " << useTreeCCbin_ << endl;
  cout << "Use DeepCSV is " << useDeepCSV_ << endl;
  cout << "The apply b-tag scale factors flag is " << applyBTagSF_ << endl;
  cout << "The apply pileup weight flag is " << applyPuWeight_ << endl;
  cout << "The custom pileup weight flag is " << customPuWeight_ << endl;


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
  if (deltaPhi_ == "ldp" && isSkim_) key += "ldp";
  if (!runBlock_.empty()) key += runBlock_;

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

  std::vector<TString> trigger;
  try {trigger = triggerMap_.at(sample);}
  catch (const std::out_of_range& oor) {trigger.clear();}
  if (trigger.empty()) {
    trigCuts_ = "1";
  } else {
    int Ntrig = trigger.size();
    if (Ntrig > 1) trigCuts_ = TString("(");
    for (auto theTrigger : trigger)
      trigCuts_ += TString("(TriggerPass[")+theTrigger+TString("]==1) + ");
    trigCuts_.Replace(trigCuts_.Length()-3, 3, "");
    if (Ntrig > 1) trigCuts_ += TString(")");
  }

  // // commonCuts_ = "(JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0 && BadPFMuonFilter && PFCaloMETRatio < 5)";  // Troy revision+
  // if (trigger.empty()) {
  //   commonCuts_ = "JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";  // Troy revision-;  moved JetID to !isSkim_
  // } else {
  //   commonCuts_ = "JetID==1&& globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";  // Troy revision-;  moved JetID to !isSkim_
  // }

  // When it's available, replace globalTightHalo2016Filter with globalSuperTightHalo2016Filter
  commonCuts_ = "globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";
  if (ntupleVersion_ != "V12") commonCuts_ += " && ecalBadCalibFilter==1";  // Added for 94X
  if (!isMC_) commonCuts_ += " &&  eeBadScFilter==1";
  if (!isSkim_) {
    commonCuts_ += " && JetIDclean";
    commonCuts_ += " && PFCaloMETRatio < 5";
    commonCuts_ += " && HT5clean/HTclean <= 2";
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming, single lepton
    // commonCuts_ += " && noFakeJet";  // Defined in loop, applied in skimming FastSim
  }

  massCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyMassCut_)
    massCut_ = "@ZCandidates.size()==1 && ZCandidates.M()>=76.188 && ZCandidates.M()<=106.188";

  ptCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyPtCut_)
    ptCut_ = "@ZCandidates.size()==1 && ZCandidates.Pt()>=200.";

  photonDeltaRcut_ = "1";
  if (sampleKey == "photon") {
    if (applyPtCut_) ptCut_ = "@Photons.size()==1 && Photons.Pt()>=200.";
    if (isMC_ && applyMinDeltaRCut_) photonDeltaRcut_ = "madMinPhotonDeltaR>=0.4";
  }
				       
    // 	if(extraCuts!=None):
    //     cuts+=extraCuts

  if (isSkim_) {
    HTcut_ = std::string("HT>=") + std::to_string(kinThresholds_[0][1]);
    NJetscut_ = std::string("NJets>=") + std::to_string(nJetThresholds_[0]);
  } else {
    HTcut_ = std::string("HTclean>=") + std::to_string(kinThresholds_[0][1]);
    NJetscut_ = std::string("NJetsclean>=") + std::to_string(nJetThresholds_[0]);
  }
  MHTcut_ = MHTCutMap_.at(deltaPhi_);
  objcut_ = objCutMap_.at(sampleKey);
  minDphicut_ = minDphiCutMap_.at(deltaPhi_);

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
  cuts += minDphicut_;
  cuts += massCut_;
  cuts += ptCut_;
  cuts += photonDeltaRcut_;
  cuts += commonCuts_;
  cuts += trigCuts_;

    // 	  if(applySF):
    //     cuts*=bJetCutsSF[bJetBin]
    // 	else:
    //     cuts+=bJetCuts[bJetBin]
  if(applyPuWeight_ && !customPuWeight_) cuts *= "puWeight*(1)";
    // 	  if(type(extraWeight) is str):
    //     extraWeight+="*(1)"
    //     cuts*=extraWeight

  return cuts;
 
}  // ======================================================================================

void
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample, std::vector<histConfig*>& histograms, TCut baselineCuts) {
  //
  // Define N - 1 (or N - multiple) cuts, book histograms.  Traverse the chain and fill.
  //
  if (verbosity_ >= 1) cout << endl << "baseline = " << endl << baselineCuts << endl << endl;
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

  cutHistos cutHistFiller(chain, forNotify);  // for cutFlow histograms
  for (auto & hg : histograms) {
    if (hg->NbinsY > 0) {
      hg->hist = new TH2F(hg->name, hg->title, hg->NbinsX, hg->rangeX.first, hg->rangeX.second,
			  hg->NbinsY, hg->rangeY.first, hg->rangeY.second);
      hg->hist->SetOption("colz");
    } else {
      hg->hist = new TH1F(hg->name, hg->title, hg->NbinsX, hg->rangeX.first, hg->rangeX.second);
      hg->hist->SetOption("hist");
      hg->hist->SetMarkerSize(0);
    }
    hg->hist->GetXaxis()->SetTitle(hg->axisTitles.first);
    hg->hist->GetYaxis()->SetTitle(hg->axisTitles.second);
    if (hg->name.Contains(TString("hFilterCuts"))) {
      hg->hist->GetXaxis()->SetBinLabel(1, "None");
      hg->hist->GetXaxis()->SetBinLabel(2, "TightHalo");
      hg->hist->GetXaxis()->SetBinLabel(3, "SupTgtHalo");
      hg->hist->GetXaxis()->SetBinLabel(4, "HBENoise");
      hg->hist->GetXaxis()->SetBinLabel(5, "HBEIsoNoise");
      hg->hist->GetXaxis()->SetBinLabel(6, "EcalDeadCell");
      hg->hist->GetXaxis()->SetBinLabel(7, "BadChCand");
      hg->hist->GetXaxis()->SetBinLabel(8, "BadPFMu");
      hg->hist->GetXaxis()->SetBinLabel(9, "NVtx");
      hg->hist->GetXaxis()->SetBinLabel(10, "eeBadSc");
      hg->hist->GetXaxis()->SetBinLabel(11, "ecalBadCal");
      hg->hist->GetXaxis()->SetBinLabel(12, "");
      hg->hist->GetXaxis()->SetBinLabel(13, "JetID");
      hg->hist->GetXaxis()->SetBinLabel(14, "PFCaloMETR");
      hg->hist->GetXaxis()->SetBinLabel(15, "HT5/HT");
      hg->hist->GetXaxis()->LabelsOption("vu");
    }
    if (hg->name.Contains(TString("hCut"))) {
      hg->NminusOneCuts = "1";
      hg->hist->GetXaxis()->SetBinLabel(1, "None");
      hg->hist->GetXaxis()->SetBinLabel(2, "HT");
      hg->hist->GetXaxis()->SetBinLabel(3, "MHT");
      hg->hist->GetXaxis()->SetBinLabel(4, "NJets");
      hg->hist->GetXaxis()->SetBinLabel(5, "mnDphi");
      hg->hist->GetXaxis()->SetBinLabel(6, "objects");
      hg->hist->GetXaxis()->SetBinLabel(7, "Zpt");
      hg->hist->GetXaxis()->SetBinLabel(8, "Zmass");
      hg->hist->GetXaxis()->SetBinLabel(9, "Trigger");
      hg->hist->GetXaxis()->SetBinLabel(10, "Filters");
      hg->hist->GetXaxis()->LabelsOption("vu");
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

  Long64_t Nentries = chain->GetEntries();
  if (verbosity_ >= 1) cout << "Nentries in tree = " << Nentries << endl;
  int count = 0;
  int countInFile = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (verbosity_ >= 1 && count % 100000 == 0) cout << "Entry number " << count << endl;

    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) {
    	if (btagcorr_) btagcorr_->SetEffs(thisFile);
	if (verbosity_ >= 1) cout << "Current file in chain: " << thisFile->GetName() << endl;
      }
      chain->GetEntry(entry);  // Pull in tree variables for reinitialization
      checkActiveTrigPrescales(sample);
      if (isMC_ && verbosity_ >= 1) cout << "MC weight for this file is " << Weight << endl;
    }
    chain->GetEntry(entry);
    cleanVars();  // If unskimmed input, copy <var>clean to <var>
    Int_t BTagsOrig = setBTags();
    // if (countInFile <= 100) cout << "BTagsOrig, BTags, BTagsDeepCSV = " << BTagsOrig << ", " << BTags << ", " << BTagsDeepCSV << endl;

    if (ZCandidates->size() > 1 && verbosity_ >= 2) cout << ZCandidates->size() << " Z candidates found" << endl;
    // baselineFormula->GetNdata();
    // double baselineWt = baselineFormula->EvalInstance(0);

    // Compute event weight factors
    Double_t eventWt = 1;
    Double_t PUweight = 1;
    if (applyPuWeight_ && customPuWeight_) {
      // This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    if (isMC_) {
      eventWt = 1000*intLumi_*Weight*PUweight;
      if (eventWt < 0) eventWt *= -1;
    }

    for (auto & hg : histograms) {
      if (hg->name.Contains(TString("hCut"))) {
	cutHistFiller.fill((TH1F*) hg->hist, eventWt);
	continue;
      }
      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt != 0) {
	if (hg->dvalue != nullptr) {
	  hg->hist->Fill(*(hg->dvalue), selWt*eventWt);
	}
	else if (hg->ivalue != nullptr) {
	  hg->hist->Fill(Double_t(*(hg->ivalue)), selWt*eventWt);
	}
	else if (hg->filler1D != nullptr) (this->*(hg->filler1D))((TH1F*) hg->hist, selWt*eventWt);
	else if (hg->filler2D != nullptr) (this->*(hg->filler2D))((TH2F*) hg->hist, selWt*eventWt);
	else cerr << "No method to fill histogram provided for " << hg->name << endl;

      }
    }  // loop over histograms
  }  // loop over entries

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

  histConfig hCC;
  hCC.name = TString("hCC_") + sample;  hCC.title = "Cut & count analysis";
  Int_t MaxBins = toCCbin_.size();
  cout << "MaxBins = " << MaxBins << endl;
  hCC.NbinsX = MaxBins;  hCC.rangeX.first = 0.5;  hCC.rangeX.second = MaxBins+0.5;
  hCC.axisTitles.first = "Njets, Nb, (HT, MHT)";  hCC.axisTitles.second = "Events / bin";
  hCC.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCC);

  histConfig hCCjb;
  hCCjb.name = TString("hCCjb_") + sample;  hCCjb.title = "Cut & count, kinematics integrated";
  Int_t MaxBinsjb = toCCbinjb_.size();
  cout << "MaxBinsjb = " << MaxBinsjb << endl;
  hCCjb.NbinsX = MaxBinsjb;  hCCjb.rangeX.first = 0.5;  hCCjb.rangeX.second = MaxBinsjb+0.5;
  hCCjb.axisTitles.first = "Njets, Nb";  hCCjb.axisTitles.second = "Events / bin";
  hCCjb.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCjb);

  histConfig hCCspl;
  hCCspl.name = TString("hCCspl_") + sample;  hCCspl.title = "Cut & count, NJet split";
  Int_t MaxBinsSpl = toCCbinSpl_.size();
  cout << "MaxBinsSpl = " << MaxBinsSpl << endl;
  hCCspl.NbinsX = MaxBinsSpl;  hCCspl.rangeX.first = 0.5;  hCCspl.rangeX.second = MaxBinsSpl+0.5;
  hCCspl.axisTitles.first = "Njets, Nb, (HT, MHT)";  hCCspl.axisTitles.second = "Events / bin";
  hCCspl.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCspl);

  histConfig hCCJb;
  hCCJb.name = TString("hCCJb_") + sample;  hCCJb.title = "Cut & count, kinematics integrated, JNet split";
  Int_t MaxBinsJb = toCCbinJb_.size();
  cout << "MaxBinsJb = " << MaxBinsJb << endl;
  hCCJb.NbinsX = MaxBinsJb;  hCCJb.rangeX.first = 0.5;  hCCJb.rangeX.second = MaxBinsJb+0.5;
  hCCJb.axisTitles.first = "Njets, Nb";  hCCJb.axisTitles.second = "Events / bin";
  hCCJb.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCJb);

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

  histConfig hnZcand;
  hnZcand.name = TString("hnZcand_") + sample;  hnZcand.title = "Number of Z candidates";
  hnZcand.NbinsX = 10;  hnZcand.rangeX.first = 0;  hnZcand.rangeX.second = 10;
  hnZcand.axisTitles.first = "N(Z candidates)";  hnZcand.axisTitles.second = "Events / bin";
  hnZcand.filler1D = &RA2bZinvAnalysis::fillnZcand;  hnZcand.omitCuts.push_back(&massCut_);
  histograms.push_back(&hnZcand);

  histConfig hZmass;
  hZmass.name = TString("hZmass_") + sample;  hZmass.title = "Z mass";
  hZmass.NbinsX = 30;  hZmass.rangeX.first = 60;  hZmass.rangeX.second = 120;
  hZmass.axisTitles.first = "M(Z) [GeV]";  hZmass.axisTitles.second = "Events / 2 GeV";
  hZmass.filler1D = &RA2bZinvAnalysis::fillZmass;  hZmass.omitCuts.push_back(&massCut_);
  histograms.push_back(&hZmass);

  histConfig hZpt;
  hZpt.name = TString("hZpt_") + sample;  hZpt.title = "Z Pt";
  hZpt.NbinsX = 60;  hZpt.rangeX.first = 0;  hZpt.rangeX.second = 3000;
  hZpt.axisTitles.first = "Pt(Z) [GeV]";  hZpt.axisTitles.second = "Events / 50 GeV";
  hZpt.filler1D = &RA2bZinvAnalysis::fillZpt;  hZpt.omitCuts.push_back(&ptCut_);  hZpt.omitCuts.push_back(&MHTcut_);
  histograms.push_back(&hZpt);

  histConfig hCutFlow;
  hCutFlow.name = TString("hCutFlow_") + sample;  hCutFlow.title = "Cut flow";
  hCutFlow.NbinsX = 10;  hCutFlow.rangeX.first = 0;  hCutFlow.rangeX.second = 10;
  hCutFlow.axisTitles.first = "";  hCutFlow.axisTitles.second = "Events surviving";
  histograms.push_back(&hCutFlow);

  histConfig hCuts;
  hCuts.name = TString("hCuts_") + sample;  hCuts.title = "Cuts passed";
  hCuts.NbinsX = 10;  hCuts.rangeX.first = 0;  hCuts.rangeX.second = 10;
  hCuts.axisTitles.first = "";  hCuts.axisTitles.second = "Events passing";
  histograms.push_back(&hCuts);

  histConfig hFilterCuts;
  hFilterCuts.name = TString("hFilterCuts_") + sample;  hFilterCuts.title = "Filter cuts";
  hFilterCuts.NbinsX = 15;  hFilterCuts.rangeX.first = 0;  hFilterCuts.rangeX.second = 15;
  hFilterCuts.axisTitles.first = "";  hFilterCuts.axisTitles.second = "Events failing";
  hFilterCuts.filler1D = &RA2bZinvAnalysis::fillFilterCuts;
  hFilterCuts.omitCuts.push_back(&commonCuts_);
  histograms.push_back(&hFilterCuts);

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
  histograms.push_back(&hTrueNumInt);

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

  // Z mass in Njet, Nb bins
  histConfig hZmass_2j0b(hZmass);
  hZmass_2j0b.name = TString("hZmass_2j0b_") + sample;  hZmass_2j0b.title = "Z mass, 2 jets & 0 b jets";
  hZmass_2j0b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_2j0b);

  histConfig hZmass_2j1b(hZmass);
  hZmass_2j1b.name = TString("hZmass_2j1b_") + sample;  hZmass_2j1b.title = "Z mass, 2 jets & 1 b jet";
  hZmass_2j1b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_2j1b);

  histConfig hZmass_2j2b(hZmass);
  hZmass_2j2b.name = TString("hZmass_2j2b_") + sample;  hZmass_2j2b.title = "Z mass, 2 jets & >=2 b jets";
  hZmass_2j2b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_2j2b);
  //
  histConfig hZmass_3j0b(hZmass);
  hZmass_3j0b.name = TString("hZmass_3j0b_") + sample;  hZmass_3j0b.title = "Z mass, 3-4 jets & 0 b jets";
  hZmass_3j0b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_3j0b);

  histConfig hZmass_3j1b(hZmass);
  hZmass_3j1b.name = TString("hZmass_3j1b_") + sample;  hZmass_3j1b.title = "Z mass, 3-4 jets & 1 b jet";
  hZmass_3j1b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_3j1b);

  histConfig hZmass_3j2b(hZmass);
  hZmass_3j2b.name = TString("hZmass_3j2b_") + sample;  hZmass_3j2b.title = "Z mass, 3-4 jets & >=2 b jets";
  hZmass_3j2b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_3j2b);
  //
  histConfig hZmass_5j0b(hZmass);
  hZmass_5j0b.name = TString("hZmass_5j0b_") + sample;  hZmass_5j0b.title = "Z mass, >=5 jets & 0 B jets";
  hZmass_5j0b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_5j0b);

  histConfig hZmass_5j1b(hZmass);
  hZmass_5j1b.name = TString("hZmass_5j1b_") + sample;  hZmass_5j1b.title = "Z mass, >=5 jets & 1 B jet";
  hZmass_5j1b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_5j1b);

  histConfig hZmass_5j2b(hZmass);
  hZmass_5j2b.name = TString("hZmass_5j2b_") + sample;  hZmass_5j2b.title = "Z mass, >=5 jets & >=2 B jets";
  hZmass_5j2b.filler1D = &RA2bZinvAnalysis::fillZmassjb;
  histograms.push_back(&hZmass_5j2b);

  bookAndFillHistograms(sample, histograms, baselineCuts);

  std::vector<TH1*> theHists;
  for (auto & thisHist : histograms) theHists.push_back(thisHist->hist);

  return theHists;

}  // ======================================================================================

void
RA2bZinvAnalysis::fillCC(TH1F* h, double wt) {

  // Filler for the analysis-bin jbk and jb histograms
  TString hName(h->GetName());
  UInt_t binCC = 0, binCCjb = 0;

  if (useTreeCCbin_) {
    binCC = RA2bin;
    if (binCC > 0) {
      if (hName.Contains("jb") || hName.Contains("Jb") || hName.Contains("spl")) {
      // 	// Calculate binCCjb
      // 	int binNjets = nJetThresholds_.size()-1;
      // 	while (NJets < nJetThresholds_[binNjets]) binNjets--;
      // 	int binNb = nbThresholds_.size()-1;
      // 	while (BTags < nbThresholds_[binNb]) binNb--;
      // 	std::vector<int> jb = {binNjets, binNb};
      // 	try {binCCjb = toCCbinjb_.at(jb);}
      // 	catch (const std::out_of_range& oor) {return;}
      // 	h->Fill(Double_t(binCCjb), wt);
      // } else {
	h->Fill(Double_t(binCC), wt);
      }
    }
  } else {
    // Calculate binCC
    int binKin = kinBin(HT, MHT);
    if (binKin < 0) return;
    int binNjets;
    if (hName.Contains("spl") || hName.Contains("Jb")) {
      binNjets = nJet1Thresholds_.size()-1;
      while (NJets < nJet1Thresholds_[binNjets]) binNjets--;
    } else {
      binNjets = nJetThresholds_.size()-1;
      while (NJets < nJetThresholds_[binNjets]) binNjets--;
    }
    if (!applyBTagSF_) {
      int binNb = nbThresholds_.size()-1;
      while (BTags < nbThresholds_[binNb]) binNb--;
      std::vector<int> jbk = {binNjets, binNb, binKin};
      if (hName.Contains("spl") || hName.Contains("Jb")) {
	try {binCC = toCCbinSpl_.at(jbk);}
	catch (const std::out_of_range& oor) {return;}
      } else {
	try {binCC = toCCbin_.at(jbk);}
	catch (const std::out_of_range& oor) {return;}  // Omitted bins j = 3,4, k = 0,3
      }
      if (verbosity_ >= 4) cout << "j = " << binNjets << ", b = " << binNb << ", k = " << binKin << ", binCC = "
				<< binCC << ", RA2bin = " << RA2bin << endl;
      if (hName.Contains("jb") || hName.Contains("Jb")) {
	// Above test on jbk needed even here, to exclude j = 3,4, k = 0,3
	std::vector<int> jb = {binNjets, binNb};
	try {binCCjb = hName.Contains("jb") ? toCCbinjb_.at(jb) : toCCbinJb_.at(jb);}
      	catch (const std::out_of_range& oor) {return;}
	h->Fill(Double_t(binCCjb), wt);
      } else {
	h->Fill(Double_t(binCC), wt);
      }
    } else {
      // apply BTagSF to all Nb bins
      if (verbosity_ >= 4) cout << "Size of input Jets = " << Jets->size()
				<< ", Jets_hadronFlavor = " << Jets_hadronFlavor->size()
				<< " Jets_HTMask = " << Jets_HTMask->size() << endl;
      vector<double> probNb = btagcorr_->GetCorrections(Jets, Jets_hadronFlavor, Jets_HTMask);
      for (int binNb = 0; binNb < (int) nbThresholds_.size(); ++binNb) {
	std::vector<int> jbk = {binNjets, binNb, binKin};
	if (hName.Contains("spl")) {
	  try {binCC = toCCbinSpl_.at(jbk);}
	  catch (const std::out_of_range& oor) {return;}
	} else {
	  try {binCC = toCCbin_.at(jbk);}
	  catch (const std::out_of_range& oor) {return;}
	}
	if (verbosity_ >= 4) cout << "j = " << binNjets << ", NbTags = " << BTags << ", b = " << binNb
				  << ", b wt = " << probNb[binNb] << ", k = " << binKin << ", binCC = " << binCC << endl;
	if (hName.Contains("jb")) {
	  std::vector<int> jb = {binNjets, binNb};
	  try {binCCjb = toCCbinjb_.at(jb);}
	  catch (const std::out_of_range& oor) {return;}
	  h->Fill(Double_t(binCCjb), wt*probNb[binNb]);
	} else {
	  h->Fill(Double_t(binCC), wt*probNb[binNb]);
	}
      }
    }  // if apply BTagSF
  }  // if useTreeCCbin

}  // ======================================================================================

void
RA2bZinvAnalysis::fillZmassjb(TH1F* h, double wt) {
  if (ZCandidates->size() == 0) return;
  Int_t j, b;
  TString hName(h->GetName());
  j = hName(hName.First('j')-1) - '0';
  b = hName(hName.First('b')-1) - '0';
  if ((j == 2 && NJets == 2) || (j == 3 && (NJets >=3 && NJets <= 4)) || (j == 5 && NJets >=5)) {
    if ((b < 2 && BTags == b) || (b == 2 && BTags >= 2)) {
      h->Fill(ZCandidates->at(0).M(), wt);
    }
  }
  
}  // ======================================================================================

void
RA2bZinvAnalysis::fillFilterCuts(TH1F* h, double wt) {
  h->Fill(0.5, wt);
  if (!(globalTightHalo2016Filter==1)) h->Fill(1.5, wt);
  // if (!(globalSuperTightHalo2016Filter==1)) h->Fill(2.5, wt);
  if (!(HBHENoiseFilter==1)) h->Fill(3.5, wt);
  if (!(HBHEIsoNoiseFilter==1)) h->Fill(4.5, wt);
  if (!(EcalDeadCellTriggerPrimitiveFilter==1)) h->Fill(5.5, wt);
  if (!BadChargedCandidateFilter) h->Fill(6.5, wt);
  if (!BadPFMuonFilter) h->Fill(7.5, wt);
  if (!(NVtx > 0)) h->Fill(8.5, wt);
  if (!(eeBadScFilter==1)) h->Fill(9.5, wt);
  if (!(ecalBadCalibFilter==1)) h->Fill(10.5, wt);
  if (!(JetID)) h->Fill(12.5, wt);
  if (!(PFCaloMETRatio < 5)) h->Fill(13.5, wt);
  if (!(HT5/HT <= 2)) h->Fill(14.5, wt);

}  // ======================================================================================

void
RA2bZinvAnalysis::fillZGmass(TH1F* h, double wt) {
  for (auto & theZ : *ZCandidates) {
    for (auto & aPhoton : *Photons) {
      TLorentzVector Zg(theZ);  Zg += aPhoton;
      h->Fill(Zg.M(), wt);
    }
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGJdR(TH1F* h, double wt) {
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
RA2bZinvAnalysis::fillGLdRnoPixelSeed(TH1F* h, double wt) {
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
RA2bZinvAnalysis::fillGLdRpixelSeed(TH1F* h, double wt) {
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

int
RA2bZinvAnalysis::kinBin(double& ht, double& mht) {
  int theBin = -1;
  int NmhtBins = kinThresholds_.size() - 1;
  if (deltaPhi_ != TString("nominal")) {
    // ldp or hdp
    if (mht < kinThresholds_[NmhtBins][0] || ht < kinThresholds_[NmhtBins][1]) return theBin;
    if (mht < kinThresholds_[0][0]) {
      theBin += 10;
      for (unsigned j = 1; j < kinThresholds_[NmhtBins].size(); ++j) {
	theBin++;
	if (ht >= kinThresholds_[NmhtBins][j] && (j == kinThresholds_[NmhtBins].size()-1 || ht < kinThresholds_[NmhtBins][j+1])) return theBin;
      }
    }
  }
  if (mht < kinThresholds_[0][0] || ht < kinThresholds_[0][1]) return theBin;
  for (int i = 0; i < NmhtBins; ++i) {
    if (mht >= kinThresholds_[i][0] && (i == NmhtBins-1 || mht < kinThresholds_[i+1][0])) {
      for (unsigned j = 1; j < kinThresholds_[i].size(); ++j) {
	theBin++;
	if (ht >= kinThresholds_[i][j] && (j == kinThresholds_[i].size()-1 || ht < kinThresholds_[i][j+1])) return theBin;
      }
    } else {
      theBin += kinThresholds_[i].size() - 1;
    }
  }
  return -2;  // Outside binned area (e.g., MHT > HT), though within preselection
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

    } else if (ntupleVersion_ == "V15") {

      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
      objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0";  // Troy mod+
      // objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
      objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracks==0 && isoPionTracks==0) || (NMuons==0 && NElectrons==2 && isoMuonTracks==0 && isoPionTracks==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
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
    } else if (ntupleVersion_ == "V15") {
      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      // objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoMuonTracksclean==0";
      objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0";
      // objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoElectronTracksclean==0";
      objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0) || (NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
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

  if (ntupleVersion_ == "V12") {
    triggerMap_["zmm"] = {"18", "20", "22", "23", "29"};
    triggerMap_["zee"] = {"3", "4", "6", "7", "11", "12"};
    triggerMap_["photon"] = {"52"};  // re-miniAOD; 51 for ReReco/PromptReco
    triggerMap_["sig"] = {"42", "43", "44", "46", "47", "48"};
  } else if (ntupleVersion_ == "V15") {
    triggerMap_["zmm"] = {"48", "50", "52", "53", "55", "63"};  // 48 prescaled in late 2017 --Owen; add 50
    triggerMap_["zee"] = {"21", "23", "28", "35", "40", "41"};
    triggerMap_["photon"] = {"141"};
    triggerMap_["sig"] = {"108", "112", "124", "128"};
  }
  triggerMap_["zll"].reserve(triggerMap_["zmm"].size() + triggerMap_["zee"].size());
  triggerMap_["zll"] = triggerMap_["zmm"];
  triggerMap_["zll"].insert(triggerMap_["zll"].end(), triggerMap_["zee"].begin(), triggerMap_["zee"].end());
  triggerMap_["sle"] = triggerMap_["sig"];
  triggerMap_["slm"] = triggerMap_["sig"];

}  // ======================================================================================

Int_t
RA2bZinvAnalysis::setBTags() {
  Int_t BTagsOrig = BTags;
  if (useDeepCSV_) {
    BTags = BTagsDeepCSV;
  } else if (ntupleVersion_ == "V15" && runBlock_.find("2016")==std::string::npos) {
    // Recompute BTags with a different discriminator threshold
    BTags = 0;
    for (size_t j = 0; j < Jets->size(); ++j) {
      if(!Jets_HTMask->at(j)) continue;
      if (Jets_bDiscriminatorCSV->at(j) > csvMthreshold_) BTags++;
    }
  }
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
  trigcutf_ = new TTreeFormula("trigcut", trigCuts_, chain);  forNotify->Add(trigcutf_);
  commoncutf_ = new TTreeFormula("commoncut", commonCuts_, chain);  forNotify->Add(commoncutf_);
}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::fill(TH1F* hcf, Double_t wt) {
  HTcutf_->GetNdata();
  MHTcutf_->GetNdata();
  NJetscutf_->GetNdata();
  minDphicutf_->GetNdata();
  objcutf_->GetNdata();
  ptcutf_->GetNdata();
  masscutf_->GetNdata();
  trigcutf_->GetNdata();
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
		  if (trigcutf_->EvalInstance(0)) {
		    hcf->Fill(8.5, wt);
		    if (commoncutf_->EvalInstance(0)) {
		      hcf->Fill(9.5, wt);
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
    if (trigcutf_->EvalInstance(0)) hcf->Fill(8.5, wt);
    if (commoncutf_->EvalInstance(0)) hcf->Fill(9.5, wt);
  }
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
      countInFile = 0;
    }
    chain->GetEntry(entry);
    countInFile++;
    if (countInFile == 1) {
      int trigNo = 0;
      for (auto & theTrigPrescale : *TriggerPrescales) {
	if (theTrigPrescale != 1) cout << trigNo << ":  " << theTrigPrescale << endl;
    	++trigNo;
      }
    }
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::checkActiveTrigPrescales(const char* sample) {
  std::vector<TString> trigger;
  try {trigger = triggerMap_.at(sample);}
  catch (const std::out_of_range& oor) {trigger.clear();}
  for (auto theTrigger : trigger) {
  stringstream ss;
  ss << theTrigger;
    string temp;
    int trgIndex;
    while (!ss.eof()) {
      ss >> temp;
      if (stringstream(temp) >> trgIndex) {
	Int_t prescale = TriggerPrescales->at(trgIndex);
	if (verbosity_ >= 1 && prescale != 1) cout << "Trigger " << trgIndex << " prescaled to " << prescale << endl;
      }
      temp = "";
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
  TTreeFormula* trigcutf_ = new TTreeFormula("trigcut", trigCuts_, chain);  forNotify->Add(trigcutf_);
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
    trigcutf_->GetNdata();  cout << "trigcut = " << trigcutf_->EvalInstance(0) << endl;
    commoncutf_->GetNdata();  cout << "commoncut = " << commoncutf_->EvalInstance(0) << endl;
    cout << "globalTightHalo2016Filter = " << globalTightHalo2016Filter << endl;
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
