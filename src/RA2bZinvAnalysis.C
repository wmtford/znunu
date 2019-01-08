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
    ("apply Z/gamma Pt cut", po::value<bool>(&applyPtCut_))
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

  // For now we have only 2016 correction files
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
  cout << "Apply photon min Delta R cut is " << applyMinDeltaRCut_ << endl;
  cout << "Use analysis bin from tree is " << useTreeCCbin_ << endl;
  cout << "Use DeepCSV is " << useDeepCSV_ << endl;
  cout << "Apply b-tag scale factors is " << applyBTagSF_ << endl;
  cout << "Apply pileup weight is " << applyPuWeight_ << endl;
  cout << "The custom pileup weight flag is " << customPuWeight_ << endl;
  
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
  else if (theSample.Contains("gjets")) key = TString("gjets");
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

  // // commonCuts_ = "(JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0 && BadPFMuonFilter && PFCaloMETRatio < 5)";  // Troy revision+
  // if (trigger.empty()) {
  //   commonCuts_ = "JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";  // Troy revision-;  moved JetID to !isSkim_
  // } else {
  //   commonCuts_ = "JetID==1&& globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";  // Troy revision-;  moved JetID to !isSkim_
  // }

  if (ntupleVersion_ == "V12" || ntupleVersion_ == "V15")
    commonCuts_ = "globalTightHalo2016Filter==1";
  else
    commonCuts_ = "globalSuperTightHalo2016Filter==1";
  commonCuts_ += " && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";
  if (ntupleVersion_ != "V12") commonCuts_ += " && ecalBadCalibFilter==1";  // Added for 94X
  if (!isMC_) commonCuts_ += " && eeBadScFilter==1";
  // Kevin Pedro, re V16:  Please note that I have included only the JetID "event cleaning" cut in these skims. 
  // The PFCaloMETRatio, HT5/HT, and noMuonJet cuts are left out,
  // so the impact of these cuts can be tested and refined if necessary.
  if (isSkim_ && ntupleVersion_ == "V16") {
    commonCuts_ += " && PFCaloMETRatio < 5";
    commonCuts_ += " && HT5/HT <= 2";
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming (xV16), single lepton
  }
  if (!isSkim_) {
    commonCuts_ += " && JetIDclean";
    commonCuts_ += " && PFCaloMETRatio < 5";
    commonCuts_ += " && HT5clean/HTclean <= 2";
    // commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming (xV16), single lepton
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
  cuts += photonDeltaRcut_;
  cuts += commonCuts_;

  return cuts;
 
}  // ======================================================================================

void
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample, std::vector<histConfig*>& histograms, TCut baselineCuts) {
  //
  // Define N - 1 (or N - multiple) cuts, book histograms.  Traverse the chain and fill.
  //
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

  effPurCorr_.getHistos(sample);  // For purity, Fdir, trigger eff, reco eff
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
    if (hg->name.Contains(TString("hCut"))) {
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
  Long64_t Nentries = chain->GetEntries();
  if (verbosity_ >= 1) cout << "Nentries in tree = " << Nentries << endl;
  int count = 0;
  int countInFile = 0;
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
	if (verbosity_ >= 1) cout << "Current file in chain: " << thisFile->GetName() << endl;
      }
      chain->GetEntry(entry);  // Pull in tree variables for reinitialization
      setTriggerIndexList(sample);
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
    if (applyPuWeight_) {
      if (customPuWeight_ && puHist_ != nullptr) {
	// This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
	PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
      } else
	PUweight = puWeight;  // Take puWeight directly from the tree
    }

    if (isMC_) {
      eventWt = 1000*intLumi_*Weight*PUweight;
      if (eventWt < 0) eventWt *= -1;
    }

    // Trigger requirements
    bool passTrg = true;
    if (!isMC_) {
      passTrg = false;
      for (auto trgIndex : triggerIndexList_)
	if (TriggerPass->at(trgIndex)) passTrg = true;
    }

    for (auto & hg : histograms) {
      if (hg->name.Contains(TString("hCut"))) {
	cutHistFiller.fill((TH1D*) hg->hist, eventWt, passTrg);
	continue;
      }
      if (!passTrg) break;

      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt == 0) continue;
      double eventWt0 = eventWt;

      if (hg->name.Contains(TString("_DR"))) {
	int CCbin = CCbins_->jbk(CCbins_->jbin(NJets), CCbins_->bbin(NJets, BTags), CCbins_->kinBin(HT, MHT));
      	if (CCbin <= 0 || BTags > 0) continue;
      	// For double ratio, apply weights for purity, Fdir, trigger eff, reco eff.
	effWt_ = effPurCorr_.weight(CCbins_, NJets, BTags, MHT, HT, *ZCandidates, *Photons, *Photons_isEB);
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
  bool isPhoton = (sampleKey == "photon");

  histConfig hCutFlow;
  hCutFlow.name = TString("hCutFlow_") + sample;  hCutFlow.title = "Cut flow";
  hCutFlow.NbinsX = 12;  hCutFlow.rangeX.first = 0;  hCutFlow.rangeX.second = 12;
  hCutFlow.axisTitles.first = "";  hCutFlow.axisTitles.second = "Events surviving";
  histograms.push_back(&hCutFlow);

  histConfig hCuts;
  hCuts.name = TString("hCuts_") + sample;  hCuts.title = "Cuts passed";
  hCuts.NbinsX = 12;  hCuts.rangeX.first = 0;  hCuts.rangeX.second = 12;
  hCuts.axisTitles.first = "";  hCuts.axisTitles.second = "Events passing";
  histograms.push_back(&hCuts);
  // These hCut* histograms are filled before the trigger cut.

  histConfig hFilterCuts;
  hFilterCuts.name = TString("hFilterCuts_") + sample;  hFilterCuts.title = "Filter cuts";
  hFilterCuts.NbinsX = 15;  hFilterCuts.rangeX.first = 0;  hFilterCuts.rangeX.second = 15;
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
  if (isPhoton) histograms.push_back(&hNJets_DRCC);
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

  histConfig hnZcand;
  hnZcand.name = TString("hnZcand_") + sample;  hnZcand.title = "Number of Z candidates";
  hnZcand.NbinsX = 10;  hnZcand.rangeX.first = 0;  hnZcand.rangeX.second = 10;
  hnZcand.axisTitles.first = "N(Z candidates)";  hnZcand.axisTitles.second = "Events / bin";
  hnZcand.filler1D = &RA2bZinvAnalysis::fillnZcand;  hnZcand.omitCuts.push_back(&massCut_);
  histograms.push_back(&hnZcand);

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
  if (!(HT5/HT <= 2)) h->Fill(14.5, wt);

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
      // photon skim cuts:  "Sum$(Photons_fullID)==1 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0"
      objCutMap_["photon"] = "@Photons.size()==1 && Sum$(Photons_nonPrompt)==0 && Photons_hasPixelSeed==0";  // Andrew recommendation, no skim cuts
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Photons.at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && NMuons==0 && NElectrons==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
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
    hcf->GetXaxis()->SetBinLabel(14, "14 PFCaloMETR");
    hcf->GetXaxis()->SetBinLabel(15, "15 HT5/HT");
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
    hcf->GetXaxis()->SetBinLabel(11, "11 Filters");
    hcf->GetXaxis()->LabelsOption("vu");
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::fill(TH1D* hcf, Double_t wt, bool passTrg) {
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
		      if (commoncutf_->EvalInstance(0)) {
			hcf->Fill(10.5, wt);
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
    if (commoncutf_->EvalInstance(0)) hcf->Fill(10.5, wt);
  }
}  // ======================================================================================

void 
RA2bZinvAnalysis::efficiencyAndPurity::openFiles() {
  TString plotDir("../plots/histograms/");
  TString effs = "effHists.root";
  purityTrigEffFile_ = new TFile((plotDir+effs).Data(), "read");

}  // ======================================================================================

void
RA2bZinvAnalysis::efficiencyAndPurity::getHistos(const char* sample) {
  // For purity, Fdir, trigger eff, reco eff
  theSample_ = TString(sample);
  hSFeff_ = nullptr;
  FdirHist_ = nullptr;
  hPurity_.clear();
  hTrigEff_.clear();
  if (theSample_.Contains("zmm") && !theSample_.Contains("tt")) {
    hPurity_.push_back((TH1F*) purityTrigEffFile_->Get("h_pur_m"));
    if (hPurity_.back() == nullptr) cout << "***** Histogram h_pur_m not found *****" << endl;
  } else if (theSample_.Contains("zee") && !theSample_.Contains("tt")) {
    hPurity_.push_back((TH1F*) purityTrigEffFile_->Get("h_pur_e"));
    if (hPurity_.back() == nullptr) cout << "***** Histogram h_pur_e not found *****" << endl;
  } else if (theSample_.Contains("photon")) {
    hPurity_.push_back((TH1F*) purityTrigEffFile_->Get("h_pur_eb"));
    if (hPurity_.back() == nullptr) cout << "***** Histogram h_pur_eb not found *****" << endl;
    hPurity_.push_back((TH1F*) purityTrigEffFile_->Get("h_pur_ec"));
    if (hPurity_.back() == nullptr) cout << "***** Histogram h_pur_ec not found *****" << endl;
    FdirHist_ = (TH1D*) purityTrigEffFile_->Get("h_bin46_NJets8910");
    if (FdirHist_ == nullptr) cout << "***** Histogram h_bin46_NJets8910 not found *****" << endl;
    // FdirGraph_ = (TGraphErrors*) purityTrigEffFile_->Get("bin46_f");
    // if (FdirGraph_ == nullptr) cout << "***** Histogram bin46_f not found *****" << endl;
  } else if (theSample_.Contains("dymm")) {
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_->Get("h_trig_m"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_m not found *****" << endl;
    hSFeff_ = (TH1F*) purityTrigEffFile_->Get("h_SFm_MHT");
    if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;
  } else if (theSample_.Contains("dyee")) {
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_->Get("h_trig_e"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_e not found *****" << endl;
    hSFeff_ = (TH1F*) purityTrigEffFile_->Get("h_SFe_MHT");  // Maybe this should be h_NJets
    if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;
  } else if (theSample_.Contains("gjets")) {
    // Troy-style trigger efficiency histograms:
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_->Get("h_trig_eb"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_eb not found *****" << endl;
    hTrigEff_.push_back((TH1F*) purityTrigEffFile_->Get("h_trig_ec"));
    if (hTrigEff_.back() == nullptr) cout << "***** Histogram h_trig_ec not found *****" << endl;
    // Sam-style trigger efficiency functions:
    fTrigEff_.push_back((TF1*) purityTrigEffFile_->Get("f_trig_eb"));
    if (fTrigEff_.back() == nullptr) cout << "***** Histogram f_trig_eb not found *****" << endl;
    fTrigEff_.push_back((TF1*) purityTrigEffFile_->Get("f_trig_ec"));
    if (fTrigEff_.back() == nullptr) cout << "***** Histogram f_trig_ec not found *****" << endl;
    hSFeff_ = (TH1F*) purityTrigEffFile_->Get("h_SFg_MHT");  // Maybe this should be h_NJets
    if (hSFeff_ == nullptr) cout << "***** Histogram h_MHT not found *****" << endl;
  }

}  // ======================================================================================

double
RA2bZinvAnalysis::efficiencyAndPurity::weight(CCbinning* CCbins, Int_t NJets, Int_t BTags,
					      Double_t MHT, Double_t HT,
					      vector<TLorentzVector> ZCandidates,
					      vector<TLorentzVector> Photons,
					      vector<double> EBphoton) {
  // For double ratio, apply weights for purity, Fdir, trigger eff, reco eff.
  double effWt = 1;

  if ((theSample_.Contains("zmm") || theSample_.Contains("zee")) && !theSample_.Contains("tt")) {
    int binCCjb = CCbins->jb(CCbins->jbin(NJets), CCbins->bbin(NJets, BTags));
    if (hPurity_[0] != nullptr) effWt *= hPurity_[0]->GetBinContent(binCCjb);

  } else if (theSample_.Contains("photon")) {
    if(EBphoton.at(0) == 1 && hPurity_[0] != nullptr) {
      int bin = hPurity_[0]->GetNbinsX();  while (MHT < hPurity_[0]->GetBinLowEdge(bin)) bin--;
      effWt *= hPurity_[0]->GetBinContent(bin);
    }
    if(EBphoton.at(0) == 0 && hPurity_[1] != nullptr) {
      int bin = hPurity_[1]->GetNbinsX();  while (MHT < hPurity_[1]->GetBinLowEdge(bin)) bin--;
      effWt *= hPurity_[1]->GetBinContent(bin);
    }
    int CCbinjk = CCbins->jk(CCbins->jbin(NJets), CCbins->kinBin(HT, MHT));
    if (CCbinjk > 0 && FdirHist_ != nullptr && CCbinjk <= FdirHist_->GetNbinsX()) effWt *= FdirHist_->GetBinContent(CCbinjk);
    // if (CCbinjk > 0 && FdirGraph_ != nullptr && CCbinjk < FdirGraph_->GetN()) effWt *= FdirGraph_->GetY()[CCbinjk-1];
    // FIXME:  For LDP, fill extra bins with value 0.825

  } else if (theSample_.Contains("dymm") || theSample_.Contains("dyee") ||
	     theSample_.Contains("ttmm") || theSample_.Contains("ttee") ||
	     theSample_.Contains("ttzmm") || theSample_.Contains("ttzee") ||
	     theSample_.Contains("VVmm") || theSample_.Contains("VVee")) {
    if (hTrigEff_[0] != nullptr) {
      Double_t zpt = ZCandidates.at(0).Pt();  zpt = max(hTrigEff_[0]->GetBinLowEdge(1), zpt);
      int bin = hTrigEff_[0]->GetNbinsX();  while (zpt < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
      // Double_t ht = HT >= hTrigEff_[0]->GetBinLowEdge(1) ? HT : hTrigEff_[0]->GetBinLowEdge(1);
      // int bin = hTrigEff_[0]->GetNbinsX();  while (ht < hTrigEff_[0]->GetBinLowEdge(bin)) bin--;
      effWt *= hTrigEff_[0]->GetBinContent(bin);
    }
    if (hSFeff_ != nullptr) {
      // Maybe this should be NJets for dyee:
      Double_t mht = MHT >= hSFeff_->GetBinLowEdge(1) ? MHT : hSFeff_->GetBinLowEdge(1);
      int bin = hSFeff_->GetNbinsX();  while (mht < hSFeff_->GetBinLowEdge(bin)) bin--;
      effWt *= hSFeff_->GetBinContent(bin);
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
    if(EBphoton.at(0) == 1 && fTrigEff_[0] != nullptr) {
      effWt *= fTrigEff_[0]->Eval(max(double(205), Photons.at(0).Pt()));  // hard-wired cutoff
    }
    if(EBphoton.at(0) == 0 && fTrigEff_[1] != nullptr) {
      effWt *= fTrigEff_[1]->Eval(max(double(205), Photons.at(0).Pt()));  // hard-wired cutoff
    }
    // effWt /= min(HT, 900.0)*0.00009615+0.9071;  // FIXME:  hard-wired correction
    if (hSFeff_ != nullptr) {
      Double_t mht = MHT >= hSFeff_->GetBinLowEdge(1) ? MHT : hSFeff_->GetBinLowEdge(1);
      int bin = hSFeff_->GetNbinsX();  while (mht < hSFeff_->GetBinLowEdge(bin)) bin--;
      effWt *= hSFeff_->GetBinContent(bin);
    }
  }
  return effWt;
  
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
