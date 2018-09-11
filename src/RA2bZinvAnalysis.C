//
//  Zinv background prediction for RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#include "RA2bZinvAnalysis.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
// #include <TLegend.h>
// #include <TGraphAsymmErrors.h>
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

#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

ClassImp(RA2bZinvAnalysis)

// ======================================================================================

RA2bZinvAnalysis::RA2bZinvAnalysis() {
  Init();
}

RA2bZinvAnalysis::RA2bZinvAnalysis(const std::string& cfg_filename) {
  Init(cfg_filename);
}

void
RA2bZinvAnalysis::Init(const std::string& cfg_filename) {

  // Set up configuration, using boost/program_options.
  po::options_description desc("Config");
  desc.add_options()
    ("verbosity", po::value<int>(&verbosity_))
    ("era", po::value<std::string>(&era_))
    ("tree path", po::value<std::string>(&treeLoc_))
    ("root file index", po::value<std::string>(&fileListsFile_))
    ("delta phi sample", po::value<std::string>(&deltaPhi_), "nominal, hdp, ldp")
    ("integrated luminosity", po::value<double>(&intLumi_))
    ("apply Z mass cut", po::value<bool>(&applyMassCut_))
    ("apply Z Pt cut", po::value<bool>(&applyPtCut_))
    ("apply photon min Delta R cut", po::value<bool>(&applyMinDeltaRCut_))
    ("use analysis bin from tree", po::value<bool>(&useTreeCCbin_))
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

  // Needed branches
  activeBranches_.push_back("NJets");
  activeBranches_.push_back("BTags");
  activeBranches_.push_back("HT");
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
    activeBranches_.push_back("HTclean");
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
  }
  activeBranches_.push_back("Muons");
  activeBranches_.push_back("Electrons");
  activeBranches_.push_back("ZCandidates");
  activeBranches_.push_back("Photons");
  activeBranches_.push_back("Photons_nonPrompt");
  activeBranches_.push_back("Photons_fullID");
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

    kinThresholds_.push_back({300, 300, 500, 1000});  // mht threshold, {ht thresholds}
    kinThresholds_.push_back({350, 350, 500, 1000});
    kinThresholds_.push_back({500, 500, 1000});
    kinThresholds_.push_back({750, 750, 1500});
    kinThresholds_.push_back({250, 300, 500, 1000}); // QCD control bins
    nJetThresholds_ = {2, 3, 5, 7, 9};
    nbThresholds_ = {0, 1, 2, 3};

    Int_t bin = 0;
    for (unsigned j = 0; j < nJetThresholds_.size(); ++j) {
      for (unsigned b = 0; b < nbThresholds_.size(); ++b) {
	if (nbThresholds_[b] > nJetThresholds_[j]) continue;  // Exclude (Njets0, Nb3)
	unsigned mmax = deltaPhi_ == TString("nominal") ? kinThresholds_.size()-1 : kinThresholds_.size();
	int k = -1;
	for (unsigned m = 0; m < mmax; ++m) {
	  for (unsigned h = 1; h < kinThresholds_[m].size(); ++h) {
	    k++;
	    if (j > 2 && (m < 2 || m == kinThresholds_.size()-1) && h == 1) continue;   // Exclude (Njets3,4; HT0,3,(6))
	    std::vector<int> jbk = {int(j), int(b), k};
	    bin++;
	    toCCbin_[jbk] = bin;
	    // cout << "Filling toCCbin; j = " << j << ", b = "  << b << ", k = " << k << ", bin = " << bin << endl;
	  }
	}
      }
    }
  } // era 2016

  kinSize_ = 0;
  int mmax = (deltaPhi_ == TString("nominal")) ? kinThresholds_.size() -1 : kinThresholds_.size();
  for (int i = 0; i < mmax; ++i)
    kinSize_ += kinThresholds_[i].size();

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
  cout << "The apply b-tag scale factors flag is " << applyBTagSF_ << endl;
  cout << "The apply pileup weight flag is " << applyPuWeight_ << endl;
  cout << "The custom pileup weight flag is " << customPuWeight_ << endl;


}  // ======================================================================================

TChain*
RA2bZinvAnalysis::getChain(const char* sample, Int_t* fCurrent, bool makeClass) {
  TString theSample(sample);
  TString key;
  if (theSample.Contains("zinv")) key = TString("zinv");
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
  if (deltaPhi_ == "ldp") key += "ldp";

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
  for (auto theBranch : activeBranches_) chain->AddBranchToCache(theBranch, true);
  // chain->AddBranchToCache("*", true);
  chain->StopCacheLearningPhase();
  cout << "Reset size of cache for chain = " << chain->GetCacheSize() << endl;

  // tmt_->Init(chain);
  setBranchAddress(chain);

  chain->SetBranchStatus("*", 0);  // disable all branches
  for (auto theBranch : activeBranches_) chain->SetBranchStatus(theBranch, 1);

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
    files.push_back(dir+buf);
    buf = buf.Strip();  // Remove any trailing whitespace
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
    if (Ntrig > 1) trigCuts_.Replace(trigCuts_.Length()-3, 3, ")");
  }

  // commonCuts_ = "(JetID==1&& HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && NVtx > 0 && BadPFMuonFilter && PFCaloMETRatio < 5)";  // Troy revision+
  if (trigger.empty()) {
    commonCuts_ = "JetID==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && NVtx > 0";  // Troy revision-
  } else {
    commonCuts_ = "JetID==1 && globalTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0";  // Troy revision-
  }
  // if (!isSkim_) {
  //   commonCuts_ += " && PFCaloMETRatio < 5";  // Applied in skimming
  //   commonCuts_ += " && noMuonJet";  // Defined in loop, applied in skimming
  //   commonCuts_ += " && noFakeJet";  // Defined in loop, applied in skimming FastSim
  // }

  massCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyMassCut_)
    massCut_ = "@ZCandidates.size()==1 && ZCandidates.M()>=76.188 && ZCandidates.M()<=106.188";

  ptCut_ = "1";
  if ((sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll") && applyPtCut_)
    ptCut_ = "@ZCandidates.size()==1 && ZCandidates.Pt()>=200.";

  photonDeltaRcut_ = "1";
  if (sampleKey == "photon") {
    if (applyPtCut_) ptCut_ = "@Photons.size()==1 && Photons->at(0).Pt()>=200.";
    if (trigger.empty() && applyMinDeltaRCut_) photonDeltaRcut_ = "madMinPhotonDeltaR>=0.4";
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
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample, std::vector<hist1D*>& histograms) {
  //
  // Define N - 1 (or N - multiple) cuts, book histograms.  Traverse the chain and fill.
  //
  TCut baselineCuts = getCuts(sample);
  if(verbosity_ >= 1) cout << endl << "baseline = " << endl << baselineCuts << endl << endl;
  Int_t fCurrent;  // current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);
  TObjArray* forNotify = new TObjArray;

  cutHistos cutHistFiller(chain, forNotify);
  for (auto & hg : histograms) {
    hg->hist = new TH1F(hg->name, hg->title, hg->Nbins, hg->lowEdge, hg->highEdge);
    hg->hist->GetXaxis()->SetTitle(hg->axisTitles.first);
    hg->hist->GetYaxis()->SetTitle(hg->axisTitles.second);
    hg->hist->SetOption("HIST");
    hg->hist->SetMarkerSize(0);
    if (hg->name.Contains(TString("hCut"))) {
      hg->NminusOneCuts = "1";
    } else {
      hg->NminusOneCuts = baselineCuts;
      for (auto cutToOmit : hg->omitCuts) hg->NminusOneCuts(*cutToOmit) = "1";
      if (strlen(hg->addCuts) != 0) hg->NminusOneCuts += TString(" && ") + hg->addCuts;
    }
    if (verbosity_ >= 1) {
      cout << "\n For sample " << sample << ", histo " << hg->name  << ", hg->omitCuts = ";
      for (auto cutToOmit : hg->omitCuts) cout << *cutToOmit << " ";
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
      if (thisFile && verbosity_ >= 1) cout << "Current file in chain: " << thisFile->GetName() << endl;
      countInFile = 0;
    }
    chain->GetEntry(entry);
    countInFile++;
    if (countInFile == 1) {
      checkActiveTrigPrescales(sample);
      if (isMC_ && verbosity_ >= 1) cout << "MC weight for this file is " << Weight << endl;
    }
    cleanVars();  // If unskimmed input, copy <var>clean to <var>

    if (ZCandidates->size() > 1 && verbosity_ >= 2) cout << ZCandidates->size() << " Z candidates found" << endl;

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
	cutHistFiller.fill(hg->hist, eventWt);
	continue;
      }
      hg->NminusOneFormula->GetNdata();
      double selWt = hg->NminusOneFormula->EvalInstance(0);
      if (selWt != 0) {
	if (hg->dvalue != nullptr) {
	  // cout << "  " << *(hg->dvalue);
	  hg->hist->Fill(*(hg->dvalue), selWt*eventWt);
	}
	else if (hg->ivalue != nullptr) {
	  // cout << "  " << *(hg->ivalue);
	  hg->hist->Fill(Double_t(*(hg->ivalue)), selWt*eventWt);
	}
	else if (hg->filler != nullptr) (this->*(hg->filler))(hg->hist, selWt*eventWt);
	else cerr << "No method to fill histogram provided for " << hg->name << endl;
      }
    }  // loop over histograms
    // cout << endl;
  }  // loop over entries

  delete forNotify;
  delete chain->GetCurrentFile();

}  // ======================================================================================

std::vector<TH1F*>
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

  std::vector<hist1D*> histograms;

  hist1D hHT;
  hHT.name = TString("hHT_") + TString(sample);  hHT.title = "HT";
  hHT.Nbins = 60;  hHT.lowEdge = 0;  hHT.highEdge = 3000;
  hHT.axisTitles.first = "HT [GeV]";  hHT.axisTitles.second = "Events / 50 GeV";
  hHT.dvalue = &HT;  hHT.omitCuts.push_back(&HTcut_);
  histograms.push_back(&hHT);

  hist1D hMHT;
  hMHT.name = TString("hMHT_") + TString(sample);  hMHT.title = "MHT";
  hMHT.Nbins = 60;  hMHT.lowEdge = 0;  hMHT.highEdge = 3000;
  hMHT.axisTitles.first = "MHT [GeV]";  hMHT.axisTitles.second = "Events / 50 GeV";
  hMHT.dvalue = &MHT;  hMHT.omitCuts.push_back(&MHTcut_);  hMHT.omitCuts.push_back(&ptCut_);
  histograms.push_back(&hMHT);

  hist1D hNJets;
  hNJets.name = TString("hNJets_") + TString(sample);  hNJets.title = "NJets";
  hNJets.Nbins = 20;  hNJets.lowEdge = 0;  hNJets.highEdge = 20;
  hNJets.axisTitles.first = "N (jets)";  hNJets.axisTitles.second = "Events / bin";
  hNJets.ivalue = &NJets;  hNJets.omitCuts.push_back(&NJetscut_);
  histograms.push_back(&hNJets);

  hist1D hBTags;
  hBTags.name = TString("hBTags_") + TString(sample);  hBTags.title = "BTags";
  hBTags.Nbins = 20;  hBTags.lowEdge = 0;  hBTags.highEdge = 20;
  hBTags.axisTitles.first = "N (b jets)";  hBTags.axisTitles.second = "Events / bin";
  hBTags.ivalue = &BTags;
  histograms.push_back(&hBTags);

  hist1D hnZcand;
  hnZcand.name = TString("hnZcand_") + TString(sample);  hnZcand.title = "Number of Z candidates";
  hnZcand.Nbins = 10;  hnZcand.lowEdge = 0;  hnZcand.highEdge = 10;
  hnZcand.axisTitles.first = "N(Z candidates)";  hnZcand.axisTitles.second = "Events / bin";
  hnZcand.filler = &RA2bZinvAnalysis::fillnZcand;  hnZcand.omitCuts.push_back(&massCut_);
  histograms.push_back(&hnZcand);

  hist1D hZmass;
  hZmass.name = TString("hZmass_") + TString(sample);  hZmass.title = "Z mass";
  hZmass.Nbins = 30;  hZmass.lowEdge = 60;  hZmass.highEdge = 120;
  hZmass.axisTitles.first = "M(Z) [GeV]";  hZmass.axisTitles.second = "Events / 2 GeV";
  hZmass.filler = &RA2bZinvAnalysis::fillZmass;  hZmass.omitCuts.push_back(&massCut_);
  histograms.push_back(&hZmass);

  hist1D hZpt;
  hZpt.name = TString("hZpt_") + TString(sample);  hZpt.title = "Z Pt";
  hZpt.Nbins = 60;  hZpt.lowEdge = 0;  hZpt.highEdge = 3000;
  hZpt.axisTitles.first = "Pt(Z) [GeV]";  hZpt.axisTitles.second = "Events / 50 GeV";
  hZpt.filler = &RA2bZinvAnalysis::fillZpt;  hZpt.omitCuts.push_back(&ptCut_);  hZpt.omitCuts.push_back(&MHTcut_);
  histograms.push_back(&hZpt);

  hist1D hCutFlow;
  hCutFlow.name = TString("hCutFlow_") + TString(sample);  hCutFlow.title = "Cut flow";
  hCutFlow.Nbins = 10;  hCutFlow.lowEdge = 0;  hCutFlow.highEdge = 10;
  hCutFlow.axisTitles.first = "";  hCutFlow.axisTitles.second = "Events surviving";
  histograms.push_back(&hCutFlow);

  hist1D hCuts;
  hCuts.name = TString("hCuts_") + TString(sample);  hCuts.title = "Cuts passed";
  hCuts.Nbins = 10;  hCuts.lowEdge = 0;  hCuts.highEdge = 10;
  hCuts.axisTitles.first = "";  hCuts.axisTitles.second = "Events passing";
  histograms.push_back(&hCuts);

  hist1D hVertices;
  hVertices.name = TString("hVertices_") + TString(sample);  hVertices.title = "Number of reco vertices";
  hVertices.Nbins = 100;  hVertices.lowEdge = 0;  hVertices.highEdge = 100;
  hVertices.axisTitles.first = "No. of vertices";  hVertices.axisTitles.second = "Events / bin";
  hVertices.ivalue = &nAllVertices;
  histograms.push_back(&hVertices);

  hist1D hTrueNumInt;
  hTrueNumInt.name = TString("hTrueNumInt_") + TString(sample);  hTrueNumInt.title = "Number of generated interactions";
  hTrueNumInt.Nbins = 100;  hTrueNumInt.lowEdge = 0;  hTrueNumInt.highEdge = 100;
  hTrueNumInt.axisTitles.first = "No. of interactions";  hTrueNumInt.axisTitles.second = "Events / bin";
  hTrueNumInt.dvalue = &TrueNumInteractions;
  histograms.push_back(&hTrueNumInt);

  // Z mass in Njet, Nb bins
  hist1D hZmass_2j0b(hZmass);
  hZmass_2j0b.name = TString("hZmass_2j0b_") + TString(sample);  hZmass_2j0b.title = "Z mass, 2 jets & 0 b jets";
  hZmass_2j0b.addCuts = isSkim_ ? "NJets==2 && BTags==0" : "NJetsclean==2 && BTagsclean==0";
  histograms.push_back(&hZmass_2j0b);

  hist1D hZmass_2j1b(hZmass);
  hZmass_2j1b.name = TString("hZmass_2j1b_") + TString(sample);  hZmass_2j1b.title = "Z mass, 2 jets & 1 b jet";
  hZmass_2j1b.addCuts = isSkim_ ? "NJets==2 && BTags==1" : "NJetsclean==2 && BTagsclean==1";
  histograms.push_back(&hZmass_2j1b);

  hist1D hZmass_2j2b(hZmass);
  hZmass_2j2b.name = TString("hZmass_2j2b_") + TString(sample);  hZmass_2j2b.title = "Z mass, 2 jets & >=2 b jets";
  hZmass_2j2b.addCuts = isSkim_ ? "NJets==2 && BTags>=2" : "NJetsclean==2 && BTagsclean>=2";
  histograms.push_back(&hZmass_2j2b);
  //
  hist1D hZmass_3j0b(hZmass);
  hZmass_3j0b.name = TString("hZmass_3j0b_") + TString(sample);  hZmass_3j0b.title = "Z mass, 3-4 jets & 0 b jets";
  hZmass_3j0b.addCuts = isSkim_ ? "NJets>=3 && NJets<=4 && BTags==0" : "NJetsclean>=3 && NJetsclean<=4 && BTagsclean==0";
  histograms.push_back(&hZmass_3j0b);

  hist1D hZmass_3j1b(hZmass);
  hZmass_3j1b.name = TString("hZmass_3j1b_") + TString(sample);  hZmass_3j1b.title = "Z mass, 3-4 jets & 1 b jet";
  hZmass_3j1b.addCuts = isSkim_ ? "NJets>=3 && NJets<=4 && BTags==1" : "NJetsclean>=3 && NJetsclean<=4 && BTagsclean==1";
  histograms.push_back(&hZmass_3j1b);

  hist1D hZmass_3j2b(hZmass);
  hZmass_3j2b.name = TString("hZmass_3j2b_") + TString(sample);  hZmass_3j2b.title = "Z mass, 3-4 jets & >=2 b jets";
  hZmass_3j2b.addCuts = isSkim_ ? "NJets>=3 && NJets<=4 && BTags>=2" : "NJetsclean>=3 && NJetsclean<=4 && BTagsclean>=2";
  histograms.push_back(&hZmass_3j2b);
  //
  hist1D hZmass_5j0b(hZmass);
  hZmass_5j0b.name = TString("hZmass_5j0b_") + TString(sample);  hZmass_5j0b.title = "Z mass, >=5 jets & 0 B jets";
  hZmass_5j0b.addCuts = isSkim_ ? "NJets>=5 && BTags==0" : "NJetsclean>=5 && BTagsclean==0";
  histograms.push_back(&hZmass_5j0b);

  hist1D hZmass_5j1b(hZmass);
  hZmass_5j1b.name = TString("hZmass_5j1b_") + TString(sample);  hZmass_5j1b.title = "Z mass, >=5 jets & 1 B jet";
  hZmass_5j1b.addCuts = isSkim_ ? "NJets>=5 && BTags==1" : "NJetsclean>=5 && BTagsclean==1";
  histograms.push_back(&hZmass_5j1b);

  hist1D hZmass_5j2b(hZmass);
  hZmass_5j2b.name = TString("hZmass_5j2b_") + TString(sample);  hZmass_5j2b.title = "Z mass, >=5 jets & >=2 B jets";
  hZmass_5j2b.addCuts = isSkim_ ? "NJets>=5 && BTags>=2" : "NJetsclean>=5 && BTagsclean>=2";
  histograms.push_back(&hZmass_5j2b);
  bookAndFillHistograms(sample, histograms);

  std::vector<TH1F*> theHists;
  theHists.push_back(hHT.hist);
  theHists.push_back(hMHT.hist);
  theHists.push_back(hNJets.hist);
  theHists.push_back(hBTags.hist);
  theHists.push_back(hnZcand.hist);
  theHists.push_back(hZmass.hist);
  theHists.push_back(hZpt.hist);
  theHists.push_back(hVertices.hist);
  theHists.push_back(hTrueNumInt.hist);
  theHists.push_back(hZmass_2j0b.hist);
  theHists.push_back(hZmass_2j1b.hist);
  theHists.push_back(hZmass_2j2b.hist);
  theHists.push_back(hZmass_3j0b.hist);
  theHists.push_back(hZmass_3j1b.hist);
  theHists.push_back(hZmass_3j2b.hist);
  theHists.push_back(hZmass_5j0b.hist);
  theHists.push_back(hZmass_5j1b.hist);
  theHists.push_back(hZmass_5j2b.hist);

  return theHists;

}  // ======================================================================================

TH1F*
RA2bZinvAnalysis::makeCChist(const char* sample) {
  //
  // Book the analysis-binned histogram, fill, and return it.
  //
  Int_t MaxBins = toCCbin_.size();
  cout << "MaxBins = " << MaxBins << endl;

  TH1F* hCCbins = new TH1F("hCCbins", "Zinv background estimate", MaxBins, 0.5, MaxBins+0.5);

  TObjArray* forNotify = new TObjArray;  // Allow for more than one TObject to notify of a new file
  
  Int_t fCurrent; //!current Tree number in a TChain
  TChain* chain = getChain(sample, &fCurrent);

  // chain->Print();

  BTagCorrector* btagcorr = nullptr;
  if (applyBTagSF_) {
    btagcorr = new BTagCorrector;
    btagcorr->SetCalib(BTagSFfile_);
  }

  // Get the baseline cuts, make a TreeFormula, and add it to the list
  // of objects to be notified as new files in the chain are encountered
  TCut baselineCuts = getCuts(sample);
  if (verbosity_ >= 1) cout << "baseline = " << baselineCuts << endl;
  // See https://root-forum.cern.ch/t/how-to-evaluate-a-formula-for-a-single-tree-event/12283
  TTreeFormula* baselineTF = new TTreeFormula("baselineTF", baselineCuts, chain);
  forNotify->Add(baselineTF);
  chain->SetNotify(forNotify);

  // Traverse the events in the chain
  Long64_t Nentries = chain->GetEntries();
  if (verbosity_ >= 1) cout << "Nentries in tree = " << Nentries << endl;
  int count = 0, outCount = 0, countInFile = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (verbosity_ >= 1 && (count % 100000 == 0)) cout << "Entry number " << count << endl;
    chain->LoadTree(entry);
    if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      TFile* thisFile = chain->GetCurrentFile();
      if (thisFile) {
    	if (verbosity_ >= 1) cout << "Current file in chain: " << thisFile->GetName() << endl;
    	if (btagcorr) btagcorr->SetEffs(thisFile);
      }
      countInFile = 0;
    }
    chain->GetEntry(entry);
    countInFile++;
    if (countInFile == 1) {
      checkActiveTrigPrescales(sample);
      if (isMC_ && verbosity_ >= 1) cout << "MC weight for this file is " << Weight << endl;
    }
    cleanVars();  // If unskimmed input, copy <var>clean to <var>

    UInt_t binCC = 0;

    // Apply baseline selection
    baselineTF->GetNdata();
    double selWt = baselineTF->EvalInstance(0);
    if (selWt == 0) continue;
    outCount++;

    // Compute event weight factors
    Double_t eventWt = 1;
    Double_t PUweight = 1;
    if (applyPuWeight_ && customPuWeight_) {
      // This recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
      PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(TrueNumInteractions, puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
    }
    if (isMC_) {
      eventWt = 1000*intLumi_*Weight*selWt*PUweight;
      if (eventWt < 0) eventWt *= -1;
    }

    if (useTreeCCbin_ && !applyBTagSF_) {
      binCC = RA2bin;
      hCCbins->Fill(Double_t(binCC), eventWt);
    } else {
      // Calculate binCC
      std::vector<int> jbk;
      int binKin = kinBin(HT, MHT);
      if (binKin < 0) continue;
      int binNjets = nJetThresholds_.size()-1;
      while (NJets < nJetThresholds_[binNjets]) binNjets--;
      if (!applyBTagSF_) {
	int binNb = nbThresholds_.size()-1;
	while (BTags < nbThresholds_[binNb]) binNb--;
	jbk = {binNjets, binNb, binKin};
	try {
	  binCC = toCCbin_.at(jbk);
	}
	catch (const std::out_of_range& oor) {
          // Omitted bins j = 3,4, k = 0,3
	  // std::cerr << "jpk out of range: " << oor.what() << ": j = " << binNjets
	  // 	    << ", b = " << binNb << ", k = " << binKin << ", RA2bin = " << RA2bin << '\n';
	  continue;
	}
	// if (outCount < 100) cout << "j = " << binNjets << ", b = " << binNb << ", k = " << binKin << ", binCC = " << binCC << ", RA2bin = " << RA2bin << endl;
	hCCbins->Fill(Double_t(binCC), eventWt);
      } else {
        // apply BTagSF to all Nb bins
	// if (count < 20 || count % 10000 == 0) cout << "Size of input Jets = " << Jets->size() << ", Jets_hadronFlavor = " << Jets_hadronFlavor->size() << " Jets_HTMask = " << Jets_HTMask->size() << endl;
	vector<double> probNb = btagcorr->GetCorrections(Jets, Jets_hadronFlavor, Jets_HTMask);
	for (int binNb = 0; binNb < (int) nbThresholds_.size(); ++binNb) {
	  jbk = {binNjets, binNb, binKin};
	  try {binCC = toCCbin_.at(jbk);}
	  catch (const std::out_of_range& oor) {continue;}
	  if (verbosity_ >= 2 && count % 100000 == 0) cout << "j = " << binNjets << ", NbTags = " << BTags << ", b = " << binNb << ", b wt = " << probNb[binNb] << ", k = " << binKin << ", binCC = " << binCC << endl;
	  hCCbins->Fill(Double_t(binCC), eventWt*probNb[binNb]);
	}
      }  // if apply BTagSF
    }  // if useTreeCCbin
  }  // End loop over entries

  delete forNotify;
  delete baselineTF;
  delete chain->GetCurrentFile();
  if (btagcorr) delete btagcorr;

  return hCCbins;

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
      objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoMuonTracks==0";
      // objCutMap_["zmm"] = "@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0;  // Troy mod+
      objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0 && (@Photons.size()==0) && isoElectronTracks==0";
      // objCutMap_["zee"] = "@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0;  // Troy mod+
      objCutMap_["zll"] = "((@Muons.size()==2 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0) || (@Muons.size()==0 && @Electrons.size()==2 && isoMuonTracks==0 && isoPionTracks==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && Photons->at(0).Pt()>=200 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";
      // objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && @Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0";  // Troy mod+
      objCutMap_["ttz"] = "@Muons.size()==0 && @Electrons.size()==0 && isoElectronTracks==0 && isoMuonTracks==0 && isoPionTracks==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "@Muons.size()==1 && @Electrons.size()==0 && isoElectronTracks==0 && isoPionTracks==0";
      objCutMap_["sle"] = "@Muons.size()==0 && @Electrons.size()==1 && isoMuonTracks==0 && isoPionTracks==0";
    } else if (ntupleVersion_ == "V15") {
    }

    minDphiCutMap_["nominal"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
    minDphiCutMap_["hdp"] = "DeltaPhi1>0.5 && DeltaPhi2>0.5 && DeltaPhi3>0.3 && DeltaPhi4>0.3";
    minDphiCutMap_["ldp"] = "(DeltaPhi1<0.5 || DeltaPhi2<0.5 || DeltaPhi3<0.3 || DeltaPhi4<0.3)";

    MHTCutMap_["nominal"] = "MHT>=300";
    MHTCutMap_["hdp"] = "MHT>=250";
    MHTCutMap_["ldp"] = "MHT>=250";

  } else {

    if (ntupleVersion_ == "V12") {
    } else if (ntupleVersion_ == "V15") {
      objCutMap_["sig"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["zmm"] = "NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoMuonTracksclean==0";
      objCutMap_["zee"] = "NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@Photons.size()==0) && isoElectronTracksclean==0";
      objCutMap_["zll"] = "((NMuons==2 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0) || (NMuons==0 && NElectrons==2 && isoMuonTracksclean==0 && isoPionTracksclean==0))";
      objCutMap_["photon"] = "Sum$(Photons_nonPrompt)==0 && Sum$(Photons_fullID)==1 && (@Photons.size()==1) && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["photonqcd"] = "Sum$(Photons_nonPrompt)!=0 && Photons->at(0).Pt()>=200 && NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["ttz"] = "NMuons==0 && NElectrons==0 && isoElectronTracksclean==0 && isoMuonTracksclean==0 && isoPionTracksclean==0 && (@GenMuons.size()==0 && @GenElectrons.size()==0 && @GenTaus.size()==0)";
      objCutMap_["slm"] = "NMuons==1 && NElectrons==0 && isoElectronTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["slm"] += " && TransverseMass(METPt,METPhi,Muons->at(0).Pt(),Muons->at(0).Phi()) < 100";  // add'l skim cut
      objCutMap_["sle"] = "NMuons==0 && NElectrons==1 && isoMuonTracksclean==0 && isoPionTracksclean==0";
      objCutMap_["sle"] += " && TransverseMass(METPt,METPhi,Electrons->at(0).Pt(),Electrons->at(0).Phi()) < 100";  // add'l skim cut
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
  hcf->GetXaxis()->SetBinLabel(1, "None");
  hcf->GetXaxis()->SetBinLabel(2, "HT");
  hcf->GetXaxis()->SetBinLabel(3, "MHT");
  hcf->GetXaxis()->SetBinLabel(4, "NJets");
  hcf->GetXaxis()->SetBinLabel(5, "mnDphi");
  hcf->GetXaxis()->SetBinLabel(6, "objects");
  hcf->GetXaxis()->SetBinLabel(7, "Zpt");
  hcf->GetXaxis()->SetBinLabel(8, "Zmass");
  hcf->GetXaxis()->SetBinLabel(9, "Trigger");
  hcf->GetXaxis()->SetBinLabel(10, "Filters");
  hcf->GetXaxis()->LabelsOption("vu");

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
