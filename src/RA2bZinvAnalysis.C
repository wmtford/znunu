//
//  Make histograms for Zinv background prediction in the RA2b analysis
//  Loosely based on Troy Mulholland's python code
//

#include "RA2bZinvAnalysis.h"
#include <TStyle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTreeCache.h>
#include "../../Analysis/KCode/KMath.h"

#include <iostream>
using std::cout;
using std::endl;

#include <sstream>
using std::stringstream;

#include <stdio.h>
// using std::sprintf;

#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// ClassImp(RA2bZinvAnalysis)

// ======================================================================================

RA2bZinvAnalysis::RA2bZinvAnalysis() {
  Config();
}

RA2bZinvAnalysis::RA2bZinvAnalysis(const string& cfg_filename, const string& runBlock) {
  Config(cfg_filename);
  if (!runBlock.empty()) runBlock_ = runBlock;  // Constructor arg overrides config
  cout << "The runBlock is " << runBlock_ << endl;
}

void
RA2bZinvAnalysis::Config(const string& cfg_filename) {

  string treeName;
  string treeLoc;
  string fileListsFile;

  // Set up configuration, using boost/program_options.
  po::options_description desc("Config");
  desc.add_options()
    ("verbosity", po::value<int>(&verbosity_))
    ("root verbosity", po::value<string>(&rootVerbosity_))  // Print, Info, Warning, Error, Break, SysError, Fatal
    ("era", po::value<string>(&era_))
    ("runBlock", po::value<string>(&runBlock_))  // May be overridden by constructor
    ("tree path", po::value<string>(&treeLoc))
    ("root file index", po::value<string>(&fileListsFile))
    ("delta phi sample", po::value<string>(&deltaPhi_), "nominal, hdp, ldp, ldpnominal")
    ("integrated luminosity", po::value<double>(&intLumi_))
    ("apply Z mass cut", po::value<bool>(&applyMassCut_))
    ("apply Z/gamma Pt cut", po::value<bool>(&applyPtCut_))
    ("min Nb cut", po::value<int>(&minNbCut_))
    ("use DeepCSV", po::value<bool>(&useDeepCSV_))
    ("apply b-tag SF", po::value<bool>(&applyBTagSF_))
    ("apply pileup weight", po::value<bool>(&applyPuWeight_))
    ("use custom pileup weight", po::value<bool>(&customPuWeight_))  // Substitute Kevin P recipe
    ("apply HEM jet veto", po::value<bool>(&applyHEMjetVeto_))
    ("apply Z Pt weight", po::value<bool>(&applyZptWt_))
    ("apply double-ratio fit weight", po::value<bool>(&applyDRfitWt_))
    ("apply scale factors to MC", po::value<bool>(&applySFwtToMC_))
    ("restrict cleanVars for raw ntuple", po::value<bool>(&restrictClean_))

    ;
  po::variables_map vm;
  std::ifstream cfi_file("RA2bZinvAnalysis.cfi");
  po::store(po::parse_config_file(cfi_file , desc, false), vm);  // false means throw exception for unrecognized option
  po::notify(vm);
  if (!cfg_filename.empty()) {
    std::ifstream cfg_file(cfg_filename);
    vm = po::variables_map();  // Clear the map.
    po::store(po::parse_config_file(cfg_file , desc, false), vm);
    po::notify(vm);
  }

  if      (rootVerbosity_ == "Print") gErrorIgnoreLevel = kPrint;
  else if (rootVerbosity_ == "Info") gErrorIgnoreLevel = kInfo;
  else if (rootVerbosity_ == "Warning") gErrorIgnoreLevel = kWarning;
  else if (rootVerbosity_ == "Error") gErrorIgnoreLevel = kError;
  else if (rootVerbosity_ == "Break") gErrorIgnoreLevel = kBreak;
  else if (rootVerbosity_ == "SysError") gErrorIgnoreLevel = kSysError;
  else if (rootVerbosity_ == "Fatal") gErrorIgnoreLevel = kFatal;

  //  Extract tree signature from the path
  string s = treeLoc, delimiter = "Run2Production";
  size_t pos = s.find(delimiter);
  if (pos == string::npos) cout << "Can't find tree version; " << delimiter << " not in tree path" << endl;
  pos += delimiter.length();
  ntupleVersion_ = s.substr(pos, 3);

  isSkim_ = treeLoc.find("Skim") != string::npos;
  if (!isSkim_) {
    treeName = "TreeMaker2/PreSelection";  // For ntuple
  } else {
    treeName = "tree";  // For skims
  }
  if (ntupleVersion_ == "V12") useDeepCSV_ = false;

  treeConfig_ = new TreeConfig(era_, ntupleVersion_, isSkim_, deltaPhi_, verbosity_, treeName,
                               treeLoc, fileListsFile, runBlock_);
  CCbins_ = new CCbinning(era_, deltaPhi_);
  effPurCorr_ = new EfficWt(deltaPhi_);

  cout << "After initialization," << endl;
  cout << "The verbosity level is " << verbosity_ << endl;
  cout << "The root verbosity level is " << rootVerbosity_ << endl;
  cout << "The ntuple version is " << ntupleVersion_ << endl;
  cout << "The input-files-are-skims flag is " << isSkim_ << endl;
  cout << "The era is " << era_ << endl;
  cout << "The integrated luminosity = " << intLumi_ << endl;
  cout << "The path to input files is " << treeLoc << endl;
  cout << "The minDeltaPhi cuts are " << deltaPhi_ << endl;
  cout << "Apply Z/gamma Pt cut is " << applyPtCut_ << endl;
  cout << "The minimum Nb is " << minNbCut_ << endl;
  cout << "Use DeepCSV is " << useDeepCSV_ << endl;
  cout << "Apply HEM jet veto for 2018HEM is " << applyHEMjetVeto_ << endl;
  cout << "Restrict cleanVars for raw ntuple is " << restrictClean_ << endl;

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

  double cpu0  = get_cpu_time();

  TString key = treeConfig_->DSkey(sample, isMC_);  // Set isMC_ here
  if (isMC_) {
    cout << "For MC data set " << sample << "," << endl;
    cout << "Apply b-tag scale factors is " << applyBTagSF_ << endl;
    cout << "Apply pileup weight is " << applyPuWeight_ << endl;
    cout << "Use custom if pileup weight is " << customPuWeight_ << endl;
    cout << "Apply Z Pt weight for 2017, 18 Z MC is " << applyZptWt_ << endl;
    cout << "Apply double-ratio fit weight to gJets is " << applyDRfitWt_ << endl;
    cout << "Apply scale factors to MC for non-DR histograms is " << applySFwtToMC_ << endl;
  }

  bool activateAll = false;
  if (activateAll) {
    cout << "Activating all branches" << endl;
    Tupl = new Ntuple(treeConfig_->getChain(key));  // Add arg = false to suppress caching
  } else {
    Tupl = new Ntuple(treeConfig_->getChain(key), treeConfig_->activeBranchList());
  }

  evSelector_ = new CutManager(sample, ntupleVersion_, isSkim_, isMC_, verbosity_,
                               era_, deltaPhi_, applyMassCut_, applyPtCut_,
			       minNbCut_, restrictClean_, CCbins_);
  CutManager::string_map sampleMap = evSelector_->sampleKeyMap();
  TString sampleKey = sampleMap.count(sample) > 0 ? sampleMap.at(sample) : "none";
  bool isZll = (sampleKey == "zmm" || sampleKey == "zee" || sampleKey == "zll");
  bool isPhoton = (sampleKey == "photon" || sampleKey =="photonqcd");

  std::vector<histConfig*> histograms, cutHistograms;

  histConfig hCutFlow;
  hCutFlow.name = TString("hCutFlow_") + sample;  hCutFlow.title = "Cut flow unweighted";
  hCutFlow.NbinsX = 12;  hCutFlow.rangeX.first = 0;  hCutFlow.rangeX.second = 12;
  hCutFlow.axisTitles.first = "";  hCutFlow.axisTitles.second = "Events surviving";
  cutHistograms.push_back(&hCutFlow);

  histConfig hCutFlowWt;
  hCutFlowWt.name = TString("hCutFlowWt_") + sample;  hCutFlowWt.title = "Cut flow with weights";
  hCutFlowWt.NbinsX = 12;  hCutFlowWt.rangeX.first = 0;  hCutFlowWt.rangeX.second = 12;
  hCutFlowWt.axisTitles.first = "";  hCutFlowWt.axisTitles.second = "Weighted events surviving";
  cutHistograms.push_back(&hCutFlowWt);

  histConfig hCutFlowSkim;
  hCutFlowSkim.name = TString("hCutFlowSkim_") + sample;  hCutFlowSkim.title = "Skim cut flow unweighted";
  hCutFlowSkim.NbinsX = 12;  hCutFlowSkim.rangeX.first = 0;  hCutFlowSkim.rangeX.second = 12;
  hCutFlowSkim.axisTitles.first = "";  hCutFlowSkim.axisTitles.second = "Events surviving";
  if (!isSkim_) cutHistograms.push_back(&hCutFlowSkim);

  histConfig hCuts;
  hCuts.name = TString("hCuts_") + sample;  hCuts.title = "Cuts passed";
  hCuts.NbinsX = 12;  hCuts.rangeX.first = 0;  hCuts.rangeX.second = 12;
  hCuts.axisTitles.first = "";  hCuts.axisTitles.second = "Events passing";
  cutHistograms.push_back(&hCuts);
  // These hCut* histograms are filled before the trigger cut.

  histConfig hFilterCutsNoWt;
  hFilterCutsNoWt.name = TString("hFilterCutsNoWt_") + sample;  hFilterCutsNoWt.title = "Filter cuts, unweighted";
  hFilterCutsNoWt.NbinsX = 16;  hFilterCutsNoWt.rangeX.first = 0;  hFilterCutsNoWt.rangeX.second = 16;
  hFilterCutsNoWt.axisTitles.first = "";  hFilterCutsNoWt.axisTitles.second = "Events failing";
  hFilterCutsNoWt.filler1D = &RA2bZinvAnalysis::fillFilterCuts;
  hFilterCutsNoWt.omitCuts.push_back(&(evSelector_->commonCuts()));
  histograms.push_back(&hFilterCutsNoWt);

  histConfig hFilterCuts;
  hFilterCuts.name = TString("hFilterCuts_") + sample;  hFilterCuts.title = "Filter cuts";
  hFilterCuts.NbinsX = 16;  hFilterCuts.rangeX.first = 0;  hFilterCuts.rangeX.second = 16;
  hFilterCuts.axisTitles.first = "";  hFilterCuts.axisTitles.second = "Events failing";
  hFilterCuts.filler1D = &RA2bZinvAnalysis::fillFilterCuts;
  hFilterCuts.omitCuts.push_back(&(evSelector_->commonCuts()));
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

  Int_t MaxBinsjk = CCbins_->binsjk();  // for the next 6 histograms
  cout << "MaxBinsjk = " << MaxBinsjk << endl;
  histConfig hCCjk;
  hCCjk.name = TString("hCCjk_") + sample;  hCCjk.title = "Cut & count, Nb = 0";
  hCCjk.NbinsX = MaxBinsjk;  hCCjk.rangeX.first = 0.5;  hCCjk.rangeX.second = MaxBinsjk+0.5;
  hCCjk.axisTitles.first = "Njets, (HT, MHT)";  hCCjk.axisTitles.second = "Events / bin";
  hCCjk.filler1D = &RA2bZinvAnalysis::fillCC;
  histograms.push_back(&hCCjk);

  histConfig hCCjkcc;
  hCCjkcc.name = TString("hCCjkcc_") + sample;  hCCjkcc.title = "Cut & count, Nb > 0, central SF";
  hCCjkcc.NbinsX = MaxBinsjk;  hCCjkcc.rangeX.first = 0.5;  hCCjkcc.rangeX.second = MaxBinsjk+0.5;
  hCCjkcc.axisTitles.first = "Njets, (HT, MHT)";  hCCjkcc.axisTitles.second = "Events / bin";
  hCCjkcc.filler1D = &RA2bZinvAnalysis::fillCC;
  if (isMC_ && applyBTagSF_ && minNbCut_ > 0) histograms.push_back(&hCCjkcc);

  histConfig hCCjkuc;
  hCCjkuc.name = TString("hCCjkuc_") + sample;  hCCjkuc.title = "Cut & count, Nb > 0, b eff SF up";
  hCCjkuc.NbinsX = MaxBinsjk;  hCCjkuc.rangeX.first = 0.5;  hCCjkuc.rangeX.second = MaxBinsjk+0.5;
  hCCjkuc.axisTitles.first = "Njets, (HT, MHT)";  hCCjkuc.axisTitles.second = "Events / bin";
  hCCjkuc.filler1D = &RA2bZinvAnalysis::fillCC;
  if (isMC_ && applyBTagSF_ && minNbCut_ > 0) histograms.push_back(&hCCjkuc);

  histConfig hCCjkcu;
  hCCjkcu.name = TString("hCCjkcu_") + sample;  hCCjkcu.title = "Cut & count, Nb > 0, b mistag SF up";
  hCCjkcu.NbinsX = MaxBinsjk;  hCCjkcu.rangeX.first = 0.5;  hCCjkcu.rangeX.second = MaxBinsjk+0.5;
  hCCjkcu.axisTitles.first = "Njets, (HT, MHT)";  hCCjkcu.axisTitles.second = "Events / bin";
  hCCjkcu.filler1D = &RA2bZinvAnalysis::fillCC;
  if (isMC_ && applyBTagSF_ && minNbCut_ > 0) histograms.push_back(&hCCjkcu);

  histConfig hCCjkdc;
  hCCjkdc.name = TString("hCCjkdc_") + sample;  hCCjkdc.title = "Cut & count, Nb > 0, b eff SF down";
  hCCjkdc.NbinsX = MaxBinsjk;  hCCjkdc.rangeX.first = 0.5;  hCCjkdc.rangeX.second = MaxBinsjk+0.5;
  hCCjkdc.axisTitles.first = "Njets, (HT, MHT)";  hCCjkdc.axisTitles.second = "Events / bin";
  hCCjkdc.filler1D = &RA2bZinvAnalysis::fillCC;
  if (isMC_ && applyBTagSF_ && minNbCut_ > 0) histograms.push_back(&hCCjkdc);

  histConfig hCCjkcd;
  hCCjkcd.name = TString("hCCjkcd_") + sample;  hCCjkcd.title = "Cut & count, Nb > 0, b mistag SF down";
  hCCjkcd.NbinsX = MaxBinsjk;  hCCjkcd.rangeX.first = 0.5;  hCCjkcd.rangeX.second = MaxBinsjk+0.5;
  hCCjkcd.axisTitles.first = "Njets, (HT, MHT)";  hCCjkcd.axisTitles.second = "Events / bin";
  hCCjkcd.filler1D = &RA2bZinvAnalysis::fillCC;
  if (isMC_ && applyBTagSF_ && minNbCut_ > 0) histograms.push_back(&hCCjkcd);

  histConfig hgenHT;
  hgenHT.name = TString("hgenHT_") + sample;  hgenHT.title = "Generated HT";
  hgenHT.NbinsX = 60;  hgenHT.rangeX.first = 0;  hgenHT.rangeX.second = 3000;
  hgenHT.axisTitles.first = "HT [GeV]";  hgenHT.axisTitles.second = "Events / 50 GeV";
  hgenHT.dvalue = &(Tupl->GenHT);
  if (isMC_) histograms.push_back(&hgenHT);

  histConfig hPrefire;
  hPrefire.name = TString("hPrefire_") + sample;  hPrefire.title = "Probability to survive trigger prefire";
  hPrefire.NbinsX = 120;  hPrefire.rangeX.first = 0;  hPrefire.rangeX.second = 1.2;
  hPrefire.axisTitles.first = "NonPrefireProb";  hPrefire.axisTitles.second = "Events / 0.01";
  hPrefire.dvalue = &(Tupl->NonPrefiringProb);
  if (isMC_) histograms.push_back(&hPrefire);

  histConfig hHT;
  hHT.name = TString("hHT_") + sample;  hHT.title = "HT";
  hHT.NbinsX = 60;  hHT.rangeX.first = 0;  hHT.rangeX.second = 3000;
  hHT.axisTitles.first = "HT [GeV]";  hHT.axisTitles.second = "Events / 50 GeV";
  hHT.dvalue = &(Tupl->HT);  hHT.omitCuts.push_back(&(evSelector_->HTcut()));
  histograms.push_back(&hHT);

  histConfig hMHT;
  hMHT.name = TString("hMHT_") + sample;  hMHT.title = "MHT";
  hMHT.NbinsX = 60;  hMHT.rangeX.first = 0;  hMHT.rangeX.second = 3000;
  hMHT.axisTitles.first = "MHT [GeV]";  hMHT.axisTitles.second = "Events / 50 GeV";
  hMHT.dvalue = &(Tupl->MHT);
  hMHT.omitCuts.push_back(&(evSelector_->MHTcut()));  hMHT.omitCuts.push_back(&(evSelector_->ptCut()));
  histograms.push_back(&hMHT);

  histConfig hNJets;
  hNJets.name = TString("hNJets_") + sample;  hNJets.title = "NJets";
  hNJets.NbinsX = 20;  hNJets.rangeX.first = 0;  hNJets.rangeX.second = 20;
  hNJets.axisTitles.first = "N (jets)";  hNJets.axisTitles.second = "Events / bin";
  hNJets.ivalue = &(Tupl->NJets);  hNJets.omitCuts.push_back(&(evSelector_->NJetscut()));
  histograms.push_back(&hNJets);

  histConfig hBTags;
  hBTags.name = TString("hBTags_") + sample;  hBTags.title = "BTags";
  hBTags.NbinsX = 20;  hBTags.rangeX.first = 0;  hBTags.rangeX.second = 20;
  hBTags.axisTitles.first = "N (b jets)";  hBTags.axisTitles.second = "Events / bin";
  hBTags.ivalue = &(Tupl->BTags);  hBTags.omitCuts.push_back(&(evSelector_->Nbcut()));
  histograms.push_back(&hBTags);

  histConfig hVertices;
  hVertices.name = TString("hVertices_") + sample;  hVertices.title = "Number of reco vertices";
  hVertices.NbinsX = 100;  hVertices.rangeX.first = 0;  hVertices.rangeX.second = 100;
  hVertices.axisTitles.first = "No. of vertices";  hVertices.axisTitles.second = "Events / bin";
  hVertices.ivalue = &(Tupl->nAllVertices);
  histograms.push_back(&hVertices);

  histConfig hTrueNumInt;
  hTrueNumInt.name = TString("hTrueNumInt_") + sample;  hTrueNumInt.title = "Number of generated interactions";
  hTrueNumInt.NbinsX = 100;  hTrueNumInt.rangeX.first = 0;  hTrueNumInt.rangeX.second = 100;
  hTrueNumInt.axisTitles.first = "No. of interactions";  hTrueNumInt.axisTitles.second = "Events / bin";
  hTrueNumInt.dvalue = &(Tupl->TrueNumInteractions);
  if (isMC_) histograms.push_back(&hTrueNumInt);

  histConfig hPUwtvsNint;
  hPUwtvsNint.name = TString("hPUwtvsNint_") + sample;
  hPUwtvsNint.title = "PU wt vs Number of generated interactions";
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
  hHT_DR.dvalue = &(Tupl->HT);
  hHT_DR.omitCuts.push_back(&(evSelector_->Nbcut()));
  histograms.push_back(&hHT_DR);

  histConfig hHT_DR_xWt;  // for weighted centers of bins
  hHT_DR_xWt.name = TString("hHT_DR_xWt_") + sample;  hHT_DR_xWt.title = "HT bin values";
  hHT_DR_xWt.NbinsX = hHT_DR.NbinsX;
  hHT_DR_xWt.binsX = hHT_DR.binsX;
  hHT_DR_xWt.axisTitles.first = "HT [GeV]";  hHT_DR_xWt.axisTitles.second = "Bin value / bin";
  hHT_DR_xWt.filler1D = &RA2bZinvAnalysis::fillHT_DR_xWt;
  hHT_DR_xWt.omitCuts.push_back(&(evSelector_->Nbcut()));
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
  hHT1_DRCC.dvalue = &(Tupl->HT);  hHT1_DRCC.omitCuts.push_back(&(evSelector_->HTcut()));
  hHT1_DRCC.omitCuts.push_back(&(evSelector_->Nbcut()));
  if (isPhoton) histograms.push_back(&hHT1_DRCC);
  histConfig hHT1_DRCC_xWt;  // for weighted centers of bins
  hHT1_DRCC_xWt.name = TString("hHT1_DRCC_xWt_") + sample;  hHT1_DRCC_xWt.title = "HT CC bin values";
  hHT1_DRCC_xWt.NbinsX = hHT1_DRCC.NbinsX;
  hHT1_DRCC_xWt.binsX = hHT1_DRCC.binsX;
  hHT1_DRCC_xWt.axisTitles.first = "HT [GeV]";  hHT1_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hHT1_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillHT_DR_xWt;
  hHT1_DRCC_xWt.omitCuts.push_back(&(evSelector_->Nbcut()));
  hHT1_DRCC_xWt.omitCuts.push_back(&(evSelector_->HTcut()));
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
  hHT2_DRCC.dvalue = &(Tupl->HT);  hHT2_DRCC.omitCuts.push_back(&(evSelector_->HTcut()));
  hHT2_DRCC.omitCuts.push_back(&(evSelector_->Nbcut()));
  if (isPhoton) histograms.push_back(&hHT2_DRCC);
  histConfig hHT2_DRCC_xWt;  // for weighted centers of bins
  hHT2_DRCC_xWt.name = TString("hHT2_DRCC_xWt_") + sample;  hHT2_DRCC_xWt.title = "HT CC bin values";
  hHT2_DRCC_xWt.NbinsX = hHT2_DRCC.NbinsX;
  hHT2_DRCC_xWt.binsX = hHT2_DRCC.binsX;
  hHT2_DRCC_xWt.axisTitles.first = "HT [GeV]";  hHT2_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hHT2_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillHT_DR_xWt;
  hHT2_DRCC_xWt.omitCuts.push_back(&(evSelector_->Nbcut()));
  hHT2_DRCC_xWt.omitCuts.push_back(&(evSelector_->HTcut()));
  if (isPhoton) histograms.push_back(&hHT2_DRCC_xWt);

  Double_t hMHT_DR_bins[7] = {300., 350., 400., 450., 600., 750, 900.};  // Check if CCbinning changes
  histConfig hMHT_DR;
  hMHT_DR.name = TString("hMHT_DR_") + sample;  hMHT_DR.title = "MHT for DR fit";
  hMHT_DR.NbinsX = 6;
  hMHT_DR.binsX = hMHT_DR_bins;
  hMHT_DR.axisTitles.first = "MHT [GeV]";  hMHT_DR.axisTitles.second = "Events / bin";
  hMHT_DR.dvalue = &(Tupl->MHT);
  hMHT_DR.omitCuts.push_back(&(evSelector_->Nbcut()));
  histograms.push_back(&hMHT_DR);

  histConfig hMHT_DR_xWt;  // for weighted centers of bins
  hMHT_DR_xWt.name = TString("hMHT_DR_xWt_") + sample;  hMHT_DR_xWt.title = "MHT bin values";
  hMHT_DR_xWt.NbinsX = hMHT_DR.NbinsX;
  hMHT_DR_xWt.binsX = hMHT_DR.binsX;
  hMHT_DR_xWt.axisTitles.first = "MHT [GeV]";  hMHT_DR_xWt.axisTitles.second = "Bin value / bin";
  hMHT_DR_xWt.filler1D = &RA2bZinvAnalysis::fillMHT_DR_xWt;
  hMHT_DR_xWt.omitCuts.push_back(&(evSelector_->Nbcut()));
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
  hMHT_DRCC.dvalue = &(Tupl->MHT);
  hMHT_DRCC.omitCuts.push_back(&(evSelector_->Nbcut()));
  hMHT_DRCC.omitCuts.push_back(&(evSelector_->MHTcut()));
  hMHT_DRCC.omitCuts.push_back(&(evSelector_->ptCut()));
  if (isPhoton) histograms.push_back(&hMHT_DRCC);
  histConfig hMHT_DRCC_xWt;  // for weighted centers of bins
  hMHT_DRCC_xWt.name = TString("hMHT_DRCC_xWt_") + sample;  hMHT_DRCC_xWt.title = "MHT CC bin values";
  hMHT_DRCC_xWt.NbinsX = hMHT_DRCC.NbinsX;
  hMHT_DRCC_xWt.binsX = hMHT_DRCC.binsX;
  hMHT_DRCC_xWt.axisTitles.first = "MHT [GeV]";  hMHT_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hMHT_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillMHT_DR_xWt;
  hMHT_DRCC_xWt.omitCuts.push_back(&(evSelector_->Nbcut()));
  hMHT_DRCC_xWt.omitCuts.push_back(&(evSelector_->MHTcut()));
  hMHT_DRCC_xWt.omitCuts.push_back(&(evSelector_->ptCut()));
  if (isPhoton) histograms.push_back(&hMHT_DRCC_xWt);

  histConfig hNJets_DR;
  hNJets_DR.name = TString("hNJets_DR_") + sample;  hNJets_DR.title = "NJets";
  hNJets_DR.NbinsX = 8;  hNJets_DR.rangeX.first = 1.5;  hNJets_DR.rangeX.second = 9.5;
  hNJets_DR.axisTitles.first = "N (jets)";  hNJets_DR.axisTitles.second = "Events / bin";
  hNJets_DR.ivalue = &(Tupl->NJets);
  hNJets_DR.omitCuts.push_back(&(evSelector_->Nbcut()));
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
  hNJets_DRCC.ivalue = &(Tupl->NJets);
  hNJets_DRCC.omitCuts.push_back(&(evSelector_->Nbcut()));
  histograms.push_back(&hNJets_DRCC);
  histConfig hNJets_DRCC_xWt;  // for weighted centers of bins
  hNJets_DRCC_xWt.name = TString("hNJets_DRCC_xWt_") + sample;  hNJets_DRCC_xWt.title = "NJets CC bin values";
  hNJets_DRCC_xWt.NbinsX = hNJets_DRCC.NbinsX;
  hNJets_DRCC_xWt.binsX = hNJets_DRCC.binsX;
  hNJets_DRCC_xWt.axisTitles.first = "N (jets)";  hNJets_DRCC_xWt.axisTitles.second = "Bin value / bin";
  hNJets_DRCC_xWt.filler1D = &RA2bZinvAnalysis::fillNJets_DR_xWt;
  hNJets_DRCC_xWt.omitCuts.push_back(&(evSelector_->Nbcut()));
  if (isPhoton) histograms.push_back(&hNJets_DRCC_xWt);

  histConfig hSFwt_DR;
  hSFwt_DR.name = TString("hSFwt_DR_") + sample;  hSFwt_DR.title = "Weights for efficiency, purity";
  hSFwt_DR.NbinsX = 80;  hSFwt_DR.rangeX.first = 0.8;  hSFwt_DR.rangeX.second = 1.2;
  hSFwt_DR.axisTitles.first = "SF weight";  hSFwt_DR.axisTitles.second = "Events / bin";
  hSFwt_DR.filler1D = &RA2bZinvAnalysis::fillSFwt_DR;
  hSFwt_DR.omitCuts.push_back(&(evSelector_->Nbcut()));
  histograms.push_back(&hSFwt_DR);

  histConfig hSFsys_DR;
  hSFsys_DR.name = TString("hSFsys_DR_") + sample;  hSFsys_DR.title = "Syst errors for efficiency, purity";
  hSFsys_DR.NbinsX = 100;  hSFsys_DR.rangeX.first = 0;  hSFsys_DR.rangeX.second = 0.2;
  hSFsys_DR.axisTitles.first = "SF error";  hSFsys_DR.axisTitles.second = "Events / bin";
  hSFsys_DR.filler1D = &RA2bZinvAnalysis::fillSFsys_DR;
  hSFsys_DR.omitCuts.push_back(&(evSelector_->Nbcut()));
  histograms.push_back(&hSFsys_DR);

  // end DR histograms

  histConfig hnZcand;
  hnZcand.name = TString("hnZcand_") + sample;  hnZcand.title = "Number of Z candidates";
  hnZcand.NbinsX = 10;  hnZcand.rangeX.first = 0;  hnZcand.rangeX.second = 10;
  hnZcand.axisTitles.first = "N(Z candidates)";  hnZcand.axisTitles.second = "Events / bin";
  hnZcand.filler1D = &RA2bZinvAnalysis::fillnZcand;  hnZcand.omitCuts.push_back(&(evSelector_->massCut()));
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
  hZpt.filler1D = &RA2bZinvAnalysis::fillZpt;
  hZpt.omitCuts.push_back(&(evSelector_->ptCut()));  hZpt.omitCuts.push_back(&(evSelector_->MHTcut()));
  if (isZll) histograms.push_back(&hZpt);

  histConfig hPhotonPt;
  hPhotonPt.name = TString("hPhotonPt_") + sample;  hPhotonPt.title = "Photon Pt";
  hPhotonPt.NbinsX = 60;  hPhotonPt.rangeX.first = 0;  hPhotonPt.rangeX.second = 3000;
  hPhotonPt.axisTitles.first = "Pt(Photon) [GeV]";  hPhotonPt.axisTitles.second = "Events / 50 GeV";
  hPhotonPt.filler1D = &RA2bZinvAnalysis::fillPhotonPt;
  hPhotonPt.omitCuts.push_back(&(evSelector_->ptCut()));
  hPhotonPt.omitCuts.push_back(&(evSelector_->MHTcut()));
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

  histConfig hJetDphi;
  hJetDphi.name = TString("hJetDphi_") + sample;  hJetDphi.title = "Delta phi (jet, MHT) HEM";
  hJetDphi.NbinsX = 200;  hJetDphi.rangeX.first = -1;  hJetDphi.rangeX.second = 1;
  hJetDphi.axisTitles.first = "#Delta#phi";  hJetDphi.axisTitles.second = "Events / 0.01";
  hJetDphi.filler1D = &RA2bZinvAnalysis::fillJetDphi;  hJetDphi.addCuts = "RA2bin>0";
  // histograms.push_back(&hJetDphi);  // This was for debugging

  histConfig hZmass;
  hZmass.name = TString("hZmass_") + sample;  hZmass.title = "Z mass";
  hZmass.NbinsX = 30;  hZmass.rangeX.first = 60;  hZmass.rangeX.second = 120;
  hZmass.axisTitles.first = "M(Z) [GeV]";  hZmass.axisTitles.second = "Events / 2 GeV";
  hZmass.filler1D = &RA2bZinvAnalysis::fillZmass;  hZmass.omitCuts.push_back(&(evSelector_->massCut()));
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

  // histConfig httZpt;
  // httZpt.name = TString("httZpt_") + sample;  httZpt.title = "Z Pt, ttZ candidates";
  // httZpt.NbinsX = 60;  httZpt.rangeX.first = 0;  httZpt.rangeX.second = 3000;
  // httZpt.axisTitles.first = "Pt(Z) [GeV]";  httZpt.axisTitles.second = "Events / 50 GeV";
  // httZpt.filler1D = &RA2bZinvAnalysis::fillttZpt;
  // httZpt.omitCuts.push_back(&(evSelector_->ptCut()));
  // httZpt.omitCuts.push_back(&(evSelector_->MHTcut()));
  // httZpt.omitCuts.push_back(&(evSelector_->HTcut()));
  // if (isZll) histograms.push_back(&httZpt);

  // histConfig hZmass_sfLepTksVeto(hZmass);
  // hZmass_sfLepTksVeto.name = TString("hZmass_sfLepTksVeto_") + sample;  hZmass_sfLepTksVeto.title = "Z mass, SF lepton vetoed";
  // hZmass_sfLepTksVeto.omitCuts.push_back(&(evSelector_->isoSFlepTksCut()));
  // hZmass_sfLepTksVeto.addCuts = evSelector_->isoSFlepTksVeto().Data();
  // histograms.push_back(&hZmass_sfLepTksVeto);

  // histConfig hZmass_photonVeto(hZmass);
  // hZmass_photonVeto.name = TString("hZmass_photonVeto_") + sample;  hZmass_photonVeto.title = "Z mass, photon vetoed";
  // hZmass_photonVeto.omitCuts.push_back(&(evSelector_->photonVeto()));  hZmass_photonVeto.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hZmass_photonVeto);

  // histConfig hGpt;
  // hGpt.name = TString("hGpt_") + sample;  hGpt.title = "Photon Pt";
  // hGpt.NbinsX = 60;  hGpt.rangeX.first = 0;  hGpt.rangeX.second = 3000;
  // hGpt.axisTitles.first = "Pt(gamma) [GeV]";  hGpt.axisTitles.second = "Events / 50 GeV";
  // hGpt.filler1D = &RA2bZinvAnalysis::fillGpt;
  // hGpt.omitCuts.push_back(&(evSelector_->photonVeto()));  hGpt.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hGpt);

  // histConfig hZGmass;
  // hZGmass.name = TString("hZGmass_") + sample;  hZGmass.title = "Z-gamma mass";
  // hZGmass.NbinsX = 100;  hZGmass.rangeX.first = 0;  hZGmass.rangeX.second = 2000;
  // hZGmass.axisTitles.first = "M(Z gamma) [GeV]";  hZGmass.axisTitles.second = "Events / 20 GeV";
  // hZGmass.filler1D = &RA2bZinvAnalysis::fillZGmass;
  // hZGmass.omitCuts.push_back(&(evSelector_->photonVeto()));  hZGmass.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hZGmass);

  // histConfig hZGdRvsM;
  // hZGdRvsM.name = TString("hZGdRvsM_")+sample;
  // hZGdRvsM.title = "Min Delta R vs M(Zgamma) for Z(ll) leptons";
  // hZGdRvsM.NbinsX = 100;  hZGdRvsM.rangeX.first = 0;  hZGdRvsM.rangeX.second = 200;
  // hZGdRvsM.NbinsY = 40;  hZGdRvsM.rangeY.first = 0;  hZGdRvsM.rangeY.second = 0.02;
  // hZGdRvsM.axisTitles.first = "M(Z gamma) [GeV]";  hZGdRvsM.axisTitles.second = "DR(l gamma)";
  // hZGdRvsM.filler2D = &RA2bZinvAnalysis::fillZGdRvsM;
  // hZGdRvsM.omitCuts.push_back(&(evSelector_->photonVeto()));  hZGdRvsM.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hZGdRvsM);

  // histConfig hGJdR;
  // hGJdR.name = TString("hGJdR_") + sample;  hGJdR.title = "min DR(photon-jet)";
  // hGJdR.NbinsX = 200;  hGJdR.rangeX.first = 0;  hGJdR.rangeX.second = 4;
  // hGJdR.axisTitles.first = "Delta R";  hGJdR.axisTitles.second = "Events / 0.02";
  // hGJdR.filler1D = &RA2bZinvAnalysis::fillGJdR;
  // hGJdR.omitCuts.push_back(&(evSelector_->photonVeto()));  hGJdR.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hGJdR);

  // histConfig hGLdRnoPixelSeed;
  // hGLdRnoPixelSeed.name = TString("hGLdRnoPixelSeed_") + sample;  hGLdRnoPixelSeed.title = "min DR(photon-lepton), noPixelSeed";
  // hGLdRnoPixelSeed.NbinsX = 200;  hGLdRnoPixelSeed.rangeX.first = 0;  hGLdRnoPixelSeed.rangeX.second = 4;
  // hGLdRnoPixelSeed.axisTitles.first = "Delta R";  hGLdRnoPixelSeed.axisTitles.second = "Events / 0.02";
  // hGLdRnoPixelSeed.filler1D = &RA2bZinvAnalysis::fillGLdRnoPixelSeed;
  // hGLdRnoPixelSeed.omitCuts.push_back(&(evSelector_->photonVeto()));  hGLdRnoPixelSeed.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hGLdRnoPixelSeed);

  // histConfig hGLdRpixelSeed;
  // hGLdRpixelSeed.name = TString("hGLdRpixelSeed_") + sample;  hGLdRpixelSeed.title = "min DR(photon-lepton), pixelSeed";
  // hGLdRpixelSeed.NbinsX = 200;  hGLdRpixelSeed.rangeX.first = 0;  hGLdRpixelSeed.rangeX.second = 4;
  // hGLdRpixelSeed.axisTitles.first = "Delta R";  hGLdRpixelSeed.axisTitles.second = "Events / 0.02";
  // hGLdRpixelSeed.filler1D = &RA2bZinvAnalysis::fillGLdRpixelSeed;
  // hGLdRpixelSeed.omitCuts.push_back(&(evSelector_->photonVeto()));  hGLdRpixelSeed.addCuts = evSelector_->photonCut().Data();
  // histograms.push_back(&hGLdRpixelSeed);

  bookAndFillHistograms(sample, histograms, cutHistograms);
  double cpu1  = get_cpu_time();
  cout << endl << "CPU Time  = " << cpu1  - cpu0  << endl << endl;

  std::vector<TH1*> theHists;
  for (auto & thisHist : cutHistograms) theHists.push_back(thisHist->hist);
  for (auto & thisHist : histograms) theHists.push_back(thisHist->hist);

  delete Tupl;
  delete evSelector_;

  return theHists;

}  // ======================================================================================

void
RA2bZinvAnalysis::bookAndFillHistograms(const char* sample,
                                        vector<histConfig*>& histograms,
                                        vector<histConfig*>& cutHistograms) {
  //
  // Define N - 1 (or N - multiple) cuts, book histograms.  Traverse the chain and fill.
  //

  CutManager::string_map sampleMap = evSelector_->sampleKeyMap();
  TString sampleKey = sampleMap.count(sample) > 0 ? sampleMap.at(sample) : "none";

  TCut baselineCuts = evSelector_->baseline();
  Tupl->setTF("baseline", baselineCuts);

  // For B-tagging corrections or modified b tagging selection
  btagsf_ = new BTagSF();

  // Z Pt weights
  Double_t ptBins[297] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.,64.,65.,66.,67.,68.,69.,70.,71.,72.,73.,74.,75.,76.,77.,78.,79.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,200.,202.,204.,206.,208.,210.,212.,214.,216.,218.,220.,222.,224.,226.,228.,230.,232.,234.,236.,238.,240.,242.,244.,246.,248.,250.,252.,254.,256.,258.,260.,262.,264.,266.,268.,270.,272.,274.,276.,278.,280.,282.,284.,286.,288.,290.,292.,294.,296.,298.,300.,304.,308.,312.,316.,320.,324.,328.,332.,336.,340.,344.,348.,352.,356.,360.,364.,368.,372.,376.,380.,384.,388.,392.,396.,400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,500.,520.,540.,560.,580.,600.,650.,700.,750.,800.,900.,1000.};
  Double_t ptWgts[297] =
    {0.615, 0.626, 0.649, 0.688, 0.739, 0.796, 0.852, 0.899, 0.936, 0.969, 0.995, 1.015, 1.034, 1.050, 1.065, 1.079, 1.092, 1.103, 1.118, 1.126, 1.139, 1.146, 1.156, 1.160, 1.164, 1.164, 1.165, 1.165, 1.164, 1.164, 1.164, 1.161, 1.158, 1.156, 1.154, 1.146, 1.147, 1.140, 1.135, 1.134, 1.129, 1.128, 1.121, 1.115, 1.110, 1.109, 1.111, 1.103, 1.094, 1.092, 1.093, 1.089, 1.085, 1.084, 1.074, 1.074, 1.067, 1.068, 1.062, 1.066, 1.063, 1.055, 1.050, 1.054, 1.048, 1.044, 1.042, 1.043, 1.034, 1.032, 1.035, 1.033, 1.028, 1.030, 1.029, 1.030, 1.012, 1.012, 1.018, 1.013, 1.011, 1.000, 1.007, 1.015, 0.989, 1.001, 0.994, 0.990, 0.990, 0.986, 0.981, 0.986, 0.972, 0.971, 0.977, 0.974, 0.978, 0.966, 0.985, 0.978, 0.966, 0.972, 0.969, 0.971, 0.964, 0.962, 0.948, 0.952, 0.943, 0.973, 0.936, 0.945, 0.935, 0.944, 0.953, 0.943, 0.935, 0.926, 0.946, 0.940, 0.943, 0.933, 0.920, 0.912, 0.923, 0.904, 0.919, 0.937, 0.941, 0.929, 0.923, 0.902, 0.925, 0.927, 0.896, 0.926, 0.905, 0.920, 0.908, 0.889, 0.910, 0.913, 0.911, 0.889, 0.881, 0.888, 0.904, 0.886, 0.891, 0.915, 0.879, 0.884, 0.900, 0.874, 0.850, 0.874, 0.873, 0.872, 0.894, 0.880, 0.870, 0.871, 0.863, 0.883, 0.863, 0.892, 0.845, 0.863, 0.892, 0.851, 0.878, 0.847, 0.881, 0.809, 0.860, 0.829, 0.851, 0.871, 0.846, 0.825, 0.860, 0.853, 0.894, 0.807, 0.793, 0.863, 0.832, 0.829, 0.870, 0.850, 0.831, 0.820, 0.829, 0.839, 0.880, 0.831, 0.804, 0.836, 0.822, 0.796, 0.812, 0.831, 0.830, 0.827, 0.802, 0.781, 0.855, 0.798, 0.774, 0.790, 0.810, 0.825, 0.799, 0.790, 0.811, 0.778, 0.803, 0.779, 0.795, 0.768, 0.809, 0.762, 0.767, 0.752, 0.822, 0.777, 0.791, 0.799, 0.721, 0.776, 0.718, 0.798, 0.754, 0.749, 0.758, 0.847, 0.766, 0.775, 0.718, 0.756, 0.826, 0.792, 0.792, 0.767, 0.679, 0.694, 0.750, 0.738, 0.707, 0.709, 0.711, 0.754, 0.762, 0.717, 0.722, 0.692, 0.714, 0.726, 0.703, 0.693, 0.704, 0.704, 0.653, 0.718, 0.748, 0.713, 0.733, 0.741, 0.742, 0.652, 0.586, 0.644, 0.656, 0.702, 0.721, 0.682, 0.758, 0.662, 0.578, 0.654, 0.723, 0.634, 0.680, 0.656, 0.602, 0.590, 0.566, 0.623, 0.562, 0.589, 0.604, 0.554, 0.525, 0.552, 0.503, 0.563, 0.476};

  cutHistos cutHistFiller(Tupl, evSelector_);  // for cutFlow histograms

  if (verbosity_ >= 1) cout << "\nFor sample " << sample << ", baseline = " << endl << baselineCuts << endl << endl;
  if (verbosity_ >= 1) cout << "commonCuts = " << evSelector_->commonCuts() << endl << endl;
  bookHistograms(cutHistograms, baselineCuts, cutHistFiller);
  bookHistograms(histograms, baselineCuts, cutHistFiller);
  Tupl->setNotify();

  // Traverse the tree and fill histograms
  vector<unsigned> triggerIndexList;
  int currentYear = -1;
  double MCwtCorr = 1.;
  int count = 0, countInFile = 0, countInSel = 0, countNegWt = 0;
  Long64_t Nentries = Tupl->fChain->GetEntries();
  if (verbosity_ >= 1) cout << endl << "Nentries in tree = " << Nentries << endl << endl;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (verbosity_ >= 1 && count % 100000 == 0) cout << "Entry number " << count << endl;

    Long64_t centry = Tupl->LoadTree(entry);
    if (centry < 0) {
      cout << "LoadTree returned " << centry;
      break;
    }
    Tupl->GetEntry(entry);
    // cout << endl << "First entry" << endl;  Show();  break;
    if (Tupl->newFileInChain()) {
      // New input root file encountered
      Tupl->setNewFileInChain(false);
      countInFile = 0;
      TFile* thisFile = Tupl->fChain->GetCurrentFile();
      TString path = thisFile->GetName();
      if (verbosity_ >= 1) cout << "Current file in chain: " << path << endl;
      // cout << "First event, new file setup " << endl;  Show();
      int theYear = -1;
      if (isMC_) {
        // Look for strings in root filenames that indicate the running year
        if      (path.Contains("MC2016") || path.Contains("Summer16v3")) theYear = EfficWt::Year2016;
        else if (path.Contains("MC2017") || path.Contains("Fall17")) theYear = EfficWt::Year2017;
        else if (path.Contains("MC2018") || path.Contains("Autumn18")) theYear = EfficWt::Year2018;
      } else {
        if      (Tupl->RunNum < CutManager::Start2017) theYear = EfficWt::Year2016;
        else if (Tupl->RunNum < CutManager::Start2018) theYear = EfficWt::Year2017;
        else                                           theYear = EfficWt::Year2018;
      }
      if (theYear != currentYear) {
        currentYear = theYear;
        cout << "currentYear = " << currentYear << endl;
        // Load the year-dependent correction factors.
        effPurCorr_->getHistos(sample, currentYear);  // For purity, Fdir, trigger eff, reco eff
        if (ntupleVersion_ == "V15" && currentYear != EfficWt::Year2016) csvMthreshold_ = 0.8838;
        if (isMC_) {
          if (currentYear == EfficWt::Year2016) {
            // For now we have only 2016 pileup correction files
            if (applyPuWeight_ && customPuWeight_) {
              TFile* pufile = TFile::Open("../../Analysis/corrections/PileupHistograms_0121_69p2mb_pm4p6.root", "READ");
              puHist_ = (TH1*) pufile->Get("pu_weights_down");
            }
	  }
	  if (btagsf_) {
	    if (minNbCut_ == 0) btagsf_->SetCalib((unsigned) currentYear);
	    else btagsf_->SetCalib((unsigned) currentYear, BTagEntryS::OP_TIGHT);
	  }
	}
      }
      // Set MCwtCorr for this file
      if (isMC_ && isSkim_ &&
          (path.Contains("V16") || path.Contains("V17")) &&
          (currentYear == EfficWt::Year2017 || currentYear == EfficWt::Year2018) &&
          // (path.Contains("MC2017") || path.Contains("MC2018")) &&
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
      if (isMC_ && verbosity_ >= 1) cout << "MC weight for this file is " << Tupl->Weight
                                         << " times correction " << MCwtCorr << endl;
      if (isMC_ && applyBTagSF_) btagsf_->SetEffs(thisFile);
      evSelector_->setTriggerIndexList(sample, &triggerIndexList, Tupl);
    }  // newFileInChain

    countInFile++;
    // if (countInFile == 1) cout << "After get first entry in file, status of HT = " << Tupl->fChain->GetBranchStatus("HT")
    //                         << ", JetIDAK8 = " << Tupl->fChain->GetBranchStatus("JetIDAK8") << endl;

    if (!isSkim_ && !restrictClean_) {
      if (TString(sample).Contains("mm") && Tupl->NMuons != 2) continue;
      if (TString(sample).Contains("ee") && Tupl->NElectrons != 2) continue;
      if ((TString(sample).Contains("photon") || TString(sample).Contains("gjets")) && Tupl->Photons->size() == 0) continue;
    }
    cleanVars();  // If unskimmed input, copy <var>clean to <var>
    Int_t BTagsOrig = setBTags(currentYear);
    // if (countInFile <= 100) cout << "BTagsOrig, BTags, BTagsDeepCSV = " << BTagsOrig
    //                           << ", " << Tupl->BTags << ", " << Tupl->BTagsDeepCSV << endl;

    if (Tupl->ZCandidates->size() > 1 && verbosity_ >= 2) cout << Tupl->ZCandidates->size() << " Z candidates found" << endl;
    // double baselineWt = Tupl->TFvalue("baseline");

    // Compute event weight factors
    Double_t eventWt = 1, MCwt = 1, PUweight = 1, NoPrefireWt = 1, ZPtWt = 1;
    if (isMC_) {
      MCwt = 1000*intLumi_*Tupl->Weight*MCwtCorr;  // MCwtCorr for 2017 MC
      if (applyPuWeight_) {
        // Pileup weight for 2016
        if (customPuWeight_ && puHist_ != nullptr) {
          // This PU weight recipe from Kevin Pedro, https://twiki.cern.ch/twiki/bin/viewauth/CMS/RA2b13TeVProduction
          PUweight = puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(min(Tupl->TrueNumInteractions,
                                                                             puHist_->GetBinLowEdge(puHist_->GetNbinsX()+1))));
        } else
          PUweight = Tupl->puWeight;  // Take puWeight directly from the tree
        MCwt *= PUweight;
      }

      if (currentYear == EfficWt::Year2016 || currentYear == EfficWt::Year2017) {
        // Apply L1 prefire weight for 2016 and 2017
        if (sampleKey.Contains("ee")) {
          for (unsigned j = 0; j < Tupl->Jets->size(); ++j) {
            NoPrefireWt *= effPurCorr_->prefiring_weight_jet(Tupl->Jets, j);
          }
          double eeNoPFwt = 1;
          for (unsigned e = 0; e < Tupl->Electrons->size(); ++e) {
            double w = effPurCorr_->prefiring_weight_electron(Tupl->Electrons, e);
            if (w < eeNoPFwt) eeNoPFwt = w;
          }
          NoPrefireWt *= eeNoPFwt;
        } else {
          NoPrefireWt = Tupl->NonPrefiringProb;
        }
        MCwt *= NoPrefireWt;
      }  // 2016 or 2017

      if (applyZptWt_ && (ntupleVersion_ == "V16" || ntupleVersion_ == "V17")
	  && (TString(sample).Contains("dy") || TString(sample).Contains("zinv"))
          && (currentYear == EfficWt::Year2017 || currentYear == EfficWt::Year2018)) {
        // Apply Z Pt weight for 2017 MC
        double ptZ = getGenPtZ();
        ZPtWt = 1.0;
        if (ptZ > 0.0) {
          int iptbin;
          for (iptbin = 1; iptbin<297; iptbin++) {
            if (ptZ < ptBins[iptbin]) break;
          }
          ZPtWt *= ptWgts[iptbin-1];
        }
        MCwt *= ZPtWt;
      }  // DY or Zinv in 2017 or 2018, for V16 and V17

      if (minNbCut_ == 1 && applyBTagSF_) {
	// Compute B tag SF via "simple correction (1a)"
	double BSFwt = Tupl->BTags > 0 ?
	  btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
			  Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll)
	  : 1;
	MCwt *= BSFwt;
      }

      eventWt *= MCwt;
    }  // isMC_

    effWt_ = 1;  effSys_ = 0;
    bool effValid = false;

    // Skim cuts with no branches in unskimmed ntuples
    bool passEcalNoiseJetFilter = true, passPhotonVeto = true, DiLepton = true;

    if (!isSkim_ && ntupleVersion_ != "V12" && ntupleVersion_ != "V15" ) {
      Tupl->EcalNoiseJetFilter = !isMC_ && currentYear == EfficWt::Year2017 ?
        ecalNoiseJetFilter() : true;
      passEcalNoiseJetFilter = Tupl->EcalNoiseJetFilter;

      // Photon veto
      int NumPhotons = 0;
      for (unsigned p = 0; p < Tupl->Photons->size(); ++p) {
      	if(!(Tupl->Photons_hasPixelSeed->at(p)) && Tupl->Photons->at(p).Pt() > 100) ++NumPhotons;
      	// if(!(Tupl->Photons_hasPixelSeed->at(p)) && Tupl->Photons_fullID->at(p)
      	//    && Tupl->Photons->at(p).Pt() > 100) ++NumPhotons;
      }
      if (NumPhotons > 0) passPhotonVeto = false;

      // Dilepton requirement:  2 iso, medium ID, opposite-charge leptons (skim cut)
      if (TString(sample).Contains("zmm"))
        DiLepton = diLepton(Tupl->NMuons, Tupl->Muons, Tupl->Muons_passIso, Tupl->Muons_charge, Tupl->Muons_mediumID);
      if (TString(sample).Contains("zee"))
        DiLepton = diLepton(Tupl->NElectrons, Tupl->Electrons, Tupl->Electrons_passIso, Tupl->Electrons_charge);
    }

    // Cuts not amenable to implementation via TTreeFormulas
    // Trigger requirements
    bool passTrg = true;
    if (!isMC_) {
      passTrg = false;
      for (auto trgIndex : triggerIndexList)
        if (Tupl->TriggerPass->at(trgIndex) == 1) passTrg = true;
	// 1 = Trigger is found in the list/menu, and it fired
	// 0 = Trigger is found in the list/menu, but it did not fire
	// -1 = Trigger is not not found in the list/menu
    }

    // HEM veto for 2018HEM
    bool passHEM = true;
    if (applyHEMjetVeto_ && !passHEMjetVeto()) passHEM = false;

    // Fill cut flow histograms before imposing inline requirements
    for (auto & hg : cutHistograms) {
      double cutHistWt = 1;
      if (hg->name.Contains("Wt")) {
        cutHistWt = eventWt;
        if (applySFwtToMC_ && isMC_) {
	  pair<double, double> efficiency = effPurCorr_->weight(Tupl, CCbins_, applyDRfitWt_, currentYear);
	  effWt_ = efficiency.first;  effSys_ = efficiency.second;  effValid = true;
	  cutHistWt *= effWt_;
	}
      }
      bool zeroChg = true;
      if ((TString(sample).Contains("zmm") || TString(sample).Contains("zee")) && !(DiLepton && passPhotonVeto)) zeroChg = false;
      if (hg->name.Contains("Skim")) {
	// cout << "NJets = " << Tupl->NJets << ", cut = " << evSelector_->NJetSkimCut() << ", pass = " << Tupl->TFvalue("NJetSkimCut");
	cutHistFiller.fill((TH1D*) hg->hist, cutHistWt, DiLepton, passPhotonVeto);
      } else {
	cutHistFiller.fill((TH1D*) hg->hist, cutHistWt, zeroChg, passTrg, passHEM, passEcalNoiseJetFilter);
      }
    }

    if ((TString(sample).Contains("zmm") || TString(sample).Contains("zee")) && !(DiLepton && passPhotonVeto)) continue;
    if (!passTrg) continue;
    if (!passHEM) continue;
    if (!passEcalNoiseJetFilter) continue;
    // (For a test) select events with electron (photon) in HEM region
    // bool keep = false; for (auto & theE : *(Tupl->Electrons)) {if (!passHEMobjVeto(theE, 30, false)) keep = true;}  if (!keep) continue;
    // bool keep = false; for (auto & theG : *(Tupl->Photons)) {if (!passHEMobjVeto(theG, 30, false)) keep = true;}  if (!keep) continue;

    int CCbin = CCbins_->jbk(CCbins_->jbin(Tupl->NJets), CCbins_->bbin(Tupl->NJets, Tupl->BTags),
                             CCbins_->kinBin(Tupl->HT, Tupl->MHT));
    if (isSkim_ && minNbCut_ == 0 && (UInt_t) CCbin != Tupl->RA2bin && !(CCbin == -1 && Tupl->RA2bin == 0)) {
      cout << "CCbin = " << CCbin << ", != RA2bin = " << Tupl->RA2bin
           << ", NJets = " << Tupl->NJets << ", j = " << CCbins_->jbin(Tupl->NJets)
           << ", Nb = " << Tupl->BTags << ", b = " << CCbins_->bbin(Tupl->NJets, Tupl->BTags)
           << ", HT = " << Tupl->HT << ", MHT = " << Tupl->MHT << ", k = " << CCbins_->kinBin(Tupl->HT, Tupl->MHT)
           << endl;
    }

    if (!effValid) {
      pair<double, double> efficiency = effPurCorr_->weight(Tupl, CCbins_, applyDRfitWt_, currentYear);
      effWt_ = efficiency.first;  effSys_ = efficiency.second;  effValid = true;
    }

    for (auto & hg : histograms) {
      double selWt = Tupl->TFvalue(hg->name);
      if (selWt == 0) continue;

      if (hg->name.Contains("hCC_")) {
	countInSel++;
	if (eventWt < 0) countNegWt++;
      }

      double eventWt0 = eventWt;
      if ((applySFwtToMC_ && isMC_) || hg->name.Contains(TString("_DR"))) {
	// For MC, or double ratio, apply weights for trigger eff, reco eff.
	if (hg->name.Contains(TString("_DR"))) {
	  if (Tupl->BTags > 0) continue;
	  if (CCbin <= 0) continue;
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

  if (btagsf_) delete btagsf_;

}  // ======================================================================================

void
RA2bZinvAnalysis::bookHistograms(vector<histConfig*>& histoList, TCut baselineCuts, cutHistos cutHistFiller) {
  // Book histograms from configuration histoList
  for (auto & hg : histoList) {
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
    if (hg->name.Contains("Cut")) cutHistFiller.setAxisLabels((TH1D*) hg->hist);
    if (hg->name.Contains("hCut") || hg->name.Contains(TString("hgen"))) {
      hg->NminusOneCuts = "1";
    } else {
      hg->NminusOneCuts = baselineCuts;
      for (auto cutToOmit : hg->omitCuts) hg->NminusOneCuts(*cutToOmit) = "1";
      if (strlen(hg->addCuts) != 0) hg->NminusOneCuts += TString(" && ") + hg->addCuts;
    }
    if (verbosity_ >= 1) {
      cout << hg->name << ", omitCuts = ";
      for (auto cutToOmit : hg->omitCuts) cout << *cutToOmit << " ";
      cout << "; addCuts = " << hg->addCuts;
      cout << ", cuts = " << endl << hg->NminusOneCuts << endl << endl;
    }
    Tupl->setTF(hg->name, hg->NminusOneCuts);
  }
}  // ======================================================================================

Int_t
RA2bZinvAnalysis::setBTags(int runYear) {
  Int_t BTagsOrig = Tupl->BTags;
  if (useDeepCSV_) {
    if (minNbCut_ > 0) {  // Boosted H uses cut on Nb(tight) > 0
      Tupl->BTags = 0;
      for (size_t j = 0; j < (Tupl->Jets)->size(); ++j) {
	if(Tupl->Jets_HTMask->at(j) &&
	   Tupl->Jets_bJetTagDeepCSVBvsAll->at(j) > btagsf_->WPvalue()) Tupl->BTags++;
      }
    } else {
      Tupl->BTags = Tupl->BTagsDeepCSV;
    }
  } else if (ntupleVersion_ == "V15" && runYear != EfficWt::Year2016) {
    // Recompute BTags with a different discriminator threshold
    Tupl->BTags = 0;
    for (size_t j = 0; j < (Tupl->Jets)->size(); ++j) {
      if (Tupl->Jets_HTMask->at(j) &&
	  Tupl->Jets_bDiscriminatorCSV->at(j) > csvMthreshold_) Tupl->BTags++;
    }
  }
  // From Rishi email of 20 Feb 2019, DeepCSV 2018 WP: 0.4184 and the 2017 WP: 0.4941
  return BTagsOrig;
}  // ======================================================================================

bool
RA2bZinvAnalysis::diLepton(Int_t NLeptons, vector<TLorentzVector>* Leptons,
			   vector<bool>* Leptons_passIso, vector<int>* Leptons_charge,
			   vector<bool>* Leptons_mediumID) {
  if (NLeptons != 2) return false;
  bool DiLepton = false;
  for (size_t m1 = 0; m1 < Leptons->size(); ++m1) {
    if (m1+1 < Leptons->size() && (Leptons_mediumID == nullptr || Leptons_mediumID->at(m1)) && Leptons_passIso->at(m1)) {
      for (size_t m2 = m1+1; m2 < Leptons->size(); ++m2) {
	if ((Leptons_mediumID == nullptr || Leptons_mediumID->at(m2)) && Leptons_passIso->at(m2)
	    && Leptons_charge->at(m1) != Leptons_charge->at(m2)) {
	  DiLepton = true;
	  break;
	}
      }
    }
    if (DiLepton) break;
  }
  return DiLepton;
}  // ======================================================================================

bool
RA2bZinvAnalysis::ecalNoiseJetFilter() {
  int counter = 0;
  bool goodJet[2] = {true, true};
  for (unsigned j = 0; j < Tupl->Jets->size(); ++j){
    if (counter >= 2) break;
    const auto& Jet = Tupl->Jets->at(j);
    if (Jet.Pt() > 30 && abs(Jet.Eta()) < 5.0 && abs(Jet.Eta()) > 2.4) {
      double dphi = abs(KMath::DeltaPhi(Jet.Phi(), Tupl->MHTPhi));
      if (Jet.Pt() > 250 && (dphi > 2.6 or dphi < 0.1)) goodJet[counter] = false;
      ++counter;
    }
  }
  return goodJet[0] && goodJet[1];
}  // ======================================================================================

void
RA2bZinvAnalysis::fillCC(TH1D* h, double wt) {

  // Filler for the analysis-bin histograms and variants
  TString hName(h->GetName());
  Int_t binCC = 0, binCCjb = 0;

  int binKin = CCbins_->kinBin(Tupl->HT, Tupl->MHT);
  if (binKin < 0) return;
  int binNjets = (hName.Contains("spl") || hName.Contains("Jb")) ? CCbins_->Jbin(Tupl->NJets)
                                                                 : CCbins_->jbin(Tupl->NJets);
  if (binNjets < 0) return;
  if (isMC_ && applyBTagSF_ && minNbCut_ == 0) {
    // apply BTagSF to all Nb bins
    if (verbosity_ >= 4) cout << "Size of input Jets = " << Tupl->Jets->size()
			      << ", Jets_hadronFlavor = " << Tupl->Jets_hadronFlavor->size()
			      << " Jets_HTMask = " << Tupl->Jets_HTMask->size() << endl;
    vector<double> probNb = btagsf_->GetCorrections(Tupl->Jets, Tupl->Jets_hadronFlavor,
						      Tupl->Jets_HTMask);
    for (int b = 0; b < (int) probNb.size(); ++b) {
      int NbinsB = hName.Contains("spl") || hName.Contains("Jb") ? CCbins_-> binsB(binNjets)
	                                                         : CCbins_-> binsb(binNjets);
      int binNb = min(b, NbinsB-1);
      binCC = hName.Contains("spl") || hName.Contains("Jb") ? 
	CCbins_->Jbk(binNjets, binNb, binKin) : CCbins_->jbk(binNjets, binNb, binKin);
      if (binCC <= 0) break;
      if (verbosity_ >= 4) cout << "j = " << binNjets << ", NbTags = " << Tupl->BTags
				<< ", b = " << b << ", b wt = " << probNb[b]
				<< ", k = " << binKin << ", binCC = " << binCC << endl;
      if (hName.Contains("jb") || hName.Contains("Jb")) {
	binCCjb = hName.Contains("jb") ? CCbins_->jb(binNjets, binNb) : CCbins_->Jb(binNjets, binNb);
	if (binCCjb <= 0) break;
	h->Fill(Double_t(binCCjb), wt*probNb[b]);
      } else if (hName.Contains("jk")) {
	if (b > 0) break;
	int binCCjk = CCbins_->jk(binNjets, binKin);
	if (binCCjk > 0) h->Fill(Double_t(binCCjk), wt*probNb[b]);
      } else {
	h->Fill(Double_t(binCC), wt*probNb[b]);
      }
    }
  } else {
    int binNb = CCbins_->bbin(Tupl->NJets, Tupl->BTags);
    binCC = (hName.Contains("spl") || hName.Contains("Jb")) ?
      CCbins_->Jbk(binNjets, binNb, binKin) :
      CCbins_->jbk(binNjets, binNb, binKin);
    if (binCC <= 0) return;
    if (verbosity_ >= 4) cout << "j = " << binNjets << ", b = " << binNb << ", k = " << binKin
			      << ", binCC = " << binCC << ", RA2bin = " << Tupl->RA2bin << endl;
    if (hName.Contains("jb") || hName.Contains("Jb")) {
      // Above test on j, b, k needed even here, to exclude j = 3,4, k = 0,3
      binCCjb = hName.Contains("jb") ? CCbins_->jb(binNjets, binNb) : CCbins_->Jb(binNjets, binNb);
      if (binCCjb > 0) h->Fill(Double_t(binCCjb), wt);
    } else if (hName.Contains("jk")) {
      int binCCjk = CCbins_->jk(binNjets, binKin);
      if (binCCjk <= 0) return;
      if (hName.Contains("jk_")) {
	if (Tupl->BTags == 0) h->Fill(Double_t(binCCjk), wt);
      } else if (isMC_ && applyBTagSF_ && minNbCut_ > 0 && Tupl->BTags > 0) {
	// Histograms for minNbCut = 1, with b tag systematics
	double BSFwt = 1;
	if (hName.Contains("uc")) {
	  BSFwt = btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				  Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll, 1, 0);
	  BSFwt /= btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				   Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll);
	} else if (hName.Contains("cu")) {
	  BSFwt = btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				  Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll, 0, 1);
	  BSFwt /= btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				   Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll);
	} else if (hName.Contains("dc")) {
	  BSFwt = btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				  Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll, -1, 0);
	  BSFwt /= btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				   Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll);
	} else if (hName.Contains("cd")) {
	  BSFwt = btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				  Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll, 0, -1);
	  BSFwt /= btagsf_->weight(Tupl->Jets, Tupl->Jets_hadronFlavor,
				   Tupl->Jets_HTMask, Tupl->Jets_bJetTagDeepCSVBvsAll);
	}
	h->Fill(Double_t(binCCjk), BSFwt*wt);
      }
    } else {
      h->Fill(Double_t(binCC), wt);
    }
  }  // if apply BTagSF

}  // ======================================================================================

void
RA2bZinvAnalysis::fillZmassjb(TH1D* h, double wt) {
  if (Tupl->ZCandidates->size() == 0) return;
  Int_t j, b;
  TString hName(h->GetName());
  j = hName(hName.First('j')-1) - '0';
  b = hName(hName.First('b')-1) - '0';
  // j = 2, 3, 5 in the histogram name really means
  // NJets in the first, second, or third and higher NJets bin, respectively.
  if ((j == 2 && Tupl->NJets >= CCbins_->nJetThreshold(0) && Tupl->NJets < CCbins_->nJetThreshold(1)) ||
      (j == 3 && (Tupl->NJets >= CCbins_->nJetThreshold(1) && Tupl->NJets < CCbins_->nJetThreshold(2))) ||
      (j == 5 && Tupl->NJets >= CCbins_->nJetThreshold(2))) {
    if ((b < 2 && Tupl->BTags == b) || (b == 2 && Tupl->BTags >= 2)) {
      h->Fill(Tupl->ZCandidates->at(0).M(), wt);
    }
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::fillFilterCuts(TH1D* h, double wgt) {
  double wt = TString(h->GetName()).Contains("NoWt") ? 1 : wgt;
  h->Fill(0.5, wt);
  if (!(Tupl->globalTightHalo2016Filter==1)) h->Fill(1.5, wt);
  if (ntupleVersion_ != "V12" && ntupleVersion_ != "V15" &&
      !(Tupl->globalSuperTightHalo2016Filter==1)) h->Fill(2.5, wt);
  if (!(Tupl->HBHENoiseFilter==1)) h->Fill(3.5, wt);
  if (!(Tupl->HBHEIsoNoiseFilter==1)) h->Fill(4.5, wt);
  if (!(Tupl->EcalDeadCellTriggerPrimitiveFilter==1)) h->Fill(5.5, wt);
  if (!Tupl->BadChargedCandidateFilter) h->Fill(6.5, wt);
  if (!Tupl->BadPFMuonFilter) h->Fill(7.5, wt);
  if (!(Tupl->NVtx > 0)) h->Fill(8.5, wt);
  if (!(Tupl->eeBadScFilter==1)) h->Fill(9.5, wt);
  if (!isMC_ && !(Tupl->ecalBadCalibFilter==1)) h->Fill(10.5, wt);
  if (!(Tupl->JetID)) h->Fill(12.5, wt);
  if (!(Tupl->PFCaloMETRatio < 5)) h->Fill(13.5, wt);

  bool fillHTRatioDPhi;
  if (era_ == "2016") {
    fillHTRatioDPhi = Tupl->HT5 / Tupl->HT <= 2;
  } else {
    fillHTRatioDPhi = isSkim_ ? Tupl->HTRatioDPhiFilter : 
      Tupl->HT5 / Tupl->HT < 1.2 ? true : (Tupl->DeltaPhi1 >= 1.025 * Tupl->HT5 / Tupl->HT - 0.5875);
  }
  if (!fillHTRatioDPhi) h->Fill(14.5, wt);

  bool fillEcalNoiseJetFilter = true;
  if (!isMC_ && Tupl->RunNum >= CutManager::Start2017 && Tupl->RunNum < CutManager::Start2018) {
    fillEcalNoiseJetFilter = Tupl->EcalNoiseJetFilter;
    if (!fillEcalNoiseJetFilter) h->Fill(15.5, wt);
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::fillZGmass(TH1D* h, double wt) {
  for (auto & theZ : *(Tupl->ZCandidates)) {
    for (auto & aPhoton : *(Tupl->Photons)) {
      TLorentzVector Zg(theZ);  Zg += aPhoton;
      h->Fill(Zg.M(), wt);
    }
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGJdR(TH1D* h, double wt) {
  TLorentzVector thePhoton = Tupl->Photons->at(0);
  if (Tupl->Jets->size() > 0) {
    Double_t dR = 999.;
    for (auto & thisJet : *(Tupl->Jets))
	dR = thePhoton.DeltaR(thisJet) < dR ? thePhoton.DeltaR(thisJet) : dR;
    h->Fill(dR, wt);
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillZGdRvsM(TH2F* h, double wt) {
  if (Tupl->ZCandidates->size() == 0) return;
  TLorentzVector theZ = Tupl->ZCandidates->at(0);
  TLorentzVector thePhoton = Tupl->Photons->at(0);
  TLorentzVector Zg(theZ);  Zg += thePhoton;
  if (Tupl->Muons->size() > 0) {
    Double_t dR = 999.;
    for (auto & thisMuon : *(Tupl->Muons))
	dR = thePhoton.DeltaR(thisMuon) < dR ? thePhoton.DeltaR(thisMuon) : dR;
    h->Fill(Zg.M(), dR, wt);
  }
  if (Tupl->Electrons->size() > 0) {
    Double_t dR = 999.;
    for (auto & thisElectron : *(Tupl->Electrons))
	dR = thePhoton.DeltaR(thisElectron) < dR ? thePhoton.DeltaR(thisElectron) : dR;
    h->Fill(Zg.M(), dR, wt);
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGLdRnoPixelSeed(TH1D* h, double wt) {
  for (size_t i = 0; i < Tupl->Photons->size(); ++i) {
    // if Tupl->(Photons_fullID->size()) cout << "Photons_fullID = " << Tupl->Photons_fullID->at(i) << endl;
    // if (Tupl->Photons_passElectronVeto->size()) cout << "Photons_passElectronVeto = " << Tupl->Photons_passElectronVeto->at(i) << endl;
    if (!Tupl->Photons_fullID->at(i)) continue;
    if (Tupl->Photons_hasPixelSeed->at(i)) continue;
    if (Tupl->Muons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisMuon : *(Tupl->Muons))
	dR = Tupl->Photons->at(i).DeltaR(thisMuon) < dR ? Tupl->Photons->at(i).DeltaR(thisMuon) : dR;
      h->Fill(dR, wt);
      // if (Tupl->Photons_isEB->size() && dR < 0.02) cout << "isEB photon = " << Tupl->Photons_isEB->at(i) << endl;
      // if (Tupl->Photons_hasPixelSeed->size() && dR < 0.02) cout << "Photons_hasPixelSeed, fake = " << Tupl->Photons_hasPixelSeed->at(i) << endl;
      // if (Tupl->Photons_hasPixelSeed->size() && dR >= 0.02) cout << "Photons_hasPixelSeed, real = " << Tupl->Photons_hasPixelSeed->at(i) << endl;
    }
    if (Tupl->Electrons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisElectron : *(Tupl->Electrons))
	dR = Tupl->Photons->at(i).DeltaR(thisElectron) < dR ?Tupl-> Photons->at(i).DeltaR(thisElectron) : dR;
      h->Fill(dR, wt);
      // if (Tupl->Photons_isEB->size() && dR < 0.02) cout << "isEB photon = " << Tupl->Photons_isEB->at(i) << endl;
      // if (Tupl->Photons_hasPixelSeed->size() && dR < 0.02) cout << "Photons_hasPixelSeed, fake = " << Tupl->Photons_hasPixelSeed->at(i) << endl;
      // if (Tupl->Photons_hasPixelSeed->size() && dR >= 0.02) cout << "Photons_hasPixelSeed, real = " << Tupl->Photons_hasPixelSeed->at(i) << endl;
    }
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::fillGLdRpixelSeed(TH1D* h, double wt) {
  for (size_t i = 0; i < Tupl->Photons->size(); ++i) {
    if (!Tupl->Photons_fullID->at(i)) continue;
    if (!Tupl->Photons_hasPixelSeed->at(i)) continue;
    if (Tupl->Muons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisMuon : *(Tupl->Muons))
	dR = Tupl->Photons->at(i).DeltaR(thisMuon) < dR ? Tupl->Photons->at(i).DeltaR(thisMuon) : dR;
      h->Fill(dR, wt);
    }
    if (Tupl->Electrons->size() > 0) {
      Double_t dR = 999.;
      for (auto & thisElectron : *(Tupl->Electrons))
	dR = Tupl->Photons->at(i).DeltaR(thisElectron) < dR ? Tupl->Photons->at(i).DeltaR(thisElectron) : dR;
      h->Fill(dR, wt);
    }
  }
}  // ======================================================================================

RA2bZinvAnalysis::cutHistos::cutHistos(Ntuple* tuple, CutManager* selector)
  : tuple_(tuple), selector_(selector) {
  tuple_->setTF("HTcut", selector_->HTcut());  cutNames_.push_back("HTcut");
  tuple_->setTF("MHTcut", selector_->MHTcut());  cutNames_.push_back("MHTcut");
  tuple_->setTF("NJetscut", selector_->NJetscut());  cutNames_.push_back("NJetscut");
  tuple_->setTF("minDphicut", selector_->minDphicut());  cutNames_.push_back("minDphicut");
  tuple_->setTF("objcut", selector_->objcut());  cutNames_.push_back("objcut");
  tuple_->setTF("ptCut", selector_->ptCut());  cutNames_.push_back("ptCut");
  tuple_->setTF("massCut", selector_->massCut());  cutNames_.push_back("massCut");
  tuple_->setTF("photonDeltaRcut", selector_->photonDeltaRcut());  cutNames_.push_back("photonDeltaRcut");
  tuple_->setTF("commonCuts", selector_->commonCuts());  cutNames_.push_back("commonCuts");
  // Skim cuts
  tuple_->setTF("NJetSkimCut", selector_->NJetSkimCut());  skimCutNames_.push_back("NJetSkimCut");
  tuple_->setTF("HTSkimCut", selector_->HTSkimCut());  skimCutNames_.push_back("HTSkimCut");
  tuple_->setTF("MHTSkimCut", selector_->MHTSkimCut());  skimCutNames_.push_back("MHTSkimCut");
  tuple_->setTF("MHTHTRatioSkimCut", selector_->MHTHTRatioSkimCut());  skimCutNames_.push_back("MHTHTRatioSkimCut");
  tuple_->setTF("ElectronVetoSkimCut", selector_->ElectronVetoSkimCut());  skimCutNames_.push_back("ElectronVetoSkimCut");
  tuple_->setTF("IsoElectronTrackVetoSkimCut", selector_->IsoElectronTrackVetoSkimCut());  skimCutNames_.push_back("IsoElectronTrackVetoSkimCut");
  tuple_->setTF("IsoPionTrackVetoSkimCut", selector_->IsoPionTrackVetoSkimCut());  skimCutNames_.push_back("IsoPionTrackVetoSkimCut");
  tuple_->setTF("DeltaPhiSkimCut", selector_->DeltaPhiSkimCut());  skimCutNames_.push_back("DeltaPhiSkimCut");
  tuple_->setTF("JetIDSkimCut", selector_->JetIDSkimCut());  skimCutNames_.push_back("JetIDSkimCut");

}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::setAxisLabels(TH1D* hcf) {
  if (TString(hcf->GetName()).Contains("hFilterCuts")) {
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
  if (TString(hcf->GetName()).Contains("hCut")) {
    if (TString(hcf->GetName()).Contains("Skim")) {
      hcf->GetXaxis()->SetBinLabel(1, "1 None");
      hcf->GetXaxis()->SetBinLabel(2, "2 NJets");
      hcf->GetXaxis()->SetBinLabel(3, "3 HT");
      hcf->GetXaxis()->SetBinLabel(4, "4 MHT");
      hcf->GetXaxis()->SetBinLabel(5, "5 MHT<HT");
      hcf->GetXaxis()->SetBinLabel(6, "6 DiLepton");
      hcf->GetXaxis()->SetBinLabel(7, "7 ele veto");
      hcf->GetXaxis()->SetBinLabel(8, "8 isoEtrk");
      hcf->GetXaxis()->SetBinLabel(9, "9 isoPitrk");
      hcf->GetXaxis()->SetBinLabel(10, "10 pho veto");
      hcf->GetXaxis()->SetBinLabel(11, "11 mnDphi");
      hcf->GetXaxis()->SetBinLabel(12, "12 JetID");
    } else {
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
    }
    hcf->GetXaxis()->LabelsOption("vu");
  }
}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::fill(TH1D* hcf, Double_t wt, bool zeroChg, bool passTrg, bool passHEM, bool passEcalNoiseJetFilter) {
  // Fill cut and cut flow histograms.  If null histogram pointer passed, just print cut value.

  bool pass = false;
  bool isCutFlow = TString(hcf->GetName()).Contains("Flow");
  TString cutName;
  Size_t j = 0;
  if (hcf != nullptr) hcf->Fill(0.5, wt);
  for (Size_t i = 0; i < cutNames_.size() + 2; ++i) {
    if (i == 8) {
      cutName = "Trigger";
      pass = passTrg;
    } else if (i == 9) {
      cutName = "HEM";
      pass = passHEM;
    } else {
      cutName = cutNames_.at(j);
      pass = tuple_->TFvalue(cutName) > 0;
      if (cutName == "objcut" && !zeroChg) pass = false;
      if (cutName == "commonCuts" && !passEcalNoiseJetFilter) pass = false;
      j++;
    }
    if (hcf == nullptr) {
      cout << "Cut " << cutName << " = " << pass << endl;
      continue;
    }
    if (pass) {
      if (hcf != nullptr) hcf->Fill(1.5 + i, wt);
    } else if (isCutFlow) {
      break;
    }
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::cutHistos::fill(TH1D* hcf, Double_t wt, bool DiLepton, bool passPhotonVeto) {
  // Fill skim cut flow histograms.  If null histogram pointer passed, just print cut value.

  bool pass = false;
  bool isCutFlow = TString(hcf->GetName()).Contains("Flow");
  TString cutName;
  Size_t j = 0;
  if (hcf != nullptr) hcf->Fill(0.5, wt);
  for (Size_t i = 0; i < skimCutNames_.size() + 2; ++i) {
    if (i == 4) {
      cutName = "DiLepton";
      pass = DiLepton;
    } else if (i == 8) {
      cutName = "phoVeto";
      pass = passPhotonVeto;
    } else {
      cutName = skimCutNames_.at(j);
      pass = tuple_->TFvalue(cutName) > 0;
      j++;
    }
    if (hcf == nullptr) {
      cout << "Cut " << cutName << " = " << pass << endl;
      continue;
    }
    if (pass) {
      if (hcf != nullptr) hcf->Fill(1.5 + i, wt);
    } else if (isCutFlow) {
      break;
    }
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::checkTrigPrescales(const char* sample) {
  TString key = treeConfig_->DSkey(sample, isMC_);  // Set isMC_ here
  Tupl = new Ntuple(treeConfig_->getChain(key), treeConfig_->activeBranchList());
  Long64_t Nentries = Tupl->fChain->GetEntries();
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    Long64_t centry = Tupl->LoadTree(entry);
    if (centry < 0) {
      cout << "LoadTree returned " << centry;
      break;
    }
    Tupl->GetEntry(entry);
    // cout << endl << "First entry" << endl;  Show();  break;
    if (Tupl->newFileInChain()) {
      // New input root file encountered
      Tupl->setNewFileInChain(false);
      TFile* thisFile = Tupl->fChain->GetCurrentFile();
      if (thisFile) cout << "Current file in chain: " << thisFile->GetName() << endl;
      for (unsigned int ti = 0; ti < Tupl->TriggerNames->size(); ++ti) {
	Int_t prescale = Tupl->TriggerPrescales->at(ti);
	// if (prescale != 1)
	  cout << "Trigger " << Tupl->TriggerNames->at(ti)
	       << " (" << ti << ") prescaled by " << prescale << endl;
      }
    }
    delete Tupl;
  }

}  // ======================================================================================

void
RA2bZinvAnalysis::dumpSelEvIDs(const char* sample, const char* idFileName) {
  //
  // For debugging, where individual events need to be selected
  //
  FILE* idFile;
  idFile = fopen(idFileName, "w");

  TString key = treeConfig_->DSkey(sample, isMC_);  // Set isMC_ here
  Tupl = new Ntuple(treeConfig_->getChain(key), treeConfig_->activeBranchList());

  evSelector_ = new CutManager(sample, ntupleVersion_, isSkim_, isMC_, verbosity_,
			       era_, deltaPhi_, applyMassCut_, applyPtCut_, minNbCut_, restrictClean_, CCbins_);
  TCut baselineCuts = evSelector_->baseline();
  if (verbosity_ >= 1) cout << endl << "baseline = " << endl << baselineCuts << endl << endl;
  Tupl->setTF("baseline", baselineCuts);
  cutHistos cutHistFiller(Tupl, evSelector_);  // to dump cuts

  Tupl->setNotify();

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

  Long64_t Nentries = Tupl->fChain->GetEntries();
  if (verbosity_ >= 1) cout << "Nentries in tree = " << Nentries << endl;
  int count = 0;
  for (Long64_t entry = 0; entry < Nentries; ++entry) {
    count++;
    if (verbosity_ >= 1 && count % 100000 == 0) cout << "Entry number " << count << endl;
    Long64_t centry = Tupl->LoadTree(entry);
    if (centry < 0) {
      cout << "LoadTree returned " << centry;
      break;
    }
    Tupl->GetEntry(entry);
    // cout << endl << "First entry" << endl;  Show();  break;
    if (Tupl->newFileInChain()) {
      // New input root file encountered
      Tupl->setNewFileInChain(false);
      TFile* thisFile = Tupl->fChain->GetCurrentFile();
      if (thisFile) {
	if (verbosity_ >= 1) cout << "Current file in fChain: " << thisFile->GetName() << endl;
      }
    }

    if (std::find(EvtNums.begin(), EvtNums.end(), Tupl->EvtNum) == EvtNums.end()) continue;
    printf("%15u %15u %15llu\n", Tupl->RunNum, Tupl->LumiBlockNum, Tupl->EvtNum);
    fprintf(idFile, "%15u %15u %15llu\n", Tupl->RunNum, Tupl->LumiBlockNum, Tupl->EvtNum);

    cleanVars();  // If unskimmed input, copy <var>clean to <var>
    int currentYear = 0;  // FIXME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Int_t BTagsOrig = setBTags(currentYear);

    cout << "baseline = " << Tupl->TFvalue("baseline") << endl;
    cutHistFiller.fill(nullptr, 1, true, true, true, true);  // FIXME evaluate zeroChg, trigger, HEM, EcalNoiseJet selections
    cout << "globalTightHalo2016Filter = " << Tupl->globalTightHalo2016Filter << endl;
    cout << "globalSuperTightHalo2016Filter = " << Tupl->globalSuperTightHalo2016Filter << endl;
    cout << "HBHENoiseFilter = " << Tupl->HBHENoiseFilter << endl;
    cout << "HBHEIsoNoiseFilter = " << Tupl->HBHEIsoNoiseFilter << endl;
    cout << "eeBadScFilter = " << Tupl->eeBadScFilter << endl;
    cout << "EcalDeadCellTriggerPrimitiveFilter = " << Tupl->EcalDeadCellTriggerPrimitiveFilter << endl;
    cout << "BadChargedCandidateFilter = " << Tupl->BadChargedCandidateFilter << endl;
    cout << "BadPFMuonFilter = " << Tupl->BadPFMuonFilter << endl;
    cout << "NVtx = " << Tupl->NVtx << endl;
    cout << "PFCaloMETRatio = " << Tupl->PFCaloMETRatio << endl;
    cout << "JetID = " << Tupl->JetID << endl;
    cout << "HT = " << Tupl->HT << ", HT5 = " << Tupl->HT5 << ", HT5/HT = " << Tupl->HT5 / Tupl->HT << endl;
    // cout << "JetIDclean = " << Tupl->JetIDclean << endl;
    // cout << "HTclean = " << Tupl->HTclean << ", HT5clean = " << Tupl->HT5clean
    // 	 << ", HT5clean/HTclean = " << Tupl->HT5clean / Tupl->HTclean << endl;

    // if (Tupl->TFvalue("baseline")) {
      // if (std::find(EvtNums.begin(), EvtNums.end(), Tupl->EvtNum) != EvtNums.end()) {
	// cout << "Passes baseline" << endl;
	// printf("%15u %15u %15llu\n", Tupl->RunNum, Tupl->LumiBlockNum, Tupl->EvtNum);
	// fprintf(idFile, "%15u %15u %15llu\n", Tupl->RunNum, Tupl->LumiBlockNum, Tupl->EvtNum);
	// cout << " Tupl->PFCaloMETRatio = " << Tupl->PFCaloMETRatio << ", HT5/HT = " << Tupl->HT5 / Tupl->HT << endl;
      // }
    // }
  }
  fclose(idFile);
  delete Tupl;
  delete evSelector_;

}  // ======================================================================================
