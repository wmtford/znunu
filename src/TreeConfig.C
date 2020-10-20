//
//  Select root files and active branches, and build chain
//

#include "TreeConfig.h"

#include <fstream>
using std::ifstream;

// ClassImp(TreeConfig)

TreeConfig::TreeConfig(const string& era, const TString& ntupleVersion, const bool& isSkim,
		       const string& deltaPhi, const int& verbosity, const string& treeName,
		       const string& treeLoc, const string& fileListsFile, const string& runBlock) :
  era_(era), ntupleVersion_(ntupleVersion), isSkim_(isSkim), deltaPhi_(deltaPhi), verbosity_(verbosity),
  treeName_(treeName), treeLoc_(treeLoc), fileListsFile_(fileListsFile), runBlock_(runBlock) {

}  // ======================================================================================

TString
TreeConfig::DSkey(const char* dataSet, bool& isMC) {

  TString theSample = dataSet;
  TString key;
  isMC_ = true;  // Depends on dataSet
  if      (theSample.Contains("zinv")) key = "zinv";
  else if (theSample.Contains("ttzvv")) key = "ttzvv";
  else if (theSample.Contains("dymm")) key = "dymm";
  else if (theSample.Contains("dyee")) key = "dyee";
  else if (theSample.Contains("ttzmm")) key = isSkim_ ? "ttzmm" : "ttz";  
  else if (theSample.Contains("ttzee")) key = isSkim_ ? "ttzee" : "ttz";
  else if (theSample.Contains("VVmm")) key = "VVmm";
  else if (theSample.Contains("VVee")) key = "VVee";
  else if (theSample.Contains("ttmm")) key = "ttmm";
  else if (theSample.Contains("ttee")) key = "ttee";
  else if (theSample.Contains("zmm")) {
    key = "zmm";  // ( and not "tt")
    isMC_ = false;
  }
  else if (theSample.Contains("zee")) {
    key = "zee";  // ( and not "tt")
    isMC_ = false;
  }
  else if (theSample.Contains("photon")) {
    key = "photon";
    isMC_ = false;
  }
  else if (theSample.Contains("gjetsqcd")) key = "gjetsqcd";
  else if (theSample.Contains("gjets")) key = "gjets";  // ( and not "qcd")
  else {
    cout << "DSkey:  unknown dataSet '" << dataSet << "'" << endl;
    return nullptr;
  }
  if (deltaPhi_.find("ldp") != string::npos && isSkim_) key += "ldp";
  if (!runBlock_.empty()) key += runBlock_;  key("HEM") = "";

  isMC = isMC_;
  return key;

}  // ======================================================================================

TChain*
TreeConfig::getChain(TString key) {
  // key is the data set identifier in the file lists file, from DSkey().

  TChain* chain = new TChain(treeName_.data());
  std::vector<TString> files = fileList(key);
  for (auto file : files) {
    if (verbosity_ >= 2) cout << file << endl;
    chain->Add(file);
  }

  return chain;

}  // ======================================================================================

std::vector<TString>
TreeConfig::fileList(TString sampleKey) {

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

vector<const char*>*
TreeConfig::activeBranchList() {

  // Needed branches
  activeBranches = new vector<const char*>;
  activeBranches->push_back("NJets");
  activeBranches->push_back("BTags");
  activeBranches->push_back("HT");
  activeBranches->push_back("HT5");
  activeBranches->push_back("MHT");
  activeBranches->push_back("MHTPhi");
  activeBranches->push_back("JetID");
  activeBranches->push_back("Jets");
  activeBranches->push_back("Jets_hadronFlavor");
  activeBranches->push_back("Jets_HTMask");
  activeBranches->push_back("isoElectronTracks");
  activeBranches->push_back("isoMuonTracks");
  activeBranches->push_back("isoPionTracks");
  activeBranches->push_back("DeltaPhi1");
  activeBranches->push_back("DeltaPhi2");
  activeBranches->push_back("DeltaPhi3");
  activeBranches->push_back("DeltaPhi4");
  if (isSkim_) {
    activeBranches->push_back("RA2bin");
    activeBranches->push_back("METRatioFilter");
    activeBranches->push_back("EcalNoiseJetFilter");
    activeBranches->push_back("HTRatioDPhiFilter");
    activeBranches->push_back("HTRatioFilter");
  } else {
    activeBranches->push_back("NJetsclean");
    activeBranches->push_back("BTagsclean");
    activeBranches->push_back("BTagsDeepCSVclean");
    activeBranches->push_back("HTclean");
    activeBranches->push_back("HT5clean");
    activeBranches->push_back("MHTclean");
    activeBranches->push_back("MHTPhiclean");
    activeBranches->push_back("JetIDclean");
    activeBranches->push_back("Jetsclean");
    activeBranches->push_back("Jetsclean_hadronFlavor");
    activeBranches->push_back("Jetsclean_HTMask");
    activeBranches->push_back("isoElectronTracksclean");
    activeBranches->push_back("isoMuonTracksclean");
    activeBranches->push_back("isoPionTracksclean");
    activeBranches->push_back("DeltaPhi1clean");
    activeBranches->push_back("DeltaPhi2clean");
    activeBranches->push_back("DeltaPhi3clean");
    activeBranches->push_back("DeltaPhi4clean");
  }
  if (ntupleVersion_ != "V12") {
    activeBranches->push_back("NMuons");
    activeBranches->push_back("NElectrons");
    activeBranches->push_back("BTagsDeepCSV");
    activeBranches->push_back("ecalBadCalibFilter");
  }
  activeBranches->push_back("RunNum");
  activeBranches->push_back("LumiBlockNum");
  activeBranches->push_back("EvtNum");
  activeBranches->push_back("Jets_bDiscriminatorCSV");
  activeBranches->push_back("Muons");
  activeBranches->push_back("Muons_passIso");
  activeBranches->push_back("Muons_mediumID");
  activeBranches->push_back("Muons_charge");
  activeBranches->push_back("Electrons");
  activeBranches->push_back("Electrons_passIso");
  activeBranches->push_back("Electrons_mediumID");
  activeBranches->push_back("Electrons_charge");
  activeBranches->push_back("ZCandidates");
  activeBranches->push_back("Photons");
  activeBranches->push_back("Photons_nonPrompt");
  activeBranches->push_back("Photons_fullID");
  activeBranches->push_back("Photons_hasPixelSeed");
  activeBranches->push_back("Photons_isEB");
  activeBranches->push_back("NVtx");
  if (ntupleVersion_ != "V18") activeBranches->push_back("TriggerNames");
  activeBranches->push_back("TriggerPass");
  activeBranches->push_back("TriggerPrescales");
  activeBranches->push_back("HBHENoiseFilter");
  activeBranches->push_back("HBHEIsoNoiseFilter");
  activeBranches->push_back("eeBadScFilter");
  activeBranches->push_back("EcalDeadCellTriggerPrimitiveFilter");
  activeBranches->push_back("globalTightHalo2016Filter");
  activeBranches->push_back("globalSuperTightHalo2016Filter");
  activeBranches->push_back("BadChargedCandidateFilter");
  activeBranches->push_back("BadPFMuonFilter");
  activeBranches->push_back("PFCaloMETRatio");
  activeBranches->push_back("nAllVertices");
  if (isMC_) {
    activeBranches->push_back("puWeight");
    activeBranches->push_back("Weight");
    activeBranches->push_back("TrueNumInteractions");
    activeBranches->push_back("madMinPhotonDeltaR");
    activeBranches->push_back("GenMuons");
    activeBranches->push_back("GenElectrons");
    activeBranches->push_back("GenHT");
    activeBranches->push_back("GenTaus");
    activeBranches->push_back("GenParticles");
    activeBranches->push_back("GenParticles_PdgId");
    activeBranches->push_back("GenParticles_Status");
    if (ntupleVersion_ != "V12" && ntupleVersion_ != "V15") {
      activeBranches->push_back("NonPrefiringProb");
      // activeBranches->push_back("NonPrefiringProbUp");
      // activeBranches->push_back("NonPrefiringProbDown");
    }
  }

  return activeBranches;

}  // ======================================================================================
