//
//  Read TChain from files
//

#ifndef TREEANALYSISBASE_H
#define TREEANALYSISBASE_H

#define _CRT_SECURE_NO_WARNINGS

#include <TChain.h>
#include <TChainElement.h>
/* #include <TTreeCache.h> */
#include <TString.h>
#include <TRegexp.h>
#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <sstream>
using std::stringstream;
using std::string;

#include <stdio.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

class TreeAnalysisBase {

public:

  TreeAnalysisBase() : chain_(nullptr) {TreeConfig("none");};
  TreeAnalysisBase(const char* dataSet, const string& cfg_filename = "",
		   const string& runBlock = "") : runBlock_(runBlock) {
    TreeConfig(dataSet, cfg_filename);
    if (!runBlock.empty()) runBlock_ = runBlock;
    cout << "The runBlock is " << runBlock_ << endl;
    chain_ = getChain(dataSet);
  };
  virtual ~TreeAnalysisBase() {}

  void TreeConfig(const char* dataSet, const string& cfg_filename="") {

    // Set up configuration, using boost/program_options.
    po::options_description desc("Config");
    desc.add_options()
      ("Tree.verbosity", po::value<int>(&verbosity_))
      ("Tree.runBlock", po::value<string>(&runBlock_))  // May be overridden by constructor
      ("Tree.tree path", po::value<string>(&treeLoc_))
      ("Tree.root file index", po::value<string>(&fileListsFile_))
      ("Tree.delta phi sample", po::value<string>(&deltaPhi_), "nominal, hdp, ldp, ldpnominal")

      ;
    po::variables_map vm;
    std::ifstream cfi_file("RA2bZinvAnalysis.cfi");
    po::store(po::parse_config_file(cfi_file, desc, true), vm);  // true means ignore unregistered options
    po::notify(vm);
    if (!cfg_filename.empty()) {
      std::ifstream cfg_file(cfg_filename);
      vm = po::variables_map();  // Clear the map.
      po::store(po::parse_config_file(cfg_file, desc, true), vm);  // true means ignore unregistered options
      po::notify(vm);
    }

    //  Extract tree signature from the path
    string s = treeLoc_, delimiter = "Run2Production";
    size_t pos = s.find(delimiter);
    if (pos == string::npos) cout << "Can't find tree version; " << delimiter << " not in tree path" << endl;
    pos += delimiter.length();
    ntupleVersion_ = s.substr(pos, 3);

    isSkim_ = treeLoc_.find("Skim") != string::npos;
    if (!isSkim_) {
      treeName_ = "TreeMaker2/PreSelection";  // For ntuple
    } else {
      treeName_ = "tree";  // For skims
    }

    isMC_ = strncmp("zmm", dataSet, 3) && strncmp("zee", dataSet, 3) && strncmp("photon", dataSet, 3);

  }  // ======================================================================================

  TChain* getChain(const char* dataSet) {

    /* bool activateAllBranches = false;  // Can be set true for debugging */
    /* cout << "activateAllBranches = " << activateAllBranches << endl; */

    TString theSample = dataSet;
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
    else if (theSample.Contains("zmm")) key = TString("zmm");  // ( and not "tt")
    else if (theSample.Contains("zee")) key = TString("zee");  // ( and not "tt")
    else if (theSample.Contains("photon")) key = TString("photon");
    else if (theSample.Contains("gjetsqcd")) key = TString("gjetsqcd");
    else if (theSample.Contains("gjets")) key = TString("gjets");  // ( and not "qcd")
    else {
      cout << "getChain:  unknown dataSet '" << dataSet << "'" << endl;
      return nullptr;
    }
    if (deltaPhi_.find("ldp") != string::npos && isSkim_) key += "ldp";
    if (!runBlock_.empty()) key += runBlock_;  key("HEM") = "";

    TChain* chain = new TChain(treeName_.data());
    std::vector<TString> files = fileList(key);
    for (auto file : files) {
      if (verbosity_ >= 2) cout << file << endl;
      chain->Add(file);
    }

    return chain;

  }  // ======================================================================================

  std::vector<TString> fileList(TString sampleKey) {

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

protected:
  TString ntupleVersion_;
  bool isSkim_;
  bool isMC_;
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  int verbosity_;
  string treeName_;
  string treeLoc_;
  string fileListsFile_;
  string runBlock_;
  TChain* chain_;

};

#endif
