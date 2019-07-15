//
//  Select root files and active branches, and build chain
//

#ifndef TREECONFIG_H
#define TREECONFIG_H

#define _CRT_SECURE_NO_WARNINGS

#include <TString.h>
#include <TChain.h>
#include <TRegexp.h>

#include <iostream>
using std::cout;
using std::endl;

class TreeConfig {

public:
  TreeConfig(const string& era, const TString& ntupleVersion, const bool& isSkim,
	     const string& deltaPhi, const int& verbosity, const string& treeName, const string& treeLoc,
	     const string& fileListsFile, const string& runBlock);
  virtual ~TreeConfig() {
    delete activeBranches;
  };
  TString DSkey(const char* dataSet, bool& isMC);
  TChain* getChain(TString key);
  vector<const char*>* activeBranchList();

private:
  string era_;  // "2016", "Run2", ...
  TString ntupleVersion_;
  bool isSkim_;
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  int verbosity_;
  bool isMC_;
  string treeName_;
  string treeLoc_;
  string fileListsFile_;
  string runBlock_;
  vector<const char*>* activeBranches;

  vector<TString> fileList(TString sampleKey);

  /* ClassDef(TreeConfig, 1) // 2nd arg is ClassVersionID */
};

#endif
