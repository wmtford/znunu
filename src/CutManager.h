//
//  Define cuts for RA2b analysis
//

#ifndef CUTMANAGER_H
#define CUTMANAGER_H

#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <TString.h>

#include <iostream>
using std::cout;
using std::endl;

class CutManager {

 public:
  CutManager(const TString sample, const TString ntupleVersion, bool isSkim, bool isMC);
  virtual ~CutManager() {};
  TCut baseline() {return cuts_;};

 private:
  void fillCutMaps();
  TString ntupleVersion_;
  bool isSkim_;
  bool isMC_;
  vstring_map triggerMapByName_;
  string_map objCutMap_;
  string_map minDphiCutMap_;
  string_map MHTCutMap_;
  string_map sampleKeyMap_;
  TCut cuts_;
  TString HTcut_;
  TString MHTcut_;
  TString NJetscut_;
  TString objcut_;
  TString minDphicut_;
  TString commonCuts_;
  TString ptCut_;
  TString massCut_;
  TString photonDeltaRcut_;
  // Following 4 cuts for special studies
  TString isoSFlepTksVeto_;
  TString isoSFlepTksCut_;
  TString photonVeto_;
  TString photonCut_;

  /* ClassDef(CutManager, 1) // 2nd arg is ClassVersionID */
};

#endif
