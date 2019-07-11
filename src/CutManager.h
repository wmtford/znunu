//
//  Define cuts for RA2b analysis
//

#ifndef CUTMANAGER_H
#define CUTMANAGER_H

#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <TString.h>
#include <TCut.h>
#include "CCbinning.h"

#include <iostream>
using std::cout;
using std::endl;

class CutManager {

 public:
  CutManager(const TString sample, const TString ntupleVersion, bool isSkim, bool isMC, int verbosity,
	     std::string era, string deltaPhi, bool applyMassCut, bool applyPtCut, CCbinning* CCbins);
  virtual ~CutManager() {};

  enum yearFirstRun {Start2016 = 271036, Start2017 = 294645, Start2018 = 315252, StartHEM = 319077, Start2018C = 319313};
  typedef std::map<TString, TString> string_map;
  typedef std::map<TString, std::vector<TString> > vstring_map;

  void setTriggerIndexList(const char* sample, vector<unsigned>* triggerIndexList,
			   vector<string>* TriggerNames, vector<int>* TriggerPrescales);
  const string_map& sampleKeyMap() const {return sampleKeyMap_;};
  const vstring_map& triggerMapByName() const {return triggerMapByName_;};
  const TCut baseline() const {return cuts_;};
  const TString& HTcut() const {return HTcut_;};
  const TString& MHTcut() const {return MHTcut_;};
  const TString& NJetscut() const {return NJetscut_;};
  const TString& objcut() const {return objcut_;};
  const TString& minDphicut() const {return minDphicut_;};
  const TString& commonCuts() const {return commonCuts_;};
  const TString& ptCut() const {return ptCut_;};
  const TString& massCut() const {return massCut_;};
  const TString& photonDeltaRcut() const {return photonDeltaRcut_;};
  const TString& isoSFlepTksVeto() const {return isoSFlepTksVeto_;};
  const TString& isoSFlepTksCut() const {return isoSFlepTksCut_;};
  const TString& photonVeto() const {return photonVeto_;};
  const TString& photonCut() const {return photonCut_;};

 private:
  void fillCutMaps();
  TString ntupleVersion_;
  bool isSkim_;
  bool isMC_;
  int verbosity_;
  std::string era_;  // "2016", "Run2"
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  bool applyMassCut_;
  bool applyPtCut_;
  CCbinning* CCbins_;
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
