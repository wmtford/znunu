//
//  Define cuts for RA2b analysis
//

#ifndef CUTMANAGER_H
#define CUTMANAGER_H

#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <TString.h>
#include <TRegexp.h>
#include <TCut.h>
#include "Ntuple.h"
#include "CCbinning.h"

#include <iostream>

class CutManager {

 public:
  CutManager(const TString sample, const TString ntupleVersion, bool isSkim, bool isMC,
	     int verbosity, string era, string deltaPhi, bool applyMassCut,
	     bool applyPtCut, int minNb, bool restrictClean, CCbinning* CCbins);
  virtual ~CutManager() {};

  enum yearFirstRun {Start2016 = 271036, Start2017 = 294645, Start2018 = 315252, StartHEM = 319077, Start2018C = 319313};
  typedef map<TString, TString> string_map;
  typedef map<TString, vector<TString> > vstring_map;

  void setTriggerIndexList(const char* sample, vector<unsigned>* triggerIndexList, Ntuple* ntuple);
  const string_map& sampleKeyMap() const {return sampleKeyMap_;};
  const vstring_map& triggerMapByName() const {return triggerMapByName_;};
  const TString skimCut(const char* cutStringClean, bool restrictClean = false) const {
    TString cutClean = cutStringClean;
    if (!(restrictClean || restrictClean_)) return (cutClean);
    TString cut = cutClean;
    while (cut.Index("clean") != -1) cut("clean") = "";
    return("((@Jetsclean.size()!=0 && " + cutClean + ") || (@Jetsclean.size()==0 && " + cut + "))");
  };
  const TCut baseline() const {return cuts_;};
  const TString& HTcut() const {return HTcut_;};
  const TString& MHTcut() const {return MHTcut_;};
  const TString& NJetscut() const {return NJetscut_;};
  const TString& Nbcut() const {return Nbcut_;};
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
  const TString& NJetSkimCut() const {return NJetSkimCut_;};
  const TString& HTSkimCut() const {return HTSkimCut_;};
  const TString& MHTSkimCut() const {return MHTSkimCut_;};
  const TString& MHTHTRatioSkimCut() const {return MHTHTRatioSkimCut_;};
  /* const TString& DiMuonSkimCut() const {return DiMuonSkimCut_;}; */
  const TString& ElectronVetoSkimCut() const {return ElectronVetoSkimCut_;};
  const TString& IsoElectronTrackVetoSkimCut() const {return IsoElectronTrackVetoSkimCut_;};
  const TString& IsoPionTrackVetoSkimCut() const {return IsoPionTrackVetoSkimCut_;};
  /* const TString& PhotonSkimCut() const {return PhotonSkimCut_;}; */
  const TString& DeltaPhiSkimCut() const {return DeltaPhiSkimCut_;};
  const TString& JetIDSkimCut() const {return JetIDSkimCut_;};

 private:
  void fillCutMaps();
  TString ntupleVersion_;
  bool isSkim_;
  bool isMC_;
  int verbosity_;
  string era_;  // "2016", "Run2"
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  bool applyMassCut_;
  bool applyPtCut_;
  bool restrictClean_;
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
  TString Nbcut_;
  TString objcut_;
  TString minDphicut_;
  TString commonCuts_;
  TString ptCut_;
  TString massCut_;
  TString photonDeltaRcut_;

  // Following are cuts used when skims are made, DiMuon sample
  TString NJetSkimCut_;
  TString HTSkimCut_;
  TString MHTSkimCut_;
  TString MHTHTRatioSkimCut_;
  /* TString DiMuonSkimCut_; */
  TString ElectronVetoSkimCut_;
  TString IsoElectronTrackVetoSkimCut_;
  TString IsoPionTrackVetoSkimCut_;
  /* TString PhotonSkimCut_; */
  TString DeltaPhiSkimCut_;
  TString JetIDSkimCut_;

  // Following 4 cuts for special studies
  TString isoSFlepTksVeto_;
  TString isoSFlepTksCut_;
  TString photonVeto_;
  TString photonCut_;

  /* ClassDef(CutManager, 1) // 2nd arg is ClassVersionID */
};

#endif
