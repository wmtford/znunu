//
//  Accessor class for a tree.
//  Optimize for a user-provided vector of active branches.
//  Manage TTreeFormulas.
//

#ifndef NTUPLE_H
#define NTUPLE_H

#ifndef NtupleClass_cxx
#define NtupleClass_cxx
#endif

#define _CRT_SECURE_NO_WARNINGS

#include "NtupleClass.h"
#include <TTreeCache.h>
#include <TTreeFormula.h>

#include <iostream>
using std::cout;
using std::endl;

void NtupleClass::Loop() {}

class Ntuple : public NtupleClass {

public:
 Ntuple() : NtupleClass(0) {};
 Ntuple(TTree *tree=0) : NtupleClass(tree) {
    forNotify_ = new TObjArray;
    forNotify_->SetOwner();  // so that TreeFormulas will be deleted
  };
  virtual ~Ntuple() {
    delete forNotify_;
  };

  Bool_t Notify() override {newFileInChain_ = kTRUE;  return(kTRUE);};

  const Bool_t& newFileInChain() const {return newFileInChain_;};

  void setNewFileInChain(const Bool_t newFile) {newFileInChain_ = newFile;};

  void setTF(const TString name, const TString formula) {
    TTreeFormula* theTF = new TTreeFormula(name, formula, fChain);
    TFmap_[name] = theTF;
    forNotify_->Add(theTF);
  };

  void setNotify() {fChain->SetNotify(forNotify_);};

  double TFvalue(const TString key) {
    if (TFmap_.count(key) == 0) {
      cout << "Ntuple::TFvalue: no key matching " << key << endl;
      return -1;
    }
    TTreeFormula* theTF = TFmap_.at(key);
    theTF->GetNdata();
    return theTF->EvalInstance(0);
  };

  void optimizeTree(vector<const char*> activeBranches, const bool activateAll = false) {
  cout << "activateAllBranches = " << activateAll << endl;

  if (!activateAll) {
    fChain->SetBranchStatus("*", 0);  // disable all branches
    // cout << "Active branches:" << endl;
    for (auto theBranch : activeBranches) {
      // cout << theBranch << endl;
      fChain->SetBranchStatus(theBranch, 1);
    }
  }
  cout << "Initial size of cache for fChain = " << fChain->GetCacheSize() << endl;
  TTreeCache::SetLearnEntries(1);
  fChain->SetCacheSize(200*1024*1024);
  fChain->SetCacheEntryRange(0, fChain->GetEntries());
  if (activateAll) {
    fChain->AddBranchToCache("*", true);
  } else {
    for (auto theBranch : activeBranches) fChain->AddBranchToCache(theBranch, true);
  }
  fChain->StopCacheLearningPhase();
  cout << "Reset size of cache for fChain = " << fChain->GetCacheSize() << endl;

}  // ======================================================================================

  void cleanVars() {
    // For TreeMaker these replacements are made when skims are produced,
    // so this should be called only for raw ntuple data sets.
    NJets = NJetsclean;
    BTags = BTagsclean;
    BTagsDeepCSV = BTagsDeepCSVclean;
    HT = HTclean;
    HT5 = HT5clean;
    MHT = MHTclean;
    JetID = JetIDclean;
    Jets = Jetsclean;
    Jets_hadronFlavor = Jetsclean_hadronFlavor;
    Jets_HTMask = Jetsclean_HTMask;
    isoElectronTracks = isoElectronTracksclean;
    isoMuonTracks = isoMuonTracksclean;
    isoPionTracks = isoPionTracksclean;
    DeltaPhi1 = DeltaPhi1clean;
    DeltaPhi2 = DeltaPhi2clean;
    DeltaPhi3 = DeltaPhi3clean;
    DeltaPhi4 = DeltaPhi4clean;
  };

private:
  Bool_t newFileInChain_;
  std::map<TString, TTreeFormula*> TFmap_;
  TObjArray* forNotify_;

  /* ClassDef(Ntuple, 1) // 2nd arg is ClassVersionID */
};

#endif
