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
  Ntuple(TTree* tree = 0) : Ntuple(tree, nullptr) {};
  Ntuple(TTree* tree, Bool_t doBuildCache) : Ntuple(tree, nullptr, doBuildCache) {}
  Ntuple(TTree* tree, const vector<const char*> *activeBranches,
	 Bool_t doBuildCache = true) : NtupleClass(tree) {
    forNotify_ = new TObjArray;
    forNotify_->SetOwner();  // so that TreeFormulas will be deleted
    if (activeBranches != nullptr) setActiveBranches(activeBranches);
    if (doBuildCache) buildCache(activeBranches);
  };
  virtual ~Ntuple() {
    delete forNotify_;
  };
  // Record when a new file in the chain is encountered:
  Bool_t Notify() override {newFileInChain_ = kTRUE;  return(kTRUE);};
  // User access to new file status:
  const Bool_t& newFileInChain() const {return newFileInChain_;};
  // User reset new file status:
  void setNewFileInChain(const Bool_t newFile) {newFileInChain_ = newFile;};
  // User register a new TTreeFormula:
  void setTF(const TString name, const TString formula) {
    TTreeFormula* theTF = new TTreeFormula(name, formula, fChain);
    TFmap_[name] = theTF;
    forNotify_->Add(theTF);
  };
  // User call after all TTreeFormulas registered, so they will be notified of a new file:
  void setNotify() {fChain->SetNotify(forNotify_);};
  // User access to the value of a TTreeFormula, by name:
  double TFvalue(const TString key) {
    if (TFmap_.count(key) == 0) {
      cout << "Ntuple::TFvalue: no key matching " << key << endl;
      return -1;
    }
    TTreeFormula* theTF = TFmap_.at(key);
    theTF->GetNdata();
    return theTF->EvalInstance(0);
  };

// ====================================================================================== */

private:
  Bool_t newFileInChain_;
  std::map<TString, TTreeFormula*> TFmap_;
  TObjArray* forNotify_;

  void setActiveBranches(const vector<const char*> *activeBranches) {
    fChain->SetBranchStatus("*", 0);  // disable all branches
    for (auto theBranch : *activeBranches) {
      fChain->SetBranchStatus(theBranch, 1);
    }
  };
  void buildCache(const vector<const char*> *activeBranches) {
    cout << "Initial size of cache for fChain = " << fChain->GetCacheSize() << endl;
    TTreeCache::SetLearnEntries(1);
    fChain->SetCacheSize(200*1024*1024);
    fChain->SetCacheEntryRange(0, fChain->GetEntries());
    if (activeBranches == nullptr) {
      fChain->AddBranchToCache("*", true);
    } else {
      for (auto theBranch : *activeBranches) fChain->AddBranchToCache(theBranch, true);
    }
    fChain->StopCacheLearningPhase();
    cout << "Reset size of cache for fChain = " << fChain->GetCacheSize() << endl;
  };

  /* ClassDef(Ntuple, 1) // 2nd arg is ClassVersionID */
};

#endif
