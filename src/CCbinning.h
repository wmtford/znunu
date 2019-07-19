//
//  Binning maps for RA2b analysis
//

#ifndef CCBINNING_H
#define CCBINNING_H

#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <TString.h>

#include <iostream>
using std::cout;
using std::endl;

class CCbinning {

public:
  CCbinning(const string& era = "Run2", const string& deltaPhi = "nominal");
  virtual ~CCbinning() {};
  int kinBin(double& ht, double& mht);
  unsigned kinSize() {return kinSize_;};

  vector<int> nJetThresholds() {
    vector<int> jThresh;
    for (Size_t j = 0; j < jbThresholds_.size(); ++j) jThresh.push_back(jbThresholds_[j][0]);
    return jThresh;
  };
  int nJetThreshold(int j) {return jbThresholds_[j][0];};
  vector<int> nJet1Thresholds() {
    vector<int> jThresh;
    for (Size_t j = 0; j < JbThresholds_.size(); ++j) jThresh.push_back(JbThresholds_[j][0]);
    return jThresh;
  };
  int nJet1Threshold(int j) {return JbThresholds_[j][0];};
  vector<int> nbThresholds(int jbin) {
    vector<int> bThresh;
    for (Size_t b = 1; b < jbThresholds_[jbin].size(); ++b) bThresh.push_back(jbThresholds_[jbin][b]);
    return bThresh;
  };
  int nbThreshold(int j, int b) {return jbThresholds_[j][b];};
  vector< vector<int> > jetSubBins() {return jetSubBins_;};

  typedef map<vector<int>, Int_t> ivector_map;
  ivector_map toCCbin() {return toCCbin_;};
  ivector_map toCCbinjb() {return toCCbinjb_;};
  ivector_map toCCbinjk() {return toCCbinjk_;};
  ivector_map toCCbinSpl() {return toCCbinSpl_;};
  ivector_map toCCbinJb() {return toCCbinJb_;};
  Int_t bins() {return toCCbin_.size();};
  Int_t binsjb() {return toCCbinjb_.size();};
  Int_t binsjk() {return toCCbinjk_.size();};
  Int_t binsSpl() {return toCCbinSpl_.size();};
  Int_t binsJb() {return toCCbinJb_.size();};
  Size_t binsj() {return jbThresholds_.size();};
  Size_t binsJ() {return JbThresholds_.size();};
  Size_t binsb(int j) {return jbThresholds_[j].size()-1;};
  Size_t binsB(int J) {return JbThresholds_[J].size()-1;};
  int jbin(int nJets) {
    if (nJets < jbThresholds_[0][0]) return -1;
    int bin = jbThresholds_.size()-1;
    while (nJets < jbThresholds_[bin][0]) bin--;
    return bin;
  };
  int Jbin(int nJets) {
    if (nJets < JbThresholds_[0][0]) return -1;
    int bin = JbThresholds_.size()-1;
    while (nJets < JbThresholds_[bin][0]) bin--;
    return bin;
  };
  int bbin(int nJets, int Nb) {
    if (nJets < jbThresholds_[0][0]) return -1;
    int jbin = jbThresholds_.size()-1;
    while (nJets < jbThresholds_[jbin][0]) jbin--;
    int bbin = jbThresholds_[jbin].size()-1;
    while (Nb < jbThresholds_[jbin][bbin]) bbin--;
    return bbin - 1;
  };

  vector< vector<double> > kinThresholds() {return kinThresholds_;};
  double mhtThreshold(int m) {return kinThresholds_[m][0];};
  double htThreshold(int m, int h) {return kinThresholds_[m][1+h];};
  Size_t binsmht() {return kinThresholds_.size();};
  Size_t binsht(int m) {return kinThresholds_[m].size() - 1;};

  int jbk(int j, int b, int k) {
    if (k < 0) return -1;
    vector<int> jbk = {j, b, k};
    ivector_map::iterator iter = toCCbin_.find(jbk);
    if (iter == toCCbin_.end()) return -1;
    return iter->second;
  };
  int jb(int j, int b) {
    vector<int> jb = {j, b};
    ivector_map::iterator iter = toCCbinjb_.find(jb);
    if (iter == toCCbinjb_.end()) return -1;
    return iter->second;
  };
  int jk(int j, int k) {
    if (k < 0) return -1;
    vector<int> jk = {j, k};
    ivector_map::iterator iter = toCCbinjk_.find(jk);
    if (iter == toCCbinjk_.end()) return -1;
    return iter->second;
  };
  int Jbk(int J, int b, int k) {
    if (k < 0) return -1;
    vector<int> Jbk = {J, b, k};
    ivector_map::iterator iter = toCCbinSpl_.find(Jbk);
    if (iter == toCCbinSpl_.end()) return -1;
    return iter->second;
  };
  int Jb(int J, int b) {
    vector<int> Jb = {J, b};
    ivector_map::iterator iter = toCCbinJb_.find(Jb);
    if (iter == toCCbinJb_.end()) return -1;
    return iter->second;
  };

private:
  string era_;  // "2016", ...
  string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  vector< vector<double> > kinThresholds_;
  vector< vector<int> > jbThresholds_;
  vector< vector<int> > JbThresholds_;
  unsigned kinSize_;
  unsigned kinSizeNominal_;
  vector< vector<int> > jetSubBins_;
  ivector_map toCCbin_, toCCbinSpl_, toCCbinJb_, toCCbinjb_, toCCbinjk_;

  /* ClassDef(CCbinning, 1) // 2nd arg is ClassVersionID */
};

#endif
