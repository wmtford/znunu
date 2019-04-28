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
  CCbinning(const std::string& era = "2016", const std::string& deltaPhi = "nominal");
  virtual ~CCbinning() {};
  int kinBin(double& ht, double& mht);
  unsigned kinSize() {return kinSize_;};

  std::vector<int> nJetThresholds() {
    std::vector<int> jThresh;
    for (Size_t j = 0; j < jbThresholds_.size(); ++j) jThresh.push_back(jbThresholds_[j][0]);
    return jThresh;
  };
  int nJetThreshold(int j) {return jbThresholds_[j][0];};
  std::vector<int> nJet1Thresholds() {
    std::vector<int> jThresh;
    for (Size_t j = 0; j < JbThresholds_.size(); ++j) jThresh.push_back(JbThresholds_[j][0]);
    return jThresh;
  };
  int nJet1Threshold(int j) {return JbThresholds_[j][0];};
  std::vector<int> nbThresholds(int jbin) {
    std::vector<int> bThresh;
    for (Size_t b = 1; b < jbThresholds_[jbin].size(); ++b) bThresh.push_back(jbThresholds_[jbin][b]);
    return bThresh;
  };
  int nbThreshold(int j, int b) {return jbThresholds_[j][b];};
  std::vector< std::vector<int> > jetSubBins() {return jetSubBins_;};

  typedef std::map<std::vector<int>, Int_t> ivector_map;
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

  std::vector< std::vector<double> > kinThresholds() {return kinThresholds_;};
  double mhtThreshold(int m) {return kinThresholds_[m][0];};
  double htThreshold(int m, int h) {return kinThresholds_[m][1+h];};
  Size_t binsmht() {return kinThresholds_.size();};
  Size_t binsht(int m) {return kinThresholds_[m].size() - 1;};

  int jbk(int j, int b, int k) {
    std::vector<int> jbk = {j, b, k};
    try {toCCbin_.at(jbk);}  catch (const std::out_of_range& oor) {return -1;}
    return(toCCbin_.at(jbk));
  };
  int jb(int j, int b) {
    std::vector<int> jb = {j, b};
    try {toCCbinjb_.at(jb);}  catch (const std::out_of_range& oor) {return -1;}
    return(toCCbinjb_.at(jb));
  };
  int jk(int j, int k) {
    std::vector<int> jk = {j, k};
    try {toCCbinjk_.at(jk);}  catch (const std::out_of_range& oor) {return -1;}
    return(toCCbinjk_.at(jk));
  };
  int Jbk(int J, int b, int k) {
    std::vector<int> Jbk = {J, b, k};
    try {toCCbinSpl_.at(Jbk);}  catch (const std::out_of_range& oor) {return -1;}
    return(toCCbinSpl_.at(Jbk));
  };
  int Jb(int J, int b) {
    std::vector<int> Jb = {J, b};
    try {toCCbinJb_.at(Jb);}  catch (const std::out_of_range& oor) {return -1;}
    return(toCCbinJb_.at(Jb));
  };

private:
  std::string era_;  // "2016", ...
  std::string deltaPhi_;  // "nominal", "hdp", "ldp", "ldpnominal"
  std::vector< std::vector<double> > kinThresholds_;
  std::vector< std::vector<int> > jbThresholds_;
  std::vector< std::vector<int> > JbThresholds_;
  unsigned kinSize_;
  unsigned kinSizeNominal_;
  std::vector< std::vector<int> > jetSubBins_;
  ivector_map toCCbin_, toCCbinSpl_, toCCbinJb_, toCCbinjb_, toCCbinjk_;

  ClassDef(CCbinning, 1) // 2nd arg is ClassVersionID
};

#endif
