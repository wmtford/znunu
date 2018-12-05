//
//  Binning maps for RA2b analysis
//

#ifndef CCBINNING_H
#define CCBINNING_H

#define _CRT_SECURE_NO_WARNINGS

#include <map>
#include <TString.h>

class CCbinning {

public:
  CCbinning(const std::string& era = "2016", const std::string& deltaPhi = "nominal");
  virtual ~CCbinning() {};
  std::vector< std::vector<double> > kinThresholds() {return kinThresholds_;};
  std::vector<int> nJetThresholds() {return nJetThresholds_;};
  std::vector<int> nJet1Thresholds() {return nJet1Thresholds_;};
  std::vector<int> nbThresholds() {return nbThresholds_;};
  unsigned kinSize() {return kinSize_;};
  std::vector< std::vector<int> > jetSubBins() {return jetSubBins_;};
  typedef std::map<std::vector<int>, Int_t> ivector_map;
  ivector_map toCCbin() {return toCCbin_;};
  ivector_map toCCbinjb() {return toCCbinjb_;};
  ivector_map toCCbinjk() {return toCCbinjk_;};
  ivector_map toCCbinSpl() {return toCCbinSpl_;};
  ivector_map toCCbinJb() {return toCCbinJb_;};
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
    try {toCCbinjb_.at(jk);}  catch (const std::out_of_range& oor) {return -1;}
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
  std::string deltaPhi_;  // "nominal", "hdp", "ldp"
  std::vector< std::vector<double> > kinThresholds_;
  std::vector<int> nJetThresholds_;
  std::vector<int> nJet1Thresholds_;
  std::vector<int> nbThresholds_;
  unsigned kinSize_;
  std::vector< std::vector<int> > jetSubBins_;
  ivector_map toCCbin_, toCCbinSpl_, toCCbinJb_, toCCbinjb_, toCCbinjk_;

  ClassDef(CCbinning, 1) // 2nd arg is ClassVersionID
};

#endif
