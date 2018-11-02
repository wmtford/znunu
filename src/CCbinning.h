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
  typedef std::map<std::vector<int>, Int_t> ivector_map;
  ivector_map toCCbin() {return toCCbin_;};
  ivector_map toCCbinSpl() {return toCCbinSpl_;};
  ivector_map toCCbinjb() {return toCCbinjb_;};

private:
  std::string era_;  // "2016", ...
  std::string deltaPhi_;  // "nominal", "hdp", "ldp"
  std::vector< std::vector<double> > kinThresholds_;
  std::vector<int> nJetThresholds_;
  std::vector<int> nJet1Thresholds_;
  std::vector<int> nbThresholds_;
  unsigned kinSize_;
  ivector_map toCCbin_, toCCbinSpl_, toCCbinjb_;

  ClassDef(CCbinning, 1) // 2nd arg is ClassVersionID
};

#endif
