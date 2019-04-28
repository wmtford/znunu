//
//  Binning maps for RA2b analysis
//

#include "CCbinning.h"
#include <algorithm>  // for std::find()

ClassImp(CCbinning)

CCbinning::CCbinning(const std::string& era, const std::string& deltaPhi) :
era_(era), deltaPhi_(deltaPhi) {

  // Exclude some bins with high NJets, low HT
  std::vector< std::vector<unsigned> > exclBins;
  exclBins.push_back(std::vector<unsigned>({3, 0, 0}));  // j, m, h
  exclBins.push_back(std::vector<unsigned>({3, 1, 0}));
  exclBins.push_back(std::vector<unsigned>({3, 4, 0}));
  exclBins.push_back(std::vector<unsigned>({4, 0, 0}));
  exclBins.push_back(std::vector<unsigned>({4, 1, 0}));
  exclBins.push_back(std::vector<unsigned>({4, 4, 0}));

  if (era_ == TString("2016")) {

    kinThresholds_.push_back({300, 300, 500, 1000});  // mht threshold, {ht thresholds}
    kinThresholds_.push_back({350, 350, 500, 1000});
    kinThresholds_.push_back({500, 500, 1000});
    kinThresholds_.push_back({750, 750, 1500});
    kinThresholds_.push_back({250, 300, 500, 1000}); // QCD control bins

    jbThresholds_.push_back({2, 0, 1, 2});  // NJets threshold, {Nb thresholds}
    jbThresholds_.push_back({3, 0, 1, 2, 3});
    jbThresholds_.push_back({5, 0, 1, 2, 3});
    jbThresholds_.push_back({7, 0, 1, 2, 3});
    jbThresholds_.push_back({9, 0, 1, 2, 3});

    JbThresholds_.push_back({2, 0, 1, 2});  // for Nb/0b extrapolation
    JbThresholds_.push_back({3, 0, 1, 2, 3});
    JbThresholds_.push_back({4, 0, 1, 2, 3});
    JbThresholds_.push_back({5, 0, 1, 2, 3});
    JbThresholds_.push_back({6, 0, 1, 2, 3});
    JbThresholds_.push_back({7, 0, 1, 2, 3});
    JbThresholds_.push_back({8, 0, 1, 2, 3});
    JbThresholds_.push_back({9, 0, 1, 2, 3});

  } else if (era_ == TString("Run2")) {
    // From Alexx Perloff RA2b talk, 8 Jan 2019
    kinThresholds_.push_back({300, 300, 600, 1200});  // mht threshold, {ht thresholds}
    kinThresholds_.push_back({350, 350, 600, 1200});
    kinThresholds_.push_back({600, 600, 1200});
    kinThresholds_.push_back({850, 850, 1700});
    kinThresholds_.push_back({250, 300, 600, 1200}); // QCD control bins

    jbThresholds_.push_back({2, 0, 1, 2});  // NJets threshold, {Nb thresholds}
    jbThresholds_.push_back({4, 0, 1, 2, 3});
    jbThresholds_.push_back({6, 0, 1, 2, 3});
    jbThresholds_.push_back({8, 0, 1, 2, 3});
    jbThresholds_.push_back({10, 0, 1, 2, 3});

    JbThresholds_.push_back({2, 0, 1, 2});  // for Nb/0b extrapolation
    JbThresholds_.push_back({3, 0, 1, 2});
    JbThresholds_.push_back({4, 0, 1, 2, 3});
    JbThresholds_.push_back({5, 0, 1, 2, 3});
    JbThresholds_.push_back({6, 0, 1, 2, 3});
    JbThresholds_.push_back({7, 0, 1, 2, 3});
    JbThresholds_.push_back({8, 0, 1, 2, 3});
    JbThresholds_.push_back({9, 0, 1, 2, 3});
    JbThresholds_.push_back({10, 0, 1, 2, 3});

  } // era

  //  The following loop fills the binning maps
  //  toCCbin_ for the main analysis NJet, Nb, (HT, MHT) histogram
  //  toCCbinjb_ for kinematics-integrated NJet, Nb histogram
  //  toCCbinSpl_ for NJet (bin width 1), Nb, (HT, MHT) histogram
  //  toCCbinJb_ for kinematics-integrated NJet (bin width 1), Nb histogram
  //  toCCbinjk_ for Nb = 0 NJet, (HT, MHT) histogram
  Int_t binJbk = 0, binjbk = 0, binJb = 0, binjb = 0, binjk = 0;
  jetSubBins_.resize(jbThresholds_.size());
  unsigned j = 0;
  for (unsigned J = 0; J < JbThresholds_.size(); ++J) {
    jetSubBins_[j].push_back(J);
    for (unsigned ib = 1; ib < jbThresholds_[j].size(); ++ib) {
      unsigned b = jbThresholds_[j][ib];
      std::vector<int> Jb = {int(J), int(b)}, jb = {int(j), int(b)};
      binJb++;
      toCCbinJb_[Jb] = binJb;
      if (JbThresholds_[J][0] == jbThresholds_[j][0]) {
	binjb++;
	toCCbinjb_[jb] = binjb;
      }
      unsigned mmax = deltaPhi_.find("nominal") != std::string::npos ? kinThresholds_.size()-1 : kinThresholds_.size();
      int k = -1;
      for (unsigned m = 0; m < mmax; ++m) {
	for (unsigned h = 1; h < kinThresholds_[m].size(); ++h) {
	  k++;
	  if (std::find(exclBins.begin(), exclBins.end(), std::vector<unsigned>({j, m, h-1})) != exclBins.end()) continue;
	  std::vector<int> Jbk = {int(J), int(b), k};
	  binJbk++;
	  toCCbinSpl_[Jbk] = binJbk;
	  if (JbThresholds_[J][0] == jbThresholds_[j][0]) {
	    std::vector<int> jbk = {int(j), int(b), k};
	    binjbk++;
	    toCCbin_[jbk] = binjbk;
	    // cout << "Filling toCCbin; j = " << j << ", b = "  << b << ", k = " << k << ", bin = " << toCCbin_[jbk] << endl;
	    if (b == 0) {
	      std::vector<int> jk = {int(j), k};
	      binjk++;
	      toCCbinjk_[jk] = binjk;
	    }
	  }
	}
      }
    }
    if (J+1 < JbThresholds_.size() && j+1 < jbThresholds_.size() && JbThresholds_[J+1][0] == jbThresholds_[j+1][0]) j++;
  }

  kinSize_ = 0;  kinSizeNominal_ = 0;
  int mmaxNominal = kinThresholds_.size() - 1;
  int mmax = deltaPhi_.find("nominal") != std::string::npos ? kinThresholds_.size() - 1 :
    kinThresholds_.size(); // No. of MHT bands
  for (int i = 0; i < mmax; ++i) {
    kinSize_ += kinThresholds_[i].size() - 1;  // First component is MHT value
    if (i < mmaxNominal) kinSizeNominal_ += kinThresholds_[i].size() - 1;
  }

}  // ======================================================================================

int
CCbinning::kinBin(double& ht, double& mht) {
  int theBin = -1;
  if (era_ == "Run2" && mht > ht) return theBin;
  int NmhtBins = kinThresholds_.size() - 1;
  if (deltaPhi_.find("nominal") == std::string::npos) {
    // ldp or hdp
    if (mht < kinThresholds_[NmhtBins][0] || ht < kinThresholds_[NmhtBins][1]) return theBin;
    if (mht < kinThresholds_[0][0]) {
      theBin += kinSizeNominal_;
      for (unsigned j = 1; j < kinThresholds_[NmhtBins].size(); ++j) {
	theBin++;
	if (ht >= kinThresholds_[NmhtBins][j] &&
	    (j == kinThresholds_[NmhtBins].size()-1 || ht < kinThresholds_[NmhtBins][j+1])) return theBin;
      }
    }
  }
  if (mht < kinThresholds_[0][0] || ht < kinThresholds_[0][1]) return theBin;
  for (int i = 0; i < NmhtBins; ++i) {
    if (mht >= kinThresholds_[i][0] && (i == NmhtBins-1 || mht < kinThresholds_[i+1][0])) {
      for (unsigned j = 1; j < kinThresholds_[i].size(); ++j) {
	theBin++;
	if (ht >= kinThresholds_[i][j] && (j == kinThresholds_[i].size()-1 || ht < kinThresholds_[i][j+1])) return theBin;
      }
    } else {
      theBin += kinThresholds_[i].size() - 1;
    }
  }
  return -2;  // Outside binned area (e.g., MHT > HT), though within preselection
}  // ======================================================================================
