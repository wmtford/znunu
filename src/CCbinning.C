//
//  Binning maps for RA2b analysis
//

#include "CCbinning.h"

ClassImp(CCbinning)

CCbinning::CCbinning(const std::string& era, const std::string& deltaPhi) :
era_(era), deltaPhi_(deltaPhi) {

  if (era_ == TString("2016")) {

    kinThresholds_.push_back({300, 300, 500, 1000});  // mht threshold, {ht thresholds}
    kinThresholds_.push_back({350, 350, 500, 1000});
    kinThresholds_.push_back({500, 500, 1000});
    kinThresholds_.push_back({750, 750, 1500});
    kinThresholds_.push_back({250, 300, 500, 1000}); // QCD control bins
    nJet1Thresholds_ = {2, 3, 4, 5, 6, 7, 8, 9};  // for Nb/0b extrapolation
    nJetThresholds_ = {2, 3, 5, 7, 9};
    nbThresholds_ = {0, 1, 2, 3};

  } // era 2016

  //  The following loop fills the binning maps
  //  toCCbin_ for the main analysis NJet, Nb, (HT, MHT) histogram
  //  toCCbinSpl_ with one NJet value per NJet bin, for Nb/0b extrapolation
  //  toCCbinjb_ for kinematics-integrated NJet, Nb histogram
  Int_t binJbk = 0, binjbk = 0, binJb = 0, binjb = 0;
  jetSubBins_.resize(nJetThresholds_.size());
  unsigned j = 0;
  for (unsigned J = 0; J < nJet1Thresholds_.size(); ++J) {
    jetSubBins_[j].push_back(J);
    for (unsigned b = 0; b < nbThresholds_.size(); ++b) {
      if (nbThresholds_[b] > nJetThresholds_[j]) continue;  // Exclude Nb > NJets
      std::vector<int> Jb = {int(J), int(b)}, jb = {int(j), int(b)};
      binJb++;
      toCCbinJb_[Jb] = binJb;
      if (nJet1Thresholds_[J] == nJetThresholds_[j]) {
	binjb++;
	toCCbinjb_[jb] = binjb;
      }
      unsigned mmax = deltaPhi_ == TString("nominal") ? kinThresholds_.size()-1 : kinThresholds_.size();
      int k = -1;
      for (unsigned m = 0; m < mmax; ++m) {
	for (unsigned h = 1; h < kinThresholds_[m].size(); ++h) {
	  k++;
	  if (j > 2 && (m < 2 || m == kinThresholds_.size()-1) && h == 1) continue;   // Exclude (Njets3,4; HT0,3,(6))
	  std::vector<int> Jbk = {int(J), int(b), k};
	  binJbk++;
	  toCCbinSpl_[Jbk] = binJbk;
	  if (nJet1Thresholds_[J] == nJetThresholds_[j]) {
	    std::vector<int> jbk = {int(j), int(b), k};
	    binjbk++;
	    toCCbin_[jbk] = binjbk;
	    // cout << "Filling toCCbin; j = " << j << ", b = "  << b << ", k = " << k << ", bin = " << toCCbin_[jbk] << endl;
	  }
	}
      }
    }
    if (J+1 < nJet1Thresholds_.size() && j+1 < nJetThresholds_.size() && nJet1Thresholds_[J+1] == nJetThresholds_[j+1]) j++;
  }

  kinSize_ = 0;
  int mmax = (deltaPhi_ == TString("nominal")) ? kinThresholds_.size() -1 : kinThresholds_.size(); // No. of MHT bands
  for (int i = 0; i < mmax; ++i)
    kinSize_ += kinThresholds_[i].size() - 1;  // First component is MHT value
}
