//
//  BTag data/MC scale factors
//  (wrapper for https://github.com/kpedro88/Analysis/btag/BTagCorrector)
//

#ifndef BTAGSF_H
#define BTAGSF_H

#include <TH2F.h>
// Following assumes that ../../Analysis is a clone of
//   https://github.com/kpedro88/Analysis/blob/SUSY2018/
#include "../../Analysis/btag/BTagCorrector.h"

class BTagSF {
 public:

  enum runYear{Year2016 = 0, Year2017 = 1, Year2018 = 2, Year2018HEP = 3, Year2018HEM = 4};
  enum workingpoint{Loose = 0, Medium = 1, Tight = 2};

  BTagSF(bool useDeepCSV=true) {
    btagcorr_ = new BTagCorrector;
    // Scale factor files, from https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation/
    BTagSFfile_.push_back(useDeepCSV ? string("../datFiles/DeepCSV_2016LegacySF_WP_V1_mod.csv")  // 2016
			  : string("../../Analysis/btag/CSVv2_Moriond17_B_H_mod.csv"));  // 2016 CSVv2
    BTagSFfile_.push_back(string("../datFiles/DeepCSV_94XSF_WP_V4_B_F_mod.csv"));  // 2017
    BTagSFfile_.push_back(string("../datFiles/DeepCSV_102XSF_WP_V1_mod.csv"));  // 2018
    BTagSFfile_.push_back(string("../datFiles/DeepCSV_102XSF_WP_V1_mod.csv"));  // 2018HEP
    BTagSFfile_.push_back(string("../datFiles/DeepCSV_102XSF_WP_V1_mod.csv"));  // 2018HEM
    // DeepCSV working points
    BTagWPval_.push_back({0.2217, 0.6321, 0.8953});  // 2016 {Loose, Medium, Tight}
    BTagWPval_.push_back({0.1522, 0.4941, 0.8001});  // 2017
    BTagWPval_.push_back({0.1241, 0.4184, 0.7527});  // 2018
    BTagWPval_.push_back({0.1241, 0.4184, 0.7527});  // 2018HEP
    BTagWPval_.push_back({0.1241, 0.4184, 0.7527});  // 2018HEM
  };

  ~BTagSF() {delete btagcorr_;};

  void SetCalib(unsigned runYearIndex, unsigned wp=1) {
    btagcorr_->SetCalib(BTagSFfile_.at(runYearIndex));
    currentWP_ = (BTagWPval_.at(runYearIndex)).at(wp);
  }
  double WPvalue() {return currentWP_;};
  void SetEffs(TFile* dataFile) {btagcorr_->SetEffs(dataFile);}
  double weight(vector<TLorentzVector>* Jets, vector<int>* Jets_flavor, vector<bool>* Jets_HTMask,
		vector<double>* Jets_bDiscriminator) {
    double sf = btagcorr_->GetSimpleCorrection(Jets, Jets_flavor, Jets_HTMask, Jets_bDiscriminator, currentWP_);
    if (isnan(sf)) {
      cout << "BTagCorrector returns NaN" << endl;
      return 1.0;
    } else {
      return sf;
    }
  };
  vector<double> GetCorrections(vector<TLorentzVector>* Jets, vector<int>* Jets_flavor,
				vector<bool>* Jets_HTMask) {
    return btagcorr_->GetCorrections(Jets, Jets_flavor, Jets_HTMask);
  };

 private:
  BTagCorrector* btagcorr_;
  vector<string> BTagSFfile_;
  vector< vector<double> > BTagWPval_;
  double currentWP_;

  /* ClassDef(BTagSF, 1) // 2nd arg is ClassVersionID */
};

#endif
