//
//  Weights to correct MC for efficiencies, data for purity
//

#ifndef EFFICWT_H
#define EFFICWT_H

#define _CRT_SECURE_NO_WARNINGS

#include <TFile.h>
#include "CCbinning.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TEfficiency.h>
#include <TLorentzVector.h>

#include <TMath.h>
using namespace TMath;

class EfficWt {
 public:

  enum runYear{Year2016 = 0, Year2017 = 1, Year2018 = 2, Year2018HEP = 3, Year2018HEM = 4};

  EfficWt() : deltaPhi_("nominal") {};
  EfficWt(string deltaPhi);
  ~EfficWt() {};

  void openFiles();
  void getHistos(const char* sample, int currentYear);
  pair<double, double> weight(CCbinning* CCbins,
			      Int_t NJets, Int_t BTags, Double_t MHT, Double_t HT,
			      vector<TLorentzVector> ZCandidates,
			      vector<TLorentzVector> Photons,
			      vector<TLorentzVector> Electrons,
			      vector<TLorentzVector> Muons,
			      vector<double> EBphoton,
			      bool applyDRfitWt,
			      int currentYear);

  double prefiring_weight_photon(vector<TLorentzVector>* Photons, unsigned p){
    if (hPrefiring_photon_ == nullptr) return 1;
    return (1 - hPrefiring_photon_->GetBinContent(hPrefiring_photon_->GetXaxis()->FindBin(Photons->at(p).Eta()),
						  hPrefiring_photon_->GetYaxis()->FindBin(Photons->at(p).Pt())));
  };
  double prefiring_weight_electron(vector<TLorentzVector>* Electrons, unsigned e){
    if (hPrefiring_photon_ == nullptr) return 1;
    return (1 - hPrefiring_photon_->GetBinContent(hPrefiring_photon_->GetXaxis()->FindBin(Electrons->at(e).Eta()),
						  hPrefiring_photon_->GetYaxis()->FindBin(Electrons->at(e).Pt())));
  };
  double prefiring_weight_jet(vector<TLorentzVector>* Jets, unsigned j){
    if (hPrefiring_jet_ == nullptr) return 1;
    return (1 - hPrefiring_jet_->GetBinContent(hPrefiring_jet_->GetXaxis()->FindBin(Jets->at(j).Eta()),
					       hPrefiring_jet_->GetYaxis()->FindBin(Jets->at(j).Pt()))) ;
  };
  double quadSum(double x, double y) {return Sqrt(Power(x,2) + Power(y,2));};
  double quadSum(double x, double y, double z) {return Sqrt(Power(x,2) + Power(y,2) + Power(z,2));};
 private:
  std::vector<TFile*> purityTrigEffFile_;
  std::vector<TFile*> photonTrigEffFile_;
  std::vector<TFile*> photonSFFile_;
  std::vector<TFile*> elecIDandIsoSFFile_;
  std::vector<TFile*> elecRecoLowSFFile_;
  std::vector<TFile*> elecRecoHighSFFile_;
  std::vector<TFile*> muonIDSFFile_;
  std::vector<TFile*> muonIsoSFFile_;
  TFile* prefiringWeightFile_;
  TH2F *hPrefiring_photon_, *hPrefiring_jet_;
  TString theSample_;
  string deltaPhi_;
  std::vector<TH1F*> hPurity_, hTrigEff_;
  std::vector<TF1*> fTrigEff_;
  std::vector<TEfficiency*> eTrigEff_;
  TF1* DRfun_;
  std::vector< std::vector<Double_t> > DRpars_;

  std::vector<TH2F*> hSFeff_;
  TGraphErrors* FdirGraph_;

  /* ClassDef(EfficWt, 1) // 2nd arg is ClassVersionID */
};

#endif
