/*
  Run from the command line with
  root -l -b Nb0bLoadClasses.C Nb0bExtrap.C
*/

#include "CCbinning.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TApplication.h"
#include "TMath.h"

using namespace std;
using namespace TMath;

void Nb0bExtrap(const string& era = "Run2", const string& deltaPhi = "nominal") {
  ofstream datFile("DY_signal.dat");
  // ofstream datFile;
  ofstream latexFile("DY_signal.tex");
  TFile* openFile(const char* fileName);
  TH1D* getHist(TFile* tFile, const char* histName);
  enum runBlock {Y2016, Y2017, Y2018AB, Y2018CD, Run2};

  int doRun = Y2017;
  bool doClosure = true;
  bool useDYMC = true;
  bool useZllData = false;
  bool usePhotonData = false;
  bool useMCJfactors = false;
  bool doJfromData = false;
  cout << "case " << doRun << ", doClosure = " << doClosure << ", useDYMC = " << useDYMC
       << ", useZllData = " << useZllData << ", usePhotonData = " << usePhotonData << endl;

  CCbinning CCbins(era, deltaPhi);
  int kinSize = CCbins.kinSize();
  int NbinsNJet = CCbins.binsj();
  vector< vector<int> > jetSubBins = CCbins.jetSubBins();
  int extrapByMCthreshold = NbinsNJet - 1;  // Use MC for the last NJets bin
  // vector<int> extrapFromRange = {5, 6};
  float relErrXsec_ttz = 0.3;

  TFile *ZllData, *photonData, *ZllXMC, *photonMC, *zinvMC;

  switch(doRun) {
  case Y2016:
    if (doClosure) {
      zinvMC = openFile("../outputs/histsZjets_2016v16_noPU.root");  if (zinvMC == nullptr) return;
    }
    else {
      ZllData = openFile("../outputs/histsDY_2016v16.root");  if (ZllData == nullptr) return;
      photonData = openFile("../outputs/histsPhoton_2016v16.root");  if (photonData == nullptr) return;
    }
    ZllXMC = openFile("../outputs/histsDYMC_2016v16_noPU.root");  if (ZllXMC == nullptr) return;
    break;

  case Y2017:
    if (doClosure) {
      zinvMC = openFile("../outputs/histsZjets_2017v16_HT17wt_ZptWt.root");  if (zinvMC == nullptr) return;
    }
    else {
      ZllData = openFile("../outputs/histsDY_2017v16.root");  if (ZllData == nullptr) return;
      photonData = openFile("../outputs/histsPhoton_2017v16.root");  if (photonData == nullptr) return;
    }
    ZllXMC = openFile("../outputs/histsDYMC_2017v16_HT17wt_ZptWt.root");  if (ZllXMC == nullptr) return;
    break;

  case Y2018AB:
    if (doClosure) {
      zinvMC = openFile("../outputs/histsZjets_2018v16_HT17wt_ZptWt.root");  if (zinvMC == nullptr) return;
    }
    else {
      ZllData = openFile("../outputs/histsDY_2018ABv16.root");  if (ZllData == nullptr) return;
      photonData = openFile("../outputs/histsPhoton_2018ABv16.root");  if (photonData == nullptr) return;
    }
    ZllXMC = openFile("../outputs//histsDYMC_2018v16_HT17wt_ZptWt.root");  if (ZllXMC == nullptr) return;
    break;

  case Y2018CD:
    if (doClosure) {
      zinvMC = openFile("../outputs/histsZjets_2018HEMv16_HT17wt_ZptWt.root");  if (zinvMC == nullptr) return;
    }
    else {
      ZllData = openFile("../outputs/histsDY_2018CDv16.root");  if (ZllData == nullptr) return;
      photonData = openFile("../outputs/histsPhoton_2018CDv16.root");  if (photonData == nullptr) return;
    }
    ZllXMC = openFile("../outputs//histsDYMC_2018HEMv16_HT17wt_ZptWt.root");  if (ZllXMC == nullptr) return;
    break;

  case Run2:
    if (doClosure) {
      if (useZllData) {
	ZllData = openFile("../outputs/histsDY_Run2v16.root");	if (ZllData == nullptr) return;
      } else if (usePhotonData) {
	photonData = openFile("../outputs/histsPhoton_Run2v16.root");  if (photonData == nullptr) return;
	// photonMC = openFile("../outputs/histsGjets_Run2v16_noPU.root");  if (photonMC == nullptr) return;
      } else {
	zinvMC = openFile("../outputs/histsZjets_Run2v16_HT17wt_ZptWt_noPU.root");  if (zinvMC == nullptr) return;
      }
    } else {
      ZllData = openFile("../outputs/histsDY_Run2v16.root");  if (ZllData == nullptr) return;
      photonData = openFile("../outputs/histsPhoton_Run2v16.root");  if (photonData == nullptr) return;
    }
    ZllXMC = openFile("../outputs//histsDYMC_Run2v16_HT17wt_ZptWt_noPU.root");  if (ZllXMC == nullptr) return;
    break;
    
  default:
    cout << "No valid run year specified" << endl;
    return;
  }

  // Get the Z->ll and photon data and MC histograms
  TH1D *hCC_zmm, *hCC_zee, *hCC_zll, *hCCjb_zmm, *hCCjb_zee, *hCCjb_zll, *hCCJb_zmm, *hCCJb_zee, *hCCJb_zll;
  TH1D *hCCJb_dymm, *hCCJb_dyee, *hCCJb_ttzmm, *hCCJb_ttzee, *hCCJb_VVmm, *hCCJb_VVee, *hCCJb_ttmm, *hCCJb_ttee;
  TH1D *hCCjb_dymm, *hCCjb_dyee, *hCCjb_ttzmm, *hCCjb_ttzee, *hCCjb_VVmm, *hCCjb_VVee, *hCCjb_ttmm, *hCCjb_ttee;
  TH1D *hCCspl_photon, *hCC_photon, *hCCJb_MCall, *hCCjb_MCall;
  if (doClosure) {
    if (useZllData) {
      hCC_zmm = getHist(ZllData, "hCC_zmm");
      hCC_zee = getHist(ZllData, "hCC_zee");
      hCC_zll = getHist(ZllData, "hCC_zee");
      hCCjb_zmm = getHist(ZllData, "hCCjb_zmm");
      hCCjb_zee = getHist(ZllData, "hCCjb_zee");
      hCCjb_zll = getHist(ZllData, "hCCjb_zll");
      hCCJb_zmm = getHist(ZllData, "hCCJb_zmm");
      hCCJb_zee = getHist(ZllData, "hCCJb_zee");
      hCCJb_zll = getHist(ZllData, "hCCJb_zll");
      hCCspl_photon = getHist(ZllData, "hCCspl_zll");
      hCC_photon = getHist(ZllData, "hCC_zll");
    } else if (usePhotonData) {
      hCC_zmm = getHist(photonData, "hCC_photon");
      hCC_zee = getHist(photonData, "hCC_photon");
      hCC_zll = getHist(photonData, "hCC_photon");
      hCCjb_zmm = getHist(photonData, "hCCjb_photon");
      hCCjb_zee = getHist(photonData, "hCCjb_photon");
      hCCjb_zll = getHist(photonData, "hCCjb_photon");
      hCCJb_zmm = getHist(photonData, "hCCJb_photon");
      hCCJb_zee = getHist(photonData, "hCCJb_photon");
      hCCJb_zll = getHist(photonData, "hCCJb_photon");
      hCCspl_photon = getHist(photonData, "hCCspl_photon");
      hCC_photon = getHist(photonData, "hCC_photon");
    } else {
      hCC_zmm = getHist(ZllXMC, "hCC_dymm");
      hCC_zee = getHist(ZllXMC, "hCC_dyee");
      hCC_zll = getHist(ZllXMC, "hCC_dyee");
      hCCjb_zmm = getHist(ZllXMC, "hCCjb_dymm");
      hCCjb_zee = getHist(ZllXMC, "hCCjb_dyee");
      hCCjb_zll = getHist(ZllXMC, "hCCjb_dyll");
      hCCJb_zmm = getHist(ZllXMC, "hCCJb_dymm");
      hCCJb_zee = getHist(ZllXMC, "hCCJb_dyee");
      hCCJb_zll = getHist(ZllXMC, "hCCJb_dyll");
      if (useDYMC) {
	hCCspl_photon = getHist(ZllXMC, "hCCspl_dyll");
	hCC_photon = getHist(ZllXMC, "hCC_dyll");
      } else {  // Use Zinv MC
	hCCspl_photon = getHist(zinvMC, "hCCspl_zinv");
	hCC_photon = getHist(zinvMC, "hCC_zinv");
      }
    }
    hCC_photon->Print("all");
  } else {  // Not closure
    hCC_zmm = getHist(ZllData, "hCC_zmm");
    hCC_zee = getHist(ZllData, "hCC_zee");
    hCCjb_zmm = getHist(ZllData, "hCCjb_zmm");
    hCCjb_zee = getHist(ZllData, "hCCjb_zee");
    hCCjb_zll = getHist(ZllData, "hCCjb_zll");
    hCCJb_zmm = getHist(ZllData, "hCCJb_zmm");
    hCCJb_zee = getHist(ZllData, "hCCJb_zee");
    hCCJb_zll = getHist(ZllData, "hCCJb_zll");
    hCCspl_photon = getHist(photonData, "hCCspl_photon");
    hCC_photon = getHist(photonData, "hCC_photon");
  }

  // Fetch and combine DY MC contributions, split
  hCCJb_dymm = getHist(ZllXMC, "hCCJb_dymm");
  hCCJb_dyee = getHist(ZllXMC, "hCCJb_dyee");
  hCCJb_ttzmm = getHist(ZllXMC, "hCCJb_ttzmm");
  hCCJb_ttzee = getHist(ZllXMC, "hCCJb_ttzee");
  hCCJb_VVmm = getHist(ZllXMC, "hCCJb_VVmm");
  hCCJb_VVee = getHist(ZllXMC, "hCCJb_VVee");
  hCCJb_ttmm = getHist(ZllXMC, "hCCJb_ttmm");
  hCCJb_ttee = getHist(ZllXMC, "hCCJb_ttee");
  hCCJb_MCall = (TH1D*) hCCJb_dymm->Clone();  hCCJb_MCall->SetName("hCCJb_MCall");
  hCCJb_MCall->Add(hCCJb_dyee);
  hCCJb_MCall->Add(hCCJb_ttzmm);
  hCCJb_MCall->Add(hCCJb_ttzee);
  hCCJb_MCall->Add(hCCJb_VVmm);
  hCCJb_MCall->Add(hCCJb_VVee);
  // hCCJb_MCall->Add(hCCJb_ttmm);  // Omit non-peaking, since data are purity corrected.
  // hCCJb_MCall->Add(hCCJb_ttee);

  // Fetch and combine DY MC contributions, unsplit
  hCCjb_dymm = getHist(ZllXMC, "hCCjb_dymm");
  hCCjb_dyee = getHist(ZllXMC, "hCCjb_dyee");
  hCCjb_ttzmm = getHist(ZllXMC, "hCCjb_ttzmm");
  hCCjb_ttzee = getHist(ZllXMC, "hCCjb_ttzee");
  hCCjb_VVmm = getHist(ZllXMC, "hCCjb_VVmm");
  hCCjb_VVee = getHist(ZllXMC, "hCCjb_VVee");
  hCCjb_ttmm = getHist(ZllXMC, "hCCjb_ttmm");
  hCCjb_ttee = getHist(ZllXMC, "hCCjb_ttee");
  hCCjb_MCall = (TH1D*) hCCjb_dymm->Clone();  hCCjb_MCall->SetName("hCCjb_MCall");
  hCCjb_MCall->Add(hCCjb_dyee);
  hCCjb_MCall->Add(hCCjb_ttzmm);
  hCCjb_MCall->Add(hCCjb_ttzee);
  hCCjb_MCall->Add(hCCjb_VVmm);
  hCCjb_MCall->Add(hCCjb_VVee);
  hCCjb_MCall->Add(hCCjb_ttmm);
  hCCjb_MCall->Add(hCCjb_ttee);
  // For calculation of the ttz cross section uncertainty from the ttz fraction.
  TH1D* hCCjb_MCttzFrac = (TH1D*) hCCjb_ttzmm->Clone();
  hCCjb_MCttzFrac->SetNameTitle("hCCjb_MCttzFrac", "hCCjb_MCttzFrac");
  hCCjb_MCttzFrac->Add(hCCjb_ttzee);
  hCCjb_MCttzFrac->Divide(hCCjb_MCall);
  // hCCjb_ttzmm->Print("all");
  // hCCjb_ttzee->Print("all");
  // hCCjb_MCttzFrac->Print("all");

  // Calculate Z->ll data stat errors
  vector<float> DYstat;
  for (int j = 0; j < NbinsNJet; ++j) {
    float b0 = 0, b0err = 0, bb = 0, bberr = 0;
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      int jj = j < extrapByMCthreshold || !useMCJfactors ? j : j-1;
      int  binCCjb = CCbins.jb(jj, b);  if (binCCjb <= 0) continue;
      if (b == 0) {
	b0 = hCCjb_zll->GetBinContent(binCCjb);
	b0err = hCCjb_zll->GetBinError(binCCjb);
	DYstat.push_back(0);
      } else {
	bb = hCCjb_zll->GetBinContent(binCCjb);
	bberr = hCCjb_zll->GetBinError(binCCjb);
	if (b0 > 0 && bb > 0)
	  if (doClosure && !(useZllData || usePhotonData))  // Use sqrt(bin wt) as error
	    DYstat.push_back(Sqrt(1/b0 + 1/bb));
	  else
	    DYstat.push_back(Sqrt(Power((b0err/b0), 2) + Power((bberr/bb), 2)));
	else
	  DYstat.push_back(1);
      }
    }
  }

  // Compute the purity uncertainties
  TFile* effFile = openFile("../plots/histograms/effHists.root");  if (effFile == nullptr) return;
  // TFile* effFile = openFile("../python/effHists.root");  if (effFile == nullptr) return;
    
  TH1D* h_pur_m = getHist(effFile, "h_pur_m");
  TH1D* h_pur_e = getHist(effFile, "h_pur_e");
  vector<float> DYpurSys;
  // Declare the bin merging map for purities
  vector< vector<int> > pbins;
  pbins.push_back({1});
  pbins.push_back({2});
  pbins.push_back({3});
  pbins.push_back({4});
  pbins.push_back({5});
  pbins.push_back({6, 7});
  pbins.push_back({6, 7});
  pbins.push_back({8, 12, 16});
  pbins.push_back({9, 13, 17});
  pbins.push_back({10, 11, 14, 15, 18, 19});
  pbins.push_back({10, 11, 14, 15, 18, 19});
  pbins.push_back({8, 12, 16});
  pbins.push_back({9, 13, 17});
  pbins.push_back({10, 11, 14, 15, 18, 19});
  pbins.push_back({10, 11, 14, 15, 18, 19});
  pbins.push_back({8, 12, 16});
  pbins.push_back({9, 13, 17});
  pbins.push_back({10, 11, 14, 15, 18, 19});
  pbins.push_back({10, 11, 14, 15, 18, 19});
  // Compute the relative purity error averaged over Zmm, Zee, and bin groups
  for (int i = 0; i < (int) pbins.size(); ++i) {
    float y = 0, ype = 0;
    for (int j = 0; j < (int) pbins.at(i).size(); ++j) {
      int bin = pbins.at(i).at(j);
      y += hCCjb_zmm->GetBinContent(bin) + hCCjb_zee->GetBinContent(bin);
      ype += hCCjb_zmm->GetBinContent(bin)*h_pur_m->GetBinError(bin)/h_pur_m->GetBinContent(bin)
       	   + hCCjb_zee->GetBinContent(bin)*h_pur_e->GetBinError(bin)/h_pur_e->GetBinContent(bin);
    }
    DYpurSys.push_back(ype/y);
  }

  vector<float> systJ = {0, 0, 0, 0};
  vector<pair<float, float>> Jextrap;

  // Calculate the MC J factors to extend the Nb/0b F factors
  // Take the Nphoton average of F over NJet sub-bins at each kin bin.
  vector<pair<float, float>> jb0;
  for (int j = 0; j <= extrapByMCthreshold; ++j) {
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	int binCC = CCbins.jbk(j, b, k);  if (binCC <= 0) continue;
	float Npho = 0;
	pair<float, float> jbb, NphoF;
	for (auto J : jetSubBins.at(j)) {
	  int binCCJb = CCbins.Jb(J, b);  if (binCCJb <= 0) continue;
	  int binCCjb = CCbins.jb(j, b);  if (binCCjb <= 0) continue;
	  pair<float, float> nMC;
	  if (doJfromData) {
	    nMC.first = hCCJb_zll->GetBinContent(binCCJb);  // for Jextrap systematic
	    nMC.second = hCCJb_zll->GetBinError(binCCJb);  //  "
	  } else {
	    nMC.first = hCCJb_MCall->GetBinContent(binCCJb);
	    nMC.second = hCCJb_MCall->GetBinError(binCCJb);
	  }
	  if (J >= (int) jb0.size()) {
	    jb0.push_back(nMC);
	    // cout << "j, J, b, k, jb0Size = " << j << ", " << J << ", " << b << ", " << k << ", " << jb0.size() << endl;
	  } else {
	    jbb.first = nMC.first / jb0.at(J).first;
	    jbb.second = nMC.first > 0 ? Sqrt(Power(nMC.second/nMC.first, 2) + Power(jb0.at(J).second/jb0.at(J).first, 2)) : 0;
	    int binCCJ0k = CCbins.Jbk(J, 0, k);  if (binCCJ0k <= 0) continue;
	    float NphoJ = hCCspl_photon->GetBinContent(binCCJ0k);
	    Npho += NphoJ;
	    NphoF.first += NphoJ*jbb.first;
	    NphoF.second += NphoJ*jbb.second;
	    // cout << "j, J, b, k, jbb, NphoJ, Npho = " << j << ", " << J << ", " << b << ", " << k << ", "
	    // 	 << jbb.first << ", relErr = " << jbb.second << ", " << NphoJ << ", " << Npho << endl;
	  }
	}
	if (b == 0) {
	  Jextrap.push_back(pair<float, float>(1, 0));
	} else {
	  if (Npho > 0)
	    Jextrap.push_back(pair<float, float>(NphoF.first/Npho, NphoF.second/Npho));
	  else
	    Jextrap.push_back(pair<float, float>(jbb.first, jbb.second));
	}
	if (j >= extrapByMCthreshold) {
	  // Compute ratio of highest to next-highest NJets bins
	  float relErrjbk = Jextrap.at(binCC - 1).second;
	  float relErrjm1bk = Jextrap.at(CCbins.jbk(j-1, b, k) - 1).second;
	  Jextrap.at(binCC - 1).first /= Jextrap.at(CCbins.jbk(j-1, b, k) - 1).first;
	  Jextrap.at(binCC - 1).second = Sqrt(Power(relErrjbk, 2) + Power(relErrjm1bk, 2));
	  // cout << "j, b, k, Jextrap = " << j << ", " << b << ", " << k << ", " <<
	  //   Jextrap.at(binCC - 1).first << ", rel err = " << Jextrap.at(binCC - 1).second << endl;
	}
      }
    }
  }

  // Jextrap systematics derived from analysis of the above for data, MC
  systJ = {0, 0.10, 0.10, 0.20};

  // Calculate Nb/0b F factors
  // Take the Nphoton average of F over NJet sub-bins at each kin bin.
  vector<float> Fextrap, Fextrapjb, fb0;
  for (int j = 0; j < NbinsNJet; ++j) {
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	int binCC = CCbins.jbk(j, b, k);  if (binCC <= 0) continue;
	int binCCjb = CCbins.jb(j, b);  if (binCCjb <= 0) continue;
	if (j < extrapByMCthreshold || !useMCJfactors) {
	  float Npho = 0, NphoF = 0, fbb = 0;
	  for (auto J : jetSubBins.at(j)) {
	    int binCCJb = CCbins.Jb(J, b);  if (binCCJb <= 0) continue;
	    float Nemuxp = 0;
	    if (doClosure) {
	      Nemuxp = hCCJb_zll->GetBinContent(binCCJb);
	    } else {
	      Nemuxp = hCCJb_zmm->GetBinContent(binCCJb) * h_pur_m->GetBinContent(binCCjb)
		+ hCCJb_zee->GetBinContent(binCCJb) * h_pur_e->GetBinContent(binCCjb);
	    }
	    if (J >= (int) fb0.size()) {
	      fb0.push_back(Nemuxp);
	      // cout << "j, J, b, k, fb0Size = " << j << ", " << J << ", " << b << ", " << k << ", " << fb0.size() << endl;
	    } else {
	      fbb = Nemuxp / fb0.at(J);
	      int binCCJ0k = CCbins.Jbk(J, 0, k);  if (binCCJ0k <= 0) continue;
	      float NphoJ = hCCspl_photon->GetBinContent(binCCJ0k);
	      // cout << "j, J, b, k, fbb, NphoJ = " << j << ", " << J << ", " << b << ", " << k << ", " << fbb << ", " << NphoJ << endl;
	      Npho += NphoJ;
	      NphoF += NphoJ*fbb;
	    }
	  }
	  if (b == 0) {
	    Fextrap.push_back(1);
	  } else {
	    float f = 0;
	    if (NphoF != 0 && Npho != 0)
	    // if (b == 0)  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      f = NphoF/Npho;
	    else {
	      // Fall back to unsplit yields
	      int binCCj0 = CCbins.jb(j, 0);
	      if (doClosure) {
		f = hCCjb_zll->GetBinContent(binCCjb);
		f /= hCCjb_zll->GetBinContent(binCCj0);
	      } else {
		f = hCCjb_zmm->GetBinContent(binCCjb) * h_pur_m->GetBinContent(binCCjb)
		  + hCCjb_zee->GetBinContent(binCCjb) * h_pur_e->GetBinContent(binCCjb);
		f /= hCCjb_zmm->GetBinContent(binCCj0) * h_pur_m->GetBinContent(binCCj0)
		  + hCCjb_zee->GetBinContent(binCCj0) * h_pur_e->GetBinContent(binCCj0);
	      }
	    }
	    Fextrap.push_back(f);
	  }
	} else {
	  // Set F(binsj) = F(binsj-1) * J(binsj-1;binsj)
	  // cout << "j, b, k, Jextrap (rel. err) = " << j << ", " << b << ", " << k << ", " <<
	  //   Jextrap.at(binCC - 1).first << " (" << Jextrap.at(binCC - 1).second << ")" << endl;
	  Fextrap.push_back(Fextrap.at(CCbins.jbk(j-1, b, k) - 1) * Jextrap.at(binCC - 1).first);
	}
      }
    }
  }

  // Get the k-integrated numbers for AN
  float yj0 = 1;
  for (int j = 0; j < NbinsNJet; ++j) {
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      bool validkbin = false;
      for (int k = 0; k < kinSize; ++k) {
	if (CCbins.jbk(j, b, k) < 0) continue;
	validkbin = true;
	int binCCjb = CCbins.jb(j, b);  if (binCCjb <= 0) continue;
	if (j < extrapByMCthreshold || !useMCJfactors) {
	  float yjb = 0;
	  if (doClosure)
	    yjb = hCCjb_zll->GetBinContent(binCCjb);
	  else
	    yjb = hCCjb_zmm->GetBinContent(binCCjb) * h_pur_m->GetBinContent(binCCjb)
	      +   hCCjb_zee->GetBinContent(binCCjb) * h_pur_e->GetBinContent(binCCjb);
	  if (b == 0) {
	    yj0 = yjb;
	    Fextrapjb.push_back(1);
	  } else {
	    Fextrapjb.push_back(yjb/yj0);
	  }
	} else {
	  Fextrapjb.push_back(Fextrapjb.at(CCbins.jb(j-1, b) - 1) * Jextrap.at(CCbins.jbk(j, b, k) - 1).first);
	}
	if (validkbin) break;
      }
    }
  }

  // kin systematics derived from analysis of the closure plot
  // vector<float> systKin = {0, 0.07, 0.10, 0.20};
  vector< vector<float> > systKin = {{0, 0.15, 0.30, 0.30},
				     {0, 0.15, 0.15, 0.30},
				     {0, 0.15, 0.15, 0.30},
				     {0, 0.15, 0.15, 0.30},
				     {0, 0.15, 0.15, 0.30}};

  TFile* histoOutFile = nullptr;
  TH1D *ZinvBGpred = nullptr, *ZinvBGsysUp = nullptr, *ZinvBGsysLow = nullptr, *hMCexp = nullptr;
  if (doClosure) {
    histoOutFile = TFile::Open("hClosure.root", "RECREATE");
    ZinvBGpred = (TH1D*) hCC_zll->Clone();  ZinvBGpred->SetNameTitle("ZinvBGpred", "Predicted Zinv yield");
    ZinvBGsysUp = (TH1D*) hCC_zll->Clone();  ZinvBGsysUp->SetNameTitle("ZinvBGsysUp", "Predicted Zinv upper error");
    ZinvBGsysLow = (TH1D*) hCC_zll->Clone();  ZinvBGsysLow->SetNameTitle("ZinvBGsysLow", "Predicted Zinv lower error");
    hMCexp = (TH1D*) hCC_photon->Clone();  hMCexp->SetNameTitle("hMCexp", "Expected Zinv yield");
  }
  // Write the output dat file
  // string datFileName("DY_");
  // if (deltaPhi == "nominal") datFileName += "signal";
  // else datFileName += deltaPhi;
  // datFileName += ".dat";
  // if (!datFile.is_open()) datFile.open(datFileName);

  // File* datFile = fopen(datFileName.c_str(), "w");
  // fprintf(datFile, "%s\n",
  char linebuf[2048];
  sprintf(linebuf, "%s\n",
  	  " j b k| | Nmumu |  Nee  | Nb/0b | stat  |MC stat| ttz SF| syst+ | syst- | sysKin| sysPur"
  	  );
  datFile << linebuf;
  for (int j = 0; j < NbinsNJet; ++j) {
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      bool usedkbin = false;
      double chisq = 0, schi = 0, sdev = 0, sdevs = 0, ksum = 0;
      // cout << "(j, b) = (" << j << ", " << b << ")" << endl;
      for (int k = 0; k < kinSize; ++k) {
	int binCC = CCbins.jbk(j, b, k);  if (binCC <= 0) continue;
	int binCCjb = CCbins.jb(j, b);  if (binCCjb <= 0) continue;
	float JextrapErr, ttzErr;
	if (j < extrapByMCthreshold || !useMCJfactors || b == 0) {
	  JextrapErr = 0;
	  ttzErr = 0;
	} else {
	  JextrapErr = Jextrap.at(binCC - 1).second;
	  // Factors in the Jextrap double ratio are fully correlated
	  ttzErr = relErrXsec_ttz * (hCCjb_MCttzFrac->GetBinContent(binCCjb) -
				     hCCjb_MCttzFrac->GetBinContent(CCbins.jb(j-1, b)) -
				     hCCjb_MCttzFrac->GetBinContent(CCbins.jb(j, 0)) +
				     hCCjb_MCttzFrac->GetBinContent(CCbins.jb(j-1, 0))
				     );
	  // Deal also with negative MC weights
	  if (ttzErr < 0) ttzErr = relErrXsec_ttz;
	}
  	// fprintf(datFile, " %1d %1d %1d| |%7d|%7d|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f\n", j, b, k,
  	sprintf(linebuf, " %1d %1d %1d| |%7d|%7d|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f\n", j, b, k,
		(int) hCC_zmm->GetBinContent(binCC),
		(int) hCC_zee->GetBinContent(binCC),
		Fextrap.at(binCC - 1),
		DYstat.at(binCCjb - 1),
		JextrapErr,
		ttzErr,
		j >= extrapByMCthreshold && useMCJfactors ? systJ.at(b) : 0,
		j >= extrapByMCthreshold && useMCJfactors ? systJ.at(b) : 0,
		systKin[j][b],
		DYpurSys.at(binCCjb - 1)
		);
	datFile << linebuf;
	if (!usedkbin) {
	  // sprintf(linebuf, "%d & %d & %d & %6.3f & %3.0f & %3.0f & $\\pm%2.0f^{+%2.0f}_{-%2.0f}$ & %2.0f & %2.0f \\\\\n",
	  sprintf(linebuf, "%d & %d & %d & %6.3f & %3.0f & %3.0f & $\\pm%2.0f\\pm%2.0f$ & %2.0f & %2.0f \\\\\n",
		  binCCjb,
		  (int) hCCjb_zmm->GetBinContent(binCCjb),
		  (int) hCCjb_zee->GetBinContent(binCCjb),
		  Fextrapjb.at(binCCjb - 1),
		  100*DYstat.at(binCCjb - 1),
		  100*DYpurSys.at(binCCjb - 1),
		  100*JextrapErr,
		  // j >= extrapByMCthreshold && useMCJfactors ? 100*systJ.at(b) : 0,  // for asymmetric error
		  j >= extrapByMCthreshold && useMCJfactors ? 100*systJ.at(b) : 0,
		  100*ttzErr,
		  100*systKin[j][b]
		  );
	  latexFile << linebuf;
	  usedkbin = true;
	}
	if (doClosure) {
	  // For closure plot
	  int binCCb0 = CCbins.jbk(j, 0, k);
	  Double_t binValue = hCC_photon->GetBinContent(binCCb0) * Fextrap.at(binCC - 1);
	  ZinvBGpred->SetBinContent(binCC, binValue);
	  ZinvBGpred->SetBinError(binCC, 0.00001*binValue);
	  Double_t binError = binValue * Sqrt(Power(DYstat.at(binCCjb - 1), 2) +
					      Power(systKin[j][b], 2));
	  ZinvBGsysUp->SetBinContent(binCC, binError);
	  ZinvBGsysLow->SetBinContent(binCC, binError);
	  // For chisq-derived systematic
	  Double_t expErr = hMCexp->GetBinError(binCC);
	  if (binValue > 0 && expErr > 0) {
	    // Double_t binErrSq = Power(binValue*DYstat.at(binCCjb - 1), 2);
	    Double_t expValue = hMCexp->GetBinContent(binCC);
	    Double_t dev = expValue - binValue;
	    chisq += Power(dev/expErr, 2);
	    schi += dev / expErr;
	    sdev += dev/binValue;
	    sdevs += Power(dev/binValue, 2);
	    ksum++;
	    // cout << "k = " << k << ", pred = " << binValue << " +/- " << binValue*DYstat.at(binCCjb - 1)
	    //      << ", exp = " << expValue << " +/- " << expErr << ", dev/pred = " << dev/binValue
	    //      << ", chisq = " << Power(dev/expErr, 2) << endl;
	  }
	}
      }  // k
      if (doClosure && b > 0) {
	Double_t rmsm = 0;
	Double_t sys = 0;
	if (ksum > 0) {
	  chisq -= Power(schi, 2) / ksum;
	  rmsm = sdevs/ksum - Power(sdev/ksum, 2);
	  sys = chisq > ksum ? rmsm*(1 - ksum/chisq) : 0;
	  // cout << "ksum = " << ksum << ", sdevs = " << sdevs << ", chisq = " << chisq << ", sys^2 = " << sys << endl;
	  rmsm = Sqrt(rmsm);
	  sys = Sqrt(sys);
	}
	Double_t chisqNorm = ksum > 1 ? chisq/(ksum-1) : chisq;
	cout << "(j, b) = (" << j << ", " << b << "):  rms = " << rmsm
	     << ", chisq/DoF = " << chisqNorm << ", sysKin = " << sys << endl;
      }
    }
  }
  datFile.close();
  latexFile.close();
  // fclose(datFile);
  if (doClosure) {
    ZinvBGpred->Draw();
    hMCexp->Draw();
    histoOutFile->Write();
  }

  gApplication->Terminate(0);

}

TFile* openFile(const char* fileName) {
  TFile* tFile = new TFile(fileName);
  if (!tFile->IsOpen()) {
    cout << "Failed to open " << fileName << endl;
    return nullptr;
  } else
    return tFile;
}

TH1D* getHist(TFile* tFile, const char* histName) {
  TH1D* tH1 = (TH1D*) tFile->Get(histName);
  if (tH1 == nullptr)
    cout << "Failed to get " << histName << " from file " << tFile->GetName() << endl;
  return tH1;
}
