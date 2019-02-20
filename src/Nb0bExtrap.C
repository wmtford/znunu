/*
  Run from the command line with
  root -l -b Nb0bLoadClasses.C Nb0bExtrap.C
*/

#include "CCbinning.h"
#include <stdio.h>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TApplication.h"
#include "TMath.h"
using TMath::Sqrt; using TMath::Power; using TMath::Max;

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;

#include "TMath.h"
using TMath::Sqrt; using TMath::Power;

void Nb0bExtrap(const std::string& era = "Run2", const std::string& deltaPhi = "nominal") {
  ofstream datFile("DY_signal.dat");
  // ofstream datFile;
  ofstream latexFile("DY_signal.tex");

  CCbinning CCbins(era, deltaPhi);
  int kinSize = CCbins.kinSize();
  int NbinsNJet = CCbins.binsj();
  std::vector< std::vector<int> > jetSubBins = CCbins.jetSubBins();

  bool useMCJfactors = false;
  bool doJfromData = false;
  int extrapByMCthreshold = NbinsNJet - 1;  // Use MC for the last NJets bin
  // std::vector<int> extrapFromRange = {5, 6};
  float relErrXsec_ttz = 0.3;

  // TFile ZllData("../outputs/histsDY_2016v16.root");  if (!ZllData.IsOpen()) return;
  // TFile photonData("../outputs/histsPhoton_2016v16.root");  if (!photonData.IsOpen()) return;
  // TFile ZllXMC("../outputs/histsDYMC_2016v16_noPU.root");  if (!ZllXMC.IsOpen()) return;
  // TFile ZllData("../outputs/histsDY_2017v16.root");  if (!ZllData.IsOpen()) return;
  // TFile photonData("../outputs/histsPhoton_2017v16.root");  if (!photonData.IsOpen()) return;
  // TFile ZllXMC("../outputs/histsDYMC_2017v16_ZptWt.root");  if (!ZllXMC.IsOpen()) return;
  // TFile ZllData("../outputs/histsDY_2018ABv16.root");  if (!ZllData.IsOpen()) return;
  // TFile photonData("../outputs/histsPhoton_2018ABv16.root");  if (!photonData.IsOpen()) return;
  // TFile ZllXMC("../outputs/histsDYMC_2018v16_ZptWt.root");  if (!ZllXMC.IsOpen()) return;
  TFile ZllData("../outputs/histsDY_2018CDv16.root");  if (!ZllData.IsOpen()) return;
  TFile photonData("../outputs/histsPhoton_2018CDv16.root");  if (!photonData.IsOpen()) return;
  TFile ZllXMC("../outputs/histsDYMC_2018HEMv16_ZptWt.root");  if (!ZllXMC.IsOpen()) return;
  // TFile ZllData("../outputs/histsDY_Run2v16.root");  if (!ZllData.IsOpen()) return;
  // TFile photonData("../outputs/histsPhoton_Run2v16.root");  if (!photonData.IsOpen()) return;
  // TFile ZllXMC("../outputs/histsDYMC_Run2v16_ZptWt.root");  if (!ZllXMC.IsOpen()) return;
  
  // Get the Z->ll and photon data and MC histograms
  TH1F* hCC_zmm = (TH1F*) ZllData.Get("hCC_zmm");
  TH1F* hCC_zee = (TH1F*) ZllData.Get("hCC_zee");
  TH1F* hCCjb_zmm = (TH1F*) ZllData.Get("hCCjb_zmm");
  TH1F* hCCjb_zee = (TH1F*) ZllData.Get("hCCjb_zee");
  TH1F* hCCjb_zll = (TH1F*) ZllData.Get("hCCjb_zll");
  TH1F* hCCspl_photon = (TH1F*) photonData.Get("hCCspl_photon");
  TH1F* hCCJb_zmm = (TH1F*) ZllData.Get("hCCJb_zmm");
  TH1F* hCCJb_zee = (TH1F*) ZllData.Get("hCCJb_zee");
  TH1F* hCCJb_zll = (TH1F*) ZllData.Get("hCCJb_zll");
  TH1F* hCCJb_dymm = (TH1F*) ZllXMC.Get("hCCJb_dymm");
  TH1F* hCCJb_dyee = (TH1F*) ZllXMC.Get("hCCJb_dyee");
  TH1F* hCCJb_ttzmm = (TH1F*) ZllXMC.Get("hCCJb_ttzmm");
  TH1F* hCCJb_ttzee = (TH1F*) ZllXMC.Get("hCCJb_ttzee");
  TH1F* hCCJb_VVmm = (TH1F*) ZllXMC.Get("hCCJb_VVmm");
  TH1F* hCCJb_VVee = (TH1F*) ZllXMC.Get("hCCJb_VVee");
  TH1F* hCCJb_ttmm = (TH1F*) ZllXMC.Get("hCCJb_ttmm");
  TH1F* hCCJb_ttee = (TH1F*) ZllXMC.Get("hCCJb_ttee");
  TH1F* hCCJb_MCall = (TH1F*) hCCJb_dymm->Clone();  hCCJb_MCall->SetName("hCCJb_MCall");
  hCCJb_MCall->Add(hCCJb_dyee);
  hCCJb_MCall->Add(hCCJb_ttzmm);
  hCCJb_MCall->Add(hCCJb_ttzee);
  hCCJb_MCall->Add(hCCJb_VVmm);
  hCCJb_MCall->Add(hCCJb_VVee);
  // hCCJb_MCall->Add(hCCJb_ttmm);  // Omit non-peaking, since data are purity corrected.
  // hCCJb_MCall->Add(hCCJb_ttee);

  // For calculation of the ttz cross section uncertainty from the ttz fraction.
  TH1F* hCCjb_dymm = (TH1F*) ZllXMC.Get("hCCjb_dymm");
  TH1F* hCCjb_dyee = (TH1F*) ZllXMC.Get("hCCjb_dyee");
  TH1F* hCCjb_ttzmm = (TH1F*) ZllXMC.Get("hCCjb_ttzmm");
  TH1F* hCCjb_ttzee = (TH1F*) ZllXMC.Get("hCCjb_ttzee");
  TH1F* hCCjb_VVmm = (TH1F*) ZllXMC.Get("hCCjb_VVmm");
  TH1F* hCCjb_VVee = (TH1F*) ZllXMC.Get("hCCjb_VVee");
  TH1F* hCCjb_ttmm = (TH1F*) ZllXMC.Get("hCCjb_ttmm");
  TH1F* hCCjb_ttee = (TH1F*) ZllXMC.Get("hCCjb_ttee");
  TH1F* hCCjb_MCall = (TH1F*) hCCjb_dymm->Clone();  hCCjb_MCall->SetName("hCCjb_MCall");
  hCCjb_MCall->Add(hCCjb_dyee);
  hCCjb_MCall->Add(hCCjb_ttzmm);
  hCCjb_MCall->Add(hCCjb_ttzee);
  hCCjb_MCall->Add(hCCjb_VVmm);
  hCCjb_MCall->Add(hCCjb_VVee);
  hCCjb_MCall->Add(hCCjb_ttmm);
  hCCjb_MCall->Add(hCCjb_ttee);
  hCCjb_ttzmm->Print("all");
  hCCjb_ttzee->Print("all");
  TH1F* hCCjb_MCttzFrac = (TH1F*) hCCjb_ttzmm->Clone();
  hCCjb_MCttzFrac->SetNameTitle("hCCjb_MCttzFrac", "hCCjb_MCttzFrac");
  hCCjb_MCttzFrac->Add(hCCjb_ttzee);
  hCCjb_MCttzFrac->Divide(hCCjb_MCall);
  hCCjb_MCttzFrac->Print("all");

  // Calculate Z->ll data stat errors
  std::vector<float> DYstat;
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
	  DYstat.push_back(Sqrt(Power((b0err/b0), 2) + Power((bberr/bb), 2)));
	else
	  DYstat.push_back(1);
      }
    }
  }

  // Compute the purity uncertainties
  TFile effFile("../plots/histograms/effHists.root");  if (!effFile.IsOpen()) return;
  // TFile effFile("../python/effHists.root");  if (!effFile.IsOpen()) return;
    
  TH1F* h_pur_m = (TH1F*) effFile.Get("h_pur_m");
  TH1F* h_pur_e = (TH1F*) effFile.Get("h_pur_e");
  std::vector<float> DYpurSys;
  // Declare the bin merging map for purities
  std::vector< std::vector<int> > pbins;
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

  std::vector<float> systJ = {0, 0, 0, 0};
  std::vector<std::pair<float, float>> Jextrap;

  // Calculate the MC J factors to extend the Nb/0b F factors
  // Take the Nphoton average of F over NJet sub-bins at each kin bin.
  std::vector<std::pair<float, float>> jb0;
  for (int j = 0; j <= extrapByMCthreshold; ++j) {
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	int binCC = CCbins.jbk(j, b, k);  if (binCC <= 0) continue;
	float Npho = 0;
	std::pair<float, float> jbb, NphoF;
	for (auto J : jetSubBins.at(j)) {
	  int binCCJb = CCbins.Jb(J, b);  if (binCCJb <= 0) continue;
	  int binCCjb = CCbins.jb(j, b);  if (binCCjb <= 0) continue;
	  std::pair<float, float> nMC;
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
	  Jextrap.push_back(std::pair<float, float>(1, 0));
	} else {
	  if (Npho > 0)
	    Jextrap.push_back(std::pair<float, float>(NphoF.first/Npho, NphoF.second/Npho));
	  else
	    Jextrap.push_back(std::pair<float, float>(jbb.first, jbb.second));
	}
	if (j >= extrapByMCthreshold) {
	  // Compute ratio of highest to next-highest NJets bins
	  float relErrjbk = Jextrap.at(binCC - 1).second;
	  float relErrjm1bk = Jextrap.at(CCbins.jbk(j-1, b, k) - 1).second;
	  Jextrap.at(binCC - 1).first /= Jextrap.at(CCbins.jbk(j-1, b, k) - 1).first;
	  Jextrap.at(binCC - 1).second = Sqrt(Power(relErrjbk, 2) + Power(relErrjm1bk, 2));
	  cout << "j, b, k, Jextrap = " << j << ", " << b << ", " << k << ", " <<
	    Jextrap.at(binCC - 1).first << ", rel err = " << Jextrap.at(binCC - 1).second << endl;
	}
      }
    }
  }

  // Jextrap systematics derived from analysis of the above for data, MC
  systJ = {0, 0.10, 0.10, 0.20};

  // Calculate Nb/0b F factors
  // Take the Nphoton average of F over NJet sub-bins at each kin bin.
  std::vector<float> Fextrap, Fextrapjb, fb0;
  for (int j = 0; j < NbinsNJet; ++j) {
    for (int b = 0; b < CCbins.binsb(j); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	int binCC = CCbins.jbk(j, b, k);  if (binCC <= 0) continue;
	int binCCjb = CCbins.jb(j, b);  if (binCCjb <= 0) continue;
	if (j < extrapByMCthreshold || !useMCJfactors) {
	  float Npho = 0, NphoF = 0, fbb = 0;
	  for (auto J : jetSubBins.at(j)) {
	    int binCCJb = CCbins.Jb(J, b);  if (binCCJb <= 0) continue;
	    float Nemuxp = hCCJb_zmm->GetBinContent(binCCJb) * h_pur_m->GetBinContent(binCCjb)
	                 + hCCJb_zee->GetBinContent(binCCJb) * h_pur_e->GetBinContent(binCCjb);
	    if (J >= (int) fb0.size()) {
	      fb0.push_back(Nemuxp);
	      // cout << "j, J, b, k, fb0Size = " << j << ", " << J << ", " << b << ", " << k << ", " << fb0.size() << endl;
	    } else {
	      fbb = Nemuxp / fb0.at(J);
	      int binCCJ0k = CCbins.Jbk(J, 0, k);  if (binCCJ0k <= 0) continue;
	      float NphoJ = hCCspl_photon->GetBinContent(binCCJ0k);
	      cout << "j, J, b, k, fbb, NphoJ = " << j << ", " << J << ", " << b << ", " << k << ", " << fbb << ", " << NphoJ << endl;
	      Npho += NphoJ;
	      NphoF += NphoJ*fbb;
	    }
	  }
	  if (b == 0) {
	    Fextrap.push_back(1);
	  } else {
	    float f = 0;
	    if (NphoF != 0 && Npho != 0)
	      f = NphoF/Npho;
	    else {
	      // Fall back to unsplit yields
	      int binCCj0 = CCbins.jb(j, 0);
	      f = hCCjb_zmm->GetBinContent(binCCjb) * h_pur_m->GetBinContent(binCCjb)
		+ hCCjb_zee->GetBinContent(binCCjb) * h_pur_e->GetBinContent(binCCjb);
	      f /= hCCjb_zmm->GetBinContent(binCCj0) * h_pur_m->GetBinContent(binCCj0)
	      	+ hCCjb_zee->GetBinContent(binCCj0) * h_pur_e->GetBinContent(binCCj0);
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
	  float yjb = hCCjb_zmm->GetBinContent(binCCjb) * h_pur_m->GetBinContent(binCCjb)
	    +         hCCjb_zee->GetBinContent(binCCjb) * h_pur_e->GetBinContent(binCCjb);
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
  std::vector<float> systKin = {0, 0.07, 0.10, 0.20};

  // Write the output dat file
  // std::string datFileName("DY_");
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
		systKin.at(b),
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
		  100*systKin.at(b)
		  );
	  latexFile << linebuf;
	  usedkbin = true;
	}
      }
    }
  }
  datFile.close();
  latexFile.close();
  // fclose(datFile);

  gApplication->Terminate(0);

}
