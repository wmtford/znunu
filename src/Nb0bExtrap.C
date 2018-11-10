/*
  Run from the command line with
  root -l -b RA2bZinvLoadClasses.C Nb0bExtrap.C
*/

#include "CCbinning.h"
#include <stdio.h>
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"

#include "TMath.h"
using TMath::Sqrt; using TMath::Power;

void Nb0bExtrap(const std::string& era = "2016", const std::string& deltaPhi = "nominal") {

  CCbinning CCmaps(era, deltaPhi);
  std::vector< std::vector<double> > kinThresholds = CCmaps.kinThresholds();
  std::vector<int> nJet1Thresholds = CCmaps.nJet1Thresholds();
  std::vector<int> nJetThresholds = CCmaps.nJetThresholds();
  std::vector<int> nbThresholds = CCmaps.nbThresholds();
  unsigned kinSize = CCmaps.kinSize();
  std::vector< std::vector<int> > jetSubBins = CCmaps.jetSubBins();
  CCbinning::ivector_map toCCbin = CCmaps.toCCbin();
  CCbinning::ivector_map toCCbinjb = CCmaps.toCCbinjb();
  CCbinning::ivector_map toCCbinSpl = CCmaps.toCCbinSpl();
  CCbinning::ivector_map toCCbinJb = CCmaps.toCCbinJb();
  int extrapByMCthreshold = 4;  // Use MC for the NJets 9+ bin
  // std::vector<int> extrapFromRange = {5, 6};

  TFile ZllData("../outputs/histsDYspl_2016v15.root");
  TFile photonData("../outputs/histsPhoton_2016v15.root");
  
  // Get the Z->ll data histograms
  TH1F* hCC_zmm = (TH1F*) ZllData.Get("hCC_zmm");
  TH1F* hCC_zee = (TH1F*) ZllData.Get("hCC_zee");
  TH1F* hCCjb_zmm = (TH1F*) ZllData.Get("hCCjb_zmm");
  TH1F* hCCjb_zee = (TH1F*) ZllData.Get("hCCjb_zee");
  TH1F* hCCjb_zll = (TH1F*) hCCjb_zmm->Clone();  hCCjb_zll->SetName("hCCjb_zll");
  hCCjb_zll->Add(hCCjb_zee);
  TH1F* hCCJb_zmm = (TH1F*) ZllData.Get("hCCJb_zmm");
  TH1F* hCCJb_zee = (TH1F*) ZllData.Get("hCCJb_zee");

  // Calculate Z->ll data stat errors
  std::vector<float> DYstat;
  for (int j = 0; j < nJetThresholds.size(); ++j) {
    float b0 = 0, b0err = 0, bb = 0, bberr = 0;
    for (int b = 0; b < nbThresholds.size(); ++b) {
      int jj = j < extrapByMCthreshold ? j : j-1;
      std::vector<int> jb = {jj, b};
      int binCCjb;
      try {binCCjb = toCCbinjb.at(jb);}
      catch (const std::out_of_range& oor) {continue;}
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
  // TFile effFile("../plots/histograms/effHists.root");
  TFile effFile("../python/effHists.root");
    
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
  for (int i = 0; i < pbins.size(); ++i) {
    float y = 0, ype = 0;
    for (int j = 0; j < pbins[i].size(); ++j) {
      int bin = pbins[i][j];
      y += hCCjb_zmm->GetBinContent(bin) + hCCjb_zee->GetBinContent(bin);
      ype += hCCjb_zmm->GetBinContent(bin)*h_pur_m->GetBinError(bin)/h_pur_m->GetBinContent(bin)
	+ hCCjb_zee->GetBinContent(bin)*h_pur_e->GetBinError(bin)/h_pur_e->GetBinContent(bin);
    }
    DYpurSys.push_back(ype/y);
  }

  // Calculate Nb/0b F factors
  // Take the Nphoton average of F over NJet sub-bins at each kin bin.
  TH1F* hCCspl_photon = (TH1F*) photonData.Get("hCCspl_photon");
  std::vector<float> Fextrap, fb0;
  for (int j = 0; j < nJetThresholds.size(); ++j) {
    float fbb = 0;
    for (int b = 0; b < nbThresholds.size(); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	try {int binCCjbk = toCCbin.at((std::vector<int>) {j, b, k});}  catch (const std::out_of_range& oor) {continue;}
	if (j < extrapByMCthreshold) {
	  float Npho = 0, NphoF = 0;
	  for (auto J : jetSubBins[j]) {
	    int binCCJb, binCCjb;
	    try {binCCJb = toCCbinJb.at((std::vector<int>) {J, b});}  catch (const std::out_of_range& oor) {continue;}
	    try {binCCjb = toCCbinjb.at((std::vector<int>) {j, b});}  catch (const std::out_of_range& oor) {continue;}
	    float Nemuxp = hCCJb_zmm->GetBinContent(binCCJb) * h_pur_m->GetBinContent(binCCjb)
	                 + hCCJb_zee->GetBinContent(binCCJb) * h_pur_e->GetBinContent(binCCjb);
	    // if (b == 0 && k == 0) {
	    if (J >= fb0.size()) {
	      fb0.push_back(Nemuxp);
	      // cout << "j, J, b, k, fb0Size = " << j << ", " << J << ", " << b << ", " << k << ", " << fb0.size() << endl;
	    } else {
	      fbb = Nemuxp / fb0.at(J);
	      int binCCJ0k;
	      try {binCCJ0k = toCCbinSpl.at((std::vector<int>) {J, 0, k});}  catch (const std::out_of_range& oor) {continue;}
	      float NphoJ = hCCspl_photon->GetBinContent(binCCJ0k);
	      // cout << "j, J, b, k, fbb, NphoJ = " << j << ", " << J << ", " << b << ", " << k << ", " << fbb << ", " << NphoJ << endl;
	      Npho += NphoJ;
	      NphoF += NphoJ*fbb;
	    }
	  }
	  if (b == 0) {
	    Fextrap.push_back(1);
	  } else {
	    Fextrap.push_back(NphoF/Npho);
	  }
	} else {
	  // Copy NJets 7-8 to 9+
	  Fextrap.push_back(Fextrap[toCCbin.at((std::vector<int>) {j-1, b, k}) - 1]);
	}
      }
    }
  }

    
  // kin systematics derived using analysis of the closure plot
  std::vector<float> systKin = {0, 0.07, 0.10, 0.20};




  // Create and write the output dat file
  FILE* outFile;
  std::string outFileName("DY_");
  if (deltaPhi == "nominal") outFileName += "signal";
  else outFileName += deltaPhi;
  outFileName += ".dat";
  outFile = fopen(outFileName.data(), "w");
  fprintf(outFile, "%s\n",
	  " j b k|| Nmumu |  Nee  | Nb/0b | stat  |MC stat| ttz SF| syst+ | syst- | sysKin| sysPur"
	  );
  for (int j = 0; j < nJetThresholds.size(); ++j) {
    for (int b = 0; b < nbThresholds.size(); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	std::vector<int> jbk = {j, b, k}, jb = {j, b};
	int binCC = 0, binCCjb = 0;
	try {binCC = toCCbin.at(jbk);}	catch (const std::out_of_range& oor) {continue;}
	try {binCCjb = toCCbinjb.at(jb);}  catch (const std::out_of_range& oor) {continue;}
	fprintf(outFile, " %1d %1d %1d||%7d|%7d|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f\n", j, b, k,
		(int) hCC_zmm->GetBinContent(binCC),
		(int) hCC_zee->GetBinContent(binCC),
		Fextrap[binCC - 1],
		DYstat[binCCjb - 1],
		0.999,
		0.999,
		0.999,
		0.999,
		systKin[b],
		DYpurSys[binCCjb - 1]
		);
      }
    }
  }
  fclose(outFile);

  // Under construction ...

}
