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
  CCbinning::ivector_map toCCbin = CCmaps.toCCbin();
  CCbinning::ivector_map toCCbinSpl = CCmaps.toCCbinSpl();
  CCbinning::ivector_map toCCbinjb = CCmaps.toCCbinjb();
  int extrapByMCthreshold = 4;  // Use MC for the NJets 9+ bin


  TFile ZllData("../outputs/histsDYspl_2016v15.root");
  TFile photonData("../outputs/histsPhoton_2016v15.root");
  
  TH1F* hCC_zmm = (TH1F*) ZllData.Get("hCC_zmm");
  TH1F* hCC_zee = (TH1F*) ZllData.Get("hCC_zee");
  
  TH1F* hCCjb_zmm = (TH1F*) ZllData.Get("hCCjb_zmm");
  TH1F* hCCjb_zee = (TH1F*) ZllData.Get("hCCjb_zee");
  TH1F* hCCjb_zll = (TH1F*) hCCjb_zmm->Clone();  hCCjb_zll->SetName("hCCjb_zll");
  hCCjb_zll->Add(hCCjb_zee);
  std::vector<float> DYstat;
  float b0 = 0, b0err = 0, bb = 0, bberr = 0;
  for (int j = 0; j < nJetThresholds.size(); ++j) {
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

  FILE* outFile;
  std::string outFileName("DY_");
  if (deltaPhi == "nominal") outFileName += "signal";
  else outFileName += deltaPhi;
  outFileName += ".dat";
  outFile = fopen(outFileName.data(), "w");
  fprintf(outFile, "%s\n",
	  " j b k|| Nmumu |  Nee  | Nb/0b | stat  |MC stat| ttz SF| syst+ | syst- | sysKin | sysPur"
	  );
  for (int j = 0; j < nJetThresholds.size(); ++j) {
    for (int b = 0; b < nbThresholds.size(); ++b) {
      for (int k = 0; k < kinSize; ++k) {
	std::vector<int> jbk = {j, b, k}, jb = {j, b};
	int binCC = 0, binCCjb = 0;
	try {binCC = toCCbin.at(jbk);}
	catch (const std::out_of_range& oor) {continue;}
	try {binCCjb = toCCbinjb.at(jb);}
	catch (const std::out_of_range& oor) {continue;}
	fprintf(outFile, " %1d %1d %1d||%7d|%7d|%7.4f|%7.4f\n", j, b, k,
		(int) hCC_zmm->GetBinContent(binCC),
		(int) hCC_zee->GetBinContent(binCC),
		0.999,
		DYstat[binCCjb - 1]
		);
      }
    }
  }
  fclose(outFile);

  // Under construction ...

}
