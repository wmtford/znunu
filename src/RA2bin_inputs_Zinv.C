//
//  Zinv inputs to RA2/b SUSY search fit
//

#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TMath.h"
using TMath::Sqrt; using TMath::Power;

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <cstring>

namespace Sample {
  enum sampleChoice {Signal, HDP, LDP};
}
using namespace Sample;

typedef std::vector<Float_t> V1F;
typedef std::vector <std::vector<Float_t> > V2F;

// void RA2bin_inputs_Zinv(Int_t doSample = Signal,
void RA2bin_inputs_Zinv(sampleChoice doSample = Signal,
			std::vector< std::pair<TString, float> > gJetsFnRoot = {std::make_pair("gJets_2016", 35.9),
			    std::make_pair("gJets_2017", 41.5), std::make_pair("gJets_2018", 59.4)},
			const TString DRfnRoot = TString("DR"),
			const TString DYfnRoot = TString("DY"),
			const TString MCfileName = TString("ZinvMCttzMC174bin.root"),
			const Float_t MClumiScale = 1.0,
			const Bool_t expandDYkin = false) {
  //
  // Supply input flat files with names
  //  <gJetsFnRoot>_signal.dat for doSample = Signal
  //  <gJetsFnRoot>_htp.dat for doSample = HDP
  //  <gJetsFnRoot>_ltp.dat for doSample = LDP
  // If MC root file has lumi different from that of the data, supply 
  //  MClumiScale = dataLumi / MClumi.
  //

  gROOT->Reset();

  // use the 'plain' style for plots (white backgrounds, etc)
  cout << "...using style 'Plain'\n";

  gROOT->SetStyle("Plain");

  // use bold lines and markers
  gStyle->SetMarkerStyle(8);  // non-scalable dot
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // For the fit/function:
  gStyle->SetOptFit(1);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(4);
  gStyle->SetFuncStyle(1);
  // gStyle->SetFuncWidth(2);

  //..Get rid of X error bars
  gStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //  ===============================================================
  //
  // Prototypes for helper functions
  //
  Int_t getData_gJets(std::vector<TString> fileNames, std::vector<float> lumi,
		      V2F& Ngobs,
		      V2F& NgobsEB,
		      V2F& NgobsEE,
		      V2F& ZgR,
		      V2F& ZgRerr,
		      V2F& gEtrg,
		      V2F& gEtrgErr,
		      V2F& gEtrgSys,
		      V2F& gSF,
		      V2F& gSFerr,
		      V2F& gFdir,
		      V2F& gFdirErrUP,
		      V2F& gFdirErrLow,
		      V2F& gFdirSys,
		      V2F& gPur,
		      V2F& gPurErr,
		      V2F& ZgDR,
		      V2F& ZgDRerrUp,
		      V2F& ZgDRerrLow);
  Int_t getData_DR(const char* fileName,
		   V2F& ZgDR,
		   V2F& ZgDRerrUp,
		   V2F& ZgDRerrLow,
		   Float_t& DRscaleErr,
		   Float_t& DY0bPurErr,
		   Float_t& DYtrigEffErr,
		   Float_t& LeptonSFerr,
		   Float_t& btag_SFerr);
  Int_t getData_DR(const char* fileName,
		   V2F& ZgDR,
		   V2F& ZgDRerrUp,
		   V2F& ZgDRerrLow);
  Int_t getData_DY(const char* fileName,
		   std::vector<V2F>& Rb0,
		   std::vector<V2F>& Rb0stat,
		   std::vector<V2F>& Rb0MCstat,
		   std::vector<V2F>& Rb0sysUp,
		   std::vector<V2F>& Rb0sysLow,
		   std::vector<V2F>& Rb0sysKin,
		   std::vector<V2F>& Rb0sysPur);
  void setCorrelationLabels(TH1* histo, const unsigned correlVar, const int minCorrelBin=9999,
			    const int maxCorrelBin=0);

  Float_t ignoreBin = -999.;
  std::vector<TString> dataFile_gJets;
  V1F vlumi;
  TString dataFile_DR;
  TString dataFile_DY;
  TString output_rootFile;
  TString output_plotFile;
  Float_t histoXlabelSize;
  Float_t canvasBottomMargin;
  Int_t MaxNjets, MaxNb, MaxKin, MaxKinDY;
  if (doSample == LDP || doSample == HDP) {
    MaxNjets = 5;
    MaxNb = 4;
    MaxKin = 13;
    histoXlabelSize = 0.015;
    canvasBottomMargin = 0.15;
    if (doSample == LDP) {
      for (auto fnroot : gJetsFnRoot) {
	dataFile_gJets.push_back(fnroot.first+"_ldp.dat");
	vlumi.push_back(fnroot.second);
      }
      dataFile_DR = DRfnRoot+TString("_ldp.dat");
      dataFile_DY = DYfnRoot+TString( "_ldp.dat");
      output_rootFile = "ZinvHistos_ldp.root";
      output_plotFile = "ZinvBGpred_ldp.pdf";
    } else if (doSample == HDP) {
      for (auto fnroot : gJetsFnRoot) {
	dataFile_gJets.push_back(fnroot.first+"_hdp.dat");
	vlumi.push_back(fnroot.second);
      }
      dataFile_DR = DRfnRoot+TString("_hdp.dat");
      dataFile_DY = DYfnRoot+TString("_hdp.dat");
      output_rootFile = "ZinvHistos_hdp.root";
      output_plotFile = "ZinvBGpred_hdp.pdf";
    }
  } else {
    MaxNjets = 5;
    MaxNb = 4;
    MaxKin = 10;
    histoXlabelSize = 0.04;
    canvasBottomMargin = 0.35;
    for (auto fnroot : gJetsFnRoot) {
      dataFile_gJets.push_back(fnroot.first+"_signal.dat");
      vlumi.push_back(fnroot.second);
    }
    dataFile_DR = DRfnRoot+TString("_signal.dat");
    dataFile_DY = DYfnRoot+TString("_signal.dat");
    output_rootFile = "ZinvHistos.root";
    output_plotFile = "ZinvBGpred.pdf";
  }
  MaxKinDY = expandDYkin ? MaxKin : 1;
  Int_t MaxBins;
  if (doSample == Signal) {
    MaxBins = (MaxNjets*MaxNb-1)*MaxKin-2*2*MaxNb;  // Exclude (Njets0,Nb3) && (Njets3,4;HT0,3)
  } else {
    MaxBins = (MaxNjets*MaxNb-1)*MaxKin-2*3*MaxNb;  // Exclude (Njets0,Nb3) && (Njets3,4;HT0,3,6)
  }
  Int_t bin;  // Bin index for loops

  // N(b=0) from gamma+jets
  V2F Ngobs(MaxNjets, V1F(MaxKin, 0));
  V2F NgobsEB(MaxNjets, V1F(MaxKin, 0));
  V2F NgobsEE(MaxNjets, V1F(MaxKin, 0));
  V2F ZgR(MaxNjets, V1F(MaxKin, 0));
  V2F ZgRerr(MaxNjets, V1F(MaxKin, 0));
  V2F gEtrg(MaxNjets, V1F(MaxKin, 0));
  V2F gEtrgErr(MaxNjets, V1F(MaxKin, 0));
  V2F gEtrgErrEff(MaxNjets, V1F(MaxKin, 0));
  V2F gEtrgSys(MaxNjets, V1F(MaxKin, 0));
  V2F gSF(MaxNjets, V1F(MaxKin, 0));
  V2F gSFerr(MaxNjets, V1F(MaxKin, 0));
  V2F gFdir(MaxNjets, V1F(MaxKin, 0));
  V2F gFdirErrUp(MaxNjets, V1F(MaxKin, 0));
  V2F gFdirErrUpEff(MaxNjets, V1F(MaxKin, 0));
  V2F gFdirErrLow(MaxNjets, V1F(MaxKin, 0));
  V2F gFdirErrLowEff(MaxNjets, V1F(MaxKin, 0));
  V2F gFdirSys(MaxNjets, V1F(MaxKin, 0));
  V2F gPur(MaxNjets, V1F(MaxKin, 0));
  V2F gPurErr(MaxNjets, V1F(MaxKin, 0));
  V2F gPurErrEff(MaxNjets, V1F(MaxKin, 0));
  V2F ZgDR(MaxNjets, V1F(MaxKin, 0));
  V2F ZgDRerrUp(MaxNjets, V1F(MaxKin, 0));
  V2F ZgDRerrLow(MaxNjets, V1F(MaxKin, 0));
  V2F ZgDR_from_gJets(MaxNjets, V1F(MaxKin, 0));
  V2F ZgDRerrUp_from_gJets(MaxNjets, V1F(MaxKin, 0));
  V2F ZgDRerrLow_from_gJets(MaxNjets, V1F(MaxKin, 0));
  Float_t DRscaleErr=0, DY0bPurErr=0, DYtrigEffErr=0, LeptonSFerr=0, btagSFerr=0;

  // Ratio Nb/N0 from DY
  std::vector < V2F >
    DYvalues(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  std::vector < V2F >
    DYstat(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  std::vector < V2F >
    DYMCstat(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  std::vector < V2F >
    DYsysNjUp(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  std::vector < V2F >
    DYsysNjLow(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  std::vector < V2F >
    DYsysKin(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  std::vector < V2F >
    DYsysPur(MaxNjets, std::vector< V1F >(MaxNb, V1F(MaxKinDY, 0)));
  // Read gJets, DR, DY data from files
  if (0 != getData_gJets(dataFile_gJets, vlumi,
			 Ngobs, NgobsEB, NgobsEE, ZgR, ZgRerr, gEtrg, gEtrgErr, gEtrgSys, gSF,
			 gSFerr, gFdir, gFdirErrUp, gFdirErrLow, gFdirSys, gPur, gPurErr, ZgDR_from_gJets,
			 ZgDRerrUp_from_gJets, ZgDRerrLow_from_gJets)) {
    cout << "Failed to get gJets data." << endl;
    return;
  }
  if (0 != getData_DR(dataFile_DR, ZgDR, ZgDRerrUp, ZgDRerrLow, DRscaleErr, DY0bPurErr, DYtrigEffErr, LeptonSFerr, btagSFerr)) {
    cout << "Failed to get DR data." << endl;
    return;
  }
  if (0 != getData_DY(dataFile_DY, DYvalues, DYstat, DYMCstat, DYsysNjUp, DYsysNjLow, DYsysKin, DYsysPur)) {
    cout << "Failed to get DY data." << endl;
    return;
  }

  Float_t NobsSum = 0, gEtrgAv = 0, gEtrgErrAv = 0, gFdirAv = 0, gFdirErrUpAv = 0, gFdirErrLowAv = 0, gPurAv = 0, gPurErrAv = 0;
  cout << "DRscaleErr = " << DRscaleErr << ", DY0bPurErr = " << DY0bPurErr << ", DYtrigEffErr = " << DYtrigEffErr << ", LeptonSFerr = " <<LeptonSFerr  << ", btagSFerr = " << btagSFerr << endl;
  bin = 0;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
      if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
      if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
      bin++;
      cout << bin << "  Ngobs = " << Ngobs[ijet][ikin]
	   << "  ZgR = " << ZgR[ijet][ikin] << " (" << 100*ZgRerr[ijet][ikin]
           << "%)  gEtrg =  " << gEtrg[ijet][ikin] << " (" << 100*gEtrgErr[ijet][ikin] << "+/-" << 100*gEtrgSys[ijet][ikin]
           << "%)  gSF =  " << gSF[ijet][ikin] << " (" << 100*gSFerr[ijet][ikin]
           << "%)  gFdir = " << gFdir[ijet][ikin]
	   << " (+" << 100*gFdirErrUp[ijet][ikin] << "-" << +100*gFdirErrLow[ijet][ikin] << "+/-" << +100*gFdirSys[ijet][ikin]
           << "%)  gPur = " << gPur[ijet][ikin] << " (" << 100*gPurErr[ijet][ikin]
           << "%)  ZgDR = " << ZgDR[ijet][ikin]
	   << " (+" << 100*ZgDRerrUp[ijet][ikin] << "-" << 100*ZgDRerrLow[ijet][ikin] << "%)"
	   << endl;

      //  Treat gFdirSys as fraction of 1-Fdir, uncorrelated, and combine with stat error
      Float_t FdirSysUp = (1/gFdir[ijet][ikin] - 1) * gFdirSys[ijet][ikin];
      gFdirErrUp[ijet][ikin] = Sqrt(Power(gFdirErrUp[ijet][ikin], 2) + Power(FdirSysUp, 2));
      Float_t FdirSysLow = (1/gFdir[ijet][ikin] - 1) * gFdirSys[ijet][ikin];
      gFdirErrLow[ijet][ikin] = Sqrt(Power(gFdirErrLow[ijet][ikin], 2) + Power(FdirSysLow, 2));

      NobsSum += Ngobs[ijet][ikin];
      gEtrgAv += Ngobs[ijet][ikin]*gEtrg[ijet][ikin];
      gEtrgErrAv += Ngobs[ijet][ikin]*gEtrgErr[ijet][ikin]*gEtrg[ijet][ikin];  // Absolute error on gEtrg
      gFdirAv += Ngobs[ijet][ikin]*gFdir[ijet][ikin];
      gFdirErrUpAv += Ngobs[ijet][ikin]*gFdirErrUp[ijet][ikin]*gFdir[ijet][ikin];  // Absolute +error on gFdir
      gFdirErrLowAv += Ngobs[ijet][ikin]*gFdirErrLow[ijet][ikin]*gFdir[ijet][ikin];  // Absolute -error on gFdir
      gPurAv += Ngobs[ijet][ikin]*gPur[ijet][ikin];
      gPurErrAv += Ngobs[ijet][ikin]*gPurErr[ijet][ikin]*gPur[ijet][ikin];  // Absolute error on gPur
    }
  }
  cout << endl;
  gEtrgAv /= NobsSum;  // Average ZgR
  gEtrgErrAv /= NobsSum;  // Average gEtrg absolute uncertainty
  gFdirAv /= NobsSum;  // Average gFdir
  gFdirErrUpAv /= NobsSum;  // Average gFdir absolute uncertainty
  gFdirErrLowAv /= NobsSum;  // Average gFdir absolute uncertainty
  gPurAv /= NobsSum;  // Average gPur
  gPurErrAv /= NobsSum;  // Average gPur absolute uncertainty
  bin = 0;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
      for (Int_t ikin=0; ikin<MaxKinDY; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
	bin++;
	cout << bin << "  DYvalue = " << DYvalues[ijet][ib][ikin]
	     << "  DYstat = " << DYstat[ijet][ib][ikin]
	     << "  DYMCstat = " << DYMCstat[ijet][ib][ikin]
             << "  DYsysNjUp = " << DYsysNjUp[ijet][ib][ikin]
	     << "  DYsysNjLow = " << DYsysNjLow[ijet][ib][ikin]
             << "  DYsysKin = " << DYsysKin[ijet][ib][ikin]
	     << "  DYsysPur = " << DYsysPur[ijet][ib][ikin] << endl;
      }
    }
  }
  cout << endl << endl;

  TFile *MCfile = 0;
  if (doSample == Signal) {
    // Input file for MC histogram
    cout << "Initializing input file " << MCfileName << endl;
    MCfile = TFile::Open(MCfileName, "READ");
    if (!MCfile) {
      cout << "Didn't find MC file " << MCfileName << endl;
    }
  }

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:"); // Without these 2 lines
  rootdir->cd();                     // the histo disappears when the macro exits.
  
  TH1F* hTemplate = new TH1F("hTemplate", "template for Zinv histograms", MaxBins, 0.5, MaxBins+0.5);
  // Tune up the template histogram
  // hTemplate->Reset();
  // hTemplate->GetXaxis()->SetLabelSize(histoXlabelSize);
  // hTemplate->GetYaxis()->SetLabelSize(0.05);
  // hTemplate->GetYaxis()->SetTitleSize(0.05);
  hTemplate->SetOption("HIST");
  hTemplate->SetMarkerSize(0);

  // hCorrelTemplate is the starting point for x-axis labels that will encode the
  // correlation patterns for the integration fitter.
  TH1F* hCorrelTemplate = (TH1F*)hTemplate->Clone("hCorrelTemplate");
  if (doSample == Signal) {  // Leave x-axis labels as numbers for hdp, ldp samples
    // Set the x-axis labels for default (no) correlations
    bin = 0;
    for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
      for (Int_t ib=0; ib<MaxNb; ++ib) {
	if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
	for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
	  if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	  if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
	  bin++;
	  char label[25];
	  Int_t iMHT = ikin/3;
	  if(ikin>=8) {
	    iMHT=3;
	  }
	  Int_t iHT = ikin;
	  sprintf(label, "%5s%1d%6s%1d%4s%1d%3s%1d", "_NJets", ijet, "_BTags", ib, "_MHT", iMHT, "_HT", iHT);
	  // cout << label << endl;
	  hCorrelTemplate->GetXaxis()->SetBinLabel(bin, label);
	}
      }
    }
  }

  // Output file
  TFile *ZinvHistos = TFile::Open(output_rootFile, "RECREATE");

  gStyle->SetOptStat("n");

  // Book histograms for N(b=0) yield factors from gamma+jets
  //  NZvv_pred = Ngobs * ZgR * gPur * ZgDR
  TH1F *hzvvgJNobs = (TH1F*)hCorrelTemplate->Clone("hzvvgJNobs");
  if (doSample == Signal) setCorrelationLabels(hzvvgJNobs, 4);  
  hzvvgJNobs->SetTitle("gamma+jets yield");
  hzvvgJNobs->GetYaxis()->SetTitle("gamma+jets yield");

  TH1F *hgJstat = (TH1F*)hTemplate->Clone("hgJstat");
  hgJstat->SetTitle("gamma+jets yield stat error");
  hgJstat->GetYaxis()->SetTitle("gamma+jets yield stat error");

  TH1F *hzvvTF = (TH1F*)hCorrelTemplate->Clone("hzvvTF");
  hzvvTF->SetTitle("Zinv transfer factor");
  hzvvTF->GetYaxis()->SetTitle("Zinv transfer factor");

  TH1F *hgJZgR = (TH1F*)hTemplate->Clone("hgJZgR");
  hgJZgR->SetTitle("gamma+jets Z/gamma ratio");
  hgJZgR->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio");

  TH1F *hgJZgR0b = new TH1F("hgJZgR0b", "gamma+jets Z/gamma ratio, Nb = 0", 46, 0.5, 46.5);
  hgJZgR0b->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio");

  // TH1F *hzvvgJZgRerr = (TH1F*)hCorrelTemplate->Clone("hzvvgJZgRerr");
  // if (doSample == Signal) setCorrelationLabels(hzvvgJZgRerr, 4);  
  // hzvvgJZgRerr->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");

  TH1F *hgJZgRerr = (TH1F*)hTemplate->Clone("hgJZgRerr");
  hgJZgRerr->SetTitle("gamma+jets Z/gamma ratio error");
  hgJZgRerr->GetYaxis()->SetTitle("gamma+jets Z/gamma ratio error");

  TH1F *hgJEtrg = (TH1F*)hTemplate->Clone("hgJEtrg");
  hgJEtrg->SetTitle("gamma+jets trigger efficiency");
  hgJEtrg->GetYaxis()->SetTitle("gamma+jets trigger efficiency");

  TH1F *hzvvgJEtrgErr = (TH1F*)hCorrelTemplate->Clone("hzvvgJEtrgErr");
  if (doSample == Signal) setCorrelationLabels(hzvvgJEtrgErr, 13);  
  hzvvgJEtrgErr->SetTitle("gamma+jets trigger efficiency error");
  hzvvgJEtrgErr->GetYaxis()->SetTitle("gamma+jets trigger efficiency error");

  TH1F *hgJSF = (TH1F*)hTemplate->Clone("hgJSF");
  hgJSF->SetTitle("gamma+jets data/MC scale factor");
  hgJSF->GetYaxis()->SetTitle("gamma+jets data/MC scale factor");

  TH1F *hgJSFerr = (TH1F*)hTemplate->Clone("hgJSFerr");
  // if (doSample == Signal) setCorrelationLabels(hgJSFerr, 15);  // Canceled by DR 
  hgJSFerr->SetTitle("gamma+jets data/MC scale factor");
  hgJSFerr->GetYaxis()->SetTitle("gamma+jets data/MC scale factor");

  TH1F *hgJFdir = (TH1F*)hTemplate->Clone("hgJFdir");
  hgJFdir->SetTitle("gamma+jets fragmentation factor");
  hgJFdir->GetYaxis()->SetTitle("gamma+jets fragmentation factor");

  TH1F *hgJFdirErrUp = (TH1F*)hTemplate->Clone("hgJFdirErrUp");
  hgJFdirErrUp->SetTitle("gamma+jets fragmentation factor +error");
  hgJFdirErrUp->GetYaxis()->SetTitle("gamma+jets fragmentation factor +error");

  TH1F *hgJFdirErrLow = (TH1F*)hTemplate->Clone("hgJFdirErrLow");
  hgJFdirErrLow->SetTitle("gamma+jets fragmentation factor -error");
  hgJFdirErrLow->GetYaxis()->SetTitle("gamma+jets fragmentation factor -error");

  TH1F *hgJPur = (TH1F*)hTemplate->Clone("hgJPur");
  hgJPur->SetTitle("gamma+jets purity");
  hgJPur->GetYaxis()->SetTitle("gamma+jets purity");

  TH1F *hzvvgJPurErr = (TH1F*)hCorrelTemplate->Clone("hzvvgJPurErr");
  if (doSample == Signal) setCorrelationLabels(hzvvgJPurErr, 15);  
  hzvvgJPurErr->SetTitle("gamma+jets purity error");
  hzvvgJPurErr->GetYaxis()->SetTitle("gamma+jets purity error");

  TH1F *hZgDR = (TH1F*)hTemplate->Clone("hZgDR");
  hZgDR->SetTitle("gamma+jets double ratio");
  hZgDR->GetYaxis()->SetTitle("gamma+jets double ratio");

  TH1F *hZgDRerrUp = (TH1F*)hTemplate->Clone("hZgDRerrUp");
  hZgDRerrUp->SetTitle("gamma+jets double ratio error");
  hZgDRerrUp->GetYaxis()->SetTitle("gamma+jets double ratio error");

  TH1F *hZgDRerrLow = (TH1F*)hTemplate->Clone("hZgDRerrLow");
  hZgDRerrLow->SetTitle("gamma+jets double ratio error");
  hZgDRerrLow->GetYaxis()->SetTitle("gamma+jets double ratio error");

  TH1F *hzvvNbCorrelUp = (TH1F*)hCorrelTemplate->Clone("hzvvNbCorrelUp");
  if (doSample == Signal) setCorrelationLabels(hzvvNbCorrelUp, 4);
  hzvvNbCorrelUp->SetTitle("Errors correlated in Nb");
  hzvvNbCorrelUp->GetYaxis()->SetTitle("Errors correlated in Nb");

  TH1F *hzvvNbCorrelLow = (TH1F*)hCorrelTemplate->Clone("hzvvNbCorrelLow");
  if (doSample == Signal) setCorrelationLabels(hzvvNbCorrelLow, 4);
  hzvvNbCorrelLow->SetTitle("Errors correlated in Nb");
  hzvvNbCorrelLow->GetYaxis()->SetTitle("Errors correlated in Nb");

  TH1F *hzvvScaleErr = (TH1F*)hCorrelTemplate->Clone("hzvvScaleErr");
  if (doSample == Signal) setCorrelationLabels(hzvvScaleErr, 15);  
  hzvvScaleErr->SetTitle("Zinv global scale error");
  hzvvScaleErr->GetYaxis()->SetTitle("Zinv global scale error");

  // Book histograms for Nb/N0 ratio values and error components
  //  Nb/N0 pred = DYvalues
  TH1F *hDYvalue = (TH1F*)hTemplate->Clone("hDYvalue");
  hDYvalue->SetTitle("DY ratio to 0b value");
  hDYvalue->GetYaxis()->SetTitle("DY ratio to 0b value");

  TH1F *hzvvDYstat = (TH1F*)hCorrelTemplate->Clone("hzvvDYstat");
  // if (doSample == Signal) setCorrelationLabels(hzvvDYstat, 11, 3999);
  if (doSample == Signal) setCorrelationLabels(hzvvDYstat, 15, 3999, 100);
  hzvvDYstat->SetTitle("DY ratio to 0b stat error");
  hzvvDYstat->GetYaxis()->SetTitle("DY ratio to 0b stat error");

  TH1F *hDYMCstat = (TH1F*)hTemplate->Clone("hDYMCstat");
  hDYMCstat->SetTitle("DY ratio to 0b stat error");
  hDYMCstat->GetYaxis()->SetTitle("DY ratio to 0b stat error");

  TH1F *hDYsysNjUp = (TH1F*)hTemplate->Clone("hDYsysNjUp");
  hDYsysNjUp->SetTitle("DY ratio to 0b syst+ error Nj extrapolation");
  hDYsysNjUp->GetYaxis()->SetTitle("DY ratio to 0b syst+ error Nj extrapolation");

  TH1F *hDYsysNjLow = (TH1F*)hTemplate->Clone("hDYsysNjLow");
  hDYsysNjLow->SetTitle("DY ratio to 0b syst- error Nj extrapolation");
  hDYsysNjLow->GetYaxis()->SetTitle("DY ratio to 0b syst- error Nj extrapolation");

  TH1F *hzvvDYMCerrUp = (TH1F*)hCorrelTemplate->Clone("hzvvDYMCerrUp");
  if (doSample == Signal) setCorrelationLabels(hzvvDYMCerrUp, 15, 9999, 4100);  // Uncorrel in Njets, but only last bin != 0
  hzvvDYMCerrUp->SetTitle("DY ratio to 0b syst+ error Nj extrapolation");
  hzvvDYMCerrUp->GetYaxis()->SetTitle("DY ratio to 0b syst+ error Nj extrapolation");

  TH1F *hzvvDYMCerrLow = (TH1F*)hCorrelTemplate->Clone("hzvvDYMCerrLow");
  if (doSample == Signal) setCorrelationLabels(hzvvDYMCerrLow, 15, 9999, 4100);  
  hzvvDYMCerrLow->SetTitle("DY ratio to 0b syst- error Nj extrapolation");
  hzvvDYMCerrLow->GetYaxis()->SetTitle("DY ratio to 0b syst- error Nj extrapolation");

  TH1F *hzvvDYsysKin = (TH1F*)hCorrelTemplate->Clone("hzvvDYsysKin");
  if (doSample == Signal) setCorrelationLabels(hzvvDYsysKin, 4, 9999, 100);  
  hzvvDYsysKin->SetTitle("DY ratio to 0b syst error kinematics dependence");
  hzvvDYsysKin->GetYaxis()->SetTitle("DY ratio to 0b syst error kinematics dependence");

  TH1F *hzvvDYsysPur = (TH1F*)hCorrelTemplate->Clone("hzvvDYsysPur");
  if (doSample == Signal) setCorrelationLabels(hzvvDYsysPur, 15, 2299);  
  hzvvDYsysPur->SetTitle("DY ratio to 0b syst error purity");
  hzvvDYsysPur->GetYaxis()->SetTitle("DY ratio to 0b syst error purity");

  // For combined results
  Float_t statErr[MaxBins], sysUp[MaxBins], sysLow[MaxBins];
  TH1F *ZinvBGpred = (TH1F*)hTemplate->Clone("ZinvBGpred");
  ZinvBGpred->SetTitle("Predicted Z#rightarrow#nu#bar{#nu} background");
  ZinvBGpred->GetYaxis()->SetTitle("Predicted Z#rightarrow#nu#bar{#nu} background");
  ZinvBGpred->SetLineColor(kBlue);
  ZinvBGpred->SetMarkerColor(kBlue);
  ZinvBGpred->SetMarkerStyle(8);
  ZinvBGpred->SetMarkerSize(0.5);
  // ZinvBGpred->SetMinimum(-1);
  TH1F *ZinvBGsysUp = (TH1F*)hTemplate->Clone("ZinvBGsysUp");
  ZinvBGsysUp->GetYaxis()->SetTitle("Upper syst. error of Z#rightarrow#nu#bar{#nu} background");
  ZinvBGsysUp->SetTitle("Upper syst. error of Z#rightarrow#nu#bar{#nu} background");
  TH1F *ZinvBGsysLow = (TH1F*)hTemplate->Clone("ZinvBGsysLow");
  ZinvBGsysLow->SetTitle("Lower syst. error of Z#rightarrow#nu#bar{#nu} background");
  ZinvBGsysLow->GetYaxis()->SetTitle("Lower syst. error of Z#rightarrow#nu#bar{#nu} background");

  // For combined results, for bins with zero control-sample events
  TH1F *ZinvBG0EVpred = (TH1F*)hTemplate->Clone("ZinvBG0EVpred");
  ZinvBG0EVpred->SetTitle("Predicted Z#rightarrow#nu#bar{#nu} background, 0 events");
  ZinvBG0EVpred->GetYaxis()->SetTitle("Predicted Z#rightarrow#nu#bar{#nu} background, 0 events");
  ZinvBG0EVpred->SetLineColor(kBlue);
  ZinvBG0EVpred->SetMarkerColor(kBlue);
  ZinvBG0EVpred->SetMarkerStyle(8);
  ZinvBG0EVpred->SetMarkerSize(0.5);
  TH1F *ZinvBG0EVsysUp = (TH1F*)hTemplate->Clone("ZinvBG0EVsysUp");
  ZinvBG0EVsysUp->SetTitle("Upper syst. error of Z#rightarrow#nu#bar{#nu} background, 0 events");
  ZinvBG0EVsysUp->GetYaxis()->SetTitle("Upper syst. error of Z#rightarrow#nu#bar{#nu} background, 0 events");
  TH1F *ZinvBG0EVsysLow = (TH1F*)hTemplate->Clone("ZinvBG0EVsysLow");
  ZinvBG0EVsysLow->SetTitle("Lower syst. error of Z#rightarrow#nu#bar{#nu} background, 0 events");
  ZinvBG0EVsysLow->GetYaxis()->SetTitle("Lower syst. error of Z#rightarrow#nu#bar{#nu} background, 0 events");

  ofstream tableFile;  tableFile.open("table_for_AN.txt");
  ofstream predFile;  predFile.open("zinvPred.txt");
  char buf[256];

  // For AN table
  const char* tbLabelNj[] = {
    (char*)("\\njets 2"), (char*)("\\njets 3-4"), (char*)("\\njets 5-6"),
    (char*)("\\njets 7-8"), (char*)("\\njets $\\ge9$")
  };
  const char* tbLabelKin[] = {
    (char*)("MHT\\_HT\\_1"), (char*)("MHT\\_HT\\_2"), (char*)("MHT\\_HT\\_3"), (char*)("MHT\\_HT\\_4"), (char*)("MHT\\_HT\\_5"),
    (char*)("MHT\\_HT\\_6"), (char*)("MHT\\_HT\\_7"), (char*)("MHT\\_HT\\_8"), (char*)("MHT\\_HT\\_9"), (char*)("MHT\\_HT\\_10")
  };
  //  
  // Fill histograms
  //
  Int_t binzb = 0;
  bin = 0;
  for (Int_t ijet=0; ijet<MaxNjets; ++ijet) {
    for (Int_t ib=0; ib<MaxNb; ++ib) {
      if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
      for (Int_t ikin=0; ikin<MaxKin; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (MaxKin > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins, QCD-binned
	bin++;
	// cout << "ijet, ib, ikin, bin = " << ijet << ", " << ib << ", " << ikin << ", " << bin << endl;
	Int_t ikinDY = MaxKinDY == MaxKin ? ikin : 0;

        //
        // Partial cancellation of photon trigger eff., fragmentation, and purity syst errors in the prediction
        //
	gEtrgErrEff[ijet][ikin] = fabs( gEtrgErr[ijet][ikin] - gEtrgErrAv / gEtrgAv );  // Frac. error on gEtrg/<gEtrg>
	gFdirErrUpEff[ijet][ikin] = fabs( gFdirErrUp[ijet][ikin] - gFdirErrUpAv / gFdirAv );  // Frac. +error on gFdir/<gFdir>
	gFdirErrLowEff[ijet][ikin] = fabs( gFdirErrLow[ijet][ikin] - gFdirErrLowAv / gFdirAv );  // Frac. -error on gFdir/<gFdir>
	gPurErrEff[ijet][ikin] = fabs( gPurErr[ijet][ikin] - gPurErrAv / gPurAv );  // Frac. error on gPur/<gPur>
	if (ib == 0) {
	  binzb++;
	  hgJZgR0b->SetBinContent(binzb, ZgR[ijet][ikin]);
	  hgJZgR0b->SetBinError(binzb, ZgRerr[ijet][ikin]);
          cout << "iJet " << ijet << " Var (err nominal, eff): "
               << " gEtrg (" << gEtrgErr[ijet][ikin] << ", " << gEtrgErrEff[ijet][ikin] << ") "
	       << " gFdir (" << gFdirErrUp[ijet][ikin] << ", " << gFdirErrUpEff[ijet][ikin] << ") "
	       << " gFdir (" << gFdirErrLow[ijet][ikin] << ", " << gFdirErrLowEff[ijet][ikin] << ") "
               << " gPur (" << gPurErr[ijet][ikin] << ", " << gPurErrEff[ijet][ikin] << ") "
	       << endl;
	}

	hzvvgJNobs->SetBinContent(bin, Ngobs[ijet][ikin]);
	Ngobs[ijet][ikin] > 0 ? hgJstat->SetBinContent(bin, 1+1/Sqrt(Ngobs[ijet][ikin])) : hgJstat->SetBinContent(bin, 1+0);
	hgJZgR->SetBinContent(bin, ZgR[ijet][ikin]);
	hZgDR->SetBinContent(bin, ZgDR[ijet][ikin]);
	hgJZgRerr->SetBinContent(bin, 1+ZgRerr[ijet][ikin]);  // uncorrelated, combined into ZgDRErrUp,Low
	hgJEtrg->SetBinContent(bin, gEtrg[ijet][ikin]);
	hzvvgJEtrgErr->SetBinContent(bin, 1+gEtrgErrEff[ijet][ikin]);
	hgJSF->SetBinContent(bin, gSF[ijet][ikin]);
	hgJSFerr->SetBinContent(bin, 1+gSFerr[ijet][ikin]);  // correlated (but Eff error = 0 since error is bin-independent)
	hgJFdir->SetBinContent(bin, gFdir[ijet][ikin]);
	hgJFdirErrUp->SetBinContent(bin, 1+gFdirErrUpEff[ijet][ikin]);  // correlated in Nb
	hgJFdirErrLow->SetBinContent(bin, 1-gFdirErrLowEff[ijet][ikin]);  // correlated in Nb
	hgJPur->SetBinContent(bin, gPur[ijet][ikin]);
	hzvvgJPurErr->SetBinContent(bin, 1+gPurErrEff[ijet][ikin]);  // correlated
	hZgDRerrUp->SetBinContent(bin, 1+ZgDRerrUp[ijet][ikin]);
	hZgDRerrLow->SetBinContent(bin, 1-ZgDRerrLow[ijet][ikin]);
	// Avoid degenerate factors by combining ZgR stat, Fdir, and ZgDR errors into hzvvNbCorrelUp,Low
	hzvvNbCorrelUp->SetBinContent(bin, 1+Sqrt(Power(ZgRerr[ijet][ikin], 2)
						  + Power(gFdirErrUpEff[ijet][ikin], 2)
						  + Power(ZgDRerrUp[ijet][ikin], 2)));  // correlated in Nb
	hzvvNbCorrelLow->SetBinContent(bin, 1-Sqrt(Power(ZgRerr[ijet][ikin], 2)
						   + Power(gFdirErrLowEff[ijet][ikin], 2)
						   + Power(ZgDRerrLow[ijet][ikin], 2)));
	// Likewise combine several globally correlated errors
	hzvvScaleErr->SetBinContent(bin, 1+Sqrt(Power(DRscaleErr, 2) + Power(btagSFerr, 2)
						+ Power(DY0bPurErr, 2) + Power(DYtrigEffErr, 2)
						+ Power(LeptonSFerr, 2)));  // Correlated globally
  	hDYvalue->SetBinContent(bin, DYvalues[ijet][ib][ikinDY]);
  	hzvvDYstat->SetBinContent(bin, 1+DYstat[ijet][ib][ikinDY]);
  	hDYMCstat->SetBinContent(bin, 1+DYMCstat[ijet][ib][ikinDY]);
  	hDYsysNjUp->SetBinContent(bin, 1+DYsysNjUp[ijet][ib][ikinDY]);
  	hDYsysNjLow->SetBinContent(bin, 1-DYsysNjLow[ijet][ib][ikinDY]);
	// and combine extrapolation J-factor errors
  	hzvvDYMCerrUp->SetBinContent(bin, 1+Sqrt(Power(DYMCstat[ijet][ib][ikinDY], 2)
					       + Power(DYsysNjUp[ijet][ib][ikinDY], 2)));
  	hzvvDYMCerrLow->SetBinContent(bin, 1-Sqrt(Power(DYMCstat[ijet][ib][ikinDY], 2)
						+ Power(DYsysNjLow[ijet][ib][ikinDY], 2)));
  	hzvvDYsysKin->SetBinContent(bin, 1+DYsysKin[ijet][ib][ikinDY]);
  	hzvvDYsysPur->SetBinContent(bin, 1+DYsysPur[ijet][ib][ikinDY]);

  	// cout << bin
	//      << " " << hDYvalue->GetXaxis()->GetBinLabel(bin)
	//      << "  " << hDYvalue->GetBinContent(bin) << endl;

	Float_t transferFactor =  ZgDR[ijet][ikin] * ZgR[ijet][ikin] / gEtrg[ijet][ikin] / gSF[ijet][ikin] * gFdir[ijet][ikin] * gPur[ijet][ikin] * DYvalues[ijet][ib][ikinDY];
	// cout << Ngobs[ijet][ikin] << " " << ZgR[ijet][ikin] << " " << gFdir[ijet][ikin] << " "
	// << gPur[ijet][ikin] << " " << DYvalues[ijet][ib][ikin] << " " << ZgDR[ijet][ikin] << endl;
	Float_t thisNgobs = Ngobs[ijet][ikin] > 0 ? Ngobs[ijet][ikin] : 1.0;
	Float_t ZinvValue = thisNgobs * transferFactor;
	hzvvTF->SetBinContent(bin, transferFactor);
	Float_t wtStat = 1/thisNgobs;
	// statErr[bin-1] = ZinvValue*Sqrt(wtStat + Power(DYstat[ijet][ib][ikinDY], 2));
	statErr[bin-1] = ZinvValue*Sqrt(wtStat);
	sysUp[bin-1] = ZinvValue*Sqrt(Power(ZgRerr[ijet][ikin], 2)
				      + Power(gEtrgErrEff[ijet][ikin], 2)
				      + 0*Power(gSFerr[ijet][ikin], 2)  // This one cancels
				      + Power(gFdirErrUpEff[ijet][ikin], 2)
				      + Power(gPurErrEff[ijet][ikin], 2)
				      + Power(ZgDRerrUp[ijet][ikin], 2)
				      + Power(btagSFerr, 2)
				      + Power(DRscaleErr, 2)
				      + Power(DY0bPurErr, 2)
				      + Power(DYtrigEffErr, 2)
				      + Power(LeptonSFerr, 2)
				      + Power(DYstat[ijet][ib][ikinDY], 2)
				      + Power(DYMCstat[ijet][ib][ikinDY], 2)
				      + Power(DYsysNjUp[ijet][ib][ikinDY], 2)
				      + Power(DYsysKin[ijet][ib][ikinDY], 2)
				      + Power(DYsysPur[ijet][ib][ikinDY], 2));
	sysLow[bin-1] = ZinvValue*Sqrt(Power(ZgRerr[ijet][ikin], 2)
				       + Power(gEtrgErrEff[ijet][ikin], 2)
                  		       + 0*Power(gSFerr[ijet][ikin], 2)  // This one cancels
				       + Power(gFdirErrLowEff[ijet][ikin], 2)
				       + Power(gPurErrEff[ijet][ikin], 2)
				       + Power(ZgDRerrLow[ijet][ikin], 2)
				       + Power(btagSFerr, 2)
				       + Power(DRscaleErr, 2)
				       + Power(DY0bPurErr, 2)
				       + Power(DYtrigEffErr, 2)
				       + Power(LeptonSFerr, 2)
				       + Power(DYstat[ijet][ib][ikinDY], 2)
				       + Power(DYMCstat[ijet][ib][ikinDY], 2)
				       + Power(DYsysNjLow[ijet][ib][ikinDY], 2)
				       + Power(DYsysKin[ijet][ib][ikinDY], 2)
				       + Power(DYsysPur[ijet][ib][ikinDY], 2));
	if (bin == 1) predFile << endl;
	predFile << bin << "  " << ZinvValue << " +/- " << statErr[bin-1] << " + " << sysUp[bin-1]
		 << " - " << sysLow[bin-1]
		 << "  TF = " << hzvvTF->GetBinContent(bin)
		 << ", Ngobs = " <<   hzvvgJNobs->GetBinContent(bin) << endl;
	if (Ngobs[ijet][ikin] > 0) {
	  ZinvBGpred->SetBinContent(bin, ZinvValue);
	  ZinvBGsysUp->SetBinContent(bin, sysUp[bin-1]);
	  ZinvBGsysLow->SetBinContent(bin, sysLow[bin-1]);
	} else {
	  ZinvBG0EVpred->SetBinContent(bin, ZinvValue);
	  ZinvBG0EVsysUp->SetBinContent(bin, sysUp[bin-1]);
	  ZinvBG0EVsysLow->SetBinContent(bin, sysLow[bin-1]);
	}

	if (MaxKin <= 10 && ib == 0) {
          //  Write out LaTex for AN table
	  Float_t ZgRcorr = ZgR[ijet][ikin] / gEtrg[ijet][ikin] / gSF[ijet][ikin];
	  Float_t ZgRsys = ZgRcorr * Sqrt(Power(gEtrgErr[ijet][ikin], 2) + Power(gSFerr[ijet][ikin], 2));
	  if (ikin == 0 || (ijet > 2 && ikin == 1)) tableFile << "\\hline" << endl;
	  sprintf(buf, "%d & %5.0f & $%5.3f\\pm%5.3f$ & $%5.3f\\pm%5.3f\\pm%5.3f$ & $%5.3f^{+%5.3f}_{-%5.3f}$ & $%5.3f\\pm%5.3f^{+%5.3f}_{-%5.3f}$ & $%6.1f\\pm%4.1f^{+%4.1f}_{-%4.1f}$ \\\\\n",
		  binzb,
		  Ngobs[ijet][ikin],
		  gPur[ijet][ikin], gPur[ijet][ikin]*gPurErr[ijet][ikin],
		  ZgRcorr, ZgRcorr*ZgRerr[ijet][ikin], ZgRsys, 
		  gFdir[ijet][ikin], gFdir[ijet][ikin]*gFdirErrUp[ijet][ikin], gFdir[ijet][ikin]*gFdirErrLow[ijet][ikin],
		  ZgDR[ijet][ikin], ZgDR[ijet][ikin]*DRscaleErr, ZgDR[ijet][ikin]*ZgDRerrUp[ijet][ikin],
		  ZgDR[ijet][ikin]*ZgDRerrLow[ijet][ikin],
		  ZinvValue, statErr[bin-1], sysUp[bin-1], sysLow[bin-1]);
	  // sprintf(buf, "%d %s, %s & %5.0f & %5.0f & $%5.3f\\pm%5.3f\\pm%5.3f$ & $%5.3f\\pm%5.3f^{+%5.3f}_{-%5.3f}$ & $%6.1f\\pm%4.1f^{+%4.1f}_{-%4.1f}$ \\\\\n",
	  // 	  binzb,
	  // 	  tbLabelNj[ijet], tbLabelKin[ikin],
	  // 	  NgobsEB[ijet][ikin], NgobsEE[ijet][ikin],
	  // 	  ZgRcorr, ZgRcorr*ZgRerr[ijet][ikin], ZgRsys, 
	  // 	  ZgDR[ijet][ikin], ZgDR[ijet][ikin]*DRscaleErr, ZgDR[ijet][ikin]*ZgDRerrUp[ijet][ikin],
	  // 	  ZgDR[ijet][ikin]*ZgDRerrLow[ijet][ikin],
	  // 	  ZinvValue, statErr[bin-1], sysUp[bin-1], sysLow[bin-1]);
	  tableFile << buf;
	}
      }  // ikin
    }  // ib
  }  // ijet
  predFile.close();
  tableFile.close();

  TGraphAsymmErrors* ZinvBGsyst = new TGraphAsymmErrors(ZinvBGpred);
  ZinvBGsyst->SetLineColorAlpha(kRed, 0);
  ZinvBGsyst->SetLineWidth(-1);
  ZinvBGsyst->SetFillColor(kRed-10);
  ZinvBGsyst->SetFillStyle(1001);
  TGraphAsymmErrors* ZinvBG0EVsyst = new TGraphAsymmErrors(ZinvBG0EVpred);
  ZinvBG0EVsyst->SetLineColor(-1);
  // ZinvBG0EVsyst->SetLineStyle(1);
  // ZinvBG0EVsyst->SetLineWidth(2);
  // gStyle->SetEndErrorSize(3);
  for (Int_t bin0=0; bin0<MaxBins; ++bin0) {
    if (ZinvBGpred->GetBinContent(bin0+1) > 0) {
      ZinvBGpred->SetBinError(bin0+1, statErr[bin0]);
      ZinvBGsyst->SetPointError(bin0, 0.5, 0.5, sysLow[bin0], sysUp[bin0]);
    } else {
      ZinvBG0EVpred->SetBinError(bin0+1, statErr[bin0]);
      ZinvBG0EVsyst->SetPointError(bin0, 0.5, 0.5, sysLow[bin0], sysUp[bin0]);
    }
  }

  TH1F *hZinvMCbin = 0;
  if (doSample == Signal) {
    // Zinv MC 
    TH1F* hZinvMC = 0;
    TString MChistoname("plot_zinv_nj5_nb4_kin10_1");
    // TString MChistoname("plot_zinv_nj4_nb4_kin10_1");
    if (MCfile != 0) {
      if (!MCfile->Get(MChistoname)) {
	cout << "No MC histogram found in this file and directory" << endl;
      } else {
	hZinvMC = (TH1F*) MCfile->Get(MChistoname);
	hZinvMCbin = (TH1F*)hTemplate->Clone("hZinvMCbin");
	hZinvMCbin->SetTitle("MC Z#rightarrow#nu#bar{#nu}");
	hZinvMCbin->GetYaxis()->SetTitle("MC Z#rightarrow#nu#bar{#nu}");
	for (Int_t bin=1; bin<=MaxBins; ++bin) {
	  hZinvMCbin->SetBinContent(bin, hZinvMC->GetBinContent(bin));
	  hZinvMCbin->SetBinError(bin, hZinvMC->GetBinError(bin));
	}
	hZinvMCbin->SetLineColor(kGreen+2);
	hZinvMCbin->Scale(MClumiScale);
      }
    }
  }

  TCanvas *Cvalue = new TCanvas("Cvalue", "DY ratio to zero b", 1000, 500);
  // Cvalue->SetTopMargin(.03);
  // Cvalue->SetBottomMargin(canvasBottomMargin);
  // Cvalue->SetRightMargin(.03);
  gStyle->SetOptStat("");

  ZinvBGpred->Draw("e0");
  ZinvBGsyst->Draw("2 same");
  if (hZinvMCbin != 0) hZinvMCbin->Draw("hist e0 same");
  ZinvBGpred->Draw("e0 same");
  Cvalue->SetLogy();
  Cvalue->Update();

  TLegend* Zinvlegend;
  if (doSample == Signal)
    Zinvlegend = new TLegend(.72, .78, .92, .9, "");
  else
    Zinvlegend = new TLegend(.72, .80, .92, .9, "");
  Zinvlegend->AddEntry(ZinvBGpred, "Yields, stat. errors", "pe");
  Zinvlegend->AddEntry(ZinvBGsyst, "Syst. errors", "f");
  if (hZinvMCbin != 0) 
    Zinvlegend->AddEntry(hZinvMCbin, "Z#rightarrow#nu#bar{#nu} simulation", "le");
  Zinvlegend->Draw();
  Cvalue->Update();
  Cvalue->SaveAs(output_plotFile);

  ZinvHistos->Write();


  TCanvas *Cerrors = new TCanvas("Cerrors", "Error components", 1000, 500);
  // Cerrors->SetTopMargin(.03);
  // Cerrors->SetBottomMargin(canvasBottomMargin);
  // Cerrors->SetRightMargin(.03);
  hgJstat->SetMinimum(-0.2);  hgJstat->SetMaximum(2.2);
  hgJstat->SetLineColor(1);  hgJstat->Draw();
  hgJZgRerr->SetLineColor(2);  hgJZgRerr->Draw("same");
  hzvvgJEtrgErr->SetLineColor(36);  hzvvgJEtrgErr->Draw("same");
  hgJSFerr->SetLineColor(4);  hgJSFerr->Draw("same");
  hgJFdirErrUp->SetLineColor(47);  hgJFdirErrUp->Draw("same");
  hgJFdirErrLow->SetLineColor(47);  hgJFdirErrLow->Draw("same");
  hzvvgJPurErr->SetLineColor(7);  hzvvgJPurErr->Draw("same");
  hzvvScaleErr->SetLineColor(6);  hzvvScaleErr->Draw("same");
  hZgDRerrUp->SetLineColor(kOrange);  hZgDRerrUp->Draw("same");
  hZgDRerrLow->SetLineColor(kOrange);  hZgDRerrLow->Draw("same");
  hzvvDYstat->SetLineColor(8);  hzvvDYstat->Draw("same");
  hDYMCstat->SetLineColor(9);  hDYMCstat->Draw("same");
  hDYsysNjUp->SetLineColor(3);  hDYsysNjUp->Draw("same");
  hDYsysNjLow->SetLineColor(3);  hDYsysNjLow->Draw("same");
  hzvvDYsysKin->SetLineColor(42);  hzvvDYsysKin->Draw("same");
  hzvvDYsysPur->SetLineColor(46);  hzvvDYsysPur->Draw("same");
  TLegend* ErrorsLegend;
  if (doSample == Signal)
    ErrorsLegend = new TLegend(.12, .60, .30, .97, "");
  else
    ErrorsLegend = new TLegend(.12, .60, .30, .97, "");
  ErrorsLegend->AddEntry(hgJstat, "gJstat");
  ErrorsLegend->AddEntry(hgJZgRerr, "RZg");
  ErrorsLegend->AddEntry(hzvvgJEtrgErr, "gEtrg");
  ErrorsLegend->AddEntry(hgJSFerr, "gSF");
  ErrorsLegend->AddEntry(hgJFdirErrUp, "gFdir");
  ErrorsLegend->AddEntry(hzvvgJPurErr, "gJpur");
  ErrorsLegend->AddEntry(hzvvScaleErr, "Scale");
  ErrorsLegend->AddEntry(hZgDRerrUp, "DR");
  ErrorsLegend->AddEntry(hzvvDYstat, "DYstat");
  ErrorsLegend->AddEntry(hDYMCstat, "DYMCstat");
  ErrorsLegend->AddEntry(hDYsysNjUp, "DYNj");
  ErrorsLegend->AddEntry(hzvvDYsysKin, "DYkin");
  ErrorsLegend->AddEntry(hzvvDYsysPur, "DYpur");
  ErrorsLegend->Draw();
  Cerrors->Update();
  if (doSample == Signal)  Cerrors->SaveAs("errorComponents.pdf");


}  // ---------------------------------------------------------------

Int_t getData_gJets(std::vector<TString> fileNames, std::vector<float> lumi,
		    V2F& Ngobs,
		    V2F& NgobsEB,
		    V2F& NgobsEE,
		    V2F& ZgR,
		    V2F& ZgRerr,
		    V2F& gEtrg,
		    V2F& gEtrgErr,
		    V2F& gEtrgSys,
		    V2F& gSF,
		    V2F& gSFerr,
		    V2F& gFdir,
		    V2F& gFdirErrUp,
		    V2F& gFdirErrLow,
		    V2F& gFdirSys,
		    V2F& gPur,
		    V2F& gPurErr,
		    V2F& ZgDR,
		    V2F& ZgDRerrUp,
		    V2F& ZgDRerrLow) {
  //
  // Read data from the file, parse, and set data arrays.
  //
  /*
  Notes on correlations:
  gEtrgErr is correlated among Nj, Nb, HT bins.          1101 = 13
  gSFerr is correlated among all bins.                   1111 = 15
  gPurErr is correlated among all bins.                  1111 = 15
  ZgRerr is uncorrelated among the 46 Nb = 0 bins.       0100 = 4, combined with
  gFdir is uncorrelated among the 46 Nb = 0 bins.        0100 = 4, combined with
  ZgDRerrUp, Low are uncorrelated among the 46 bins.     0100 = 4  this one
  As of 10 Jul 2016 we combine these last two in ZgDRerrUp, Low.
   */

  void comb(Int_t mode, V1F lumi, V1F v, V1F e, Float_t& val, Float_t& err);

  Int_t Nrow = Ngobs.size();
  Int_t Ncol = Ngobs[0].size();

  // Vectors to hold values from separate input files
  int Nfiles = fileNames.size();
  std::vector<V2F> Ngobst(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> NgobsEBt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> NgobsEEt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> ZgRt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> ZgRerrt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gEtrgt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gEtrgErrt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gEtrgSyst(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gSFt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gSFerrt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gFdirt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gFdirErrUpt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gFdirErrLowt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gFdirSyst(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gPurt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> gPurErrt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> ZgDRt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> ZgDRerrUpt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));
  std::vector<V2F> ZgDRerrLowt(Nrow, V2F(Ncol, V1F(Nfiles, 0)));

  for (int ifn = 0; ifn < Nfiles; ++ifn) {
    cout << "gJets data filename = " << fileNames.at(ifn) << ", lumi = " << lumi.at(ifn) << endl;
    ifstream dataStream;
    dataStream.open((fileNames.at(ifn)).Data()); // open the data file
    if (!dataStream.good()) {
      cout << "Open failed for file " << fileNames.at(ifn) << endl;
      return 1; // exit if file not found
    }
    char buf[512];  // line buffer
    do {
      dataStream.getline(buf, 512);  // Discard (pre-)header lines
      // cout << buf << endl;
      if (dataStream.eof()) {
	cout << "EOF found prematurely while searching for gJets header line" << endl;
	return 1;
      }
    }
    while (strncmp(buf, " binIndex", 9) != 0);
    // while (strncmp(buf, "no ", 3) != 0);
    // while (strncmp(buf, "no  Nobs", 8) != 0);
    for (Int_t ijet=0; ijet<Nrow; ++ijet) {
      for (Int_t ikin=0; ikin<Ncol; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (Ncol > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
	// cout << ijet << " " << ikin << endl;
	// Read next line of the stream and extract its data
	dataStream.getline(buf, 512);
	// cout << buf << endl;
	if (dataStream.eof()) {
	  cout << "EOF found prematurely while searching for jJets data" << endl;
	  return 1;
	}
	const char* token[100] = {}; // initialize to 0
	Int_t n = 0;
	// parse the line
	token[n] = strtok(buf, ":"); // first token
	if (token[n] == 0) {
	  cout << n << endl;
	  cout << "gj Empty line encountered" << endl;
	  return 1;
	}
	n++; token[n] = strtok(0, "|");  // for file with MC columns
	n++; token[n] = strtok(0, "|");  // for file with MC EB, EE columns
	n++; token[n] = strtok(0, "|");  // for file with MC EB, EE columns
	n++; token[n] = strtok(0, "|");             // Nobs
	sscanf(token[n], "%f", &Ngobst[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");            // NgobsEB
	sscanf(token[n], "%f", &NgobsEBt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");            // NgobsEE
	sscanf(token[n], "%f", &NgobsEEt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "(");             // Etrg
	sscanf(token[n], "%f", &gEtrgt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ",");             // EtrgSys  0100 = 4
	sscanf(token[n], "%f", &gEtrgSyst[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ")");             // EtrgErr  1101 = 13
	sscanf(token[n], "%f", &gEtrgErrt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "(");             // SF
	sscanf(token[n], "%f", &gSFt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ")");             // SFerr  1111 = 15
	sscanf(token[n], "%f", &gSFerrt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "(");             // RZg
	sscanf(token[n], "%f", &ZgRt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ")");             // RZgErr  0100 = 4
	sscanf(token[n], "%f", &ZgRerrt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "(");             // Fdir
	sscanf(token[n], "%f", &gFdirt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ",");             // FdirSys  0100 = 4
	sscanf(token[n], "%f", &gFdirSyst[ijet][ikin][ifn]);
	//******************** kludge *********************************************
	// gFdirSyst[ijet][ikin][ifn] = 0.3;
	//******************** kludge *********************************************
	n++; token[n] = strtok(0, "-");             // FdirErrUp  0100 = 4
	sscanf(token[n], "%f", &gFdirErrUpt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ")");             // FdirErrLow  0100 = 4
	sscanf(token[n], "%f", &gFdirErrLowt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "(");             // Purity
	sscanf(token[n], "%f", &gPurt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ")");             // PurityErr  1111 = 15
	sscanf(token[n], "%f", &gPurErrt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, " ");
	n++; token[n] = strtok(0, "(");             // ZgDR
	sscanf(token[n], "%f", &ZgDRt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, "-");             // ZgDRErrUp  0000 = 0 (not used)
	sscanf(token[n], "%f", &ZgDRerrUpt[ijet][ikin][ifn]);
	n++; token[n] = strtok(0, ")");             // ZgDRErrLow  0000 = 0 (not used)
	sscanf(token[n], "%f", &ZgDRerrLowt[ijet][ikin][ifn]);
      }
    }
  }

  Float_t dummy;
  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ikin=0; ikin<Ncol; ++ikin) {
      if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
      if (Ncol > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
      if (Nfiles == 1) {
	Ngobs[ijet][ikin] = Ngobst[ijet][ikin][0];
	NgobsEB[ijet][ikin] = NgobsEBt[ijet][ikin][0];
	NgobsEE[ijet][ikin] = NgobsEEt[ijet][ikin][0];
	ZgR[ijet][ikin] = ZgRt[ijet][ikin][0];
	ZgRerr[ijet][ikin] = ZgRerrt[ijet][ikin][0];
	gEtrg[ijet][ikin] = gEtrgt[ijet][ikin][0];
	gEtrgErr[ijet][ikin] = gEtrgErrt[ijet][ikin][0];
	gEtrgSys[ijet][ikin] = gEtrgSyst[ijet][ikin][0];
	gSF[ijet][ikin] = gSFt[ijet][ikin][0];
	gSFerr[ijet][ikin] = gSFerrt[ijet][ikin][0];
	gFdir[ijet][ikin] = gFdirt[ijet][ikin][0];
	gFdirErrUp[ijet][ikin] = gFdirErrUpt[ijet][ikin][0];
	gFdirErrLow[ijet][ikin] = gFdirErrLowt[ijet][ikin][0];
	gFdirSys[ijet][ikin] = gFdirSyst[ijet][ikin][0];
	gPur[ijet][ikin] = gPurt[ijet][ikin][0];
	gPurErr[ijet][ikin] = gPurErrt[ijet][ikin][0];
	ZgDR[ijet][ikin] = ZgDRt[ijet][ikin][0];
	ZgDRerrUp[ijet][ikin] = ZgDRerrUpt[ijet][ikin][0];
	ZgDRerrLow[ijet][ikin] = ZgDRerrLowt[ijet][ikin][0];
      } else {
	// Combine results from the separate files
	Ngobs[ijet][ikin] = 0;
	NgobsEB[ijet][ikin] = 0;
	NgobsEE[ijet][ikin] = 0;
	for (size_t i = 0; i < Ngobst[ijet][ikin].size(); ++i) {
	  Ngobs[ijet][ikin] += Ngobst[ijet][ikin][i];
	  NgobsEB[ijet][ikin] += NgobsEBt[ijet][ikin][i];
	  NgobsEE[ijet][ikin] += NgobsEEt[ijet][ikin][i];
	}
	comb(2, lumi, ZgRt[ijet][ikin], ZgRerrt[ijet][ikin], ZgR[ijet][ikin], ZgRerr[ijet][ikin]);
	comb(0, lumi, gEtrgt[ijet][ikin], gEtrgErrt[ijet][ikin], gEtrg[ijet][ikin], gEtrgErr[ijet][ikin]);
	comb(1, lumi, gEtrgt[ijet][ikin], gEtrgSyst[ijet][ikin], dummy, gEtrgSys[ijet][ikin]);
	comb(1, lumi, gSFt[ijet][ikin], gSFerrt[ijet][ikin], gSF[ijet][ikin], gSFerr[ijet][ikin]);
	comb(1, lumi, gFdirt[ijet][ikin], gFdirErrUpt[ijet][ikin], gFdir[ijet][ikin], gFdirErrUp[ijet][ikin]);
	comb(1, lumi, gFdirt[ijet][ikin], gFdirErrLowt[ijet][ikin], dummy, gFdirErrLow[ijet][ikin]);
	comb(1, lumi, gFdirt[ijet][ikin], gFdirSyst[ijet][ikin], dummy, gFdirSys[ijet][ikin]);
	comb(1, lumi, gPurt[ijet][ikin], gPurErrt[ijet][ikin], gPur[ijet][ikin], gPurErr[ijet][ikin]);
	comb(0, lumi, ZgDRt[ijet][ikin], ZgDRerrUpt[ijet][ikin], ZgDR[ijet][ikin], ZgDRerrUp[ijet][ikin]);
	comb(0, lumi, ZgDRt[ijet][ikin], ZgDRerrLowt[ijet][ikin], dummy, ZgDRerrLow[ijet][ikin]);
      }
    }
  }

  return 0;

}  // ---------------------------------------------------------------

void comb(Int_t mode, V1F lumi, V1F v, V1F e, Float_t& val, Float_t& err) {
  // Combine values and errors from three sources, weighted by lumi
  // Errors are relative.
  // mode = 0: inputs uncorrelated
  // mode = 1: inputs correlated
  // mode = 2: inputs 1+ correlated, indep of input 0.
  std::vector<Float_t> wt;
  Float_t sumlumi = 0;
  for (auto ilumi : lumi) sumlumi += ilumi;
  for (auto ilumi : lumi) wt.push_back(ilumi/sumlumi);
  val = 0;
  for (size_t i = 0; i < wt.size(); ++i) val += wt[i]*v[i];
  if (mode == 0) {
    err = 0;
    for (size_t i = 0; i < wt.size(); ++i) err += Power(wt[i]*v[i]*e[i], 2);
    err = Sqrt(err) / val;
  } else if (mode == 1) {
    err = 0;
    for (size_t i = 0; i < wt.size(); ++i) err += wt[i]*v[i]*e[i];
    err /= val;
  } else if (mode == 2) {
    std::vector<Float_t> wtc;  // vector of weights for correlated components
    wtc.push_back(0);  // first value is uncorrelated
    Float_t sumlumic = 0, ex0 = 0;
    for (size_t i = 1; i < lumi.size(); ++i) sumlumic += lumi.at(i);
    for (size_t i = 1; i < lumi.size(); ++i) wtc.push_back(lumi.at(i)/sumlumic);
    for (size_t i = 1; i < wtc.size(); ++i) ex0 += wtc[i]*v[i]*e[i];
    err = Sqrt(Power(wt[0]*v[0]*e[0], 2) + Power(ex0, 2)) / val;
  }

}  // ---------------------------------------------------------------

Int_t getData_DR(const char* fileName,
		 std::vector <V1F >& ZgDR,
		 std::vector <V1F >& ZgDRerrUp,
		 std::vector <V1F >& ZgDRerrLow,
		 Float_t& DRscaleErr,
		 Float_t& DY0bPurErr,
		 Float_t& DYtrigEffErr,
		 Float_t& LeptonSFerr,
		 Float_t& btagSFerr) {
  //
  // Read data from the file, parse, and set data arrays.
  //
  Int_t Nrow = ZgDR.size();
  Int_t Ncol = ZgDR[0].size();
  ifstream dataStream;
  cout << fileName << endl;
  dataStream.open(fileName); // open the data file
  if (!dataStream.good()) {
  cout << "Open failed for file " << fileName << endl;
    return 1; // exit if file not found
  }
  char buf[512];  // line buffer
  do {
    dataStream.getline(buf, 512);  // Discard header line
    // cout << buf << endl;
    if (dataStream.eof()) {
      cout << "EOF found prematurely while searching for gJets header line" << endl;
      return 1;
    }
  }
  while (strncmp(buf, "DR", 2) != 0);
  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ikin=0; ikin<Ncol; ++ikin) {
      if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
      if (Ncol > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
      // cout << ijet << " " << ikin << endl;
      // Read next line of the stream and extract its data
       dataStream.getline(buf, 512);
      // cout << buf << endl;
      if (dataStream.eof()) {
	cout << "EOF found prematurely while searching for DR data" << endl;
	return 1;
      }
      const char* token[40] = {}; // initialize to 0
      Int_t n = 0;
      // parse the line
      token[n] = strtok(buf, "|");                // ZgDR
      if (token[n] == 0) {
	cout << n << endl;
	cout << "DR Empty line encountered" << endl;
	return 1;
      }
      sscanf(token[n], "%f", &ZgDR[ijet][ikin]);
      n++; token[n] = strtok(0, "|");  // Double ratio scale error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DRscaleErr);
      n++; token[n] = strtok(0, "|");             // ZgDRErrUp  0100 = 4
      sscanf(token[n], "%f", &ZgDRerrUp[ijet][ikin]);
      n++; token[n] = strtok(0, "|");             // ZgDRErrLow  0100 = 4
      sscanf(token[n], "%f", &ZgDRerrLow[ijet][ikin]);
      n++; token[n] = strtok(0, "|");             // 0b purity
      n++; token[n] = strtok(0, "|");             // 0b purity error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DY0bPurErr);
      n++; token[n] = strtok(0, "|");             // 0b trigger eff.
      n++; token[n] = strtok(0, "|");             // 0b trigger eff. error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &DYtrigEffErr);
      n++; token[n] = strtok(0, "|");             // lepton SF
      n++; token[n] = strtok(0, "|");             // lepton SF error  1111 = 15
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &LeptonSFerr);
      n++; token[n] = strtok(0, "|");             // btag SF
      n++; token[n] = strtok(0, "|");             // btag SF error  (not used)
      if (ijet == 0 && ikin == 0) sscanf(token[n], "%f", &btagSFerr);
    }
  }
  return 0;
}  // ---------------------------------------------------------------

Int_t getData_DY(const char* fileName,
		 std::vector <V2F >& Rb0,
		 std::vector <V2F >& Rb0stat,
		 std::vector <V2F >& Rb0MCstat,
		 std::vector <V2F >& Rb0sysUp,
		 std::vector <V2F >& Rb0sysLow,
		 std::vector <V2F >& Rb0sysKin,
		 std::vector <V2F >& Rb0sysPur) {
  //
  // Read data from the file, parse, and set data arrays.
  //
/*
  Notes on correlations:
  All DY nuisances except DYsysKin are correlated across all kinematic bins.
  DYstat is correlated across Njets bins for Njets >= 7; irrelevant for Nb = 0   1111 = 15, 3999, 100
  Rb0MCstat & DYsysNjUp, Low are correlated over kinematic bins;
    irrelevant for Nb < 4                                              1111 = 15, 9999, 4100
  DYsysPur is correlated across Njets bins for Njets >= 5,
    and across Nb bins for Nb >=2                                      1111 = 15, 2299        
  DYsysKin is uncorrelated (though all are zero for Nb = 0).           0100 =  4, 9999, 100
*/
  Int_t Nrow = Rb0.size();
  Int_t Ncol = Rb0[0].size();
  Int_t Nplane = Rb0[0][0].size();
  Float_t Rb0sysJup, Rb0sysJlow, Rb0systtZ;
  ifstream dataStream;
  dataStream.open(fileName); // open the data file
  if (!dataStream.good()) {
  cout << "Open failed for file " << fileName << endl;
    return 1; // exit if file not found
  }
  char buf[512];  // line buffer
  dataStream.getline(buf, 512);  // Discard header line

  for (Int_t ijet=0; ijet<Nrow; ++ijet) {
    for (Int_t ib=0; ib<Ncol; ++ib) {
      if (ijet == 0 && ib > 2) continue;  // Skip Njet=2, Nb>2 bins
      for (Int_t ikin=0; ikin<Nplane; ++ikin) {
	if (ijet > 2 && (ikin == 0 || ikin == 3)) continue;  // Skip high-Njet, low-HT bins
	if (Nplane > 10 && ijet > 2 && ikin == 6) continue;  // Skip high-Njet, low-HT bins
        //
        // Read next line of the stream and extract its data
        //
	// cout << ijet << " " << ib << endl;
	dataStream.getline(buf, 512);
	if (dataStream.eof()) {
	  cout << "EOF found prematurely while searching for DY data" << endl;
	  return 1;
	}
	const char* token[20] = {}; // initialize to 0
	Int_t n = 0;
        // parse the line
	token[n] = strtok(buf, "|"); // first token
	if (token[n] == 0) {
	  cout << "Empty line encountered" << endl;
	  return 1;
	}
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0[ijet][ib][ikin]);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0stat[ijet][ib][ikin]);  //   1111 = 15, 3999, 100
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0MCstat[ijet][ib][ikin]); //  1111 = 15, 9999, 4100 
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0systtZ);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysJup);
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysJlow);
	Rb0sysUp[ijet][ib][ikin] = Sqrt(Power(Rb0systtZ, 2) + Power(Rb0sysJup, 2));  //  1111 = 15, 9999, 4100
	Rb0sysLow[ijet][ib][ikin] = Sqrt(Power(Rb0systtZ, 2) + Power(Rb0sysJlow, 2));  // 1111 = 15, 9999, 4100
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysKin[ijet][ib][ikin]);  // 0100 =  4, 9999, 100
	n++; token[n] = strtok(0, "|");  sscanf(token[n], "%f", &Rb0sysPur[ijet][ib][ikin]);  // 1111 = 15, 2299
      }
    }
  }
  return 0;
}  // ---------------------------------------------------------------

void setCorrelationLabels(TH1* histo, const unsigned correlVar, const int minCorrelBin=9999,
			  const int maxCorrelBin=0) {
  TString hname = histo->GetName();
  cout << hname;
  hname.Remove(0,1);  // Strip off the "h" at start of histogram name
  for (Int_t bin=1; bin<=histo->GetNbinsX(); ++bin) {
    TString binLabel = histo->GetXaxis()->GetBinLabel(bin);
    // binLabel = "_NJets<n>_BTags<m>_MHT<p>_HT<q>"
    // Remove substring for variable with correlated uncertainties.
    // correlVar bit 0 HT; bit 1 MHT; bit 2 Nb; bit 3 Njets

    if (correlVar & 8) {
      Ssiz_t index = binLabel.Index("_NJets");
      int thisMinCorrelBin = minCorrelBin/1000 % 10;
      int thisMaxCorrelBin = maxCorrelBin/1000 % 10;
      if (thisMinCorrelBin == 9 && thisMaxCorrelBin == 0) {
	binLabel.Remove(index, 7);  // Correlated across all NJets bins
      } else {
	TString theBin = binLabel(index+6,1);  // Correlated across some NJets bins
	char newBin[10];
	if (thisMinCorrelBin < 9) sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), thisMinCorrelBin));
	if (thisMaxCorrelBin > 0) sprintf(newBin, "%1d", TMath::Max(theBin.Atoi(), thisMaxCorrelBin));
	binLabel.Replace(index+6, 1, newBin);
      }
    }

    if (correlVar & 4) {
      Ssiz_t index = binLabel.Index("_BTags");
      int thisMinCorrelBin = minCorrelBin/100 % 10;
      int thisMaxCorrelBin = maxCorrelBin/100 % 10;
      if (thisMinCorrelBin == 9 && thisMaxCorrelBin == 0) {
    	binLabel.Remove(index, 7);  // Correlated across all BTag bins
      } else {
	TString theBin = binLabel(index+6,1);  // Correlated across some BTag bins
	char newBin[10];
	if (thisMinCorrelBin < 9) sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), thisMinCorrelBin));
	if (thisMaxCorrelBin > 0) sprintf(newBin, "%1d", TMath::Max(theBin.Atoi(), thisMaxCorrelBin));
	binLabel.Replace(index+6, 1, newBin);
      }
    }

    if (correlVar & 2) {
      Ssiz_t index = binLabel.Index("_MHT");
      int thisMinCorrelBin = minCorrelBin/10 % 10;
      int thisMaxCorrelBin = maxCorrelBin/10 % 10;
      if (thisMinCorrelBin == 9 && thisMaxCorrelBin == 0) {
    	binLabel.Remove(index, 5);  // Correlated across all MHT bins
      } else {
	TString theBin = binLabel(index+4,1);  // Correlated across some MHT bins
	char newBin[10];
	if (thisMinCorrelBin < 9) sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), thisMinCorrelBin));
	if (thisMaxCorrelBin > 0) sprintf(newBin, "%1d", TMath::Max(theBin.Atoi(), thisMaxCorrelBin));
	binLabel.Replace(index+4, 1, newBin);
      }
    }

    if (correlVar & 1) {
      Ssiz_t index = binLabel.Index("_HT");
      int thisMinCorrelBin = minCorrelBin % 10;
      int thisMaxCorrelBin = maxCorrelBin % 10;
      if (thisMinCorrelBin == 9 && thisMaxCorrelBin == 0) {
    	binLabel.Remove(index, 4);  // Correlated across all HT bins
      } else {
	TString theBin = binLabel(index+3,1);  // Correlated across some HT bins
	char newBin[10];
	if (thisMinCorrelBin < 9) sprintf(newBin, "%1d", TMath::Min(theBin.Atoi(), thisMinCorrelBin));
	if (thisMaxCorrelBin > 0) sprintf(newBin, "%1d", TMath::Max(theBin.Atoi(), thisMaxCorrelBin));
	binLabel.Replace(index+3, 1, newBin);
      }
    }

    binLabel.Prepend(hname);
    if (bin == 135) cout << ",  binLabel(135) = " << binLabel;
    histo->GetXaxis()->SetBinLabel(bin, binLabel.Data());
  }
  cout << endl;
}  // ---------------------------------------------------------------
