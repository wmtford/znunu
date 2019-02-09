#!/usr/bin/python
"""Run this script to set the zll related efficiencies in
plots/histograms/effHists.root
Purity fits are run using RA2b.getZmassFitPlot()
current binning for purity is 
(NJ_2: Nb0, Nb1, Nb2)
(NJ_3-4:  Nb0, Nb1, Nb2+) 
(NJ_5+:  Nb0, Nb1, Nb2+) 
Trigger efficiencies are hard coded below.
Lepton scale factors are two separate files
provided by Frank."""

from array import array
import RA2b
import ROOT
import histoZmassFits

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

doPurityFits = True
doRA2bFits = False

########## trigger effs from manuel ##############
########## from Nov 29th RA2b talk  ##############
#### dimu trigger eff high and flat ##############
#### diel trigger drops at high HT  ##############
#### diel HT binning [300-1000,1000+]
# trig_m = [(0.993,0.003)]
# trig_m_bins = array('d', [300, 2300])
# trig_title = "HT"
# trig_e = [(0.985,0.002),
#           (0.943,0.002)]
# trig_e_bins = array('d', [300, 1000, 2300])
########## trigger effs from AN-18-271
trig_title = "Pt(Z)"
trig_m = [(0.980,0.007)]  # Alt. measurement [(0.992,0.003)]
trig_m_bins = array('d', [250, 2300])
trig_e = [(0.988,0.005)]  # Alt. measurement [(0.995,0.003)]
trig_e_bins = array('d', [250, 2300])

if (doPurityFits):
    ########## run fits to get purity ################
    if (doRA2bFits):
        fit_2j = []
        fit_3to4j = []
        fit_5jplus = []

        # nb=12 is nb>=2
        for nb in [0,1,12]:
            fit_2j.append(RA2b.getZmassFitPlot(bJetBin=nb,nJetBin=1,savePlot=True,plotMC=False))  # wtf
            fit_3to4j.append(RA2b.getZmassFitPlot(bJetBin=nb,nJetBin=2,savePlot=True,plotMC=False))  # wtf
            fit_5jplus.append(RA2b.getZmassFitPlot(bJetBin=nb,extraCuts='NJets>=5',savePlot=True,plotMC=False))  # wtf
            # fit_2j.append(RA2b.getZmassFitPlot(bJetBin=nb,nJetBin=1,savePlot=True))
            # fit_3to4j.append(RA2b.getZmassFitPlot(bJetBin=nb,nJetBin=2,savePlot=True))
            # fit_5jplus.append(RA2b.getZmassFitPlot(bJetBin=nb,extraCuts='NJets>=5',savePlot=True))

        fits = [fit_2j,fit_3to4j,fit_5jplus,fit_5jplus,fit_5jplus]
    else:
        fitjb = histoZmassFits.purityFits('../outputs/histsDY_Run2v16.root', '../outputs/histsDYMC_2017v16.root')  # use Z mass histograms from RA2bZinvAnalysis
        fits = [fitjb[0], fitjb[1], fitjb[2], fitjb[2], fitjb[2]]

########## get the scale factors files and extract histograms ################
# These are now taken directly from the (in princple year-dependent) files from Frank and Emily
# SFfile_m = ROOT.TFile("../plots/histograms/SFcorrections.Muons.root", "READ")
# h_SF_m = SFfile_m.Get("h_MHT")
# h_SF_m.SetName("h_SFm_MHT")
# SFfile_e = ROOT.TFile("../plots/histograms/SFcorrections.Electrons.root", "READ")
# h_SF_e = SFfile_e.Get("h_MHT")
# h_SF_e.SetName("h_SFe_MHT")

########## get the efficiency file ################
effFile = ROOT.TFile("../plots/histograms/effHists.root","UPDATE")
# effFile = ROOT.TFile("effHists.root","UPDATE")  # wtf

if (doPurityFits):
    ######### set the purities found above ############
    h_pur_m = effFile.Get("h_pur_m")
    if (not h_pur_m):
        h_pur_m = ROOT.TH1F("h_pur_m", "Zmm purities vs (Njet, Nb)", 19, .5, 19.5)
    h_pur_m.GetXaxis().SetTitle("(NJets, Nb) bin")
    h_pur_e = effFile.Get("h_pur_e")
    if (not h_pur_e):
        h_pur_e = ROOT.TH1F("h_pur_e", "Zee purities vs (Njet, Nb)", 19, .5, 19.5)
    h_pur_e.GetXaxis().SetTitle("(NJets, Nb) bin")

    Bin = 1
    for nj in range(1,6):
        for nb in range(4):
            if(nb==3):
                if(nj==1):
                    continue
                else:
                    h_pur_m.SetBinContent(Bin,fits[nj-1][nb-1][0][0])
                    h_pur_m.SetBinError(Bin,fits[nj-1][nb-1][0][1])
                    h_pur_e.SetBinContent(Bin,fits[nj-1][nb-1][1][0])
                    h_pur_e.SetBinError(Bin,fits[nj-1][nb-1][1][1])
            else:
                h_pur_m.SetBinContent(Bin,fits[nj-1][nb][0][0])
                h_pur_m.SetBinError(Bin,fits[nj-1][nb][0][1])
                h_pur_e.SetBinContent(Bin,fits[nj-1][nb][1][0])
                h_pur_e.SetBinError(Bin,fits[nj-1][nb][1][1])
            Bin+=1
    h_pur_m.Write(h_pur_m.GetName(),2)
    h_pur_e.Write(h_pur_e.GetName(),2)

    zDict = {0: '\multirow{3}{*}{\zmm}',
             1: '\multirow{3}{*}{\zee}',}
    njDict = {0: '& $2\leq\\njets\leq3$        ',
              1: '& $4\leq\\njets\leq5$',
              2: '& $\\njets\geq6$     ',}
    for lep in range(2):
        print '\hline'
        print zDict[lep],
        for nj in range(3):
            print njDict[nj],
            for nb in range(3):
                print " & $"+str(round(fits[nj][nb][lep][0],3))+'\pm'+str(round(fits[nj][nb][lep][1],3))+"$",
            print " \\\\ "

######### set the trig effs ############
h_trig_m = effFile.Get("h_trig_m")
if (not h_trig_m):
    h_trig_m = ROOT.TH1F("h_trig_m", "Zmm trigger effs vs HT", len(trig_m_bins)-1, trig_m_bins)
else:
    h_trig_m.SetBins(len(trig_m_bins)-1, trig_m_bins)
h_trig_m.GetXaxis().SetTitle(trig_title)
h_trig_e = effFile.Get("h_trig_e")
if (not h_trig_e):
    h_trig_e = ROOT.TH1F("h_trig_e", "Zee trigger effs vs HT", len(trig_e_bins)-1, trig_e_bins)
else:
    h_trig_e.SetBins(len(trig_e_bins)-1, trig_e_bins)
h_trig_e.GetXaxis().SetTitle(trig_title)

for i in range(len(trig_m)):
    h_trig_m.SetBinContent(1,trig_m[i][0])
    h_trig_m.SetBinError(1,trig_m[i][1])

for i in range(len(trig_e)):
    h_trig_e.SetBinContent(i+1,trig_e[i][0])
    h_trig_e.SetBinError(i+1,trig_e[i][1])

h_trig_m.Write(h_trig_m.GetName(),2)
h_trig_e.Write(h_trig_e.GetName(),2)

######### set the scale factors from Frank's SF files ############
# h_SF_m.Write(h_SF_m.GetName(), 2)
# h_SF_e.Write(h_SF_e.GetName(), 2)


effFile.Close()

