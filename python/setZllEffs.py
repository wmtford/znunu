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

import RA2b
import ROOT
import histoZmassFits

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

doRA2bFits = False

########## trigger effs from manuel ##############
########## from Nov 29th RA2b talk  ##############
#### dimu trigger eff high and flat ##############
#### diel trigger drops at high HT  ##############
#### diel HT binning [300-1000,1000+]
trig_m = [(0.993,0.003)]
trig_e = [(0.985,0.002),
          (0.943,0.002)]

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
    fitjb = histoZmassFits.purityFits()  # use Z mass histograms from RA2bZinvAnalysis
    fits = [fitjb[0], fitjb[1], fitjb[2], fitjb[2], fitjb[2]]

########## get the efficiency file ################
# effFile = ROOT.TFile("../plots/histograms/effHists.root","UPDATE")
effFile = ROOT.TFile("effHists.root","UPDATE")  # wtf


######### set the purities found above ############
h_pur_m = effFile.Get("h_pur_m")
h_pur_e = effFile.Get("h_pur_e")

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
njDict = {0: '& $\\njets=2$        ',
          1: '& $3\leq\\njets\leq4$',
          2: '& $\\njets\geq5$     ',}
for lep in range(2):
    print '\hline'
    print zDict[lep],
    for nj in range(3):
        print njDict[nj],
        for nb in range(3):
            print " & $"+str(round(fits[nj][nb][lep][0],3))+'\pm'+str(round(fits[nj][nb][lep][1],3))+"$",
        print " \\\\ "

######### set the trig effs from manuel ############
h_trig_m = effFile.Get("h_trig_m1")
h_trig_e = effFile.Get("h_trig_e2")

for i in range(len(trig_m)):
    h_trig_m.SetBinContent(1,trig_m[i][0])
    h_trig_m.SetBinError(1,trig_m[i][1])

for i in range(len(trig_e)):
    h_trig_e.SetBinContent(i+1,trig_e[i][0])
    h_trig_e.SetBinError(i+1,trig_e[i][1])

h_trig_m.Write(h_trig_m.GetName(),2)
h_trig_e.Write(h_trig_e.GetName(),2)

effFile.Close()

