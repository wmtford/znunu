#!/usr/bin/python
"""Run this script to set the photon related efficiencies in
plots/histograms/effHists.root
All efficiency factors are hard coded below"""

from array import array
import ROOT

########## photon inputs: (value, absolute error) #####
## photon efficiency scale factor
## apply to MC
######## obsolete, use SFcorrections.Photons.root #####
SF = (0.99, 0.01)

## photon fragmentation
## apply to data
######## obsolete, use fragmentation.root #############
frag = (0.92, 0.07)

## photon trigger efficiency 
## Andrew provided these (see Nov. 22 email)
## apply as a SF to MC (don't apply any simulated triggers to MC)
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<300, 300-350, 350-500, 500-750, 750+]
trig_bins = array('d', [0, 300, 350, 500, 750, 1000])
trig_title = "MHT"
trig_eb = [(0.969, 0.002),
           (0.983, 0.001),
           (0.985, 0.001),
           (0.984, 0.001),
           (0.979, 0.004),]

trig_ec = [(0.953, 0.004),
           (0.974, 0.003),
           (0.984, 0.001),
           (0.989, 0.003),
           (0.980, 0.019),]

## photon purity (Andrew provides these Feb 17 email)
## apply to data
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<225, 225-250, 250-300, 300-350, 350-500, 500+]
pur_bins = array('d', [0, 225, 250, 300, 350, 500, 1000])
pur_title = "MHT"
pur_eb = [(0.9580, 0.0267),
          (0.9623, 0.0028),
          (0.9729, 0.0098),
          (0.9690, 0.0196),
          (0.9622, 0.0456),
          (0.9780, 0.0264),]

pur_ec = [(0.8879, 0.0108),
          (0.8570, 0.0114),
          (0.8969, 0.0174),
          (0.9037, 0.0175),
          (0.9342, 0.0299),
          (0.9637, 0.0169),]
###################################################

########## get the efficiency file ################
effFile = ROOT.TFile("../plots/histograms/effHists.root","UPDATE")
# effFile = ROOT.TFile("effHists.root","UPDATE")  # wtf

########## update the efficiency file #############
h_SF_g1 = effFile.Get("h_SF_g1")
h_SF_g1.SetBinContent(1,SF[0])
h_SF_g1.SetBinError(1,SF[1])
h_SF_g1.Write(h_SF_g1.GetName(),2)

h_frag1 = effFile.Get("h_frag1")
h_frag1.SetBinContent(1,frag[0])
h_frag1.SetBinError(1,frag[1])
h_frag1.Write(h_frag1.GetName(),2)

h_trig_eb = effFile.Get("h_trig_eb")
h_trig_eb.SetBins(len(trig_bins)-1, trig_bins)
h_trig_eb.GetXaxis().SetTitle(trig_title)
for i in range(len(trig_eb)):
    h_trig_eb.SetBinContent(i+1,trig_eb[i][0])
    h_trig_eb.SetBinError(i+1,trig_eb[i][1])
h_trig_eb.Write(h_trig_eb.GetName(),2)

h_trig_ec = effFile.Get("h_trig_ec")
h_trig_ec.SetBins(len(trig_bins)-1, trig_bins)
h_trig_ec.GetXaxis().SetTitle(trig_title)
for i in range(len(trig_ec)):
    h_trig_ec.SetBinContent(i+1,trig_ec[i][0])
    h_trig_ec.SetBinError(i+1,trig_ec[i][1])
h_trig_ec.Write(h_trig_ec.GetName(),2)

h_pur_eb = effFile.Get("h_pur_eb")
h_pur_eb.SetBins(len(pur_bins)-1, pur_bins)
h_pur_eb.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_eb)):
    h_pur_eb.SetBinContent(i+1,pur_eb[i][0])
    h_pur_eb.SetBinError(i+1,pur_eb[i][1])
h_pur_eb.Write(h_pur_eb.GetName(),2)

h_pur_ec = effFile.Get("h_pur_ec")
h_pur_ec.SetBins(len(pur_bins)-1, pur_bins)
h_pur_ec.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_ec)):
    h_pur_ec.SetBinContent(i+1,pur_ec[i][0])
    h_pur_ec.SetBinError(i+1,pur_ec[i][1])
h_pur_ec.Write(h_pur_ec.GetName(),2)


effFile.Close()
