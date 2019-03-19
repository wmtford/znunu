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

# Get trigger efficiencies from Sam's file for Run2, 2017
trigEffFile = ROOT.TFile("../plots/histograms/triggersRa2bRun2_v1.root", "READ")
f_trig_eb_2016 = trigEffFile.Get("tf1_SinglePhotonBarrelLoose_hists_Run2016_JetHT")
f_trig_eb_2016.SetName("f_trig_eb_2016")
f_trig_ec_2016 = trigEffFile.Get("tf1_SinglePhotonEndcapLoose_hists_Run2016_JetHT")
f_trig_ec_2016.SetName("f_trig_ec_2016")
f_trig_eb_2017 = trigEffFile.Get("tf1_SinglePhotonBarrelLoose_hists_Run2017_JetHT")
f_trig_eb_2017.SetName("f_trig_eb_2017")
f_trig_ec_2017 = trigEffFile.Get("tf1_SinglePhotonEndcapLoose_hists_Run2017_JetHT")
f_trig_ec_2017.SetName("f_trig_ec_2017")
f_trig_eb_2018 = trigEffFile.Get("tf1_SinglePhotonBarrelLoose_hists_Run2018_JetHT")
f_trig_eb_2018.SetName("f_trig_eb_2018")
f_trig_ec_2018 = trigEffFile.Get("tf1_SinglePhotonEndcapLoose_hists_Run2018_JetHT")
f_trig_ec_2018.SetName("f_trig_ec_2018")

# SFfile = ROOT.TFile("../plots/histograms/SFcorrections.Photons.root", "READ")
# h_SF_g = SFfile.Get("h_MHT")
# h_SF_g.SetName("h_SFg_MHT")

FdirFile = ROOT.TFile("../plots/histograms/fragmentation.root", "READ")
h_Fdir = FdirFile.Get("bin46_NJets8910")
h_Fdir.SetName("h_bin46_NJets8910")

## photon purity (Andrew provides these)
EB = {}
EE = {}
## Feb 17, 2017 email, for era 2016 analysis
## apply to data
## split by barrel (eb) and endcap (ec)
## binning in MHT: [<225, 225-250, 250-300, 300-350, 350-500, 500+]
# pur_title = "MHT"
# pur_bins = array('d', [0, 225, 250, 300, 350, 500, 1000])
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

# Run 2, Andrew's Jan. 27, 2019 email:
pur_title = "MHT"
pur_bins = array('d', [0, 300, 350, 600, 2500])  # First bin measured from 250, last bin to 1000
# 2016:
pur_eb_2016 = []
pur_ec_2016 = []
EB['avg'] = [0.9479783904003333, 0.9612409354666666, 0.9657104593263334, 0.9819289765320001]
EB['avgErr'] = [0.031016399071333245, 0.01721366277833347, 0.015233616496666658, 0.008033783562999885]
EE['avg'] = [0.898884215471, 0.9232578112983333, 0.9296478600020001, 0.9149493319753333]
EE['avgErr'] = [0.04399141754999991, 0.04426306202466668, 0.02750138579299999, 0.047701439598333395]
for i in range(0,4):
    pur_eb_2016.append( (EB['avg'][i], EB['avgErr'][i]) )
    pur_ec_2016.append( (EE['avg'][i], EE['avgErr'][i]) )

# 2017:
pur_eb_2017 = []
pur_ec_2017 = []
# Andrew's Jan. 27, 2019 email
# EB['avg'] [0.9405091609876667, 0.9394312251323332, 0.9529458327956667, 0.9742914092513333]
# EB['avgErr'] [0.009999098545333318, 0.014585434074333214, 0.011037659601666694, 0.005863896771333255]
# EE['avg'] [0.9599167930893332, 0.9644865701173333, 0.9468864281673334, 0.8846943180493333]
# EE['avgErr'] [0.013086806293333253, 0.014466961288666758, 0.02317752860333333, 0.0495032579343333]
# Tribeni 19 March, 2019 e-mail:  photon purity (Gjets_0p4)
EB['avg'] = [0.9340764831373334, 0.9320274834113333, 0.945903291095, 0.9736938879303333]
EB['avgErr'] = [0.011099681395666638, 0.016106591866333342, 0.01255540386100007, 0.006075190538333275]
EE['avg'] = [0.9447821990946667, 0.9473740980093334, 0.931639983517, 0.9158964067776667]
EE['avgErr'] = [0.017447269969666612, 0.021160995631666624, 0.028633940094, 0.03600229982366665]
for i in range(0,4):
    pur_eb_2017.append( (EB['avg'][i], EB['avgErr'][i]) )
    pur_ec_2017.append( (EE['avg'][i], EE['avgErr'][i]) )

# 2018:
pur_eb_2018 = []
pur_ec_2018 = []
# Andrew's Jan. 27, 2019 email
# EB['avg'] [0.934561081295, 0.935941864678, 0.9528332801676668, 0.9698275681449999]
# EB['avgErr'] [0.010881257075999962, 0.015425372193999976, 0.010485711729666725, 0.0064986703679998925]
# EE['avg'] [0.955755605372, 0.9603043977869999, 0.948359821654, 0.8919489522903333]
# EE['avgErr'] [0.013496460590999959, 0.016012506070000043, 0.016055838115, 0.032147232375666634]
# Tribeni 19 March, 2019 e-mail:  photon purity (Gjets_0p4)
EB['avg'] = [0.9280759956736667, 0.9285315544210001, 0.9457477890636666, 0.9691437937430001]
EB['avgErr'] = [0.011992508267333335, 0.016934029379000126, 0.011962957321666634, 0.006785803926000034]
EE['avg'] = [0.9405511421239999, 0.9439760420550001, 0.9345166398170001, 0.9283238776746666]
EE['avgErr'] = [0.017382353302999975, 0.022191131761999938, 0.019755000849999926, 0.021305790389333334]
for i in range(0,4):
    pur_eb_2018.append( (EB['avg'][i], EB['avgErr'][i]) )
    pur_ec_2018.append( (EE['avg'][i], EE['avgErr'][i]) )


###################################################
# H_SF_g1 = effFile.Get("h_SF_g1")
# h_SF_g1.SetBinContent(1,SF[0])
# h_SF_g1.SetBinError(1,SF[1])
# h_SF_g1.Write(h_SF_g1.GetName(),2)

# h_frag1 = effFile.Get("h_frag1")
# h_frag1.SetBinContent(1,frag[0])
# h_frag1.SetBinError(1,frag[1])
# h_frag1.Write(h_frag1.GetName(),2)

########## get the efficiency file ################
effFile = ROOT.TFile("../plots/histograms/effHists.root","UPDATE")
# effFile = ROOT.TFile("effHists.root","UPDATE")  # wtf

########## update the efficiency file #############

# f_trig_eb_2016.Write(f_trig_eb_2016.GetName(),2)
# f_trig_ec_2016.Write(f_trig_ec_2016.GetName(),2)
# f_trig_eb_2017.Write(f_trig_eb_2017.GetName(),2)
# f_trig_ec_2017.Write(f_trig_ec_2017.GetName(),2)
# f_trig_eb_2018.Write(f_trig_eb_2018.GetName(),2)
# f_trig_ec_2018.Write(f_trig_ec_2018.GetName(),2)

# h_SF_g.Write(h_SF_g.GetName(), 2)
h_Fdir.Write(h_Fdir.GetName(), 2)

h_trig_eb = effFile.Get("h_trig_eb")
if (not h_trig_eb):
    h_trig_eb = ROOT.TH1F("h_trig_eb", "gJets trigger efficiencies, EB", len(trig_bins)-1, trig_bins)
else:
    h_trig_eb.SetBins(len(trig_bins)-1, trig_bins)
h_trig_eb.GetXaxis().SetTitle(trig_title)
for i in range(len(trig_eb)):
    h_trig_eb.SetBinContent(i+1,trig_eb[i][0])
    h_trig_eb.SetBinError(i+1,trig_eb[i][1])
# h_trig_eb.Write(h_trig_eb.GetName(),2)

h_trig_ec = effFile.Get("h_trig_ec")
if (not h_trig_ec):
    h_trig_ec = ROOT.TH1F("h_trig_ec", "gJets trigger efficiencies, EC", len(trig_bins)-1, trig_bins)
else:
    h_trig_ec.SetBins(len(trig_bins)-1, trig_bins)
h_trig_ec.GetXaxis().SetTitle(trig_title)
for i in range(len(trig_ec)):
    h_trig_ec.SetBinContent(i+1,trig_ec[i][0])
    h_trig_ec.SetBinError(i+1,trig_ec[i][1])
# h_trig_ec.Write(h_trig_ec.GetName(),2)

h_pur_eb = effFile.Get("h_pur_eb")
if (not h_pur_eb):
    h_pur_eb = ROOT.TH1F("h_pur_eb", "gJets purities, EB", len(pur_bins)-1, pur_bins)
else:
    h_pur_eb.SetBins(len(pur_bins)-1, pur_bins)
h_pur_eb.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_eb)):
    h_pur_eb.SetBinContent(i+1,pur_eb[i][0])
    h_pur_eb.SetBinError(i+1,pur_eb[i][1])
h_pur_eb.Write(h_pur_eb.GetName(),2)

h_pur_ec = effFile.Get("h_pur_ec")
if (not h_pur_ec):
    h_pur_ec = ROOT.TH1F("h_pur_ec", "gJets purities, EC", len(pur_bins)-1, pur_bins)
else:
    h_pur_ec.SetBins(len(pur_bins)-1, pur_bins)
h_pur_ec.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_ec)):
    h_pur_ec.SetBinContent(i+1,pur_ec[i][0])
    h_pur_ec.SetBinError(i+1,pur_ec[i][1])
h_pur_ec.Write(h_pur_ec.GetName(),2)

h_pur_eb_2016 = effFile.Get("h_pur_eb_2016")
if (not h_pur_eb_2016):
    h_pur_eb_2016 = ROOT.TH1F("h_pur_eb_2016", "gJets purities, EB", len(pur_bins)-1, pur_bins)
else:
    h_pur_eb_2016.SetBins(len(pur_bins)-1, pur_bins)
h_pur_eb_2016.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_eb_2016)):
    h_pur_eb_2016.SetBinContent(i+1,pur_eb_2016[i][0])
    h_pur_eb_2016.SetBinError(i+1,pur_eb_2016[i][1])
h_pur_eb_2016.Write(h_pur_eb_2016.GetName(),2)

h_pur_ec_2016 = effFile.Get("h_pur_ec_2016")
if (not h_pur_ec_2016):
    h_pur_ec_2016 = ROOT.TH1F("h_pur_ec_2016", "gJets purities, EC", len(pur_bins)-1, pur_bins)
else:
    h_pur_ec_2016.SetBins(len(pur_bins)-1, pur_bins)
h_pur_ec_2016.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_ec_2016)):
    h_pur_ec_2016.SetBinContent(i+1,pur_ec_2016[i][0])
    h_pur_ec_2016.SetBinError(i+1,pur_ec_2016[i][1])
h_pur_ec_2016.Write(h_pur_ec_2016.GetName(),2)

h_pur_eb_2017 = effFile.Get("h_pur_eb_2017")
if (not h_pur_eb_2017):
    h_pur_eb_2017 = ROOT.TH1F("h_pur_eb_2017", "gJets purities, EB", len(pur_bins)-1, pur_bins)
else:
    h_pur_eb_2017.SetBins(len(pur_bins)-1, pur_bins)
h_pur_eb_2017.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_eb_2017)):
    h_pur_eb_2017.SetBinContent(i+1,pur_eb_2017[i][0])
    h_pur_eb_2017.SetBinError(i+1,pur_eb_2017[i][1])
h_pur_eb_2017.Write(h_pur_eb_2017.GetName(),2)

h_pur_ec_2017 = effFile.Get("h_pur_ec_2017")
if (not h_pur_ec_2017):
    h_pur_ec_2017 = ROOT.TH1F("h_pur_ec_2017", "gJets purities, EC", len(pur_bins)-1, pur_bins)
else:
    h_pur_ec_2017.SetBins(len(pur_bins)-1, pur_bins)
h_pur_ec_2017.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_ec_2017)):
    h_pur_ec_2017.SetBinContent(i+1,pur_ec_2017[i][0])
    h_pur_ec_2017.SetBinError(i+1,pur_ec_2017[i][1])
h_pur_ec_2017.Write(h_pur_ec_2017.GetName(),2)

h_pur_eb_2018 = effFile.Get("h_pur_eb_2018")
if (not h_pur_eb_2018):
    h_pur_eb_2018 = ROOT.TH1F("h_pur_eb_2018", "gJets purities, EB", len(pur_bins)-1, pur_bins)
else:
    h_pur_eb_2018.SetBins(len(pur_bins)-1, pur_bins)
h_pur_eb_2018.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_eb_2018)):
    h_pur_eb_2018.SetBinContent(i+1,pur_eb_2018[i][0])
    h_pur_eb_2018.SetBinError(i+1,pur_eb_2018[i][1])
h_pur_eb_2018.Write(h_pur_eb_2018.GetName(),2)

h_pur_ec_2018 = effFile.Get("h_pur_ec_2018")
if (not h_pur_ec_2018):
    h_pur_ec_2018 = ROOT.TH1F("h_pur_ec_2018", "gJets purities, EC", len(pur_bins)-1, pur_bins)
else:
    h_pur_ec_2018.SetBins(len(pur_bins)-1, pur_bins)
h_pur_ec_2018.GetXaxis().SetTitle(pur_title)
for i in range(len(pur_ec_2018)):
    h_pur_ec_2018.SetBinContent(i+1,pur_ec_2018[i][0])
    h_pur_ec_2018.SetBinError(i+1,pur_ec_2018[i][1])
h_pur_ec_2018.Write(h_pur_ec_2018.GetName(),2)


effFile.Close()
