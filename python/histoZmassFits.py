#!/usr/bin/python
"""
Plot and fit Z mass histograms
"""
import ROOT
import RA2b

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

doMumu = True
doEe = True
removeDYkfactor = False
MCscaleM = 1
MCscaleE = 1
print '\n'
print "removeDYkfactor = "+str(removeDYkfactor)

# period = 5  # 2016
# lumi = 35.9
# DataFileM = ROOT.TFile('../outputs/histsDY_2016v12.root')
# DataFileE = DataFileM
# # MCfile = ROOT.TFile('../outputs/histsDYMC_2016v12.root')
# MCfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')

# period = 6  # 2017
# lumi = 41.5
# DataFileM = ROOT.TFile('../outputs/histsDYmm_2017v15.root')
# DataFileE = ROOT.TFile('../outputs/histsDYee_2017v15.root')
# MCfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
# MCscaleM = lumi/35.9
# MCscaleE = lumi/35.9

period = 8  # 2018
lumi = 14.0
DataFileM = ROOT.TFile('../outputs/histsDYmm_2018v15.root')
DataFileE = ROOT.TFile('../outputs/histsDYee_2018v15.root')
MCfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
MCscaleM = lumi/35.9
MCscaleE = 13.5/35.9

reactions = ["tt", "ttz", "VV", "dy"]

# outFile = ROOT.TFile('ZmassRooPlots.root', 'RECREATE')

d1 = {}
d1["loose"] = DataFileM.Get("hZmass_zmm")
d2 = {}
d2["loose"] = DataFileE.Get("hZmass_zee")
jbGroup = 0
for jbin in ["2j", "3j", "5j", "all"]:
  for bbin in ["0b", "1b", "2b"]:
    if (jbin == "all"):
      if (bbin == "2b"):
        histNameRoot = "hZmass_"
        plotNameRoot = "fitZmass_all_"
      else:
        continue
    else:
      histNameRoot = "hZmass_"+str(jbin)+str(bbin)+"_"
      plotNameRoot = "fitZmass_"+str(jbin)+str(bbin)+"_"

    hData = []
    d1["sig"] = DataFileM.Get(str(histNameRoot)+"zmm")
    hData.append(d1)
    d2["sig"] = DataFileE.Get(str(histNameRoot)+"zee")
    hData.append(d2)

    hsMC = []
    hsMM = ROOT.THStack("hsMM","dimuon mass [GeV]")
    for proc in reactions:
      hName = str(histNameRoot)+str(proc)+"mm"
      if (removeDYkfactor and 'dy' in hName):
        MCfile.Get(hName).Scale(MCscaleM/1.23)
      else:
        MCfile.Get(hName).Scale(MCscaleM)
      hsMM.Add(MCfile.Get(hName))
    hsMC.append(hsMM)
    hsEE = ROOT.THStack("hsEE","dielectron mass [GeV]")
    for proc in reactions:
      hName = str(histNameRoot)+str(proc)+"ee"
      if (removeDYkfactor and 'dy' in hName):
        MCfile.Get(hName).Scale(MCscaleE/1.23)
      else:
        MCfile.Get(hName).Scale(MCscaleE)
      hsEE.Add(MCfile.Get(hName))
    hsMC.append(hsEE)

    print ''
    print 'histNameRoot = '+str(histNameRoot)
    for i in range(len(hData)):
      hData[i]["loose"].Print()
      hData[i]["sig"].Print()
    for hist in hsMC:
      hist.Print()

    purityList = RA2b.getZmassFitPlot(doDiMu=doMumu, doDiEl=doEe, dataSet=hData, mcSet=hsMC, doLumi=lumi, iPeriod = period)

    for i in range(len(purityList)):
      print "purity = "+str(purityList[i][0])+" +/- "+str(purityList[i][1])

    for i in range(1,5):
      cIter = jbGroup+i
      canv = ROOT.gROOT.FindObject("canvas"+str(cIter))
      if (type(canv)==ROOT.TCanvas):
        if (i == 2):
          canv.SaveAs(str(plotNameRoot)+"mm.pdf")
        elif (i == 4):
          canv.SaveAs(str(plotNameRoot)+"ee.pdf")
        elif (cIter == 1):
          canv.SaveAs("fitZmass_allloose_mm.pdf")
        elif (cIter == 3):
          canv.SaveAs("fitZmass_allloose_ee.pdf")
    jbGroup += 4

# outFile.Write()

# def getZmassFitPlot(fitFunc=None, dataSet=None, mcSet=None, plotMC=None, doDiMu=None, doDiEl=None, doDiLep=None, getShapeFromLoose=None, nBins=None, distRange=None, nJetBin=None, bJetBin=None, kinBin=None, doVarBinning=None, binning=None, extraCuts=None, dphiCut=None, doLumi=None, do20=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, keepCanvas=None, drawText=None, text=None, textCoords=None, drawText2=None, text2=None, textCoords2=None):
