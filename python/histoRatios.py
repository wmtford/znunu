#!/usr/bin/python
"""
Make various plots of histograms and ratios
"""
import ROOT
import RA2b

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

doMumu = True
doEe = True
fScale = 1
MZmmMax = 0
MZeeMax = 2000
# MZmmMax = 700
# MZeeMax = 600
ratioMin = None
ratioMax = None

hists = []

muhists = {}
ehists = {}

Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
Nfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
if (doMumu):
  muhists['N'] = ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm", "hCutFlow_zmm", "hCuts_zmm", "hVertices_zmm"]
  muhists['D'] = ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hCutFlow_dymm", "hCuts_dymm", "hVertices_dymm"]
  hists.append(muhists)
if (doEe):
  ehists['N'] = ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee", "hCutFlow_zee", "hCuts_zee", "hVertices_zee"]
  ehists['D'] = ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hCutFlow_dyee", "hCuts_dyee", "hVertices_dyee"]
  hists.append(ehists)
legList = ['2016 DY data', '2016 DY MC']

# Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12.root')
# Nfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
# ratioMin = 0.85
# ratioMax = 1.15
# if (doMumu):
#   hists.append(["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hCutFlow_dymm", "hCuts_dymm"])
# if (doEe):
#   hists.append(["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hCutFlow_dyee", "hCuts_dyee"])
# legList = ['2016 DY MC, pileup wt','2016 DY MC']

# Dfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
# if (doMumu):
#   Nfile = ROOT.TFile('../outputs/histsDYmm_2018v15.root')
#   hists.append("hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm"])
#   fScale = 14.0/35.9
# if (doEe):
#   Nfile = ROOT.TFile('../outputs/histsDYee_2018v15.root')
#   hists.append(["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee"])
#   fScale = 13.5/35.9
# legList = ['2018 data','2016 data, scaled']

# zhists = {}
# Dfile = ROOT.TFile('~/cms/src/root/zinvData_2017Feb24/ZinvMCttzMC174binV12.root')
# Nfile = ROOT.TFile('~/cms/src/root/zinvData_2017Feb24/ZinvHistos.root')
# zhists['D'] = ["hCCzinv"]
# zhists['N'] = ["ZinvBGpred"]
# hists.append(zhists)

for NDpairs in hists:
  for i in range(len(NDpairs['N'])):
    nhName = NDpairs['N'][i]
    print str(nhName)
    hNumer = Nfile.Get(nhName)
    hNumer.SetName(str(nhName)+"N")
    dhName = NDpairs['D'][i]
    print str(dhName)
    hDen   = Dfile.Get(dhName)
    hDen.SetName(str(dhName)+"D")
    hDen.Scale(fScale)
    if ("hZmass" in nhName):
      doLogy = False
      if ("mm" in nhName and MZmmMax != 0):
        hNumer.SetMaximum(MZmmMax)
      elif ("ee" in nhName and MZeeMax != 0):
        hNumer.SetMaximum(MZeeMax)
    else:
      doLogy = True
    canvas = RA2b.getPlotAndRatio(
      numHists=hNumer, denomHists=hDen, doRatio=True,
      doLogy=doLogy, doCMSlumi=True, iPeriod=8, drawHorizontalLine=True,
      xTitle=hNumer.GetXaxis().GetTitle(), yTitle=hNumer.GetYaxis().GetTitle(),
      ratioMin=ratioMin, ratioMax=ratioMax,
      legList = legList
      )
    # For 174-bin plot
    # canvas = RA2b.getPlotAndRatio(
    #   numHists=hNumer, denomHists=hDen, doRatio=True,
    #   doLogy=doLogy, doCMSlumi=True, iPeriod=8, drawHorizontalLine=True,
    #   ratioMin=ratioMin, ratioMax=ratioMax, ratioTitle="#frac{Direct}{Prediction} ",
    #   # legList = legList
    #   # doClosureStyle=True
    #   errorBandFillStyle=3144,
    #   # errorBandColor = ROOT.kRed-10
    #   errorBandColor = 632-10,
    #   ratioTitle = "#frac{Expectation}{Prediction} ",
    #   extraText = "Preliminary",
    #   xTitle = "Search region bin number",
    #   yTitle = "Events",
    #   #legCoords = [.65, .95, .54, .79],
    #   legHeader = "Z #rightarrow #nu#bar{#nu} background",
    #   ratioMax=2.15,
    #   ratioMin=0.001,
    #   nDivRatio=505,
    #   markerSize=1.1,
    #   legList = ['Expectation from simulation','Prediction from data']
    #   )
    canvas.SaveAs(str(nhName)+".pdf")

# def getPlotAndRatio(numHists, denomHists=None, bottomPlots=None, doStack=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, ratioTitle=None, ratioMin=None, ratioMax=None, doLogy=None, doFlip=None, doDiff=None, doPull=None, makeLeg=None, legList=None, legCoords=None, textCoords=None, canvasSize=None, canvasName=None, numColors=None, denomColor=None, numMarkers=None, denomMarker=None, markerSize=None, lineWidth=None, numDrawStyles=None, denomDrawStyle=None, drawErrorBand=None, stackColors=None, axisTitleSize=None, drawVerticalLines=None, drawHorizontalLine=None, statBox=None, drawText=None, text=None, setMax=None, setMin=None, doClosureStyle=None,errorBandColor=None,errorBandFillStyle=None,legHeader=None,nDivRatio=None,doNumFill=None, hLineVal=None, hLineColors=None,nDivX=None,ratioGridx=None,ratioGridy=None,topGridx=None,topGridy=None,doRatio=None,numFillStyles=None,numFillColors=None)

