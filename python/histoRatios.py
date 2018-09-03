#!/usr/bin/python
"""
Make various plots of histograms and ratios
"""
import ROOT
import RA2b

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

doMumu = True
doEe = False
iPeriod = 5
fScale = 1
MZmmMax = 0
MZeeMax = 0
ratioMin = None
ratioMax = None

hists = []
filehists = {}

# legList = ['2016 DY data', '2016 DY MC']
# Nfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
# Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
# if (doMumu):
#   filehists['N'] = (Nfile, ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm", "hCutFlow_zmm", "hCuts_zmm", "hVertices_zmm"])
#   filehists['D'] = (Dfile, ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hCutFlow_dymm", "hCuts_dymm", "hVertices_dymm"])
#   hists.append(filehists)
# if (doEe):
#   filehists['N'] = (Nfile, ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee", "hCutFlow_zee", "hCuts_zee", "hVertices_zee"])
#   filehists['D'] = (Dfile, ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hCutFlow_dyee", "hCuts_dyee", "hVertices_dyee"])
#   hists.append(filehists)

# legList = ['2016 DY MC, pileup wt','2016 DY MC']
# Nfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
# Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12.root')
# ratioMin = 0.85
# ratioMax = 1.15
# if (doMumu):
#   histnames = ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hCutFlow_dymm", "hCuts_dymm"]
#   filehists['N'] = (Nfile, histnames)
#   filehists['D'] = (Dfile, histnames)
#   hists.append(filehists)
# if (doEe):
#   histnames = ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hCutFlow_dyee", "hCuts_dyee"]
#   filehists['N'] = (Nfile, histnames)
#   filehists['D'] = (Dfile, histnames)
#   hists.append(filehists)

iPeriod = 6
legList = ['2017 data','2016 data, scaled']
# legList = ['2017 data, no trigger 50','2017 data']
if (doMumu):
  # Nfile = ROOT.TFile('../outputs/histsDYmm_2017v15_noTrig50.root')
  Nfile = ROOT.TFile('../outputs/histsDYmm_2017v15.root')
  Dfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
  histnames = ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm"]
  filehists['N'] = (Nfile, histnames)
  filehists['D'] = (Dfile, histnames)
  hists.append(filehists)
  fScale = 41.5/35.9
if (doEe):
  Nfile = ROOT.TFile('../outputs/histsDYee_2017v15.root')
  Dfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
  histnames = ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee"]
  filehists['N'] = (Nfile, histnames)
  filehists['D'] = (Dfile, histnames)
  hists.append(filehists)
  fScale = 41.5/35.9

# iPeriod = 8
# legList = ['2018 data','2016 data, scaled']
# if (doMumu):
#   Nfile = ROOT.TFile('../outputs/histsDYmm_2018v15.root')
#   Dfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
#   histnames = ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm"]
#   filehists['N'] = (Nfile, histnames)
#   filehists['D'] = (Dfile, histnames)
#   hists.append(filehists)
#   fScale = 14.0/35.9
#   MZmmMax = 700
# if (doEe):
#   Nfile = ROOT.TFile('../outputs/histsDYee_2018v15.root')
#   Dfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
#   histnames = ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee"]
#   filehists['N'] = (Nfile, histnames)
#   filehists['D'] = (Dfile, histnames)
#   hists.append(filehists)
#   fScale = 13.5/35.9
#   MZeeMax = 600

# zhists = {}
# Dfile = ROOT.TFile('~/cms/src/root/zinvData_2017Feb24/ZinvMCttzMC174binV12.root')
# Nfile = ROOT.TFile('~/cms/src/root/zinvData_2017Feb24/ZinvHistos.root')
# zhists['D'] = (Dfile, ["hCCzinv"])
# zhists['N'] = (Nfile, ["ZinvBGpred"])
# hists.append(zhists)

for samples in hists:
  Nfile = samples['N'][0]
  Dfile = samples['D'][0]
  for i in range(len(samples['N'][1])):
    nhName = samples['N'][1][i]
    print str(nhName)
    hNumer = Nfile.Get(nhName)
    hNumer.SetName(str(nhName)+"N")
    dhName = samples['D'][1][i]
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
      doLogy=doLogy, doCMSlumi=True, iPeriod=iPeriod, drawHorizontalLine=True,
      xTitle=hNumer.GetXaxis().GetTitle(), yTitle=hNumer.GetYaxis().GetTitle(),
      ratioMin=ratioMin, ratioMax=ratioMax,
      legList = legList
      )
    # For 174-bin plot
    # canvas = RA2b.getPlotAndRatio(
    #   numHists=hNumer, denomHists=hDen, doRatio=True,
    #   doLogy=doLogy, doCMSlumi=True, iPeriod=iPeriod, drawHorizontalLine=True,
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

