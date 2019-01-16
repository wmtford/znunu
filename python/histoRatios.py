#!/usr/bin/python
"""
Make various plots of histograms and ratios
"""
import ROOT
import RA2b

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

singleOutFile = True
doMumu = True
doEe = True
doPhoton = True
norm2016 = False
fScaleM = 1
fScaleE = 1
fScaleP = 1
MZmmMax = 0
MZeeMax = 0
ratioMin = None
ratioMax = None

legList = []
reaction = -1

hists = []
filehistsM = {}
filehistsE = {}
filehistsP = {}

# iPeriod = 5
# Nfile = ROOT.TFile('../outputs/histsDY_2016v12.root')
# Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
# if (doMumu):
  # reaction += 1
  # legList.append(['Z#mu#mu 2016 data','        2016 DY MC'])
#   filehistsM['N'] = (Nfile, ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm", "hCutFlow_zmm", "hCuts_zmm", "hVertices_zmm"])
#   filehistsM['D'] = (Dfile, ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hCutFlow_dymm", "hCuts_dymm", "hVertices_dymm"])
#   hists.append(filehistsM)
# if (doEe):
  # reaction += 1
  # legList.append(['Zee 2016 data','       2016 DY MC'])
#   filehistsE['N'] = (Nfile, ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee", "hCutFlow_zee", "hCuts_zee", "hVertices_zee"])
#   filehistsE['D'] = (Dfile, ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hCutFlow_dyee", "hCuts_dyee", "hVertices_dyee"])
#   hists.append(filehistsE)

#  ========================================================================================

# iPeriod = 5
# Nfile = ROOT.TFile('../outputs/histsDYMC_2016v12_puWt.root')
# Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12.root')
# ratioMin = 0.85
# ratioMax = 1.15
# if (doMumu):
  # reaction += 1
  # legList.append(['Z#mu#mu 2016 MC, pileup wt','        2016 DY MC'])
#   histnames = ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hCutFlow_dymm", "hCuts_dymm"]
#   filehistsM['N'] = (Nfile, histnames)
#   filehistsM['D'] = (Dfile, fScaleM, histnames)
#   hists.append(filehistsM)
# if (doEe):
  # reaction += 1
  # legList.append(['Zee 2016 MC, pileup wt','       2016 DY MC'])
#   histnames = ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hCutFlow_dyee", "hCuts_dyee"]
#   filehistsE['N'] = (Nfile, histnames)
#   filehistsE['D'] = (Dfile, fScaleE, histnames)
#   hists.append(filehistsE)

#  ========================================================================================

# iPeriod = 5
# Nfile = ROOT.TFile('../outputs/histsDY_2016v12_DR0b.root')
# # Nfile = ROOT.TFile('../outputs/histsDY_2016v15_DeepCSV.root')
# # Dfile = ROOT.TFile('../outputs/histsDY_2016v15.root')
# Dfile = ROOT.TFile('../outputs/histsDYMC_2016v12_DR0b.root')
# if (doMumu):
  # reaction += 1
  # legList.append(['Z#mu#mu 2016 data V12','        2016 MC  V12'])
#   histnamesN = ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm", "hVertices_zmm", "hCC_zmm"]
#   histnamesD = ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hVertices_dymm", "hCC_dymm"]
#   filehistsM['N'] = (Nfile, histnamesN)
#   filehistsM['D'] = (Dfile, fScaleM, histnamesD)
#   hists.append(filehistsM)
# if (doEe):
  # reaction += 1
  # legList.append(['Zee 2016 data V12','       2016 MC  V12'])
#   histnamesN = ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee", "hVertices_dyee", "hCC_zee"]
#   histnamesD = ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hVertices_dyee", "hCC_dyee"]
#   filehistsE['N'] = (Nfile, histnamesN)
#   filehistsE['D'] = (Dfile, fScaleE, histnamesD)
#   hists.append(filehistsE)

#  ========================================================================================

iPeriod = 6
NfileZll = ROOT.TFile('../outputs/histsDY_2017v16.root')
NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
if (norm2016):
  # legList = ['2017 0b data','2016 0b data, scaled']
  DfileZll = ROOT.TFile('../outputs/histsDY_2016v15_DR0b.root')
  DfilePhoton = ROOT.TFile('../outputs/histsPhoton_2016v15_DR0b.root')
  fScaleM = 41.5/35.9
  fScaleE = fScaleM
  fScaleP = fScaleM
  ratioMin = 0.85
  ratioMax = 1.25
else:
  DfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16.root')
  DfilePhoton = ROOT.TFile('../outputs/histsGjets_2017v16.root')
if (doMumu):
  reaction += 1
  legList.append(['Z#mu#mu 2017 data V16','        2017 MC V16'])
  histnamesN = ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm", "hMuonEta_zmm", "hVertices_zmm", "hCC_zmm"]
  # histnamesD = histnamesN
  histnamesD = ["hHT_dymm", "hMHT_dymm", "hNJets_dymm", "hBTags_dymm", "hZmass_dymm", "hZpt_dymm", "hMuonEta_dymm", "hVertices_dymm", "hCC_dymm"]
  # histnames = ["hHT_DR_zmm", "hMHT_DR_zmm", "hNJets_DR_zmm"]
  filehistsM['N'] = (NfileZll, histnamesN)
  filehistsM['D'] = (DfileZll, fScaleM, histnamesD)
  hists.append(filehistsM)
if (doEe):
  reaction += 1
  legList.append(['Zee 2017 data V16','       2017 MC V16'])
  histnamesN = ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee", "hElectronEta_zee", "hVertices_zee", "hCC_zee"]
  # histnamesD = histnamesN
  histnamesD = ["hHT_dyee", "hMHT_dyee", "hNJets_dyee", "hBTags_dyee", "hZmass_dyee", "hZpt_dyee", "hElectronEta_dyee", "hVertices_dyee", "hCC_dyee"]
  # histnames = ["hHT_DR_zee", "hMHT_DR_zee", "hNJets_DR_zee"]
  filehistsE['N'] = (NfileZll, histnamesN)
  filehistsE['D'] = (DfileZll, fScaleE, histnamesD)
  hists.append(filehistsE)
if (doPhoton):
  reaction += 1
  legList.append(['#gamma+jets 2017 data V16','           2017 MC V16'])
  histnamesN = ["hHT_photon", "hMHT_photon", "hNJets_photon", "hBTags_photon", "hPhotonPt_photon", "hPhotonEta_photon", "hVertices_photon", "hCC_photon"]
  histnamesD = ["hHT_gjets", "hMHT_gjets", "hNJets_gjets", "hBTags_gjets", "hPhotonPt_gjets", "hPhotonEta_gjets", "hVertices_gjets", "hCC_gjets"]
  # histnames = ["hHT_DR_photon", "hMHT_DR_photon", "hNJets_DR_photon"]
  filehistsP['N'] = (NfilePhoton, histnamesN)
  filehistsP['D'] = (DfilePhoton, fScaleP, histnamesD)
  hists.append(filehistsP)

#  ========================================================================================

# iPeriod = 8
# Nfile = ROOT.TFile('../outputs/histsDY_2018_14ifb_v15_DeepCSV.root')
# Dfile = ROOT.TFile('../outputs/histsDY_2016v15_DeepCSV.root')
# fScaleM = 14.0/35.9
# fScaleE = 13.5/35.9
# if (doMumu):
  # reaction += 1
  # legList.append(['Z#mu#mu 2018 data','        2016 data scaled'])
#   # Nfile = ROOT.TFile('../outputs/histsDYmm_2018v15.root')
#   histnames = ["hHT_zmm", "hMHT_zmm", "hNJets_zmm", "hBTags_zmm", "hZmass_zmm", "hZpt_zmm", "hVertices_zmm", "hCC_zmm"]
#   filehistsM['N'] = (Nfile, histnames)
#   filehistsM['D'] = (Dfile, fScaleM, histnames)
#   hists.append(filehistsM)
#   MZmmMax = 750
# if (doEe):
  # reaction += 1
  # legList.append(['Zee 2018 data','       2016 data scaled'])
#   # Nfile = ROOT.TFile('../outputs/histsDYee_2018v15.root')
#   histnames = ["hHT_zee", "hMHT_zee", "hNJets_zee", "hBTags_zee", "hZmass_zee", "hZpt_zee", "hVertices_zee", "hCC_zee"]
#   filehistsE['N'] = (Nfile, histnames)
#   filehistsE['D'] = (Dfile, fScaleE, histnames)
#   hists.append(filehistsE)
#   # MZeeMax = 600

# zhists = {}
# Dfile = ROOT.TFile('~/cms/src/root/zinvData_2017Feb24/ZinvMCttzMC174binV12.root')
# Nfile = ROOT.TFile('~/cms/src/root/zinvData_2017Feb24/ZinvHistos.root')
# zhists['D'] = (Dfile, 1, ["hCCzinv"])
# zhists['N'] = (Nfile, ["ZinvBGpred"])
# hists.append(zhists)

canvas = ROOT.TCanvas()
if (singleOutFile):
  canvas.Print("ratioPlots.pdf[")

reaction = -1
for samples in hists:
  reaction += 1
  Nfile = samples['N'][0]
  Dfile = samples['D'][0]
  for i in range(len(samples['N'][1])):
    nhName = samples['N'][1][i]
    print str(nhName)
    hNumer = Nfile.Get(nhName)
    hNumer.SetName(str(nhName)+"N")
    dhName = samples['D'][2][i]  # scale is stored in samples['D'][1]
    print str(dhName)
    hDen   = Dfile.Get(dhName)
    hDen.SetName(str(dhName)+"D")
    hDen.Scale(samples['D'][1])
    if ("hZmass" in nhName):
      doLogy = False
      if ("mm" in nhName and MZmmMax != 0):
        hNumer.SetMaximum(MZmmMax)
      elif ("ee" in nhName and MZeeMax != 0):
        hNumer.SetMaximum(MZeeMax)
    else:
      doLogy = True
    canvTuple = RA2b.getPlotAndRatio(
      numHists=hNumer, denomHists=hDen, doRatio=True,
      doLogy=doLogy, doCMSlumi=True, iPeriod=iPeriod, drawHorizontalLine=True,
      xTitle=hNumer.GetXaxis().GetTitle(), yTitle=hNumer.GetYaxis().GetTitle(),
      ratioMin=ratioMin, ratioMax=ratioMax,
      legList = legList[reaction]
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

    if (singleOutFile):
      canvTuple[0].Print("ratioPlots.pdf")
    else:
      canvTuple[0].SaveAs(str(nhName)+".pdf")

if (singleOutFile):
  canvas.Print("ratioPlots.pdf]")



# def getPlotAndRatio(numHists, denomHists=None, bottomPlots=None, doStack=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, ratioTitle=None, ratioMin=None, ratioMax=None, doLogy=None, doFlip=None, doDiff=None, doPull=None, makeLeg=None, legList=None, legCoords=None, textCoords=None, canvasSize=None, canvasName=None, numColors=None, denomColor=None, numMarkers=None, denomMarker=None, markerSize=None, lineWidth=None, numDrawStyles=None, denomDrawStyle=None, drawErrorBand=None, stackColors=None, axisTitleSize=None, drawVerticalLines=None, drawHorizontalLine=None, statBox=None, drawText=None, text=None, setMax=None, setMin=None, doClosureStyle=None,errorBandColor=None,errorBandFillStyle=None,legHeader=None,nDivRatio=None,doNumFill=None, hLineVal=None, hLineColors=None,nDivX=None,ratioGridx=None,ratioGridy=None,topGridx=None,topGridy=None,doRatio=None,numFillStyles=None,numFillColors=None)

