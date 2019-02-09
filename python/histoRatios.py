#!/usr/bin/python
"""
Make various plots of histograms and ratios
"""
import ROOT
import RA2b

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(1)

singleOutFile = True
doDatavMC = True
doDatavData = False
doMCvMC = False
runYear = "2017"
doMumu = True
doEe = True
doLl = True
doPhoton = False
fScaleM = 1
fScaleE = 1
fScaleL = 1
fScaleP = 1
MZmmMax = 0
MZeeMax = 0
setMin = None
setMax = None
ratioMin = [None, None, None, None]
ratioMax = [None, None, None, None]

legList = []
# legCoordsDefault = [0.6,0.6,0.89,0.89]
legCoordsDefault = [0.68,0.6,0.89,0.89]
legCoords = None
drawText = False
textCoordsDefault = [0.35,0.75,.55,.88]
textCoords = None

def histNamesBuild(namesRootList, suffix):
  histNames = []
  for nameRoot in namesRootList:
    histNames.append(nameRoot+"_"+suffix)
  return histNames

reaction = -1

hists = []
filehistsM = {}
filehistsE = {}
filehistsL = {}
filehistsP = {}

#  ========================================================================================

if (doDatavMC):
  # Compare data with MC
  varNames = [['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hMuonEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hElectronEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hPhotonPt', 'hPhotonEta', 'hVertices', 'hCC', 'hCCjk']]
  sampleSuffixN = ['zmm', 'zee', 'zll', 'photon']
  sampleSuffixD = [['ttmm', 'VVmm', 'ttzmm', 'dymm'],
                   ['ttee', 'VVee', 'ttzee', 'dyee'],
                   ['ttll', 'VVll', 'ttzll', 'dyll'],
                   ['gjetsqcd', 'gjets']]
  legendsD = [['        ttbar', '        VV', '        ttZ', '        DY'],
              ['        ttbar', '        VV', '        ttZ', '        DY'],
              ['        ttbar', '        VV', '        ttZ', '        DY'],
              ['QCD', '#gamma+jets']] 
  if (runYear is "2016"):
    iPeriod = 5
    NfileZll = ROOT.TFile('../outputs/histsDY_2016v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2016v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDYMC_2016v16.root')
    DfilePhoton = ROOT.TFile('../outputs/histsGjets_2016v16.root')
    legendsN = ['Z#mu#mu'+' 2016 data', 'Zee 2016 data', 'Zll 2016 data', 'photon 2016 data']
  elif (runYear is "2017"):
    iPeriod = 6
    NfileZll = ROOT.TFile('../outputs/histsDY_2017v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16.root')
    DfilePhoton = ROOT.TFile('../outputs/histsGjets_2017v16.root')
    legendsN = ['Z#mu#mu'' 2017 data', 'Zee 2017 data', 'Zll 2017 data', 'photon 2017 data']
  elif (runYear is "Run2"):
    iPeriod = 9
    NfileZll = ROOT.TFile('../outputs/histsDY_Run2v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_Run2v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDYMC_Run2v16.root')
    DfilePhoton = ROOT.TFile('../outputs/histsGjets_Run2v16.root')
    legendsN = ['Z#mu#mu'+' Run2 data', 'Zee Run2 data', 'Zll Run2 data', 'photon Run2 data']

elif (doDatavData):
  # Data vs data
  varNames = [['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hMuonEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hElectronEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hPhotonPt', 'hPhotonEta', 'hVertices', 'hCC']]
  sampleSuffixN = ['zmm', 'zee', 'zll', 'photon']
  sampleSuffixD = [['zmm'], ['zee'], ['zll'], ['photon']]
  # sampleSuffixD = sampleSuffixN
  if (runYear is "2017"):
    iPeriod = 6
    NfileZll = ROOT.TFile('../outputs/histsDY_2017v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDY_2016v16.root')
    DfilePhoton = ROOT.TFile('../outputs/histsPhoton_2016v16.root')
    legendsN = ['Z#mu#mu 2017 data', 'Zee 2017 data', 'Zll 2017 data', 'photon 2017 data']
    legendsD = [[' 2016 scaled'], [' 2016 scaled'], [' 2016 scaled'], [' 2016 scaled']]
    fScaleM = 41.5/35.9
  elif (runYear is "2018"):
    iPeriod = 8
    NfileZll = ROOT.TFile('../outputs/histsDY_2018v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2018v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDY_2017v16.root')
    DfilePhoton = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
    legendsN = ['Z#mu#mu 2018 data', 'Zee 2018 data', 'Zll 2018 data', 'photon 2018 data']
    legendsD = [[' 2017 scaled'], [' 2017 scaled'], [' 2017 scaled'], [' 2017 scaled']]
    fScaleM = 59.4/41.5

    fScaleE = fScaleM
    fScaleL = fScaleM
    fScaleP = fScaleM

elif (doMCvMC):
  # Compare MC
  # varNames = [["hHT", "hMHT", "hNJets", "hBTags", "hCC", 'hCCjk', 'hCCjb'],
  #             ["hHT", "hMHT", "hNJets", "hBTags", "hCC", 'hCCjk', 'hCCjb'],
  #             ['hCC', 'hCCjk', 'hCCjb'],
  #             ["hHT", "hMHT", "hNJets", "hBTags", "hCC", 'hCCjk', 'hCCjb']]
  varNames = [['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hMuonEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hElectronEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hPhotonPt', 'hPhotonEta', 'hVertices', 'hCC', 'hCCjk']]
  iPeriod = 6

  # NfileZll = ROOT.TFile('../outputs/histsZjets_2017v16.root')
  # NfilePhoton = ROOT.TFile('../outputs/histsZjets_2017v16.root')
  # sampleSuffixN = ['zinv', 'zinv', 'zinv', 'zinv']
  # legendsN = ['Z#nu#nu 2017 MC', 'Z#nu#nu 2017 MC', 'Z#nu#nu 2017 MC','Z#nu#nu 2017 MC']
  NfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16.root')
  NfilePhoton = ROOT.TFile('../outputs/histsGjets_2017v16.root')
  sampleSuffixN = ["dymm", "dyee", 'dyll', "gjets"]
  legendsN = ['DY#mu#mu 2017 MC', 'DYee 2017 MC', 'DYll 2017 MC', '#gamma+jets incl. 2017']

  DfileZll = ROOT.TFile('../outputs/histsDYMC_2016v16.root')
  DfilePhoton = ROOT.TFile('../outputs/histsGjets_2016v16.root')
  sampleSuffixD = [["dymm"], ["dyee"], ['dyll'], ["gjets"]]
  # legendsD = [['DY#mu#mu 2018 MC'], ['DYee 2018 MC'], ['DYll 2018 MC'], ['#gamma+jets incl. 2018']] 
  # ratioMin = [6, 6, 6, 0]
  # ratioMax = [13, 13, 13, 1.1]
  # sampleSuffixD = sampleSuffixN
  legendsD = [['        2016'], ['        2016'], ['        2016'], ['        2016']]

  fScaleM = 41.5/35.9
  # fScaleM = 41.5/59.4
  fScaleE = fScaleM
  fScaleL = fScaleM
  fScaleP = fScaleM
  # ratioMin = [0.8, 0.8, 0.8, 0]
  # ratioMax = [1.2, 1.2, 1.2, 2]

# =============================================================================

if (doMumu):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[0], sampleSuffixN[0])
  legend.append(legendsN[0])
  histnamesD = []
  for hnd in sampleSuffixD[0]:
    histnamesD.append(histNamesBuild(varNames[0], hnd))
  for ld in legendsD[0]:
    legend.append(ld)
  filehistsM['N'] = (NfileZll, histnamesN)
  filehistsM['D'] = (DfileZll, fScaleM, histnamesD)
  hists.append(filehistsM)
  legList.append(legend)
if (doEe):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[1], sampleSuffixN[1])
  legend.append(legendsN[1])
  histnamesD = []
  for hnd in sampleSuffixD[1]:
    histnamesD.append(histNamesBuild(varNames[1], hnd))
  for ld in legendsD[1]:
    legend.append(ld)
  filehistsE['N'] = (NfileZll, histnamesN)
  filehistsE['D'] = (DfileZll, fScaleE, histnamesD)
  hists.append(filehistsE)
  legList.append(legend)
if (doLl):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[2], sampleSuffixN[2])
  legend.append(legendsN[2])
  histnamesD = []
  for hnd in sampleSuffixD[2]:
    histnamesD.append(histNamesBuild(varNames[2], hnd))
  for ld in legendsD[2]:
    legend.append(ld)
  filehistsL['N'] = (NfileZll, histnamesN)
  filehistsL['D'] = (DfileZll, fScaleL, histnamesD)
  hists.append(filehistsL)
  legList.append(legend)
if (doPhoton):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[3], sampleSuffixN[3])
  legend.append(legendsN[3])
  histnamesD = []
  for hnd in sampleSuffixD[3]:
    histnamesD.append(histNamesBuild(varNames[3], hnd))
  for ld in legendsD[3]:
    legend.append(ld)
  filehistsP['N'] = (NfilePhoton, histnamesN)
  filehistsP['D'] = (DfilePhoton, fScaleP, histnamesD)
  hists.append(filehistsP)
  legList.append(legend)

#  ========================================================================================

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
    print "nhName = "+str(nhName)
    hNumer = Nfile.Get(nhName)
    hNumer.SetName(str(nhName)+"N")
    hDenList = []
    for c in range(len(samples['D'][2])):  # Denominator may have multiple components
      dhName = samples['D'][2][c][i]  # scale is stored in samples['D'][1]
      print "dhName = "+str(dhName)
      hDen = Dfile.Get(dhName)
      hDen.SetName(str(dhName)+"D")
      hDen.Scale(samples['D'][1])
      hDenList.append(hDen)
    if ("hZmass" in nhName):
      doLogy = False
      if ("mm" in nhName and MZmmMax != 0):
        hNumer.SetMaximum(MZmmMax)
      elif ("ee" in nhName and MZeeMax != 0):
        hNumer.SetMaximum(MZeeMax)
    else:
      doLogy = True
    if ("Eta" in nhName and not "Photon" in nhName):
      setMin = 1
      # legCoords = legCoordsDefault
      # textCoords = textCoordsDefault
      # legCoords = [0.6,0.6,0.89,0.89]
      # drawText = True
      # textCoords = [0.35,0.75,.55,.88]
    else:
      setMin = None
      # drawText = True
      # textCoords = textCoordsDefault
      # legCoords = legCoordsDefault
    # print "legCoords = "+str(legCoords[0])+", "+str(legCoords[1])+", "+str(legCoords[2])+", "+str(legCoords[3])
    canvTuple = RA2b.getPlotAndRatio(
      numHists=hNumer, denomHists=hDenList, doRatio=True,
      doLogy=doLogy, doCMSlumi=True, iPeriod=iPeriod, drawHorizontalLine=True,
      xTitle=hNumer.GetXaxis().GetTitle(), yTitle=hNumer.GetYaxis().GetTitle(),
      ratioMin=ratioMin[reaction], ratioMax=ratioMax[reaction], setMin=setMin, setMax=setMax,
      legList = legList[reaction],
      # legCoords = legCoords,
      drawText = drawText, textCoords = textCoords
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

# =================================================================================


# def getPlotAndRatio(numHists, denomHists=None, bottomPlots=None, doStack=None, Title=None, xTitle=None, yTitle=None, doCMSlumi=None, iPos=None, iPeriod=None, extraText=None, ratioTitle=None, ratioMin=None, ratioMax=None, doLogy=None, doFlip=None, doDiff=None, doPull=None, makeLeg=None, legList=None, legCoords=None, textCoords=None, canvasSize=None, canvasName=None, numColors=None, denomColor=None, numMarkers=None, denomMarker=None, markerSize=None, lineWidth=None, numDrawStyles=None, denomDrawStyle=None, drawErrorBand=None, stackColors=None, axisTitleSize=None, drawVerticalLines=None, drawHorizontalLine=None, statBox=None, drawText=None, text=None, setMax=None, setMin=None, doClosureStyle=None,errorBandColor=None,errorBandFillStyle=None,legHeader=None,nDivRatio=None,doNumFill=None, hLineVal=None, hLineColors=None,nDivX=None,ratioGridx=None,ratioGridy=None,topGridx=None,topGridy=None,doRatio=None,numFillStyles=None,numFillColors=None)

