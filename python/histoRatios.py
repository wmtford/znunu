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
runYear = "Run2"
doMumu = True
doEe = True
doLl = True
doPhoton = True
nScale = [1, 1, 1, 1]
dScale = [1, 1, 1, 1]
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
    # Consider mulitplying DYMC by new 2016 k factor instead of old:  scale = 1.257/1.23
    iPeriod = 5
    NfileZll = ROOT.TFile('../outputs/histsDY_2016v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2016v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDYMC_2016v16_noPU.root')
    DfilePhoton = ROOT.TFile('../outputs/histsGjets_2016v16_DR2016wt_noPU.root')
    legendsN = ['Z#mu#mu'+' 2016 data', 'Zee 2016 data', 'Zll 2016 data', 'photon 2016 data']
  elif (runYear is "2017"):
    # Consider mulitplying DYMC by 2017 k factor instead of 2016:  scale = 1.165/1.23
    iPeriod = 6
    NfileZll = ROOT.TFile('../outputs/histsDY_2017v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
    # DfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16_noZptWt.root')
    DfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16_ZptWt.root')
    DfilePhoton = ROOT.TFile('../outputs/histsGjets_2017v16_DR2017wt.root')
    legendsN = ['Z#mu#mu'' 2017 data', 'Zee 2017 data', 'Zll 2017 data', 'photon 2017 data']
  elif (runYear is "Run2"):
    iPeriod = 9
    NfileZll = ROOT.TFile('../outputs/histsDY_Run2v161617.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_Run2v161617.root')
    DfileZll = ROOT.TFile('../outputs/histsDYMC_Run2v161617.root')
    DfilePhoton = ROOT.TFile('../outputs/histsGjets_Run2v161617_DRr2wt.root')
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
    dScale = [41.5/35.9, 41.5/35.9, 41.5/35.9, 41.5/35.9]
  elif (runYear is "2018"):
    iPeriod = 8
    NfileZll = ROOT.TFile('../outputs/histsDY_2018v16.root')
    NfilePhoton = ROOT.TFile('../outputs/histsPhoton_2018v16.root')
    DfileZll = ROOT.TFile('../outputs/histsDY_2017v16.root')
    DfilePhoton = ROOT.TFile('../outputs/histsPhoton_2017v16.root')
    legendsN = ['Z#mu#mu 2018 data', 'Zee 2018 data', 'Zll 2018 data', 'photon 2018 data']
    legendsD = [[' 2017 scaled'], [' 2017 scaled'], [' 2017 scaled'], [' 2017 scaled']]
    dScale = [59.5/41.5, 59.5/41.5, 59.5/41.5, 59.5/41.5]

elif (doMCvMC):
  # Compare MC
  # varNames = [["hHT", "hMHT", "hNJets", "hBTags", "hCC", 'hCCjk', 'hCCjb'],
  #             ["hHT", "hMHT", "hNJets", "hBTags", "hCC", 'hCCjk', 'hCCjb'],
  #             ['hCC', 'hCCjk', 'hCCjb'],
  #             ["hHT", "hMHT", "hNJets", "hBTags", "hCC", 'hCCjk', 'hCCjb']]
  varNames = [['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hgenZpt', 'hMuonEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hZmass', 'hZpt', 'hgenZpt', 'hElectronEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb'],
              ['hCC', 'hCCjk', 'hCCjb'],
              ['hHT', 'hMHT', 'hNJets', 'hBTags', 'hPhotonPt', 'hPhotonEta', 'hVertices', 'hCC', 'hCCjk', 'hCCjb']]

  if (runYear is "2017"):
    iPeriod = 6
    # NfileZll = ROOT.TFile('../outputs/histsZjets_2017v16.root')
    # NfilePhoton = ROOT.TFile('../outputs/histsZjets_2017v16.root')
    # sampleSuffixN = ['zinv', 'zinv', 'zinv', 'zinv']
    # legendsN = ['Z#nu#nu 2017 MC', 'Z#nu#nu 2017 MC', 'Z#nu#nu 2017 MC','Z#nu#nu 2017 MC']
    # NfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16_noZptWt.root')
    NfileZll = ROOT.TFile('../outputs/histsDYMC_2017v16_HT17wt_ZptWt.root')
    NfilePhoton = ROOT.TFile('../outputs/histsGjets_2017v16.root')
    legendsN = ['DY#mu#mu 2017 MC', 'DYee 2017 MC', 'DYll 2017 MC', '#gamma+jets incl. 2017']
    dScale = [41.5/35.9, 41.5/35.9, 41.5/35.9, 41.5/35.9]
  elif (runYear is "2018"):
    iPeriod = 8
    NfileZll = ROOT.TFile('../outputs/histsDYMC_2018v17.root')
    NfilePhoton = ROOT.TFile('../outputs/histsGjets_2018v17.root')
    legendsN = ['DY#mu#mu 2018 MC', 'DYee 2018 MC', 'DYll 2018 MC', '#gamma+jets incl. 2018']
    dScale = [59.6/35.9, 59.6/35.9, 59.6/35.9, 59.6/35.9]

  sampleSuffixN = ["dymm", "dyee", 'dyll', "gjets"]
  DfileZll = ROOT.TFile('../outputs/histsDYMC_2016v16_noPU.root')
  DfilePhoton = ROOT.TFile('../outputs/histsGjets_2016v16_noPU.root')
  sampleSuffixD = [["dymm"], ["dyee"], ['dyll'], ["gjets"]]
  # legendsD = [['DY#mu#mu 2018 MC'], ['DYee 2018 MC'], ['DYll 2018 MC'], ['#gamma+jets incl. 2018']] 
  # ratioMin = [6, 6, 6, 0]
  # ratioMax = [13, 13, 13, 1.1]
  # sampleSuffixD = sampleSuffixN
  legendsD = [['2016 MC scaled'], ['2016 MC scaled'], ['2016 MC scaled'], ['2016 MC scaled']]
  # ratioMin = [0.8, 0.8, 0.8, 0]
  # ratioMax = [1.2, 1.2, 1.2, 2]

  # R_ZZ, R_Zgamma, Run 2
  # iPeriod = 8
  # NfileZll = ROOT.TFile('../outputs/histsZjets_Run2v16_HT17wt_ZptWt.root')
  # NfilePhoton = ROOT.TFile('../outputs/histsZjets_Run2v16_HT17wt_ZptWt.root')
  # DfileZll = ROOT.TFile('../outputs/histsDYMC_Run2v16_HT17wt_ZptWt_noPU.root')
  # DfilePhoton = ROOT.TFile('../outputs/histsGjets_Run2v16_noPU.root')
  # doMumu = False
  # doEe = False
  # varNames = [['hCC', 'hCCjk', 'hCCjb'], ['hCC', 'hCCjk', 'hCCjb']]
  # sampleSuffixN = ["zinv", "zinv"]
  # legendsN = ['Z+jets incl. Run 2', 'Z+jets incl. Run 2']
  # sampleSuffixD = [['dyll'], ['gjets']]
  # legendsD = [['DYll Run 2 MC'], ['#gamma+jets * k Run 2']]
  # nScale = [1, 1/1.23]
  # ratioMax = [11, 1.1]

# =============================================================================
reaction = -1
if (doMumu):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
  legend.append(legendsN[reaction])
  histnamesD = []
  for hnd in sampleSuffixD[reaction]:
    histnamesD.append(histNamesBuild(varNames[reaction], hnd))
  for ld in legendsD[reaction]:
    legend.append(ld)
  filehistsM['N'] = (NfileZll, nScale[reaction], histnamesN)
  filehistsM['D'] = (DfileZll, dScale[reaction], histnamesD)
  hists.append(filehistsM)
  legList.append(legend)
if (doEe):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
  legend.append(legendsN[reaction])
  histnamesD = []
  for hnd in sampleSuffixD[reaction]:
    histnamesD.append(histNamesBuild(varNames[reaction], hnd))
  for ld in legendsD[reaction]:
    legend.append(ld)
  filehistsE['N'] = (NfileZll, nScale[reaction], histnamesN)
  filehistsE['D'] = (DfileZll, dScale[reaction], histnamesD)
  hists.append(filehistsE)
  legList.append(legend)
if (doLl):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
  legend.append(legendsN[reaction])
  histnamesD = []
  for hnd in sampleSuffixD[reaction]:
    histnamesD.append(histNamesBuild(varNames[reaction], hnd))
  for ld in legendsD[reaction]:
    legend.append(ld)
  filehistsL['N'] = (NfileZll, nScale[reaction], histnamesN)
  filehistsL['D'] = (DfileZll, dScale[reaction], histnamesD)
  hists.append(filehistsL)
  legList.append(legend)
if (doPhoton):
  reaction += 1
  legend = []
  histnamesN = histNamesBuild(varNames[reaction], sampleSuffixN[reaction])
  legend.append(legendsN[reaction])
  histnamesD = []
  for hnd in sampleSuffixD[reaction]:
    histnamesD.append(histNamesBuild(varNames[reaction], hnd))
  for ld in legendsD[reaction]:
    legend.append(ld)
  filehistsP['N'] = (NfilePhoton, nScale[reaction], histnamesN)
  filehistsP['D'] = (DfilePhoton, dScale[reaction], histnamesD)
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
  for i in range(len(samples['N'][2])):
    nhName = samples['N'][2][i]
    print "nhName = "+str(nhName)+" in file "+str(Nfile)
    hNumer = Nfile.Get(nhName)
    hNumer.SetName(str(nhName)+"N")
    hNumer.Scale(samples['N'][1])
    hDenList = []
    for c in range(len(samples['D'][2])):  # Denominator may have multiple components
      dhName = samples['D'][2][c][i]  # scale is stored in samples['D'][1]
      print "dhName = "+str(dhName)+" in file "+str(Dfile)
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
    doClosureStyle = False
    # if ("hCC" in nhName):
    #   doClosureStyle = True
    canvTuple = RA2b.getPlotAndRatio(
      numHists=hNumer, denomHists=hDenList, doRatio=True,
      doLogy=doLogy, doCMSlumi=True, iPeriod=iPeriod, drawHorizontalLine=True,
      xTitle=hNumer.GetXaxis().GetTitle(), yTitle=hNumer.GetYaxis().GetTitle(),
      ratioMin=ratioMin[reaction], ratioMax=ratioMax[reaction], setMin=setMin, setMax=setMax,
      legList = legList[reaction],
      # legCoords = legCoords,
      drawText = drawText, textCoords = textCoords,
      doClosureStyle = doClosureStyle, markerSize=0.8
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

