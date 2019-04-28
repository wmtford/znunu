#!/usr/bin/python
"""Run this script to set the photon related efficiencies in
plots/histograms/effHists.root
All efficiency factors are hard coded below"""

from array import array
import ROOT

def writeEffHists():
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
    FdirldpFile = ROOT.TFile("../plots/histograms/fragmentation_ldp.root", "READ")
    h_Fdirldpnominal = FdirldpFile.Get("bin46_NJets8910")
    h_Fdirldpnominal.SetName("h_bin46_NJets8910_ldpnominal")
    h_Fdirldp = FdirldpFile.Get("bin59_NJets8910")
    h_Fdirldp.SetName("h_bin59_NJets8910_ldp")

    ## photon purity (Andrew provides these)
    EB = {}
    EE = {}
    EB_ldp = {}
    EE_ldp = {}
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
    pur_eb_ldp_2016 = []
    pur_ec_ldp_2016 = []
    EB['avg'] = [0.9479783904003333, 0.9612409354666666, 0.9657104593263334, 0.9819289765320001]
    EB['avgErr'] = [0.031016399071333245, 0.01721366277833347, 0.015233616496666658, 0.008033783562999885]
    EE['avg'] = [0.898884215471, 0.9232578112983333, 0.9296478600020001, 0.9149493319753333]
    EE['avgErr'] = [0.04399141754999991, 0.04426306202466668, 0.02750138579299999, 0.047701439598333395]
    # Tribeni's April 26, 2910 email:
    EB_ldp['avg'] =  [0.9129327697090001, 0.8734955341873333, 0.8667444292550001, 0.909717842707]
    EB_ldp['avgErr'] = [0.059203833219000024, 0.10124510427333333, 0.08421265228800001, 0.057047576944999956]
    EE_ldp['avg'] =  [0.8120440475803333, 0.7330674656063333, 0.7199778339983333, 0.6674219311026667]
    EE_ldp['avgErr'] = [0.06759096441266665, 0.09573324213933332, 0.07272030768166671, 0.08260617818033333]
    for i in range(0,4):
        pur_eb_2016.append( (EB['avg'][i], EB['avgErr'][i]) )
        pur_ec_2016.append( (EE['avg'][i], EE['avgErr'][i]) )
        pur_eb_ldp_2016.append( (EB_ldp['avg'][i], EB_ldp['avgErr'][i]) )
        pur_ec_ldp_2016.append( (EE_ldp['avg'][i], EE_ldp['avgErr'][i]) )

    # 2017:
    pur_eb_2017 = []
    pur_ec_2017 = []
    pur_eb_ldp_2017 = []
    pur_ec_ldp_2017 = []
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
    # Tribeni's April 26, 2910 email:
    EB['avg'] = [0.867368480889, 0.8540455754413333, 0.8483960770856666, 0.9045874008393334]
    EB['avgErr'] = [0.056006291401, 0.06455268990933327, 0.05125898994866651, 0.03271043206433344]
    EE['avg'] = [0.8886791026499999, 0.897029311727, 0.8028611527306667, 0.7788664830493334]
    EE['avgErr']  = [0.045075913452000016, 0.02533815284499996, 0.05005716822833328, 0.05555676379266661]
    for i in range(0,4):
        pur_eb_2017.append( (EB['avg'][i], EB['avgErr'][i]) )
        pur_ec_2017.append( (EE['avg'][i], EE['avgErr'][i]) )
        pur_eb_ldp_2017.append( (EB_ldp['avg'][i], EB_ldp['avgErr'][i]) )
        pur_ec_ldp_2017.append( (EE_ldp['avg'][i], EE_ldp['avgErr'][i]) )

    # 2018:
    pur_eb_2018 = []
    pur_ec_2018 = []
    pur_eb_ldp_2018 = []
    pur_ec_ldp_2018 = []
    # Andrew's Jan. 27, 2019 email
    # EB['avg'] [0.934561081295, 0.935941864678, 0.9528332801676668, 0.9698275681449999]
    # EB['avgErr'] [0.010881257075999962, 0.015425372193999976, 0.010485711729666725, 0.0064986703679998925]
    # EE['avg'] [0.955755605372, 0.9603043977869999, 0.948359821654, 0.8919489522903333]
    # EE['avgErr'] [0.013496460590999959, 0.016012506070000043, 0.016055838115, 0.032147232375666634]
    # Tribeni 19 March, 2019 e-mail:  photon purity (Gjets_0p4)
    # EB['avg'] = [0.9280759956736667, 0.9285315544210001, 0.9457477890636666, 0.9691437937430001]
    # EB['avgErr'] = [0.011992508267333335, 0.016934029379000126, 0.011962957321666634, 0.006785803926000034]
    # EE['avg'] = [0.9405511421239999, 0.9439760420550001, 0.9345166398170001, 0.9283238776746666]
    # EE['avgErr'] = [0.017382353302999975, 0.022191131761999938, 0.019755000849999926, 0.021305790389333334]
    # Tribeni 6 April, 2019 e-mail:  photon purity (Gjets_0p4, 2018 MC)
    EB['avg'] = [0.9300185028443333, 0.9343408542576667, 0.9492861356343334, 0.9718661961480001]
    EB['avgErr'] = [0.017377262343666655, 0.01575348880933336, 0.011945768857666628, 0.006624156997999897]
    EE['avg'] = [0.942969000792, 0.9498724011163334, 0.939690691016, 0.933428825701]
    EE['avgErr'] = [0.02342499170900003, 0.022395310887333375, 0.024629399997999957, 0.026930692982999993]
    # Tribeni's April 26, 2910 email:
    EB['avg'] = [0.8675806532079999, 0.858558924094, 0.8755681166220001, 0.9187331436863334]
    EB['avgErr'] = [0.0605641676619999, 0.04332221426000005, 0.054710200390000074, 0.035716791730333375]
    EE['avg'] = [0.9055162794053334, 0.876280041638, 0.809691137856, 0.7472074954379999]
    EE['avgErr'] = [0.026188641671666568, 0.032801027306999964, 0.05357446361700002, 0.07041979225200012]
    for i in range(0,4):
        pur_eb_2018.append( (EB['avg'][i], EB['avgErr'][i]) )
        pur_ec_2018.append( (EE['avg'][i], EE['avgErr'][i]) )
        pur_eb_ldp_2018.append( (EB_ldp['avg'][i], EB_ldp['avgErr'][i]) )
        pur_ec_ldp_2018.append( (EE_ldp['avg'][i], EE_ldp['avgErr'][i]) )

    ########## get the efficiency file ################
    # effFile = ROOT.TFile("../plots/histograms/effHists.root", "UPDATE")
    effFile = ROOT.TFile("effHists.root", "UPDATE")  # wtf
    # effFile.Delete("h_pur_eb;1")
    # effFile.Delete("h_pur_ec;1")


    ###################################################
    # H_SF_g1 = effFile.Get("h_SF_g1")
    # h_SF_g1.SetBinContent(1,SF[0])
    # h_SF_g1.SetBinError(1,SF[1])
    # h_SF_g1.Write(h_SF_g1.GetName(),2)

    # h_frag1 = effFile.Get("h_frag1")
    # h_frag1.SetBinContent(1,frag[0])
    # h_frag1.SetBinError(1,frag[1])
    # h_frag1.Write(h_frag1.GetName(),2)

    ########## update the efficiency file #############

    # f_trig_eb_2016.Write(f_trig_eb_2016.GetName(),2)
    # f_trig_ec_2016.Write(f_trig_ec_2016.GetName(),2)
    # f_trig_eb_2017.Write(f_trig_eb_2017.GetName(),2)
    # f_trig_ec_2017.Write(f_trig_ec_2017.GetName(),2)
    # f_trig_eb_2018.Write(f_trig_eb_2018.GetName(),2)
    # f_trig_ec_2018.Write(f_trig_ec_2018.GetName(),2)

    # h_SF_g.Write(h_SF_g.GetName(), 2)
    h_Fdir.Write(h_Fdir.GetName(), 2)
    h_Fdirldpnominal.Write(h_Fdirldpnominal.GetName(), 2)
    h_Fdirldp.Write(h_Fdirldp.GetName(), 2)

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

    def writeHist(ebsample, runyear, content):
        hName = "h_pur_"+ebsample+"_"+runyear
        hist = effFile.Get(hName)
        if (not hist):
            hist = ROOT.TH1F(hName, "gJets purities, "+ebsample, len(pur_bins)-1, pur_bins)
        else:
            hist.SetBins(len(pur_bins)-1, pur_bins)
        hist.GetXaxis().SetTitle(pur_title)
        for i in range(len(content)):
            hist.SetBinContent(i+1, content[i][0])
            hist.SetBinError(i+1, content[i][1])
        hist.Write(hist.GetName(),2)

    writeHist("eb", "2016", pur_eb_2016)
    writeHist("ec", "2016", pur_ec_2016)
    writeHist("eb", "2017", pur_eb_2017)
    writeHist("ec", "2017", pur_ec_2017)
    writeHist("eb", "2018", pur_eb_2018)
    writeHist("ec", "2018", pur_ec_2018)
    writeHist("eb_ldp", "2016", pur_eb_ldp_2016)
    writeHist("ec_ldp", "2016", pur_ec_ldp_2016)
    writeHist("eb_ldp", "2017", pur_eb_ldp_2017)
    writeHist("ec_ldp", "2017", pur_ec_ldp_2017)
    writeHist("eb_ldp", "2018", pur_eb_ldp_2018)
    writeHist("ec_ldp", "2018", pur_ec_ldp_2018)


    effFile.Close()

def main():
    ROOT.gROOT.Reset()
    ROOT.gROOT.SetBatch(1)

    writeEffHists()
  
if __name__ == "__main__":
    main()
