import ROOT
from collections import OrderedDict as dict

infilenames_eta = {
    "Pt100": "trackValHistograms_trackValTree_Pt100_maxCand_5_maxChi2_2000_nSigma_3.root",
    "Pt10":  "trackValHistograms_trackValTree_Pt10_maxCand_5_maxChi2_2000_nSigma_3.root",
#    "Pt100": "trackValHistogramsPt100.root"
}

infilenames_eta_gsf = {
    "Pt100":   "trackValHistograms_trackValTree_Pt100_maxCand_5_maxChi2_2000_nSigma_3.root",
    "Pt10" :   "trackValHistograms_trackValTree_Pt10_maxCand_5_maxChi2_2000_nSigma_3.root",
#    "Pt10": "trackValHistogramsPt10GSF.root",
#    "Pt100": "trackValHistogramsPt100GSF.root"
}

infilenames_pt = {
#    "FlatPt": "trackValHistogramsFlatPt.root"
    "FlatPt": "trackValHistograms_trackValTree_FlatPt_maxCand_5_maxChi2_100_nSigma_3.root"
}

infilenames_pt_gsf = {
    "FlatPt": "trackValHistograms_trackValTree_FlatPt_maxCand_5_maxChi2_100_nSigma_3.root"
#    "FlatPt": "trackValHistogramsFlatPtGSF.root"
}

def load_input_files(indir, infilenames):
    """
    load input files. give dir path and dictionary of filenames as input
    """
    infiles = dict()

    for region, filename in infilenames.iteritems():
        inpath = indir + filename
        infiles[region] = ROOT.TFile(inpath)

    return infiles

def get_hit_efficiency_hist(infiles, var, quality, nrhits = 10):
    """
    read a dictionary of input root files in different regions ("Pt1", "Pt10", "Pt100", "FlatPt") and return a dictionary of histograms for each hit
    var -- eta or pt
    quality -- 65, 75, 85, 95 - corresponding to the pt efficiency required at each hit
    """

    hit_dir = "VariablesBySimhit/"
    hit_eff_hists = {}

    for region, histfile in infiles.iteritems():
        if region=="FlatPt":
            hit_eff_hists[region + "_barrel"] = []
            hit_eff_hists[region + "_trans"] = []
            hit_eff_hists[region + "_endcap"] = []
        else:
            hit_eff_hists[region] = [] # at each region save a vector of efficiencies for each hit

        for ihit in range(0, nrhits):
            if region=="FlatPt": #separate different eta regions
                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_barrel_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region+"_barrel"].append(histfile.Get(hist_name))

                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_trans_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region+"_trans"].append(histfile.Get(hist_name))

                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_endcap_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region+"_endcap"].append(histfile.Get(hist_name))
            else:
                hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_quality0" + str(quality)
                print "Opening histogram: " + hist_name
                hit_eff_hists[region].append( histfile.Get(hist_name) )



    return hit_eff_hists


colors = [ROOT.kBlack, ROOT.kRed, ROOT.kYellow, ROOT.kYellow-3, ROOT.kGreen, ROOT.kGreen+3, ROOT.kCyan, ROOT.kCyan+1, ROOT.kCyan+2, ROOT.kCyan+3, ROOT.kCyan+4, ROOT.kBlue, ROOT.kViolet-3, ROOT.kViolet+3]

def draw_efficiency_histograms(hists, xtitle = "none", ytitle = "none", ymax =  1., region="none", style="", logy=False):
    """
    hists -- dictionary of histograms
    plot a list of efficiency histogras
    """
    n = 0
    for hist in hists:
        hist.SetLineColor(colors[n])
        hist.SetMarkerColor(colors[n])
        hist.SetMarkerStyle(20);
        hist.SetMarkerSize(1);
        if n==0:
            hist.SetMaximum(ymax)
            hist.SetMinimum(0.)
            if logy: 
                hist.SetMinimum(0.001)
            hist.GetXaxis().SetTitleOffset(1.3)
            hist.GetYaxis().SetTitleOffset(1.4)
            hist.SetTitle("blabla")
            hist.SetStats(False)

            if( region[:3]=="Pt1"):
                hist.GetXaxis().SetTitle("#eta")
            elif( xtitle=="p_{T}"):
                hist.SetAxisRange(1., 200, 'x')

            if not xtitle == "none":
                hist.GetXaxis().SetTitle(xtitle)
            if not ytitle == "none":
                hist.GetYaxis().SetTitle(ytitle)

            print "style = " + style
            if style == "noerr":
                hist.Draw("hist")
            else:
                hist.Draw("ep")
        else:
            if style == "noerr":
                hist.Draw("histsame")
            else:
                hist.Draw("epsame")
        n = n + 1

def draw_legend(hists, pos = "down_right"):
    """
    hist - dictionary of process names and histograms
    """
    if pos == "down_right":
        leg = ROOT.TLegend(0.55,0.5,0.95,0.29);

    if pos == "up_right":
        leg = ROOT.TLegend(0.55, 0.7, 0.95, 0.89);

    if pos == "up_left":
        leg = ROOT.TLegend(0.2, 0.7, 0.45, 0.89);

    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    for process,hist in hists.iteritems():
        leg.AddEntry(hist, process)

    return leg

def add_text_box(text=''):
    """ Add a CMS blurb to a plot """
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.SetTextAlign(11)
    return latex.DrawLatex(0.17, 0.96, text)

def draw_and_save_eff(hists, var, eff_fake, is_gsf, label = "", leg_pos = "up_right", title = "", ymax_res=0.5, style=""):
    """
    hists - dictionary of process names and histograms
    var - xaxis variable
    eff_fake - "eff", "eff_seed", "eff_wrt_seed", "fake"
    label - "extra label"
    """

    c = ROOT.TCanvas("c","c", 600, 600)
    c.SetGrid()

    logy=False
    xtitle = ""
    if var == "pt":
        c.SetLogx()
        xtitle = "p_{T}"
    if var == "eta":
        xtitle = "#eta"
    if var[:6] == "nrhits":
        xtitle = "Number of sim. hits"

    ytitle = ""
    ymax = 1
    if eff_fake == "eff":
        ytitle = "Efficiency "
    if eff_fake == "eff_seed":
        ytitle = "Seeding efficiency"
    if eff_fake == "eff_wrt_seed":
        ytitle = "Reco wrt seeding efficiency"
    if eff_fake[:4] == "fake":
        ytitle = "Fake rate"
        ymax = 0.5
#        if var == "pt":
#            ymax = ymax_res

    if eff_fake[:4] == "pull":
        ytitle = "pull"
    if eff_fake[:4] == "res":
        ytitle = label
        ymax = ymax_res
#        logy=True


    if len(title)>0:
        ytitle=ytitle + " (" + title + ")"

#    if len(label) > 0:
#        ytitle = ytitle + " (" + label + ")"


    draw_efficiency_histograms(hists.values(), xtitle, ytitle, ymax, style=style, logy=logy)
    leg = draw_legend(hists, pos = leg_pos)
    leg.Draw()

    GSFstr = ""
    if(is_gsf):
        GSFstr = "_GSF"
#    if log:
#        c.SetLogy()

#    c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/" + eff_fake + "_" + var + "_" + label + GSFstr + ".pdf")
    c.SaveAs("$WORKING_DIR/plot/out_plots_paramScans/13TeV_011014/" + eff_fake + "_" + var + "_" + label + GSFstr + ".png")
    c.Close()

    
def draw_resolution(res_hist_2d, res_hist_name):
    """
    res_hist_2d -- histogram of eta/pt vs resolution variable
    res_hist_name -- name to be given to the output resolution histogram
    """

    nbinsx = res_hist_2d.GetNbinsX()
    nbinsy = res_hist_2d.GetNbinsY()

    temp_res = ROOT.TH1F("temp_res", "temp_res", nbinsy, res_hist_2d.GetBinLowEdge(1), res_hist_2d.GetBinLowEdge(nbinsx+1) )

    res_1d = ROOT.TH1F(res_hist_name,res_hist_name, nbinsx-1, res_hist_2d.GetXaxis().GetXmin(), res_hist_2d.GetXaxis().GetXmax() )

    for i in range(1, nbinsx+1 ):
        temp_res = res_hist_2d.ProjectionY("proj",i, i+1)
        temp_gaus = ROOT.TF1("temp_gaus","gaus", temp_res.GetMean()-1*temp_res.GetRMS(), temp_res.GetMean()+1*temp_res.GetRMS() )#, -0.05, 0.05)              

        r = temp_res.Fit(temp_gaus, "SRML Q")
        mean = r.Parameter(1)
        sigma = r.Parameter(2)
        
        res_1d.SetBinContent(i, sigma)

    return res_1d

