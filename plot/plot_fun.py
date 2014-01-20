import ROOT

infilenames_eta = {
#    "Pt1": "trackValHistogramsPt1.root",
    "Pt10": "trackValHistogramsPt10.root",
    "Pt100": "trackValHistogramsPt100.root"
}

infilenames_eta_gsf = {
#    "Pt1": "trackValHistogramsPt1.root",
    "Pt10": "trackValHistogramsPt10GSF.root",
    "Pt100": "trackValHistogramsPt100GSF.root"
}

infilenames_pt = {
    "FlatPt": "trackValHistogramsFlatPt.root"
}

infilenames_pt_gsf = {
    "FlatPt": "trackValHistogramsFlatPtGSF.root"
}

def load_input_files(indir, infilenames):
    """
    load input files. give dir path and dictionary of filenames as input
    """
    infiles = {}

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

def draw_efficiency_histograms(hists, xtitle = "none", ytitle = "none", ymax =  1., region="none"):
    """
    plot a list of efficiency histogras
    """
    n = 0
    for hist in hists:
        hist.SetLineColor(colors[n])
        if n==0:
            hist.SetMaximum(ymax)
            hist.SetMinimum(0.)
            hist.GetXaxis().SetTitleOffset(1.3)
            hist.GetYaxis().SetTitleOffset(1.4)
            hist.SetTitle("blabla")

            if( region[:3]=="Pt1"):
                hist.GetXaxis().SetTitle("#eta")
            elif( region[:6]=="FlatPt"):
                hist.GetXaxis().SetTitle("p_{T}")

                hist.SetAxisRange(1., 200, 'x')
            if not xtitle == "none":
                hist.GetXaxis().SetTitle(xtitle)
            if not ytitle == "none":
                hist.GetYaxis().SetTitle(ytitle)
            hist.Draw()
        else:
            hist.Draw("same")
        n = n + 1

def draw_legend(hists, pos = "down_right"):
    """
    hist - dictionary of process names and histograms
    """
    if pos == "down_right":
        leg = ROOT.TLegend(0.7,0.5,0.89,0.29);
    
    if pos == "up_right":
        leg = ROOT.TLegend(0.7,0.7,0.89,0.89);

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

def draw_and_save_eff(hists, var, eff_fake, is_gsf, leg_pos = "up_right", title = ""):
    """
    hists - dictionary of process names and histograms
    var - xaxis variable
    eff_fake - "eff", "eff_seed", "fake"
    """

    c = ROOT.TCanvas("c","c")
    c.SetGrid()

    xtitle = ""
    if var == "pt":
        c.SetLogx()
        xtitle = "p_{T}"
    if var == "eta":
        xtitle = "#eta"
    if var == "nrhits":
        xtitle = "Number of sim. hits"
        
    ytitle = ""
    ymax = 1
    if eff_fake == "eff":
        ytitle = "Efficiency"
    if eff_fake == "eff_seed":
        ytitle = "Seeding efficiency"
    if eff_fake == "eff_wrt_seed":
        ytitle = "Reco wrt seeding efficiency"
    
    if eff_fake[:4] == "fake":
        ytitle = "fake rate"
        ymax = 0.3

    draw_efficiency_histograms(hists.values(), xtitle, ytitle, ymax)
    leg = draw_legend(hists, pos = leg_pos)
    leg.Draw()

    GSFstr = ""
    if(is_gsf):
        GSFstr = "_GSF"

    c.SaveAs("$WORKING_DIR/plot/out_plots/" + eff_fake + "_" + var + GSFstr + ".pdf")
    c.SaveAs("$WORKING_DIR/plot/out_plots/" + eff_fake + "_" + var + GSFstr + ".png")
    c.Close()
