import ROOT

infilenames_eta = {
#    "Pt1": "trackValHistogramsPt1.root",
    "Pt10": "trackValHistogramsPt10.root",
    "Pt100": "trackValHistogramsPt100.root"
}

infilenames_pt = {
    "FlatPt": "trackValHistogramsFlatPt.root"
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

def draw_efficiency_histograms(hists, xtitle = "none", region="none"):
    """
    plot a list of efficiency histogras
    """
    n = 0
    for hist in hists:
        hist.SetLineColor(colors[n])
        if n==0:
            hist.SetMaximum(1.)
            hist.SetMinimum(0.)
            if( region[:3]=="Pt1"):
                hist.GetXaxis().SetTitle("#eta")
            elif( region[:6]=="FlatPt"):
                hist.GetXaxis().SetTitle("p_{T}")
                hist.SetAxisRange(1., 200, 'x')
            if not xtitle == "none":
                hist.GetXaxis().SetTitle(xtitle)
            hist.Draw()
        else:
            hist.Draw("same")
        n = n + 1

def draw_legend(hists):
    """
    hist - dictionary of process names and histograms
    """
    leg = ROOT.TLegend(0.4,0.6,0.89,0.89);
    for process,hist in hists.iteritems():
        leg.AddEntry(hist, process)
    leg.Draw()

#        eff_hist[region]
