import ROOT

infilenames_eta = {
#    "Pt1": "trackValHistogramsPt1.root",
#    "Pt10": "trackValHistogramsPt10.root",
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
        hit_eff_hists[region] = [] # at each region save a vector of efficiencies for each hit
        for ihit in range(0, nrhits):
            hist_name = hit_dir + "eff_" + var + "_athit_" + str(ihit) + "_quality0" + str(quality) 
            print "Opening histogram: " + hist_name
            hit_eff_hists[region].append( histfile.Get(hist_name) )


    
    return hit_eff_hists
        
        
#        eff_hist[region]
