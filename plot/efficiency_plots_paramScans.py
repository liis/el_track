import sys, os
import ROOT

import tdrstyle
tdrstyle.tdrstyle()

from plot_fun import infilenames_eta, infilenames_pt, infilenames_eta_gsf, infilenames_pt_gsf, load_input_files, draw_efficiency_histograms, draw_legend, draw_and_save_eff

from plot_fun_paramScans import get_infilenames_by_params

indir = "$WORKING_DIR/tree_to_histo/histograms/" #location of input histograms -->clean up!

print "Opening input files from " + indir

is_gsf = 0
#--------- run over all combinations of parameter values in list of input files, have multiple values for only 1 parameter for comparisons ----------

parameters_maxCand = {
    "maxChi2": [2000],
    "nSigma": [3],
    "maxCand": [3, 5, 6, 7],
    } #consider all combinations of parameters

parameters_maxChi2 = {
    "maxChi2": [50, 100, 300, 2000],
    "nSigma": [3],
    "maxCand": [5]
    }

parameters_nSigma = {
    "maxChi2": [2000],
    "nSigma": [2, 3, 4, 5],
    "maxCand": [5],
    }

parameter_sets = [
#    parameters_maxCand, 
    parameters_maxChi2, 
#    parameters_nSigma
    ]

for parameters in parameter_sets:

    infilenames_pt10 = get_infilenames_by_params(parameters, "Pt10")
    infilenames_pt100 = get_infilenames_by_params(parameters, "Pt100")
    infilenames_flatPt = get_infilenames_by_params(parameters, "FlatPt")

    infiles_pt10 = load_input_files(indir, infilenames_pt10)
    infiles_pt100 = load_input_files(indir, infilenames_pt100)
    infiles_flatPt = load_input_files(indir, infilenames_flatPt)

##-------------- Group histograms wrt. Pt regions (Pt10, Pt100)--------------
    from collections import OrderedDict as dict
    fake_eta = {}
    fake_eta["Pt10"] = {}
    fake_eta["Pt100"] = {}
    fake_eta["FlatPt"] = {}

    eff_eta = {}
    eff_eta["Pt10"] = {}
    eff_eta["Pt100"] = {}
    eff_eta["FlatPt"] = {}

    eff_pt = {}
    fake_pt = {}
    eta_regions = ["barrel", "endcap", "trans"]

    for eta_region in eta_regions:
        eff_pt["FlatPt_" + eta_region] = {}
        fake_pt["FlatPt_" + eta_region] = {}
    
    for cutstring in infiles_pt10:
        eff_eta["Pt10"][cutstring] = infiles_pt10[cutstring].Get("eff_eta")
        fake_eta["Pt10"][cutstring] = infiles_pt10[cutstring].Get("fake_rate_eta")

        eff_eta["Pt100"][cutstring] = infiles_pt100[cutstring].Get("eff_eta")
        fake_eta["Pt100"][cutstring] = infiles_pt100[cutstring].Get("fake_rate_eta")

        for eta_region in eta_regions:
            eff_pt["FlatPt_" + eta_region][cutstring] = infiles_flatPt[cutstring].Get("eff_pt_" + eta_region)
            fake_pt["FlatPt_" + eta_region][cutstring] = infiles_flatPt[cutstring].Get("fake_pt_" + eta_region)


    print "Plotting efficiencies and fake rates"

    print eff_eta["Pt10"].values()

    for scan in parameters:
        if len(parameters[scan]) > 1:
            sel_str = scan + "Scan"

    draw_and_save_eff(eff_eta["Pt10"], "eta", "eff", is_gsf=False, label=sel_str+"_Pt10", leg_pos="down_right", title="el. p_{T} = 10")
    draw_and_save_eff(fake_eta["Pt10"], "eta", "fake", is_gsf=False, label=sel_str+"_Pt10", leg_pos="up_right", title="el. p_{T} = 10")

    draw_and_save_eff(eff_eta["Pt100"], "eta", "eff", is_gsf=False, label=sel_str+"_Pt100", leg_pos="down_right", title="el. p_{T} = 100")
    draw_and_save_eff(fake_eta["Pt100"], "eta", "fake", is_gsf=False, label=sel_str+"_Pt100", leg_pos="up_right", title="el. p_{T} = 100")

    for eta_region in eta_regions:
        draw_and_save_eff(eff_pt["FlatPt_" + eta_region], "pt", "eff", is_gsf=False, label=sel_str+"_FlatPt_"+eta_region, leg_pos="down_right", title="barrel el.")
        draw_and_save_eff(fake_pt["FlatPt_" + eta_region], "pt", "fake", is_gsf=False, label=sel_str+"_FlatPt_"+eta_region, leg_pos="up_right", title="barrel el.")

#draw_efficiency_histograms(eff_eta["Pt10"].values(), xtitle="#eta", ytitle="Efficiency", ymax = 1.)

#---------------pt----------------------
#draw_and_save_eff(eff_pt, "pt", "eff", is_gsf = is_gsf, leg_pos = "down_right", title="Total efficiency")
#draw_and_save_eff(eff_eta, "eta", "eff", is_gsf = is_gsf, leg_pos="down_right", title= "Total efficiency")

#----------------- fake --------------
#draw_and_save_eff(fake_eta, "eta", "fake", is_gsf = is_gsf, leg_pos = "up_right")
#draw_and_save_eff(fake_pt, "pt", "fake", is_gsf = is_gsf, leg_pos = "up_right")

