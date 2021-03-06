import sys, os
import ROOT
ROOT.gSystem.Load('RooDoubleCB/RooDoubleCB')
from ROOT import RooDoubleCB

import tdrstyle
tdrstyle.tdrstyle()
ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)

from plot_fun import infilenames_eta, infilenames_pt, infilenames_eta_gsf, infilenames_pt_gsf, load_input_files, draw_efficiency_histograms, draw_legend, draw_and_save_res, draw_resolution

from plot_fun_paramScans import get_infilenames_by_params

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--indir', dest='indir', required=True) #directory for input histograms
args = parser.parse_args()
indir = args.indir #"../tree_to_histo/histograms/" #location of input histograms -->clean up!

print "Opening input files from " + indir

is_gsf = 1
outdir="./out_plots_paramScans/13TeV_011014_res"

#--------- run over all combinations of parameter values in list of input files, have multiple values for only 1 parameter for comparisons ----------

parameters_maxCand = {
    "maxChi2": [2000],
    "nSigma": [3],
 #   "maxCand": [1,2, 3, 4, 5, 6, 7],
    "maxCand": [1,2]
    } #consider all combinations of parameters

parameters_maxChi2 = {
    #    "maxChi2": [10, 30, 50, 100, 300, 2000],
    "maxChi2": [2, 5, 8, 10, 30, 50, 100, 300, 2000],
    "nSigma": [3],
    "maxCand": [5]
    }

parameters_nSigma = {
    "maxChi2": [2000],
    #    "nSigma": [5, 6],
    "nSigma": [1, 2, 3, 4, 5, 6],
    #    "nSigma": [ 2, 3, 4, 5, 6],
    "maxCand": [5],
    }

parameters_test = {
  "maxChi2": [30, 50, 100, 300, 2000],
  "nSigma": [3],
  "maxCand": [5],
}

parameter_sets = [
    parameters_maxCand,
#    parameters_maxChi2,
#    parameters_nSigma,
#    parameters_test,
    ]

input_files = [
#    "Pt10",
#    "Pt100",
#    "FlatPt",
    "Zee"
]

infilenames = {} # filenames
infiles = {} # ROOT files
eta_regions = ["barrel", "endcap", "trans"]

ii = 0
for parameters in parameter_sets:

    from collections import OrderedDict as dict
    res_eta_dxy_68 = {}
    res_eta_dz_68 = {}
    res_eta_cotth_68 = {}
    res_eta_cotth_68 = {}
    res_eta_pt_68 = {}
    res_eta_phi_68 = {}
    res_eta_dxy_95 = {}
    res_eta_dz_95 = {}
    res_eta_cotth_95 = {}
    res_eta_cotth_95 = {}
    res_eta_pt_95 = {}
    res_eta_phi_95 = {}
                              
    res_pt_dxy_68 = {}
    res_pt_dz_68 = {}
    res_pt_cotth_68 = {}
    res_pt_cotth_68 = {}
    res_pt_pt_68 = {}
    res_pt_phi_68 = {}
    res_pt_dxy_95 = {}
    res_pt_dz_95 = {}
    res_pt_cotth_95 = {}
    res_pt_cotth_95 = {}
    res_pt_pt_95 = {}
    res_pt_phi_95 = {}
    
    for input_file in input_files:
        infilenames[input_file] = get_infilenames_by_params("efficiencyHistograms", parameters, input_file, isGsf = is_gsf)
        infiles[input_file] = load_input_files(indir, infilenames[input_file])

        res_eta_dxy_68[input_file] = dict()
        res_eta_dxy_95[input_file] = dict()
        res_eta_dz_68[input_file] = dict()
        res_eta_dz_95[input_file] = dict()
        res_eta_cotth_68[input_file] = dict()
        res_eta_cotth_95[input_file] = dict()
        res_eta_pt_68[input_file] = dict()
        res_eta_pt_95[input_file] = dict()
        res_eta_phi_68[input_file] = dict()
        res_eta_phi_95[input_file] = dict()

        for eta_region in eta_regions:
            res_pt_dxy_68[input_file + "_" + eta_region] = dict()
            res_pt_dxy_95[input_file + "_" + eta_region] = dict()
            res_pt_dz_68[input_file + "_" + eta_region] = dict()
            res_pt_dz_95[input_file + "_" + eta_region] = dict()
            res_pt_cotth_68[input_file + "_" + eta_region] = dict()
            res_pt_cotth_95[input_file + "_" + eta_region] = dict()
            res_pt_pt_68[input_file + "_" + eta_region] = dict()
            res_pt_pt_95[input_file + "_" + eta_region] = dict()
            res_pt_phi_68[input_file + "_" + eta_region] = dict()
            res_pt_phi_95[input_file + "_" + eta_region] = dict()

        for cutstring in infiles[input_file]:
            if ii == 0: # for control fit plots show only for the first entry (otherwise too many)
                isFirst = True
            else:
                isFirst = False
            ii+=1
            # print "isFirst = " +str(isFirst)

            if input_file != "FlatPt": #skip eta histograms for flatPt
                test = infiles[input_file][cutstring].Get("resolutions_eta/res_dxy_vs_eta")
                #                test.Draw() ## for debug -- check if necessary res histograms are present

#                res_eta_dxy_68[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_dxy_vs_eta"), "res_dxy_vs_eta_"+input_file+"_"+cutstring + "_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
#                res_eta_dxy_95[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_dxy_vs_eta"), "res_dxy_vs_eta_"+input_file+"_"+cutstring + "_95", outdir=outdir, do_control_fit_plots=isFirst, mode="95")

#                res_eta_dz_68[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_dz_vs_eta"), "res_dz_vs_eta"+input_file+"_"+cutstring + "_68", outdir=outdir, do_control_fit_plots=isFirst, mode = "68")
#                res_eta_dz_95[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_dz_vs_eta"), "res_dz_vs_eta"+input_file+"_"+cutstring + "_95", outdir=outdir, do_control_fit_plots=isFirst, mode = "95")
#                res_eta_cotth_68[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_cotth_vs_eta"), "res_cotth_vs_eta"+input_file+"_"+cutstring + "_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
#                res_eta_cotth_95[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_cotth_vs_eta"), "res_cotth_vs_eta"+input_file+"_"+cutstring + "_95", outdir=outdir, do_control_fit_plots=isFirst, mode="95")
#                res_eta_pt_68[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_pt_vs_eta"), "res_pt_vs_eta"+input_file+"_"+cutstring + "_68", outdir=outdir, do_control_fit_plots=isFirst, mode = "68")
#                res_eta_pt_95[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_pt_vs_eta"), "res_pt_vs_eta"+input_file+"_"+cutstring + "_95", outdir=outdir, do_control_fit_plots=isFirst, mode = "95")
#                res_eta_phi_68[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_phi_vs_eta"), "res_phi_vs_eta"+input_file+"_"+cutstring + "_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
#                res_eta_phi_95[input_file][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_eta/res_phi_vs_eta"), "res_phi_vs_eta"+input_file+"_"+cutstring + "_95", outdir=outdir, do_control_fit_plots=isFirst, mode="95")


            if input_file != "Pt10" and input_file != "Pt100": # skip pt histograms for fixed pt samples
                print res_pt_dxy_68
                for eta_region in eta_regions:
                    for cutstring in infiles[input_file]:
                        #                       res_pt_dxy_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_dxy_vs_pt"), "res_dxy_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        #                       res_pt_dxy_95[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_dxy_vs_pt"), "res_dxy_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_95", outdir=outdir, do_control_fit_plots=isFirst, mode="95")
                        #                        res_pt_dz_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_dz_vs_pt"), "res_dz_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        #                        res_pt_dz_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_dz_vs_pt"), "res_dz_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        #                        res_pt_cotth_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_cotth_vs_pt"), "res_cotth_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        #                        res_pt_cotth_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_cotth_vs_pt"), "res_cotth_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        #                        res_pt_phi_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_phi_vs_pt"), "res_phi_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        #                        res_pt_phi_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_phi_vs_pt"), "res_phi_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        res_pt_pt_68[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_pt_vs_pt"), "res_pt_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_68", outdir=outdir, do_control_fit_plots=isFirst, mode="68")
                        res_pt_pt_95[input_file + "_" + eta_region][cutstring] = draw_resolution( infiles[input_file][cutstring].Get("resolutions_pt/res_pt_vs_pt"), "res_pt_vs_pt_"+input_file+"_"+eta_region+"_"+cutstring+"_95", outdir=outdir, do_control_fit_plots=isFirst, mode="95")



#----------------------------------- plot and save output ---------------------------------------------------

    print "Plotting efficiencies and fake rates"

    for scan in parameters:
        if len(parameters[scan]) > 1:
            sel_str = scan + "Scan"
        elif len(parameters[scan]) == 1:
            "Insert at least two entries in the scan!"


    for input_file in input_files:
        print "Drawing histogram from file: " + str(input_file)

        if input_file != "flatPt":
            print "saving eta resolution"
            
            # ------------------- eta resolutions -----------------------
            #draw_and_save_res(res_eta_dxy_68[input_file], res_eta_dxy_95[input_file], "eta", "dxy", sel_str, is_gsf, outdir=outdir, logy=True)
            #draw_and_save_res(res_eta_dz_68[input_file], res_eta_dz_95[input_file], "eta", "dz", sel_str, is_gsf, outdir=outdir, logy=True)
            #draw_and_save_res(res_eta_cotth_68[input_file], res_eta_cotth_95[input_file], "eta", "cotth", sel_str, is_gsf, outdir=outdir, logy=True)
            #draw_and_save_res(res_eta_pt_68[input_file], res_eta_pt_95[input_file], "eta", "pt", sel_str, is_gsf, outdir=outdir, logy=True)
            #draw_and_save_res(res_eta_phi_68[input_file], res_eta_phi_95[input_file], "eta", "phi", sel_str, is_gsf, outdir=outdir, logy=True)

            # --------------------- pt resolutions ------------------------

        if input_file != "Pt10" and input_file != "Pt100":
            for eta_region in eta_regions:
#                draw_and_save_res(res_pt_dxy_68[input_file + "_" + eta_region], res_pt_dxy_95[input_file + "_" + eta_region], "pt", "dxy", sel_str, is_gsf, outdir=outdir, logy=True)
#                draw_and_save_res(res_pt_dz_68[input_file + "_" + eta_region], res_pt_dz_95[input_file + "_" + eta_region], "pt", "dz", sel_str, is_gsf, outdir=outdir, logy=True)
#                draw_and_save_res(res_pt_cotth_68[input_file + "_" + eta_region], res_pt_cotth_95[input_file + "_" + eta_region], "pt", "cotth", sel_str, is_gsf, outdir=outdir, logy=True)
#                draw_and_save_res(res_pt_phi_68[input_file + "_" + eta_region], res_pt_phi_95[input_file + "_" + eta_region], "pt", "phi", sel_str, is_gsf, outdir=outdir, logy=True)
                draw_and_save_res(res_pt_pt_68[input_file + "_" + eta_region], res_pt_pt_95[input_file + "_" + eta_region], "pt", "pt", sel_str, is_gsf, outdir=outdir, logy=True)

##end
