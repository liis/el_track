import sys, os
import ROOT

import tdrstyle
tdrstyle.tdrstyle()

from plot_fun import infilenames_eta, infilenames_pt, infilenames_eta_gsf, infilenames_pt_gsf, load_input_files, draw_efficiency_histograms, draw_legend, draw_and_save_eff

from plot_fun_paramScans import get_infilenames_by_params

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--indir', dest='indir', required=True) #directory for input histograms
args = parser.parse_args()
indir = args.indir #"../tree_to_histo/histograms/" #location of input histograms -->clean up!

print "Opening input files from " + indir

is_gsf = 0
#--------- run over all combinations of parameter values in list of input files, have multiple values for only 1 parameter for comparisons ----------

parameters_maxCand = {
    "maxChi2": [2000],
    "nSigma": [3],
    "maxCand": [3, 4, 5, 6, 7],
    } #consider all combinations of parameters

parameters_maxChi2 = {
    "maxChi2": [10, 30, 50, 100, 300, 2000],
    "nSigma": [3],
    "maxCand": [5]
     }

parameters_nSigma = {
    "maxChi2": [2000],
#    "nSigma": [1, 2, 3, 4, 5],
    "nSigma": [2, 3, 4, 5, 6],
    "maxCand": [5],
    }

parameters_test = {
  "maxChi2": [ 50, 100, 300, 2000],
  "nSigma": [3],
  "maxCand": [5],
}

parameter_sets = [
    parameters_maxCand,
    parameters_maxChi2,
    parameters_nSigma,
#  parameters_test,
    ]

input_files = [
    "Pt10",
    "Pt100",
#    "FlatPt",
    "Zee"
]

infilenames = {} # filenames
infiles = {} # ROOT files

for parameters in parameter_sets:

    from collections import OrderedDict as dict
    fake_eta = {}
    eff_eta = {}
    eff_seed_eta = {}
    eff_wrt_seed_eta = {}
    eff_eta_sim = {}

    eff_pt = {}
    eff_seed_pt = {}
    eff_wrt_seed_pt = {}
    fake_pt = {}
    eta_regions = ["barrel", "endcap", "trans"]


    for input_file in input_files:
        infilenames[input_file] = get_infilenames_by_params("efficiencyHistograms", parameters, input_file)
        infiles[input_file] = load_input_files(indir, infilenames[input_file])

        if input_file != "FlatPt":
            fake_eta[input_file] = dict()
            eff_eta[input_file] = dict()
            eff_wrt_seed_eta[input_file] = dict()
            eff_seed_eta[input_file] = dict()
#            eff_eta_sim[input_file] = dict()

            for cutstring in infiles[input_file]:
                eff_eta[input_file][cutstring] = infiles[input_file][cutstring].Get("efficiencies/eff_eta")
                eff_wrt_seed_eta[input_file][cutstring] = infiles[input_file][cutstring].Get("efficiencies/eff_wrt_seed_eta")
                eff_seed_eta[input_file][cutstring] = infiles[input_file][cutstring].Get("efficiencies/eff_seed_eta")
#                eff_eta_sim[input_file][cutstring] = infiles[input_file][cutstring].Get("eff_eta_simMatch_sel")
                fake_eta[input_file][cutstring] = infiles[input_file][cutstring].Get("efficiencies/fake_rate_eta")

        if input_file != "Pt10" and input_file != "Pt100":
            for eta_region in eta_regions:
                eff_pt[input_file + "_" + eta_region] = dict()
                eff_wrt_seed_pt[input_file + "_" + eta_region] = dict()
                eff_seed_pt[input_file + "_" + eta_region] = dict()
#                eff_pt_sim[input_file + "_" + eta_region] = dict()
                fake_pt[input_file + "_" + eta_region] = dict()

                for cutstring in infiles[input_file]:
                    eff_pt[input_file + "_" + eta_region][cutstring] = infiles[input_file][cutstring].Get("efficiencies/eff_pt_" + eta_region)
                    eff_wrt_seed_pt[input_file + "_" + eta_region][cutstring] = infiles[input_file][cutstring].Get("efficiencies/eff_wrt_seed_pt_" + eta_region)
                    eff_seed_pt[input_file + "_" + eta_region][cutstring] = infiles[input_file][cutstring].Get("efficiencies/eff_seed_pt_" + eta_region)

                    fake_pt[input_file + "_" + eta_region][cutstring] = infiles[input_file][cutstring].Get("efficiencies/fake_rate_pt_" + eta_region)


#----------------------------------- plot and save output ---------------------------------------------------

    print "Plotting efficiencies and fake rates"

    for scan in parameters:
        if len(parameters[scan]) > 1:
            sel_str = scan + "Scan"


    for input_file in input_files:
        print "Drawing histogram from file: " + str(input_file)
        print str(eff_eta[input_file])
        if input_file != "FlatPt":
            draw_and_save_eff(eff_eta[input_file], "eta", "eff", is_gsf=False, label=sel_str+"_" + input_file, leg_pos="down_right", title=input_file)
            draw_and_save_eff(eff_seed_eta[input_file],"eta", "eff_seed", is_gsf=False, label=sel_str+"_"+input_file, leg_pos="down_right",title=input_file)
            draw_and_save_eff(eff_wrt_seed_eta[input_file],"eta", "eff_wrt_seed", is_gsf=False, label=sel_str+"_"+input_file, leg_pos="down_right",title=input_file)
#            draw_and_save_eff(eff_eta_sim[input_file], "eta", "eff", is_gsf=False, label="sim_" + sel_str + "_" + input_file, leg_pos="down_right", title="sim, " + input_file)
            draw_and_save_eff(fake_eta[input_file], "eta", "fake", is_gsf=False, label=sel_str+"_" + input_file, leg_pos="up_right", title=input_file)

        if input_file != "Pt10" and input_file != "Pt100":
            for eta_region in eta_regions:
                draw_and_save_eff(eff_pt[input_file + "_" + eta_region], "pt", "eff", is_gsf=False, label=sel_str+"_" + input_file + "_" + eta_region, leg_pos="down_right", title=eta_region + " el. , " + input_file)
                draw_and_save_eff(eff_seed_pt[input_file + "_" + eta_region], "pt", "eff_seed", is_gsf=False, label=sel_str+"_" + input_file + "_" + eta_region, leg_pos="down_right", title=eta_region + " el. , " + input_file)
                draw_and_save_eff(eff_wrt_seed_pt[input_file + "_" + eta_region], "pt", "eff_wrt_seed", is_gsf=False, label=sel_str+"_" + input_file + "_" + eta_region, leg_pos="down_right", title=eta_region + " el. , " + input_file)

                draw_and_save_eff(fake_pt[input_file + "_" + eta_region], "pt", "fake", is_gsf=False, label=sel_str+"_" + input_file + "_" + eta_region, leg_pos="up_right", title=eta_region + " el." + input_file)

##end                
