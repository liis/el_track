import re

from odict import OrderedDict as dict
tracking_cfg_parameters = dict()

from sub_fun import read_template, write_outfile, create_cmssw_cfg_from_template, create_batch_submission_script, create_varstrings
from infilelists import infilelist_Zee

# Chi2 scan
tracking_cfg_parameters["maxCand"] = [5, 5, 5, 5, 5, 5 ]
tracking_cfg_parameters["maxChi2"] = [10, 30, 50, 100, 300, 2000]
tracking_cfg_parameters["nSigma"] = [3, 3, 3, 3, 3, 3 ]

# Max cand scan
tracking_cfg_parameters["maxCand"] += [3, 4, 6, 7]
tracking_cfg_parameters["maxChi2"] += [2000, 2000, 2000, 2000]
tracking_cfg_parameters["nSigma"] += [3, 3, 3, 3]

# nSigma scan
tracking_cfg_parameters["maxCand"] += [5, 5, 5]
tracking_cfg_parameters["maxChi2"] += [2000, 2000, 2000]
tracking_cfg_parameters["nSigma"] += [2, 4, 5]

print "Creating " + str(len(tracking_cfg_parameters["maxCand"]) )+ " jobs"

varstrings = create_varstrings(tracking_cfg_parameters, iter = len(tracking_cfg_parameters["maxCand"]), skip_default = False)

filelists = {
#    "Flat_Pt": "",
#    "Pt10": "",
#    "Pt100": "",
    "Zee": infilelist_Zee,
    }


for varstr in varstrings:
    for dataset in filelists:
        cmssw_cfg_name = create_cmssw_cfg_from_template("./templates/makeTrackValTree_reTrk_template.py", varstr, outdir = "./batch_jobs/input_batch", mode = "batch", dataset = dataset, infiles = filelists[dataset])

        create_batch_submission_script( cmssw_cfg_name )
    
       
