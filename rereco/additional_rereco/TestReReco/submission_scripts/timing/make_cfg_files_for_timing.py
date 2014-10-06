import re
import sys
sys.path.insert(0, '../../submission_scripts')

from odict import OrderedDict as dict
tracking_cfg_parameters = dict()

from sub_fun import read_template, write_outfile, create_cmssw_cfg_from_template, create_check_timing_script, create_varstrings
from infilelists import infilelist_sam

# Chi2 scan
tracking_cfg_parameters["maxCand"] = [5, 5, 5]
tracking_cfg_parameters["maxChi2"] = [5, 30, 2000]
tracking_cfg_parameters["nSigma"] = [3, 3, 3]

# Max cand scan
tracking_cfg_parameters["maxCand"] += [1, 3, 5, 7]
tracking_cfg_parameters["maxChi2"] += [2000, 2000, 2000, 2000]
tracking_cfg_parameters["nSigma"] += [3, 3, 3, 3]

# nSigma scan
tracking_cfg_parameters["maxCand"] += [5, 5, 5]
tracking_cfg_parameters["maxChi2"] += [2000, 2000, 2000]
tracking_cfg_parameters["nSigma"] += [1, 3, 5]

print "Creating " + str(len(tracking_cfg_parameters["maxCand"]) )+ " jobs"

varstrings = create_varstrings(tracking_cfg_parameters, iter = len(tracking_cfg_parameters["maxCand"]), skip_default = False)

filelists = {
#    "Flat_Pt": "",
#    "Pt10": "",
#    "Pt100": "",
    "Zee": infilelist_sam,
    }

for dataset in filelists:
    for varstr in varstrings:
        if dataset == "Zee":
            isSinglePart=False
        else:
            isSinglePart=True

        cmssw_cfg_name = create_cmssw_cfg_from_template("../templates/makeTrackValTree_reTracking_timing_template.py", varstr, isSinglePart, outdir = "./input_timing", mode = "batch", dataset = dataset, infiles = filelists[dataset])

        create_check_timing_script(cmssw_cfg_name)
    
       
