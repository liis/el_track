import re
import sys

sys.path.insert(0, '/afs/cern.ch/user/l/liis/scratch0/tracking/el_track/rereco/CMSSW/src/TestReReco/submission_scripts')
from sub_fun import read_template, write_outfile, create_cmssw_cfg_from_template, create_crab_cfg_from_template, create_varstrings

from collections import OrderedDict as dict
tracking_cfg_parameters = dict()

# Chi2 scan
tracking_cfg_parameters["maxCand"] = [5, 5, 5, 5, 5, 5 ]
tracking_cfg_parameters["maxChi2"] = [10, 30, 50, 100, 300, 2000]
tracking_cfg_parameters["nSigma"] = [3, 3, 3, 3, 3, 3 ]

# Max cand scan
tracking_cfg_parameters["maxCand"] += [3, 4, 6, 7]
tracking_cfg_parameters["maxChi2"] += [2000, 2000, 2000, 2000]
tracking_cfg_parameters["nSigma"] += [3, 3, 3, 3]

# nSigma scan
tracking_cfg_parameters["maxCand"] += [5, 5, 5, 5]
tracking_cfg_parameters["maxChi2"] += [2000, 2000, 2000, 2000]
tracking_cfg_parameters["nSigma"] += [2, 4, 5, 6]

print "Creating " + str(len(tracking_cfg_parameters["maxCand"]) )+ " jobs"

varstrings = create_varstrings(tracking_cfg_parameters, iter = len(tracking_cfg_parameters["maxCand"]), skip_default = False)


datasets = { 
#    "Flat_Pt": "",
#    "Pt10": "/SingleElectronPt10_GENRAW/liis-SingleElectronPt10_RECO-b0d97c8144eaac9090ed0cd9df0f13de/USER",
#    "Pt100": "/SingleElectronPt100_GENRAW/liis-SingleElectronPt100_RECO-b0d97c8144eaac9090ed0cd9df0f13de/USER",
#    "Zee": "/DYJetsToLL_M-50_13TeV-pythia6/phys_egamma-EGM711_PU40bx25_POSTLS171_V11_RECODEBUG-v1-ffac44eb0cb582bdcc6ecfb3c5f327a8/USER",
    "FlatPt": "/SingleElectronFlatPt_GENRAW/liis-SingleElectronFlatPt_RECO-b0d97c8144eaac9090ed0cd9df0f13de/USER",
    }



outdir = "./input_crab"
for varstr in varstrings:
    #--- create cmssw cfg files for each set of parameters
    create_cmssw_cfg_from_template("../templates/makeTrackValTree_reTracking_crab_template.py", varstr, outdir=outdir, mode = "crab")
    for dataset in datasets:
        #--- create crab.cfg files
        create_crab_cfg_from_template("../templates/crab_template.cfg", varstr, dataset, datasets[dataset], outdir = outdir)


        
