import re

from odict import OrderedDict as dict
tracking_cfg_parameters = dict()

tracking_cfg_parameters["maxCand"] = [5, 5, 5, 6  ]
tracking_cfg_parameters["maxChi2"] = [2000, 1000, 300, 1000 ]
tracking_cfg_parameters["nSigma"] = [3, 3, 3, 3 ]

datasetnames = {"FlatPt": "/SingleElMinusFlatLogPt_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-8de1ffbdb519f9edbafc5606a1926f13/USER",
                "Pt10": "/SingleElMinusPt10_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG-v2/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER",
                "Pt100": "/SingleElMinusPt100_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG-v2/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER"
    }

def read_template(filename):
    f = open(filename)
    s = f.read()
    f.close()
    return s

def write_outfile(input, outname):
    print "Saving modified file as: " + outname
    out_file = open(outname, "w")
    out_file.write(input)
    out_file.close()
                    
def create_crab_cfg_from_template(template, varstr, dataset, outdir = ""):
    """
    dotaset = Pt100, Pt10, FlatPt
    """
    input = read_template(template)
    input = input.replace("CUTVALS", varstr)    
    input = input.replace("DATASET", dataset)
    input = input.replace("DSNAME", datasetnames[dataset])

    if len(outdir):
        outdir = outdir + "/"

    out_name = outdir + "crab_" + varstr + "_" + dataset + ".cfg"
    write_outfile(input, out_name)

def create_cmssw_cfg_from_template(template, varstr, outdir = ""):
    """
    template -- template file name
    varstr -- string of variables and their values, separated by _
    """
    input = read_template(template)

    input = input.replace("VARSTR", varstr)
    input = input.replace("MAXCAND", varstr.rsplit("maxCand_")[1].rsplit("_")[0] )
    input = input.replace("MAXCHI2", varstr.rsplit("maxChi2_")[1].rsplit("_")[0] )
    input = input.replace("NSIGMA", varstr.rsplit("nSigma_")[1].rsplit("_")[0] )

    if len(outdir):
        outdir = outdir + "/"
    out_name = outdir + "makeTrackValTree_reTrk_" + varstr + ".py"
    write_outfile(input, out_name)

#def make_varstrings( cfg_params ):


def create_varstrings(tracking_cfg_parameters, iter = 0, skip_default = True):
    """
    tracking_cfg_parameters -- dictionary of varnames and lists of possible values
    iter -- number of parametersets to run 1 ... b
    skip_default -- skip running on the default set of parameters (under index 0)
    """
    varstrings = []
    for i in range(iter):
        varstr = ""
        if skip_default and i==0:
            continue
        for varname, varvalue in tracking_cfg_parameters.iteritems():
            varstr = varstr + varname + "_" + str(varvalue[i]) + "_"    
        varstr = varstr[:-1]
        varstrings.append(varstr)
        
    return varstrings


varstrings = create_varstrings(tracking_cfg_parameters, iter = 4)

for varstr in varstrings:
    create_cmssw_cfg_from_template("./templates/makeTrackValTree_reTrk_template.py", varstr, outdir = "input_crab") # run both retracking and tree production
    for datasetname in datasetnames:
        create_crab_cfg_from_template("./templates/crab_template.cfg", varstr, dataset = datasetname, outdir = "input_crab")
    
       
