import re, ntpath
import time, datetime


from collections import OrderedDict as dict
tracking_cfg_parameters = dict()

tracking_cfg_parameters["maxCand"] = [5, 5, 5, 5, 5, 5 ]
tracking_cfg_parameters["maxChi2"] = [10, 30, 50, 100, 300, 2000]
tracking_cfg_parameters["nSigma"] = [3, 3, 3, 3, 3, 3 ]

datasetnames = {
#    "FlatPt": "/SingleElMinusFlatLogPt_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-8de1ffbdb519f9edbafc5606a1926f13/USER",
#    "Pt10": "/SingleElMinusPt10_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG-v2/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER",
#    "Pt100": "/SingleElMinusPt100_CMSSW_4_2_8-START42_V12_GEN-SIM-DIGI-RAW-HLTDEBUG-v2/mangano-CMSSW_4_2_8-START42_V12_GEN-SIM-RECO-v3-7bc796286602c18e9ed77a7f93a692b8/USER",
    "Zee": "/RelValZEE/CMSSW_4_2_9_HLT1_patch1-START42_V14B_RelVal_ZEErv_20Jun2013-v1/GEN-SIM-DIGI-RAW-HLTDEBUG",
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
                    
def create_crab_cfg_from_template(template, varstr, dataset, datasetname, leadingVertexOnly, isGsf=True, outdir = ""):
    """
    dataset = Pt100, Pt10, FlatPt
    datasetname --> corresponding full datasetname
    leadingVertexOnly -- consider only tracks from the leading vertex for fake rate studies
    """
    input = read_template(template)
    input = input.replace("CUTVALS", varstr)    
    input = input.replace("DATASET", dataset)
    input = input.replace("DSNAME", datasetname)

    input = input.replace("OUTDIR", outdir)
    if len(outdir):
        outdir = outdir + "/"

    if leadingVertexOnly:
        singlePartLabel = "_leadVtx"
        input = input.replace("EVTS_PER_JOB", "1500")
    else:
        singlePartLabel = ""
        input = input.replace("EVTS_PER_JOB", "4000")

    if isGsf:
        gsfLabel = "_gsf"
    else:
        gsfLabel = ""

    input = input.replace("SINGLEPART", singlePartLabel)
    input = input.replace("ISGSF", gsfLabel)

    out_name = outdir + "crab_" + varstr + "_" + dataset +  gsfLabel + ".cfg"
    write_outfile(input, out_name)

def create_cmssw_cfg_from_template(template, varstr, leadingVertexOnly, isGSF=True, outdir = "", mode = "batch", isGsf=True, dataset = "", infiles = [], infiles_sec = []):
    """
    template -- template file name
    varstr -- string of variables and their values, separated by _
    """
    input = read_template(template)

    if isGsf:
        gsfLabel = "_gsf"
    else:
        gsfLabel = ""

    input = input.replace("VARSTR", varstr)
    input = input.replace("MAXCAND", varstr.rsplit("maxCand_")[1].rsplit("_")[0] )
    input = input.replace("MAXCHI2", varstr.rsplit("maxChi2_")[1].rsplit("_")[0] )
    input = input.replace("NSIGMA", varstr.rsplit("nSigma_")[1].rsplit("_")[0] )
    input = input.replace("ISGSF", str(isGSF))
    input = input.replace("GSFLABEL", gsfLabel)
    input = input.replace("LEADINGVERTEXONLY", str(leadingVertexOnly))

    #    if mode == "crab":
    #        input = input.replace("OUTFILENAME", " 'trackValTree_reTrk.root' ")

    if not mode == "crab":
        input = input.replace("OUTFILENAME", " ' " + "../output_batch/trackValTree_" + dataset + "_" + varstr + ".root" + " ' ")

        input = input.replace("INFILELIST", str(infiles))
        input = input.replace("SECFILELIST", str(infiles_sec))

    if leadingVertexOnly:
        varstr = varstr+"_leadVtx"

    if isGsf:
        varstr = varstr+"_gsf"
    
    if len(outdir):
        outdir = outdir + "/"
    out_name = outdir + "makeTrackValTree_reTrk_" + varstr + ".py"
    write_outfile(input, out_name)
    return out_name


def create_batch_submission_script(cmssw_cfg_file): # create .sh file to submit cmsRun job to batch
    outname = cmssw_cfg_file.split(".py")[0] +".sh"
    out_file = open(outname, "w")
    
    out_file.write("#!/bin/bash\n\n")

    out_file.write('cd ${CMSSW_BASE}/src/TestReReco/submission_scripts/batch_jobs/input_batch \n')
    out_file.write('eval `scramv1 runtime -sh`\n')
    out_file.write("cmsRun " + ntpath.basename(cmssw_cfg_file) + "\n")
    
    out_file.close()

    print "Saved submission script as: " + outname
                        

def create_check_timing_script(cmssw_cfg_file):
    outname = cmssw_cfg_file.split(".py")[0] +".sh"
    cutstr = (cmssw_cfg_file.split(".py")[0]).split("makeTrackValTree_reTrk_")[1]

    ts = time.time()
    timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d%H%M%S')
    
    timinglogfile = "DetailedInfo_" + cutstr + ".txt"
    out_file = open(outname, "w")

    out_file.write("#!/bin/bash\n\n")
#    out_file.write('cd ./input_timing \n')
    out_file.write('cd ${CMSSW_BASE}/src/TestReReco/submission_scripts/timing/input_timing \n')

    out_file.write('timingdir=timing_' + cutstr + '_' +'"$(date +"%s")" \n') #create the timestamp at execution

    out_file.write('mkdir -p $timingdir \n')
    out_file.write('eval `scramv1 runtime -sh`\n')

#    out_file.write("cmsRun " + ntpath.basename(cmssw_cfg_file) + " 2>" + timinglogfile +   "\n")
    out_file.write("cmsRun " + ntpath.basename(cmssw_cfg_file) + "\n")

    out_file.write("echo 'Saving timing information...' " + "\n") 
    out_file.write("mv " + timinglogfile + " $timingdir \n")
    out_file.write("cd $timingdir \n")
    out_file.write("grep TimeModule\> " + timinglogfile +  " > TimingInfo.txt" + "\n")
    out_file.write("echo 'running timing.cpp' \n")
    out_file.write("g++ ${CMSSW_BASE}/src/TestReReco/timing.cpp \n")
    out_file.write("./a.out>final_timing_output.txt \n")
    out_file.write("grep CkfTrackCandidateMaker final_timing_output.txt  >> final_timing_eltrk_summary.txt \n")
    out_file.write("grep Gsf final_timing_output.txt  >> final_timing_eltrk_summary.txt \n")
    out_file.write("grep SeedGeneratorFromRegionHitsEDProducer final_timing_output.txt  >> final_timing_eltrk_summary.txt \n")
    out_file.write("echo '...done'")
    out_file.close()
    print "Saved submission script as: " + outname


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
    
       
