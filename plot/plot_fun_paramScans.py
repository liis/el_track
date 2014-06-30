import ROOT
from collections import OrderedDict as dict

def get_infilenames_by_params(filename, parameters, ptmode):
    print "Opening input files for " + ptmode + " samples"

    infiles = dict()

    nr_maxChi2 = len(parameters["maxChi2"])
    nr_nSigma = len(parameters["nSigma"])
    nr_maxCand = len(parameters["maxCand"])

    for it_maxChi2 in range(0, nr_maxChi2):
        for it_nSigma in range(0, nr_nSigma):
            for it_maxCand in range(0, nr_maxCand):
                cutstring = "maxCand_" + str(parameters["maxCand"][it_maxCand]) + "_maxChi2_" + str(parameters["maxChi2"][it_maxChi2]) + "_nSigma_" + str(parameters["nSigma"][it_nSigma])

                infile = filename + "_" + ptmode + "_" + cutstring + ".root"
                print "Adding: " + infile + " to the analysis "
                infiles[cutstring] = infile
    return infiles
