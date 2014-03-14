import ROOT

def get_infilenames_by_params(parameters, ptmode):
    print "Opening input files for " + ptmode + " samples"
    infiles = {}

    nr_maxChi2 = len(parameters["maxChi2"])
    nr_nSigma = len(parameters["nSigma"])
    nr_maxCand = len(parameters["maxCand"])

    for it_maxChi2 in range(0, nr_maxChi2):
        for it_nSigma in range(0, nr_nSigma):
            for it_maxCand in range(0, nr_maxCand):
                cutstring = "maxCand_" + str(parameters["maxCand"][it_maxCand]) + "_maxChi2_" + str(parameters["maxChi2"][it_maxChi2]) + "_nSigma_" + str(parameters["nSigma"][it_nSigma])

                infile = "trackValHistograms_" + ptmode + "_" + cutstring + ".root"
                print "Adding: " + infile + " to the analysis "
                infiles[cutstring] = infile
    return infiles
