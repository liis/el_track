from array import *
maxpart = 100

int_list = ["np_reco", "np_gen", "np_fake", "np_gen_toReco"]
int_array_list = [] #["reco_genMatched", "gen_recoMatched"]
double_list = []
double_array_list = ["reco_pt", "reco_eta", "reco_phi", "gen_pt", "gen_eta", "gen_phi", "fake_pt", "fake_eta", "fake_phi", "gen_matched_eta", "gen_matched_phi", "gen_matched_pt"] 

var_list = int_list + int_array_list + double_list + double_array_list

def var_type(var_name):
    """
    return the variable type for tree entry
    """

    if var_name in int_list: return array('i',[0])
    if var_name in int_array_list: return array('i',[0]*maxpart)
    if var_name in double_list: return array('d', [0])
    if var_name in double_array_list: return array('d', [0]*maxpart)
