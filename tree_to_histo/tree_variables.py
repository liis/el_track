from array import *
import numpy as np
maxpart = 100
maxhit = 30

int_list = ["np_reco", "np_gen", "np_fake", "np_gen_toReco", "gen_matchedSeedOkCharge"]
int_array_list = ["gen_nrMatchedSeedHits", "gen_nrUniqueSimHits", "gen_reco_matched", "is_ecalDrivenSeed", "is_trackerDrivenSeed", "gen_matchedSeedOkCharge"] 

int_2D_array_list = ["gen_hit_subdetector", "gen_hit_layer"]
float_2D_array_list = ["gen_hit_pt"]

double_list = []
double_array_list_hits = []
double_array_list_res_pull = ["pt_pull", "qoverp_pull", "theta_pull", "phi_pull", "d0_pull", "z0_pull", "gen_matched_rec_pt", "gen_matched_cotth", "gen_matched_rec_cotth", "gen_matched_rec_qoverp", "gen_matched_qoverp", "gen_matched_z0", "gen_matched_rec_z0", "gen_matched_rec_phi", "gen_matched_d0", "gen_matched_rec_d0"] #, "gen_matched_theta", "gen_matched_rec_theta"]
double_array_list = ["reco_pt", "reco_eta", "reco_phi", "gen_pt", "gen_eta", "gen_phi", "fake_pt", "fake_eta", "fake_phi", "gen_matched_eta", "gen_matched_phi", "gen_matched_pt", "gen_dxy", "gen_dz", "gen_bremFraction"] + double_array_list_hits + double_array_list_res_pull

var_list = int_list + int_array_list + double_list + double_array_list + int_2D_array_list + float_2D_array_list

def var_type(var_name):
    """
    return the variable type for tree entry
    """

    if var_name in int_list: return array('i',[0])
    if var_name in double_list: return array('d', [0])
    if var_name in int_array_list: return array('i',[0]*maxpart)
    if var_name in double_array_list: return array('d', [0]*maxpart)
    if var_name in int_2D_array_list: return np.array( np.zeros( (maxpart, maxhit) ), dtype=np.int32) 
    if var_name in float_2D_array_list: return np.array( np.zeros( (maxpart, maxhit) ), dtype=np.float32 ) 

#if var_name in int_2D_array_list return array('i', [0]*maxpart) #np.array(np.zeros(maxpart,maxhit), dtype=int) )
