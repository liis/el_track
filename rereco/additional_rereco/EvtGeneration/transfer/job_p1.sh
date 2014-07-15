#!/bin/bash

export SCRAM_ARCH="slc5_amd64_gcc462"
source $VO_CMS_SW_DIR/cmsset_default.sh
eval `scramv1 runtime -sh`

data_replica.py --discovery --to-site T3_CH_PSI filelist_test.txt /store/user/liis/GSF_tracking_samples/test_Zee_62_RECO_2