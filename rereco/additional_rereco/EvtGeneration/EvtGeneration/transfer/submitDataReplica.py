#!/usr/bin/env python


import commands
import re
import os

def submit( list ):

    input = open( list+'.txt' )
    lines = input.readlines()
    
    counter = 0
    for line in lines:
        counter += 1
        
        newinputname = 'job_'+list+'_p'+str(counter)+'.txt'
        newinput = open( newinputname ,'w')
        newinput.write( line )
        
        scriptName = 'job_p'+str(counter)+'.sh'
        f = open( scriptName,'w')
        
        f.write('#!/bin/bash\n\n')
#        f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/test/\n')
        f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
        f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
        f.write('eval `scramv1 runtime -sh`\n')
        f.write('\n')
#        f.write('data_replica.py --from-site T2_EE_Estonia --to-site T3_CH_PSI '+newinputname+'  /store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V4\n')

        f.write('data_replica.py --discovery --to-site T3_CH_PSI filelist_test.txt /store/user/liis/GSF_tracking_samples/test_Zee_62_RECO_2')
        f.close()
        os.system('chmod +x '+scriptName)
        
        submitToQueue = 'qsub -V -cwd -l h_vmem=1G -q all.q -N job_p'+str(counter)+' '+scriptName
        print submitToQueue
        os.system(submitToQueue)
        
        input.close()
        print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
        #############################
        
submit( 'filelist_test' )
