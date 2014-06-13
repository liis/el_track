#!/bin/bash

#abort on error
set -e

#print all lines
set -x

export SCRAM_ARCH=slc6_amd64_gcc481
scram project -n CMSSW CMSSW CMSSW_7_1_0_pre8

cd CMSSW/src

#run cmsenv
echo "running cmsenv in  "`pwd`
eval `scramv1 runtime -sh`

#https://twiki.cern.ch/twiki/bin/view/Sandbox/ASmallGitCMSSWTutorial
echo "Getting additional packages"
git cms-addpkg TrackingTools/KalmanUpdators
git cms-addpkg  RecoTracker/CkfPattern

scramv1 b -j 9

cd $CMSSW_BASE/../..
echo "Setting up CMSSW submodule"
