#!/bin/bash

#abort on error
set -e

#print all lines
set -x

#export SCRAM_ARCH=slc6_amd64_gcc481 # do when running at EE
scram project -n CMSSW CMSSW CMSSW_7_1_0_pre7

cd CMSSW/src

#run cmsenv
echo "running cmsenv in  "`pwd`
eval `scramv1 runtime -sh`

#https://twiki.cern.ch/twiki/bin/view/Sandbox/ASmallGitCMSSWTutorial
echo "Getting additional packages"
git cms-addpkg TrackingTools/KalmanUpdators
git cms-addpkg  RecoTracker/CkfPattern
cp -r ../../additional_rereco . #copy additional packages

scramv1 b -j 9

cd $CMSSW_BASE/../..
echo "Setting up CMSSW submodule"
