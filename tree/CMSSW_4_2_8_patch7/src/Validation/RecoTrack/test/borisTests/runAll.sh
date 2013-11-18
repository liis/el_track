#!/bin/bash

if [ "$1" == "1" ] 
then 
# KF python files
cat MTV.template.py|sed 's/#CFG_TpBarrEff#//g' | sed 's/#CFG_TkEff#//g' > elFlatPtBarrelForEff.py ; 
cat MTV.template.py|sed 's/#CFG_TpTraEff#//g'  | sed 's/#CFG_TkEff#//g' > elFlatPtTransitionForEff.py ;
cat MTV.template.py|sed 's/#CFG_TpEndEff#//g'  | sed 's/#CFG_TkEff#//g' > elFlatPtEndcapForEff.py ;

cat MTV.template.py|sed 's/#CFG_TpFakes#//g' | sed 's/#CFG_TkBarrFake#//g' > elFlatPtBarrelForFake.py ; 
cat MTV.template.py|sed 's/#CFG_TpFakes#//g' | sed 's/#CFG_TkTraFake#//g'  > elFlatPtTransitionForFake.py ;
cat MTV.template.py|sed 's/#CFG_TpFakes#//g' | sed 's/#CFG_TkEndFake#//g'  > elFlatPtEndcapForFake.py ;

cat MTV.template.py|sed 's/#CFG_TpVsEta#//g' | sed 's/#CFG_TkVsEta#//g' > elVsEta.py ;


# GSF python files
cat elFlatPtBarrelForEff.py |sed 's/#CFG_GSF#//g' > elFlatPtBarrelForEffMode.py ;
cat elFlatPtTransitionForEff.py |sed 's/#CFG_GSF#//g' > elFlatPtTransitionForEffMode.py ;
cat elFlatPtEndcapForEff.py |sed 's/#CFG_GSF#//g' > elFlatPtEndcapForEffMode.py ;

cat elFlatPtBarrelForFake.py |sed 's/#CFG_GSF#//g' > elFlatPtBarrelForFakeMode.py ;
cat elFlatPtTransitionForFake.py |sed 's/#CFG_GSF#//g' > elFlatPtTransitionForFakeMode.py ;
cat elFlatPtEndcapForFake.py |sed 's/#CFG_GSF#//g' > elFlatPtEndcapForFakeMode.py ;

cat elVsEta.py | sed 's/#CFG_GSF#//g' > elVsEtaMode.py ; 


#
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtBarrelForEff.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtBarrelForEff#' > crab.elFlatPtBarrelForEff.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtTransitionForEff.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtTransitionForEff#' > crab.elFlatPtTransitionForEff.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtEndcapForEff.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtEndcapForEff#' > crab.elFlatPtEndcapForEff.cfg

cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtBarrelForFake.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtBarrelForFake#' > crab.elFlatPtBarrelForFake.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtTransitionForFake.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtTransitionForFake#' > crab.elFlatPtTransitionForFake.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtEndcapForFake.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtEndcapForFake#' > crab.elFlatPtEndcapForFake.cfg


cat crab.template.cfg | sed 's#CFG_MTV#elVsEta.py#' | sed 's/#CFG_KF_Pt1#//' | \
    sed 's#CFG_NAME#crab_0_elPt1Tk#' > crab.elPt1Tk.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elVsEta.py#' | sed 's/#CFG_KF_Pt10#//' | \
    sed 's#CFG_NAME#crab_0_elPt10Tk#' > crab.elPt10Tk.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elVsEta.py#' | sed 's/#CFG_KF_Pt100#//' | \
    sed 's#CFG_NAME#crab_0_elPt100Tk#' > crab.elPt100Tk.cfg



cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtBarrelForEffMode.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtBarrelForEffMode#' > crab.elFlatPtBarrelForEffMode.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtTransitionForEffMode.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtTransitionForEffMode#' > crab.elFlatPtTransitionForEffMode.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtEndcapForEffMode.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtEndcapForEffMode#' > crab.elFlatPtEndcapForEffMode.cfg

cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtBarrelForFakeMode.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtBarrelForFakeMode#' > crab.elFlatPtBarrelForFakeMode.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtTransitionForFakeMode.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtTransitionForFakeMode#' > crab.elFlatPtTransitionForFakeMode.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elFlatPtEndcapForFakeMode.py#' | sed 's/#CFG_KF_FlatPt#//' | \
    sed 's#CFG_NAME#crab_0_elFlatPtEndcapForFakeMode#' > crab.elFlatPtEndcapForFakeMode.cfg


cat crab.template.cfg | sed 's#CFG_MTV#elVsEtaMode.py#' | sed 's/#CFG_KF_Pt1#//' | \
    sed 's#CFG_NAME#crab_0_elPt1Mode#' > crab.elPt1Mode.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elVsEtaMode.py#' | sed 's/#CFG_KF_Pt10#//' | \
    sed 's#CFG_NAME#crab_0_elPt10Mode#' > crab.elPt10Mode.cfg
cat crab.template.cfg | sed 's#CFG_MTV#elVsEtaMode.py#' | sed 's/#CFG_KF_Pt100#//' | \
    sed 's#CFG_NAME#crab_0_elPt100Mode#' > crab.elPt100Mode.cfg


fi


if [ "$1" == "2" ] 
then 

crab -create -cfg crab.elPt1Tk.cfg &
crab -create -cfg crab.elPt10Tk.cfg &
crab -create -cfg crab.elPt100Tk.cfg &

crab -create -cfg crab.elFlatPtBarrelForEff.cfg &
crab -create -cfg crab.elFlatPtTransitionForEff.cfg &
crab -create -cfg crab.elFlatPtEndcapForEff.cfg &

crab -create -cfg crab.elFlatPtBarrelForFake.cfg &
crab -create -cfg crab.elFlatPtTransitionForFake.cfg &
crab -create -cfg crab.elFlatPtEndcapForFake.cfg &

crab -create -cfg crab.elPt1Mode.cfg &
crab -create -cfg crab.elPt10Mode.cfg &
crab -create -cfg crab.elPt100Mode.cfg &

crab -create -cfg crab.elFlatPtBarrelForEffMode.cfg &
crab -create -cfg crab.elFlatPtTransitionForEffMode.cfg &
crab -create -cfg crab.elFlatPtEndcapForEffMode.cfg &

crab -create -cfg crab.elFlatPtBarrelForFakeMode.cfg &
crab -create -cfg crab.elFlatPtTransitionForFakeMode.cfg &
crab -create -cfg crab.elFlatPtEndcapForFakeMode.cfg &


fi

if [ "$1" == "3" ] 
then 

crab -submit -c crab_0_elPt1Tk &
crab -submit -c crab_0_elPt10Tk &
crab -submit -c crab_0_elPt100Tk &

crab -submit -c crab_0_elFlatPtBarrelForEff &
crab -submit -c crab_0_elFlatPtTransitionForEff &
crab -submit -c crab_0_elFlatPtEndcapForEff &

crab -submit -c crab_0_elFlatPtBarrelForFake &
crab -submit -c crab_0_elFlatPtTransitionForFake &
crab -submit -c crab_0_elFlatPtEndcapForFake &

crab -submit -c crab_0_elPt1Mode &
crab -submit -c crab_0_elPt10Mode &
crab -submit -c crab_0_elPt100Mode &

crab -submit -c crab_0_elFlatPtBarrelForEffMode &
crab -submit -c crab_0_elFlatPtTransitionForEffMode &
crab -submit -c crab_0_elFlatPtEndcapForEffMode &

crab -submit -c crab_0_elFlatPtBarrelForFakeMode &
crab -submit -c crab_0_elFlatPtTransitionForFakeMode &
crab -submit -c crab_0_elFlatPtEndcapForFakeMode &


fi



if [ "$1" == "4" ] 
then 
rm crab.elPt1Tk.cfg
rm crab.elPt10Tk.cfg
rm crab.elPt100Tk.cfg
rm crab.elFlatPtTransitionForFakeMode.cfg
rm crab.elFlatPtTransitionForFake.cfg
rm crab.elFlatPtTransitionForEffMode.cfg
rm crab.elFlatPtTransitionForEff.cfg
rm crab.elFlatPtEndcapForFakeMode.cfg
rm crab.elFlatPtEndcapForFake.cfg
rm crab.elFlatPtEndcapForEffMode.cfg
rm crab.elFlatPtEndcapForEff.cfg
rm crab.elFlatPtBarrelForFakeMode.cfg
rm crab.elFlatPtBarrelForFake.cfg
rm crab.elFlatPtBarrelForEffMode.cfg
rm crab.elFlatPtBarrelForEff.cfg
rm crab.elPt1Mode.cfg
rm crab.elPt10Mode.cfg
rm crab.elPt100Mode.cfg


rm elVsEtaMode.py
rm elVsEta.py
rm elFlatPtTransitionForFakeMode.py
rm elFlatPtTransitionForFake.py
rm elFlatPtTransitionForEffMode.py
rm elFlatPtTransitionForEff.py
rm elFlatPtEndcapForFakeMode.py
rm elFlatPtEndcapForFake.py
rm elFlatPtEndcapForEffMode.py
rm elFlatPtEndcapForEff.py
rm elFlatPtBarrelForFakeMode.py
rm elFlatPtBarrelForFake.py
rm elFlatPtBarrelForEffMode.py
rm elFlatPtBarrelForEff.py

rm elVsEtaMode.pyc
rm elVsEta.pyc
rm elFlatPtTransitionForFakeMode.pyc
rm elFlatPtTransitionForFake.pyc
rm elFlatPtTransitionForEffMode.pyc
rm elFlatPtTransitionForEff.pyc
rm elFlatPtEndcapForFakeMode.pyc
rm elFlatPtEndcapForFake.pyc
rm elFlatPtEndcapForEffMode.pyc
rm elFlatPtEndcapForEff.pyc
rm elFlatPtBarrelForFakeMode.pyc
rm elFlatPtBarrelForFake.pyc
rm elFlatPtBarrelForEffMode.pyc
rm elFlatPtBarrelForEff.pyc

fi

if [ "$1" == "5" ] 
then 
for x in crab_* ; do crab -status -get -c $x & done;
fi

if [ "$1" == "6" ] 
then 
for x in crab_0_el*; do cp $x/res/multitrackvalidator_1_*.root puppa/${x/crab_0_/}.root ; done;
tar -czvf puppa.tgz puppa/

fi