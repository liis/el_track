#!/bin/bash

INDIR=$1 # the dir, where the input trees are located
#SUBMIT_BATCH=1

cd $INDIR
TREES=`ls trackValTree_Zee_maxCand_*.root`
echo processing list of input trees $TREES
cd -
for TREE in $TREES
do
  INFILE=$INDIR$TREE
  echo Processing file $INFILE
  sh make_histograms.sh $INFILE
#  qsub make_histograms.sh $INFILE
done