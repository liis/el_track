#!/bin/bash

INDIR=$1 # the dir, where the input trees are located
#OUTDIR=$2 # directory where to store output histograms
OUTDIR="./histograms/13TeV_v1"

SUBMIT_BATCH=1

cd $INDIR
TREES=`ls trackValTree_*.root`
echo processing list of input trees $TREES
cd -

for TREE in $TREES
do
  INFILE=$INDIR$TREE
  echo Processing file $INFILE

  echo submit_batch = $SUBMIT_BATCH
  if [ "$SUBMIT_BATCH" = 1 ] 
      then echo submitting to batch: $INFILE
      qsub make_histograms.sh $INFILE $OUTDIR
  else 
      echo "Run locally: " $INFILE
      sh make_histograms.sh $INFILE $OUTDIR
  fi
    
done