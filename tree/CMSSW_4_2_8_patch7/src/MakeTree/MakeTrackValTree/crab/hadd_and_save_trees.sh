#!/bin/bash                                                                                                                         
INDIR=$1
OUTDIR=output_crab #specify directory, where to put the crab output *.root files
if [ -z "$INDIR" ]; then echo "Usage: $0 INDIR"; exit 1; fi
cd $INDIR
INDIR=`pwd`

CRABDIRNAME=crab_0_
CRABDIRS=`find $INDIR -name "$CRABDIRNAME*"`

for CRABDIR in $CRABDIRS; do
echo "Getting histograms from:" $CRABDIR
#cd $CRABDIR/res
CRABDIRNAME=`basename $CRABDIR`
OUTFILENAME="trackValTree_"${CRABDIRNAME:7}".root" #omit the first part of the crabdir name
OUTPATH=$OUTDIR"/"$OUTFILENAME
#echo $OUTPATH
hadd $OUTPATH $CRABDIR/res/*.root
echo "Wrote histograms to:"$OUTFILENAME
#cd $INDIR
done
