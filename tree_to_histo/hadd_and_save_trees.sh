#!/bin/bash                                                                                                                         

if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT to use hadd!"
 exit
fi

INDIR=$1 #directory of the crab output directories with root-files
OUTDIR=$2 #specify directory, where to put the crab output *.root files
#OUTDIR=output_crab 

OVERWRITE=1
GET_FROM_STORAGE=1

if [ -z "$INDIR" ]; then echo "Usage: INDIR, where the crabdirs are has to be given as a first argument"; exit 1; fi
if [ -z "$OUTDIR" ]; then echo "Specify OUTDIR, where the output root trees are saved as a second argument after INDIR"; exit 1; fi

cd $INDIR
INDIR=`pwd`

echo "looking for crab output directories in "$INDIR
CRABDIRNAME=FlatPt
CRABDIRS=`find $INDIR -name "*$CRABDIRNAME*"`

cd -

for CRABDIR in $CRABDIRS; do
    echo "Getting histograms from:" $CRABDIR
    CRABDIRNAME=`basename $CRABDIR`
    OUTFILENAME="trackValTree_"${CRABDIRNAME}".root"

    OUTPATH=$OUTDIR"/"$OUTFILENAME

    if [ GET_FROM_STORAGE == 0 ]; then
	INPATH=$CRABDIR/res
    else
	INPATH=$CRABDIR
    fi
    
    if [ -e $OUTPATH ] && [ ! $OVERWRITE = 1 ]; then
	echo "File exists, skipping"
    else
	echo "Writing merged histograms to:"$OUTPATH 
	hadd -f $OUTPATH $INPATH/*.root
   fi

done
