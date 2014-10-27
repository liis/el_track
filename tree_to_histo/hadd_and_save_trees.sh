#!/bin/bash

if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT to use hadd!"
 exit
fi

INDIR=$1 #directory of the crab output directories with root-files
OUTDIR=$2 #specify directory, where to put the crab output *.root files
#OUTDIR=output_crab 

IS_PSI=0 #need to copy locally first

OVERWRITE=1 
GET_FROM_STORAGE=1 #look for crab_0 directories (0) or not (1)

if [ -z "$INDIR" ]; then echo "Usage: INDIR, where the crabdirs are has to be given as a first argument"; exit 1; fi
if [ -z "$OUTDIR" ]; then echo "Specify OUTDIR, where the output root trees are saved as a second argument after INDIR"; exit 1; fi

cd $INDIR
INDIR=`pwd`

echo "looking for crab output directories in "$INDIR
CRABDIRNAME=Zee_full #part of directory name to search
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
    elif [ $IS_PSI == 1 ]; then # this is impossibly slow!
	echo "Working at PSI: copy files and hadd locally"

	TMPDIR=$OUTDIR/`basename $CRABDIR`
	if [ ! -e $TMPDIR ]; then
	    mkdir $TMPDIR
	fi
	for ROOTFILE in `ls $CRABDIR`; do
	#	    echo "srmcp -2 "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="$CRABDIR/$ROOTFILE file:///$TMPDIR/`basename $ROOTFILE`"
	    srmcp -2 "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="$CRABDIR/$ROOTFILE file:///$TMPDIR/`basename $ROOTFILE`
	done
	hadd -f $OUTPATH $TMPDIR/*.root
	rm -r $TMPDIR
    else
	echo "Writing merged histograms to:"$OUTPATH 
	hadd -f $OUTPATH $INPATH/*.root
    fi

done
