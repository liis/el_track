
#!/bin/bash                                                                                                                         

if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT to use hadd!"
 exit
fi

INDIR=$1 #directory, where the crab_dirs are
OUTDIR=$2 #specify directory, where to put the crab output *.root files
OUTDIR=output_crab 

GET_FROM_STORAGE=0

if [ -z "$INDIR" ]; then echo "Usage: INDIR, where the crabdirs are has to be given as a first argument"; exit 1; fi
cd $INDIR
INDIR=`pwd`

echo $INDIR

echo "looking for crab_0_* directories in "$INDIR
CRABDIRNAME=crab_0_
CRABDIRS=`find $INDIR -name "$CRABDIRNAME*"`

echo $CRABDIRS
cd -

for CRABDIR in $CRABDIRS; do
    echo "Getting histograms from:" $CRABDIR
    CRABDIRNAME=`basename $CRABDIR`
    OUTFILENAME="trackValTree_"${CRABDIRNAME:7}".root" #omit the first part of the crabdir name

    OUTPATH=$OUTDIR"/"$OUTFILENAME

    if [ $GET_FROM_STORAGE == 0 ]; then
	INPATH=$CRABDIR/res
    else
	INPATH=$CRABDIR
    fi
    hadd -f $OUTPATH $INPATH/*.root
    echo "Wrote merged histograms to:"$OUTPATH
done

cp $OUTDIR/TrackValTree_*.root ~/tracking/el_track/tree_to_histo/input_trees/rereco_trees/