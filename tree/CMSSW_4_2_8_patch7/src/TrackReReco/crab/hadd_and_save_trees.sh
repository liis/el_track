
#!/bin/bash                                                                                                                         
INDIR=$1 #directory, where the crab_dirs are
OUTDIR=$2
OUTDIR=output_crab #specify directory, where to put the crab output *.root files
if [ -z "$INDIR" ]; then echo "Usage: $0 INDIR"; exit 1; fi
cd $INDIR
INDIR=`pwd`

CRABDIRNAME=crab_0_
CRABDIRS=`find $INDIR -name "$CRABDIRNAME*"`

cd -

for CRABDIR in $CRABDIRS; do
echo "Getting histograms from:" $CRABDIR
CRABDIRNAME=`basename $CRABDIR`
OUTFILENAME="trackValTree_"${CRABDIRNAME:7}".root" #omit the first part of the crabdir name

OUTPATH=$OUTDIR"/"$OUTFILENAME
echo hadd -f $OUTPATH $CRABDIR/res/*.root
hadd -f $OUTPATH $CRABDIR/res/*.root
echo "Wrote histograms to:"$OUTFILENAME
done

