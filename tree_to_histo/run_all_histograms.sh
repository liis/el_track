#! /usr/bin/Python

INDIR=$1
cd $INDIR
TREES=`ls trackValTree_*.root`
echo processing list of input trees $TREES
cd -
for TREE in $TREES
do
  python make_histograms.py --infile=$TREE 
done