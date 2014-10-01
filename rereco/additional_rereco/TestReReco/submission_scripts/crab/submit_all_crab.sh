INDIR=$1

CRABFILES=`ls $INDIR/crab_maxCand_*.cfg`

for CRABFILE in $CRABFILES
do
    crab -create -submit -cfg $CRABFILE
done