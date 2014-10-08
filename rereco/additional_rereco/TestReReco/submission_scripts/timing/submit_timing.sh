INDIR=$1
cd $INDIR

JOBS=`ls makeTrackValTree*.sh`
for JOB in $JOBS
do
  echo submitting $JOB
  substring="bsub -q 1nd -J $JOB < $JOB" #for lxplus
  #  substring="qsub -V -cwd -q all.q $JOB"

  echo $substring
  bsub -q 1nd -J $JOB < $JOB

done
