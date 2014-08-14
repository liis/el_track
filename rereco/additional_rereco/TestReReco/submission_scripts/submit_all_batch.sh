INDIR="./batch_jobs/input_batch/"
cd $INDIR

JOBS=`ls makeTrackValTree*.sh`

for JOB in $JOBS
do
  echo submitting $JOB
  qsub -V -cwd -q all.q $JOB
done