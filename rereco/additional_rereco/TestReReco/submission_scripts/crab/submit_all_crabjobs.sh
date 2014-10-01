#source /opt/software/CRAB_2_8_6_slurmlocal/crab.sh 
voms-proxy-init -voms cms

INDIR="input_crab" #directory with all crab.cfg files
cd $INDIR

CRABJOBS=`ls crab_*.cfg`
echo $CRABJOBS

for CRABJOB in $CRABJOBS
do
  echo submitting crabjob $CRABJOB
  crab -create -submit -cfg $CRABJOB
done