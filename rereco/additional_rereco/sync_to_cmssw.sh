echo Sync TestReReco...
#rm -r TestReReco
rsync -r --exclude submission_scripts/crab/input_crab --exclude submission_scripts/crab/output_crab TestReReco ../CMSSW/src/
echo Sync MakeTree...
#rm -r MakeTree
rsync -r MakeTree ../CMSSW/src/
echo Sync EvtGeneration...
#rm -r EvtGeneration
rsync -r EvtGeneration ../CMSSW/src/
echo Sync RawTests...
#rm -r RawTests
rsync -r RawTests ../CMSSW/src/-

echo ...done