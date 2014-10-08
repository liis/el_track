echo Sync TestReReco...
rsync -r --exclude submission_scripts/crab/input_crab --exclude submission_scripts/crab/output_crab ../CMSSW/src/TestReReco/ TestReReco
echo Sync MakeTree...
rsync -r ../CMSSW/src/MakeTree/ MakeTree
echo Sync EvtGeneration...
rsync -r ../CMSSW/src/EvtGeneration/ EvtGeneration
echo Sync RawTests...
rsync -r ../CMSSW/src/RawTests/ RawTests
echo ...done