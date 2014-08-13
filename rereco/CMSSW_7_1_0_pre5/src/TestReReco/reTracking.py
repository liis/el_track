import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger",
                                         default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
                                    )

# source
readFiles = cms.untracked.vstring()

secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( ['file:/hdfs/cms/store/user/liis/El_GSF_studies/test/0E8FC0BF-A9BC-E311-AABA-02163E00E7DF.root'])

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

### conditions
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#process.GlobalTag.globaltag = 'PRE_MC_71_V2'
### standard includes

## Geometry and Detector Conditions (needed for a few patTuple production steps)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
#process.GlobalTag.globaltag = 'PRE_STA71_V2'


process.load('Configuration/StandardSequences/Services_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("CondCore.DBCommon.CondDBSetup_cfi")

#process.options = cms.untracked.PSet(
#    allowUnscheduled = cms.untracked.bool( True )
#)

maxCandDefault = 5
maxChi2Default = 2000
nSigmaDefault = 3.0

maxCand = 6
maxChi2 = 100
nSigma = 4

########################################################################
# to change parameters  as in slides from A.Tropiano

process.TrajectoryBuilderForElectrons.maxCand = cms.int32( maxCand )
process.ElectronChi2.MaxChi2 = cms.double( maxChi2 )
process.ElectronChi2.nSigma = cms.double( nSigma )

########################################################################

#process.GlobalPaositionSource = cms.ESSource("PoolDBESSource",
#                                            process.CondDBSetup,
#                                            # Reading from oracle (instead of Frontier) needs the following shell variable setting (in zsh):
#                                            # export CORAL_AUTH_PATH=/afs/cern.ch/cms/DB/conddb
#                                            # string connect = "oracle://cms_orcoff_int2r/CMS_COND_ALIGNMENT"
#                                            # untracked uint32 authenticationMethod = 1
#                                            toGet = cms.VPSet(cms.PSet(
#    record = cms.string('GlobalPositionRcd'),
#    tag = cms.string('IdealGeometry')
#    )),
#                                            connect = cms.string('sqlite_file:output.db')
#                                            )


# sequence for re-running gsfTracking over RECO
process.myGsfReco = cms.Sequence(
#    process.localreco
    process.siPixelRecHits*
    process.siStripMatchedRecHits* #make local hits
    process.MeasurementTrackerEvent* #ADD
    process.iterTracking*
    process.electronSeedsSeq* #ADD
#    process.trackCollectionMerging  #REMOVE
#    process.newCombinedSeeds* # contained in electronSeedSeq 
    process.electronSeeds*    #produced merged collection of TkDriven and Ecaldriven seeds
    process.electronCkfTrackCandidates*process.electronGsfTracks #run electron tracking
    )
###############
outdir = "out_tests/"
outfilename = outdir + "reGsfTracking_maxCand_" + str(maxCand) + "_MaxChi2_" + str(maxChi2) + "_nSigma_" + str(nSigma) + ".root"
print "Writing output to file: " + outfilename

process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
                                           splitLevel = cms.untracked.int32(0),
                                           eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                           outputCommands = process.RECODEBUGEventContent.outputCommands,
                                           fileName = cms.untracked.string( outfilename ),
                                           dataset = cms.untracked.PSet(
    filterName = cms.untracked.string(''),
    dataTier = cms.untracked.string('GEN-SIM-RECODEBUG')
    )
                                           )

process.RECODEBUGoutput.outputCommands.extend(cms.untracked.vstring('keep *_elGsfTracksWithQuality_*_*'))

#---------------- high purity selection of reco::Tracks---------------
process.load("PhysicsTools.RecoAlgos.recoTrackSelector_cfi")
process.cutsRecoTracksHp = process.recoTrackSelector.clone()
process.cutsRecoTracksHp.quality = cms.vstring("highPurity")
process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)


process.ValidationSelectors = cms.Sequence(
        process.cutsRecoTracksHp
            )

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# paths
process.p = cms.Path(
    process.myGsfReco *
    process.ValidationSelectors 
#    process.printEventContent 
    )

process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)

process.schedule = cms.Schedule(
    process.p
#    process.RECODEBUGoutput_step
    )
