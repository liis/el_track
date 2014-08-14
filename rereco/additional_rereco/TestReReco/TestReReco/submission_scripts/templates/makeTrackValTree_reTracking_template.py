import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger = cms.Service("MessageLogger",
                                         default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
                                    )

# source
# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles,
                     secondaryFileNames=secFiles)

infilelist = INFILELIST
secfilelist = SECFILELIST
readFiles.extend( infilelist )
secFiles.extend( secfilelist )

process.source = source


process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')
process.GlobalTag.globaltag = 'PRE_STA71_V2::All'

process.load('Configuration/StandardSequences/Services_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff")

process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Reconstruction_cff")



from SimGeneral.MixingModule.trackingTruthProducer_cfi import *

#process.TrajectoryBuilderForElectrons.estimator = cms.string('Chi4A')

maxCandDefault = 5
maxChi2Default = 2000
nSigmaDefault = 3.0

maxCand = MAXCAND
maxChi2 = MAXCHI2
nSigma = NSIGMA

########################################################################
# to change parameters  as in slides from A.Tropiano


process.TrajectoryBuilderForElectrons.maxCand = cms.int32( maxCand )
process.ElectronChi2.MaxChi2 = cms.double( maxChi2 )
process.ElectronChi2.nSigma = cms.double( nSigma )

########################################################################

process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")
# sequence for re-running gsfTracking over RECO
process.myGsfReco = cms.Sequence(
    process.siPixelRecHits
    *process.siStripMatchedRecHits #make local hits
    *process.MeasurementTrackerEvent #ADD
    *process.siPixelClusterShapeCache # needed to add when moving from CMSSW_7_1_0_pre5 to pre7
    *process.iterTracking
    
    *process.electronSeedsSeq #ADD
    *process.electronSeeds    #produced merged collection of TkDriven and Ecaldriven seeds
    *process.electronCkfTrackCandidates*process.electronGsfTracks #run electron tracking
    #    *process.printEventContent
)

outdir = "out_tests/"
outfilename = OUTFILENAME
print "Writing output to file: " + outfilename

process.TFileService = cms.Service("TFileService", # if save
                                   fileName = cms.string(outfilename)
                                   )



process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string("reco") #Quality_SimToReco = shared hits/#reco(or #sim)
process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.75)
#process.TrackAssociatorByHits.AbsoluteNumberOfHits = cms.bool(True)
#---------------- high purity selection of reco::Tracks---------------
process.load("PhysicsTools.RecoAlgos.recoTrackSelector_cfi")
process.cutsRecoTracksHp = process.recoTrackSelector.clone()
process.cutsRecoTracksHp.quality = cms.vstring("highPurity")
process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)


process.ValidationSelectors = cms.Sequence(
    process.cutsRecoTracksHp
    )

#--------------------------- tree maker --------------------------
process.load("MakeTree.MakeTrackValTree.maketrackvaltree_cfi") # for writing output to a flat tree

process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.preValidation = cms.Sequence(
    process.simHitTPAssocProducer
    )

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# paths
process.p = cms.Path(
    process.myGsfReco 
    *process.ValidationSelectors
    *process.preValidation
#    *process.printEventContent 
    *process.trackValTreeMaker
    )


process.schedule = cms.Schedule(
    process.p
    )
