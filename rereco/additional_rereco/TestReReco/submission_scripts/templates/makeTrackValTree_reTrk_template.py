import FWCore.ParameterSet.Config as cms

process = cms.Process("reGsfTracking")

runDigi = 1

# message logger
process.MessageLogger = cms.Service("MessageLogger",
     default = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)

#if not runDigi:
#    readFiles.extend( ['file:/hdfs/cms/store/user/liis/El_GSF_studies/Pt10/step2_9_1_DD1.root'])
#else:

infilelist = INFILELIST
readFiles.extend( infilelist )


process.source = source

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V11::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 50

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

process.load("SimG4Core.Application.g4SimHits_cfi")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# Reconstruction geometry service
#process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
# geometry (Only Ecal)
#process.load("Geometry.EcalCommonData.EcalOnly_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.EcalTrigTowerConstituents_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")

# use trivial ESProducer for tests
#process.load("CalibCalorimetry.EcalTrivialCondModules.EcalTrivialCondRetriever_cfi")

# ECAL 'slow' digitization sequence
process.load("SimCalorimetry.Configuration.ecalDigiSequence_cff")

# sequence for re-running gsfTracking over RECO
process.myGsfReco_base = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits* #make local hits
    process.iterTracking*process.trackCollectionMerging*  #CTF iterative tracking
    process.newCombinedSeeds* #together with EG Cluasters, input for ecalDriven seeds
    process.electronSeeds*    #produced merged collection of TkDriven and Ecaldriven seeds
    process.electronCkfTrackCandidates*process.electronGsfTracks #run electron tracking
    )
###############


process.additional_digi = cms.Sequence( # if running on digi input
        process.siPixelDigis*
        process.siStripDigis*
        process.ecalDigis*
        process.ecalPreshowerDigis*
        process.hcalDigis*
        process.muonDTDigis*
        process.muonCSCDigis*
        process.muonRPCDigis*
        process.castorDigis
        )


process.myGsfReco_forDigi = cms.Sequence(
        process.additional_digi *
        process.localreco *
        process.globalreco *
#        process.highlevelreco *
        process.myGsfReco_base 
    )


if runDigi:
    process.myGsfReco = process.myGsfReco_forDigi
else:
    process.myGsfReco = process.myGsfReco_base


outfilename = OUTFILENAME
process.TFileService = cms.Service("TFileService", # if save 
                                   fileName = cms.string(outfilename)
                                   )

#####################################

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string("reco") #Quality_SimToReco = shared hits/#reco(or #sim)
process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.75)
#process.TrackAssociatorByHits.AbsoluteNumberOfHits = cms.bool(True)
process.load("MakeTree.MakeTrackValTree.maketrackvaltree_cfi") # for writing output to a flat tree

#---------------- high purity selection of reco::Tracks---------------
process.load("PhysicsTools.RecoAlgos.recoTrackSelector_cfi")
process.cutsRecoTracksHp = process.recoTrackSelector.clone()
process.cutsRecoTracksHp.quality = cms.vstring("highPurity")
process.cutsRecoTracksHp.minAbsEta = cms.double(0.0)
process.cutsRecoTracksHp.maxAbsEta = cms.double(2.5)

process.ValidationSelectors = cms.Sequence(
    process.cutsRecoTracksHp
    )

process.elGsfTracksWithQuality = cms.EDProducer("AnalyticalGsfTrackSelector",
                                                max_d0 = cms.double(100.0),
                                                minNumber3DLayers = cms.uint32(3),
                                                applyAbsCutsIfNoPV = cms.bool(False),
                                                qualityBit = cms.string('highPurity'),
                                                minNumberLayers = cms.uint32(3),
                                                chi2n_par = cms.double(0.7),
                                                nSigmaZ = cms.double(3.0),
                                                dz_par2 = cms.vdouble(0.4, 4.0),
                                                applyAdaptedPVCuts = cms.bool(True),
                                                dz_par1 = cms.vdouble(0.35, 4.0),
                                                copyTrajectories = cms.untracked.bool(True),
                                                vtxNumber = cms.int32(-1),
                                                keepAllTracks = cms.bool(True),
                                                maxNumberLostLayers = cms.uint32(2),
                                                beamspot = cms.InputTag("offlineBeamSpot"),
                                                copyExtras = cms.untracked.bool(True),
                                                vertexCut = cms.string('ndof>=2&!isFake'),
                                                max_z0 = cms.double(100.0),
                                                src = cms.InputTag("zeroStepWithTightQuality"),
                                                vertices = cms.InputTag("pixelVertices"),
                                                d0_par2 = cms.vdouble(0.4, 4.0),
                                                d0_par1 = cms.vdouble(0.3, 4.0),
                                                res_par = cms.vdouble(0.003, 0.001)
                                                )
process.elGsfTracksWithQuality.src = cms.InputTag("electronGsfTracks")

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# paths
process.p = cms.Path(
    process.myGsfReco *
    process.ValidationSelectors *
    process.elGsfTracksWithQuality *
#    process.printEventContent *
    process.TrackValTreeMaker
)

process.schedule = cms.Schedule(
      process.p
      )


