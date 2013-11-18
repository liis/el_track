import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:aod.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('trackValTree.root')
                                   )


process.load("MakeTree.MakeTrackValTree.maketrackvaltree_cfi")

#process.demo = cms.EDAnalyzer('MakeTrackValTree'
#)



process.p = cms.Path(
    process.demo
    )

