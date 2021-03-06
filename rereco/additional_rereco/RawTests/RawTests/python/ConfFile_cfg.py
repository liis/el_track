import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    '/store/user/liis/GSF_tracking_samples/test_Zee_62_RECO_2/044C97C0-B675-E311-B49C-0025901D4D76.root'
#    '/store/user/liis/GSF_tracking_samples/Zee_sam_GenSimRaw/044C97C0-B675-E311-B49C-0025901D4D76.root'
    '/store/user/liis/GSF_tracking_samples/RelValZll_HLTDEBUG_PU50/00E5A57A-3D13-E411-99D5-0025905A6076.root'
    )
)

process.demo = cms.EDAnalyzer('RawTests'
)


process.p = cms.Path(process.demo)
