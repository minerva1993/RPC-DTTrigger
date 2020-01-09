import FWCore.ParameterSet.Config as cms

process = cms.Process("DTRPC")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/Mu_FlatPt2to100-pythia8-gun__PU200_106X_upgrade2023_realistic_v3-v2__FFCFF986-ED0B-B74F-B253-C511D19B8249.root'
    )
)

process.DTRPCTimingUpdate = cms.EDProducer('DTRPCTimingUpdate',
  src = cms.InputTag('simMuonRPCDigis')

)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring('drop *',
    'keep *_simMuonRPCDigis__*'
  )
)


process.p = cms.Path(process.DTRPCTimingUpdate)

process.e = cms.EndPath(process.out)
