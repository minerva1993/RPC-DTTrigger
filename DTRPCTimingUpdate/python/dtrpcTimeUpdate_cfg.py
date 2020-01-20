import FWCore.ParameterSet.Config as cms

process = cms.Process("DTRPC")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D38Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D38_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# run the local reco on the flight 
process.load('RecoLocalMuon.Configuration.RecoLocalMuon_cff')
process.rpcRecHits.rpcDigiLabel = cms.InputTag('simMuonRPCDigis')
process.dt1DRecHits.dtDigiLabel = "simMuonDTDigis"
process.p = cms.Path(process.rpcRecHits * process.dt1DRecHits * process.dt4DSegments)

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/Mu_FlatPt2to100-pythia8-gun__PU200_106X_upgrade2023_realistic_v3-v2__FFCFF986-ED0B-B74F-B253-C511D19B8249.root'
    )
)

process.DTRPCTimingUpdate = cms.EDProducer('DTRPCTimingUpdate',
  src = cms.InputTag('simMuonRPCDigis'),
  dt4DSegments = cms.InputTag('dt4DSegments','','')
)

process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile.root'),
  outputCommands = cms.untracked.vstring('drop *',
    'keep *_simMuonRPCDigis__*'
  )
)


process.p = cms.Path(process.DTRPCTimingUpdate)

process.e = cms.EndPath(process.out)
