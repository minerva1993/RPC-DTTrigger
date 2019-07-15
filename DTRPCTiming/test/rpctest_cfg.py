import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
#options.register('Labels', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "1: 1DIGI, 4:4DIGI")
options.parseArguments()

process = cms.Process("rpcNtupler")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits

process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       'file:/afs/cern.ch/user/j/jipark/work/public/rpcDTTrigger/CMSSW_10_6_1_patch2/src/RPC-DTTrigger/RPCRecHitDTProducer/test.root'
    )
)

process.rpcntupler = cms.EDAnalyzer("DTRPCTiming",
#  rpcRecHits = cms.InputTag('rpcRecHits'),
#  simMuonRPCDigis = cms.InputTag('simMuonRPCDigis'),
#  simHitLl = cms.untracked.InputTag('g4SimHits','MuonDTHits'),
  simMuonRPCDigis = cms.InputTag('rpcRecHits'),
  dt4DSegments = cms.InputTag('dt4DSegments'),
#  label = cms.untracked.int32(options.Labels)
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('output.root')
 )

process.p = cms.Path(process.rpcntupler)
