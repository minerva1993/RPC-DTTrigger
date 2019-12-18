import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register('Labels', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "1: 1DIGI, 4:4DIGI")
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
       #'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/CMSSW_10_6_1_patch2/src/RPC-DTTrigger/RPCRecHitDTProducer/test_1.root',
       #'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/CMSSW_10_6_1_patch2/src/RPC-DTTrigger/RPCRecHitDTProducer/test_2.root',
       #'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/CMSSW_10_6_1_patch2/src/RPC-DTTrigger/RPCRecHitDTProducer/test_3.root',
       #'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/CMSSW_10_6_1_patch2/src/RPC-DTTrigger/RPCRecHitDTProducer/test_4.root',
       #'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/CMSSW_10_6_1_patch2/src/RPC-DTTrigger/RPCRecHitDTProducer/test_5.root',
      'file:/afs/cern.ch/work/j/jipark/public/rpcDTTrigger/Mu_FlatPt2to100-pythia8-gun__PU200_106X_upgrade2023_realistic_v3-v2__FFCFF986-ED0B-B74F-B253-C511D19B8249.root'
    )
)

# run the local reco on the flight 
process.load('RecoLocalMuon.Configuration.RecoLocalMuon_cff')
process.rpcRecHits.rpcDigiLabel = cms.InputTag('simMuonRPCDigis')
process.dt1DRecHits.dtDigiLabel = "simMuonDTDigis"

process.p = cms.Path(process.rpcRecHits * process.dt1DRecHits * process.dt4DSegments)

#process.tmpOut = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('localMuonReco.root'),
#    outputCommands = cms.untracked.vstring(#'keep *',
#                              'drop *',
#                              'keep *_dt4DSegments_*_*')
#)
#process.tmpOutPath = cms.EndPath(process.tmpOut)
#
#process.rpcObjectTask = cms.Task()
#process.schedule = cms.Schedule(
#    process.p, process.tmpOutPath,
#    tasks = process.rpcObjectTask
#)

process.rpcntupler = cms.EDAnalyzer("DTRPCTiming",
  DTsimHitLabel = cms.untracked.InputTag('g4SimHits','MuonDTHits'),
  simMuonRPCDigis = cms.InputTag('rpcRecHits'),
  rpcSimLinkLabel = cms.InputTag('simMuonRPCDigis','RPCDigiSimLink'),
  dt4DSegments = cms.InputTag('dt4DSegments','',''),
  label = cms.untracked.int32(options.Labels)#Needs for sim val. calculation
)

process.rpcntupler.dt4DSegments = cms.InputTag('dt4DSegments','','rpcNtupler')

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('output.root')
 )

process.p += process.rpcntupler
