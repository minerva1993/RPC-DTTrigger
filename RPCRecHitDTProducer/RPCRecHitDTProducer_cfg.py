import FWCore.ParameterSet.Config as cms

from FWCore.PythonUtilities.LumiList import LumiList
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")

process = cms.Process("REClUSTERIZATION")

process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#for 93X sample
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "93X_upgrade2023_realistic_v5"

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

#RPCRecHit
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import rpcRecHits
process.rpcRecHits = rpcRecHits
process.rpcRecHits.rpcDigiLabel = cms.InputTag("simMuonRPCDigis")

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/eos/cms/store/group/dpg_dt/comm_dt/TriggerSimulation/SamplesReco/SingleMu_FlatPt-2to100/Version_10_5_0/SimRECO_1.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.p = cms.Path( process.rpcRecHits )

# Output
process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        "drop *",
        #, "keep *_simMuonRPCDigis__*",
        "keep *_g4SimHits_MuonDTHits_*",
        "keep *_g4SimHits_MuonRPCHits_*",
        "keep *_rpcRecHits__*",
        "keep *_simMuonRPCDigis_RPCDigiSimLink_*",
        "keep *_dt4DSegments_*_*",
        "keep *_simDtTriggerPrimitiveDigis__*"),
    fileName = cms.untracked.string("test.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
)


process.e = cms.EndPath(process.out)
