from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = False
config.General.transferOutputs = True
config.General.requestName = 'DTRPC_TriggerDev'

config.section_('JobType')
config.JobType.psetName = 'RPCRecHitDTProducer_cfg.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['DTRPC_TriggerDev.root']

config.section_('Data')
#config.Data.inputDataset = '/SingleMu_FlatPt-2to100/PhaseIIFall17D-L1TPU200_93X_upgrade2023_realistic_v5-v1/GEN-SIM-DIGI-RAW'
config.Data.inputDataset = '/Mu_FlatPt2to100-pythia8-gun/PhaseIITDRSpring19DR-PU200_106X_upgrade2023_realistic_v3-v2/GEN-SIM-DIGI-RAW'
config.Data.publication = False
config.Data.unitsPerJob = 180

config.section_('Site')
#config.Site.whitelist = ['T2_IT_Bari']
#config.Site.storageSite = 'T2_KR_KNU'
config.Site.storageSite = 'T3_KR_KISTI'

