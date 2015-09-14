from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'GGF120_8TeV_v25'
#config.General.requestName = 'GGF125_8TeV_v25'
#config.General.requestName = 'VBF125_8TeV_v25'
#config.General.requestName = 'DY_8TeV_v25'
#config.General.requestName = 'SUSYGGF120_8TeV_v25'
config.General.requestName = 'SUSYGGF130_8TeV_v25'
config.General.workArea = 'HToTauTau'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'diTauPairVertexStudy_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/GluGluToHToTauTau_M-120_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
#config.Data.inputDataset = '/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
#config.Data.inputDataset = '/VBF_HToTauTau_M-125_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
#config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-120_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-130_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1300
config.Data.outLFNDirBase = '/store/user/anayak/HToTauTau/TauTauFullyHad/VertexStudy/05June2015'
config.Data.publication = False
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_DE_DESY'
