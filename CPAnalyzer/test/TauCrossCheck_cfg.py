###############################################################################
## Pat-tuple for taul+tau analysis
##
##
##
###############################################################################

from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")

postfix     = "PFlow"
runOnMC     = True

#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START53_V23::All')
else:
    process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All')
    
    
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #"/store/mc/Summer12_DR53X/GluGluToHToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/52ACA0C7-53E9-E111-A786-002618943919.root"
    #"/store/mc/Summer12_DR53X/VBF_HToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/7CF1CED7-A6ED-E111-BB07-002481ACF3B0.root"
    "/store/mc/Summer12_DR53X/VBF_HToTauTau_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/8AE4E054-85ED-E111-81F5-1CC1DE0570A0.root"
    #"file:pickevents_merged.root"
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
    inputCommands=cms.untracked.vstring(
    'keep *'
    #'drop *PFTau*_*_*_*'
    )
 )

#process.add_(cms.Service("PrintLoadingPlugins"))

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:182771'
#    )
################### event content ##################

process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

################### filters log  ####################
                
################### gen listing  ####################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer(
        "ParticleListDrawer",
            src = cms.InputTag("genParticles"),
            maxEventsToPrint  = cms.untracked.int32(1)
            )

from PhysicsTools.PatAlgos.tools.coreTools import *
removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=None)
removeCleaning(process,
               outputInProcess=None)

restrictInputToAOD(process, ['All'])

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.tauMatch.maxDeltaR                = 0.15
process.tauMatch.resolveAmbiguities       = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR          = 0.15
process.tauGenJetMatch.maxDPtRel          = 999

#######################Run vertex fit for tau-pair and store tree ############
process.xcAnalyzer = cms.EDAnalyzer("TauIDCrossCheck",
                                    tauTag       = cms.InputTag("selectedPatTaus")
                                    )
####################### final sequences ##############################

process.commonOfflineSequence = cms.Sequence(
    #process.atLeastOneGoodVertexSequence* #This sequence crashes, I don't understand the reason 
    process.PFTau*
    process.patDefaultSequence
    )

process.treeTau = cms.Sequence(
    process.commonOfflineSequence*
    process.xcAnalyzer
    +process.printTree1
    )

#process.p = cms.Path(process.printEventContent+process.skim)
process.pTau = cms.Path(process.treeTau)

####################################################################
del process.outpath

############## tfileservice ##############################
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName = cms.string('xchecktau.root')


