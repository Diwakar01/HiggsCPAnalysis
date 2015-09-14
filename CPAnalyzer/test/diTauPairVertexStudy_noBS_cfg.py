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
runOnEmbed  = False
embedType   = "RhEmbedMuTau" #"PfEmbed" or "RhEmbedMuTau" or "RhEmbedElTau"

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

################### jet sequence ####################

################### jet sequence ####################
process.load('RecoJets.Configuration.RecoPFJets_cff')

################### bTag ##############################

process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
if runOnEmbed:
    #process.load('RecoBTag.Configuration.RecoBTag_cff')
    #process.load('RecoJets.JetAssociationProducers.ak5JTA_cff')
    process.ak5JetTracksAssociatorAtVertex.jets   = cms.InputTag("ak5PFJets")
    if "PfEmbed" in embedType:
        process.ak5JetTracksAssociatorAtVertex.tracks = cms.InputTag("tmfTracks")
        
## Plus, add this to your path:
#process.ak5JetTracksAssociatorAtVertex*process.btagging
        
################### vertex sequence ####################
        
process.selectedPrimaryVertices = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid && ndof >= 4 && z > -24 && z < 24 && position.Rho < 2."),
    filter = cms.bool(False)
    )

process.primaryVertexCounter = cms.EDFilter(
    "VertexCountFilter",
    src = cms.InputTag('selectedPrimaryVertices'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

################### pat specific ####################
################### pat specific ####################

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

from PhysicsTools.PatAlgos.tools.coreTools import *
if not runOnMC:
    removeMCMatching(process,["All"])

removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=None)
removeCleaning(process,
               outputInProcess=None)

restrictInputToAOD(process, ['All'])

#from LLRAnalysis.Utilities.customizePAT  import *
#addSelectedPFlowParticle(process)

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, postfix)

from PhysicsTools.PatAlgos.tools.jetTools import *

switchJetCollection(process,cms.InputTag('ak5PFJets'),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5PF', ['L2Relative', 'L3Absolute',]),
                    doType1MET   = False,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True,
                    jetIdLabel   = 'ak5'
                    )

JEClevels = cms.vstring(['L2Relative', 'L3Absolute'])
if runOnMC:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
    
process.patJetCorrFactors.levels = JEClevels
process.patJetCorrFactors.rho    = cms.InputTag('kt6PFJets','rho')
process.patJetCorrFactors.useRho = True
process.patJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors"))
'''
if runOnMC:
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.Configuration.RecoGenJets_cff")
    process.genJetsNoNu = cms.Sequence(process.genParticlesForJetsNoNu*
                                       process.ak5GenJetsNoNu)
    process.patDefaultSequence.replace(process.patJetGenJetMatch,
                                       process.genJetsNoNu*
                                       process.patJetGenJetMatch)
    process.patJetGenJetMatch.matched = cms.InputTag("ak5GenJetsNoNu")

'''
#################### tau sequence #######################

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.tauMatch.maxDeltaR                = 0.15
process.tauMatch.resolveAmbiguities       = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR          = 0.15
process.tauGenJetMatch.maxDPtRel          = 999

if "PfEmbed" in embedType:
    process.hpsPFTauPrimaryVertexProducer.TrackCollectionTag = cms.InputTag("tmfTracks")

######################## pat::jet ################################

getattr(process,"selectedPatJets").cut = cms.string('pt>10 && abs(eta)<5.0')


######################## pat::trigger ############################

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'
if runOnEmbed:
    process.patTriggerEvent.processName = 'HLT'
    
if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'
    if runOnEmbed:
        process.patTrigger.processName = 'HLT'

######################## Tau-Pair Selectiong ###############################

process.tauPtEta  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTaus"),
    cut = cms.string("pt>10 && abs(eta)<2.3"),
    filter = cms.bool(False)
    )

process.tauPtEtaCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEta"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaID  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTaus"),
    cut = cms.string(process.tauPtEta.cut.value()+
                     " && tauID('decayModeFinding')>0.5"),
                     #" && (tauID('decayModeFindingNewDMs')>0.5 | tauID('decayModeFindingOldDMs')>0.5)"), ## IN newTauID
                     #" && abs(userFloat('dzWrtPV'))<0.2"),
    filter = cms.bool(False)
    )

process.tauPtEtaIDCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEtaID"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaIDAgMu  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTaus"),
    cut = cms.string(process.tauPtEtaID.cut.value()+
                     " && tauID('againstMuonTight')>0.5"),
    filter = cms.bool(False)
    )

process.tauPtEtaIDAgMuCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEtaIDAgMu"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    )
process.tauPtEtaIDAgMuAgElec  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTaus"),
    cut = cms.string(process.tauPtEtaIDAgMu.cut.value()+
                     " && tauID('againstElectronLoose')>0.5"),
    filter = cms.bool(False)
    )

process.tauPtEtaIDAgMuAgElecCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEtaIDAgMuAgElec"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    )

process.tauPtEtaIDAgMuAgElecRelIso  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTaus"),
    cut = cms.string(process.tauPtEtaIDAgMuAgElec.cut.value()+
                     " && tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') > 0.5 "),
    filter = cms.bool(False)
    )

process.tauPtEtaIDAgMuAgElecRelIsoCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("tauPtEtaIDAgMuAgElecRelIso"),
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999),
    )

#######################Run vertex fit for tau-pair ##################
#from MyRootMaker.MyRootMaker.runVertexReFitForTauPairs_cfi import *
#getRefittedVerticesForTauPair(process, "selectedPrimaryVertices", "tauPtEtaIDAgMuAgElecRelIso")
#
from RecoTauTag.RecoTau.PFTauPrimaryVertexProducer_cfi import *
process.pfTauPrimaryVertices = PFTauPrimaryVertexProducer.clone()
process.pfTauPrimaryVertices.PVTag = cms.InputTag("offlinePrimaryVerticesWithBS")
process.pfTauPrimaryVertices.beamSpot = cms.InputTag("offlineBeamSpot")
process.pfTauPrimaryVertices.useSelectedTaus = cms.bool(True)
process.pfTauPrimaryVertices.discriminators = cms.VPSet(cms.PSet(discriminator = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),selectionCut = cms.double(0.5)),cms.PSet(discriminator = cms.InputTag('hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits'),selectionCut = cms.double(0.5)),cms.PSet(discriminator = cms.InputTag('hpsPFTauDiscriminationByTightElectronRejection'),selectionCut = cms.double(0.5)),cms.PSet(discriminator = cms.InputTag('hpsPFTauDiscriminationByTightMuonRejection'),selectionCut = cms.double(0.5)))
#######################Run vertex fit for tau-pair and store tree ############
process.vertexAnalyzer = cms.EDAnalyzer("VertexAnalyzer",
                                        tauTag       = cms.InputTag("tauPtEtaIDAgMuAgElecRelIso"),
                                        #PVTag        = cms.InputTag("selectedPrimaryVertices"),
                                        PVTag        = cms.InputTag("offlinePrimaryVertices"),
                                        #PVTag        = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                        beamSpot     = cms.InputTag("offlineBeamSpot"),
                                        TrackCollectionTag = cms.InputTag("generalTracks"),
                                        JetTag       = cms.InputTag("patJets"), 
                                        PFJetId = cms.PSet( version = cms.string("FIRSTDATA"), quality = cms.string("LOOSE") ),
                                        genParticles = cms.InputTag("genParticles"),
                                        useBeamSpot = cms.bool(False),
                                        tauPVTag    = cms.InputTag("pfTauPrimaryVertices","PFTauPrimaryVertices")
                                        )
####################### final sequences ##############################

process.atLeastOneGoodVertexSequence = cms.Sequence(
    process.selectedPrimaryVertices*
    process.primaryVertexCounter
    )

process.tauLegSequence = cms.Sequence(
    (process.tauPtEta*process.tauPtEtaCounter) *
    (process.tauPtEtaID*process.tauPtEtaIDCounter) *
    (process.tauPtEtaIDAgMu*process.tauPtEtaIDAgMuCounter)*
    (process.tauPtEtaIDAgMuAgElec*process.tauPtEtaIDAgMuAgElecCounter)*
    process.tauPtEtaIDAgMuAgElecRelIso*process.tauPtEtaIDAgMuAgElecRelIsoCounter
    )

########################## paths ###############################
process.commonOfflineSequence = cms.Sequence(
    #process.atLeastOneGoodVertexSequence* #This sequence crashes, I don't understand the reason 
    process.PFTau*process.pfTauPrimaryVertices*
    (process.ak5JetTracksAssociatorAtVertex*process.btagging)*
    process.patDefaultSequence
    )

process.treeDiTau = cms.Sequence(
    process.commonOfflineSequence*
    process.tauLegSequence*
    process.vertexAnalyzer
    ##process.runVtxReFitSequence
    #+process.printTree1
    )

#process.p = cms.Path(process.printEventContent+process.skim)
process.pDiTau = cms.Path(process.treeDiTau)

########################## output ###############################
#
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.out.outputCommands = cms.untracked.vstring('drop *',
#                                                   *patEventContentNoCleaning )
#process.out.outputCommands.extend( cms.vstring(
#
#    'keep *_TriggerResults_*_*',
#    'keep *_hltTriggerSummaryAOD_*_*',
#    'keep recoGenParticles_genParticles*_*_*',
#    'keep *_patTriggerEvent_*_*',
#    'keep *_patTrigger_*_*',
#    'keep *_offlinePrimaryVertices_*_*',
#    'keep *_selectedPrimaryVertices_*_*',
#    'keep *_*_*_PAT'
#    )
# )
#
#process.out.SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('pDiTau')
#    )
#
#process.out.fileName = cms.untracked.string('patTuples_DiTauStream.root')
#
#process.outpath = cms.EndPath(process.out)
####################################################################
del process.outpath

############## tfileservice ##############################
process.load("PhysicsTools.UtilAlgos.TFileService_cfi")
process.TFileService.fileName = cms.string('vtxanalysis.root')


