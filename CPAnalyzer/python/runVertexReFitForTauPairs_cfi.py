import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet

def getRefittedVerticesForTauPair(process, vertexColl, tauColl, postfix="") :

    print "getRefittedVerticesForTauPair:"

    
    process.diTauPrimaryVertices = cms.EDProducer("DiTauPrimaryVertexProducer",
                                                  PFTauTagLeg1 = cms.InputTag(tauColl),
                                                  PFTauTagLeg2 = cms.InputTag(tauColl),
                                                  PVTag        = cms.InputTag(vertexColl),
                                                  beamSpot     = cms.InputTag("offlineBeamSpot"),
                                                  TrackCollectionTag = cms.InputTag("generalTracks"),
                                                  useBeamSpot = cms.bool(False)
                                                  )

    #Make pair of taus and re-fit vertex for each pair
    diTauInputList = []
    runVtxReFitSequence = cms.Sequence()
    for idxLeg1 in range(5):
        for idxLeg2 in range(5):
            
            srcLeg1 = None
            if idxLeg2 <= idxLeg1:
                continue
            moduleNameLeg1 = "%sLeg1comb%i%i%s" % (tauColl, idxLeg1, idxLeg2, postfix)
            moduleLeg1 = cms.EDProducer("SinglePatTauPicker",
                                        src = cms.InputTag(tauColl),
                                        itemNumber = cms.uint32(idxLeg1),
                                        verbose = cms.untracked.bool(False)
                                        )
            setattr(process, moduleNameLeg1, moduleLeg1)
            runVtxReFitSequence += moduleLeg1
            srcLeg1 = moduleNameLeg1
            
            moduleNameLeg2 = "%sLeg2comb%i%i%s" % (tauColl, idxLeg1, idxLeg2, postfix)
            moduleLeg2 = cms.EDProducer("SinglePatTauPicker",
                                        src = cms.InputTag(tauColl),
                                        itemNumber = cms.uint32(idxLeg2),
                                        verbose = cms.untracked.bool(False)
                                        )
            setattr(process, moduleNameLeg2, moduleLeg2)
            runVtxReFitSequence += moduleLeg2
            srcLeg2 = moduleNameLeg2

            moduleNameVtxFit = "diTauPrimaryVertexComb%i%i%s" % (idxLeg1, idxLeg2 ,postfix)
            moduleVtxReFit = process.diTauPrimaryVertices.clone(
                PFTauTagLeg1 = cms.InputTag(srcLeg1),
                PFTauTagLeg2 = cms.InputTag(srcLeg2)
                )
            setattr(process, moduleNameVtxFit, moduleVtxReFit)
            runVtxReFitSequence += moduleVtxReFit
    runVtxReFitSequenceName = "runVtxReFitSequence%s" % postfix
    setattr(process, runVtxReFitSequenceName, runVtxReFitSequence)
    return runVtxReFitSequence
            
