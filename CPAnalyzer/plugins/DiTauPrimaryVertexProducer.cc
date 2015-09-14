#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h" 
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h" 
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/RefProd.h"

#include <memory>

using namespace reco;
using namespace edm;
using namespace std;

class DiTauPrimaryVertexProducer : public EDProducer {
public:

  explicit DiTauPrimaryVertexProducer(const edm::ParameterSet& iConfig);
  ~DiTauPrimaryVertexProducer();
  virtual void produce(edm::Event&,const edm::EventSetup&);
private:
  edm::InputTag PFTauTagLeg1_;
  edm::InputTag PFTauTagLeg2_;
  edm::InputTag PVTag_;
  edm::InputTag beamSpotTag_;
  edm::InputTag TrackCollectionTag_;
  edm::InputTag TracksTag_;
  bool useBeamSpot_;
};

DiTauPrimaryVertexProducer::DiTauPrimaryVertexProducer(const edm::ParameterSet& iConfig):
  PFTauTagLeg1_(iConfig.getParameter<edm::InputTag>("PFTauTagLeg1")),
  PFTauTagLeg2_(iConfig.getParameter<edm::InputTag>("PFTauTagLeg2")),
  PVTag_(iConfig.getParameter<edm::InputTag>("PVTag")),
  beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpot")),
  TrackCollectionTag_(iConfig.getParameter<edm::InputTag>("TrackCollectionTag")),
  useBeamSpot_(iConfig.getParameter<bool>("useBeamSpot"))
{
  produces<VertexCollection>(""); 
}

DiTauPrimaryVertexProducer::~DiTauPrimaryVertexProducer(){

}

void DiTauPrimaryVertexProducer::produce(edm::Event& iEvent,const edm::EventSetup& iSetup){
  // Obtain Collections
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  
  edm::Handle<std::vector<reco::PFTau> > TauLeg1;
  iEvent.getByLabel(PFTauTagLeg1_,TauLeg1);

  edm::Handle<std::vector<reco::PFTau> > TauLeg2;
  iEvent.getByLabel(PFTauTagLeg2_,TauLeg2);

  edm::Handle<reco::VertexCollection > PV;
  iEvent.getByLabel(PVTag_,PV);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(beamSpotTag_,beamSpot);

  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent.getByLabel(TrackCollectionTag_,trackCollection);

  std::auto_ptr<VertexCollection>  VertexCollection_out= std::auto_ptr<VertexCollection>(new VertexCollection);

  //Get signal tracks from Taus
  std::vector<reco::TrackBaseRef> SignalTracks;
  if(TauLeg1.isValid()){
    for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < TauLeg1->size(); iPFTau++) {
      reco::PFTauRef RefPFTau(TauLeg1, iPFTau);
      ///////////////////////////////////////////////////////////////////////////////////////////////
      // Get tracks form PFTau daugthers
      const reco::PFCandidateRefVector & cands = RefPFTau->signalPFChargedHadrCands(); 
      for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
	if(iter->get()->trackRef().isNonnull()) SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
	else if(iter->get()->gsfTrackRef().isNonnull()){SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
      }
    }
  }

  if(TauLeg2.isValid()){
    for(reco::PFTauCollection::size_type iPFTau = 0; iPFTau < TauLeg2->size(); iPFTau++) {
      reco::PFTauRef RefPFTau(TauLeg2, iPFTau);
      ///////////////////////////////////////////////////////////////////////////////////////////////
      // Get tracks form PFTau daugthers
      const reco::PFCandidateRefVector & cands = RefPFTau->signalPFChargedHadrCands();
      for (reco::PFCandidateRefVector::const_iterator iter = cands.begin(); iter!=cands.end(); iter++){
        if(iter->get()->trackRef().isNonnull()) SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
        else if(iter->get()->gsfTrackRef().isNonnull()){SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
      }
    }
  }

  // Get Primary vertex
  reco::Vertex thePV;
  if(PV->size() > 0)
    thePV=PV->front();

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Get Non-Tau tracks 
  reco::TrackCollection nonTauTracks;
  if (trackCollection.isValid()) {
    // remove tau tracks and only tracks associated with the vertex
    unsigned int idx = 0;
    for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
      reco::TrackRef tmpRef(trackCollection, idx);
      reco::TrackRef tmpRefForBase=tmpRef;
      bool isSigTrk = false;
      bool fromVertex=false;
      for (unsigned int sigTrk = 0; sigTrk < SignalTracks.size(); sigTrk++) {
	if (reco::TrackBaseRef(tmpRefForBase)==SignalTracks.at(sigTrk)){isSigTrk = true; break;}
      }
      for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=thePV.tracks_begin();vtxTrkRef<thePV.tracks_end();vtxTrkRef++){
	if(thePV.trackWeight(*vtxTrkRef)>0 ){
	  if((*vtxTrkRef)==reco::TrackBaseRef(tmpRefForBase)){fromVertex=true; break;}
	}
      }
      if (!isSigTrk && fromVertex) nonTauTracks.push_back(*iTrk);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Refit the vertex
  TransientVertex transVtx;
  std::vector<reco::TransientTrack> transTracks;
  for (reco::TrackCollection::iterator iter=nonTauTracks.begin(); iter!=nonTauTracks.end(); ++iter){
    transTracks.push_back(transTrackBuilder->build(*iter));
  }
  bool FitOk(true);
  AdaptiveVertexFitter avf;
  avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception
  try{
    if(!useBeamSpot_){transVtx = avf.vertex(transTracks);}
    else{transVtx = avf.vertex(transTracks,*beamSpot);}
  }catch(...){
    FitOk=false;
  }
  reco::Vertex primaryVertexReFit;
  if(FitOk)primaryVertexReFit=transVtx;

  VertexCollection_out->push_back(primaryVertexReFit);
  iEvent.put(VertexCollection_out);
}

DEFINE_FWK_MODULE(DiTauPrimaryVertexProducer);
