#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include <memory>
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <utility>
#include <map>
#include "TMath.h"
#include <iostream>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <Math/SMatrixDfwd.h>

using namespace reco;
using namespace edm;
using namespace std;

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point3D; 
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Vector3D;
//typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > SymMatrix33;
//ROOT::Math::SMatrixSym3D
typedef ROOT::Math::SMatrixSym3D SymMatrix33;

struct DiTauInfo 
{ 
  DiTauInfo(){}; 
  int diTauCharge_; 
  double sumPt_; 
  int index1_;
  int index2_;
}; 

struct SortDiTauPairs 
{ 
  bool operator() (const DiTauInfo t1, const DiTauInfo t2) 
  { 
    // 1st criterion: OS 
    if ( t1.diTauCharge_ < t2.diTauCharge_ ) return true; 
    if ( t1.diTauCharge_ > t2.diTauCharge_ ) return false; 
    // 2nd criterion: sumPt of diTau pair 
    return (t1.sumPt_ > t2.sumPt_);  
  } 
}; 

class VertexAnalyzer : public edm::EDAnalyzer {
public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit VertexAnalyzer(const edm::ParameterSet& iConfig);
  ~VertexAnalyzer();
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void beginJob() ;
  void endJob() ;
  
private:
  edm::InputTag TauTag_;
  edm::InputTag PVTag_;
  edm::InputTag beamSpotTag_;
  edm::InputTag tauPVTag_;
  edm::InputTag TrackCollectionTag_;
  edm::InputTag jetTag_;
  edm::InputTag genParticlesTag_;
  bool useBeamSpot_;

  PFJetIDSelectionFunctor pfjetIDFunctor_;

  TTree* tree_;

  unsigned long run_,event_,lumi_;
  int index_;
  
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauLegsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genDiTauLegsP4_; 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTausP4_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* diTauLegsLchP4_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genVP4_;
  int genVPid_;
  std::vector< int > *genTausPid_;
  std::vector< int > *genTausCharge_;
  std::vector< int > *genTausStatus_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTauPSonsP4_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTauNSonsP4_;
  std::vector< int > *genTauPSonsPid_;
  std::vector< int > *genTauPSonsCharge_;
  std::vector< int > *genTauPSonsStatus_;
  std::vector< int > *genTauNSonsPid_;
  std::vector< int > *genTauNSonsCharge_;
  std::vector< int > *genTauNSonsStatus_;
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* JetsP4_;

  int tightestHPSDB3HWPLeg1_;
  int tightestHPSDB3HWPLeg2_;
  float hpsDB3HLeg1_;
  float hpsDB3HLeg2_;
  float diTauCharge_;
  float chargeLeg1_;
  int decayModeLeg1_;
  int decayModeLeg2_;
  int genDecayModeLeg1_;
  int genDecayModeLeg2_;

  float VtxX_;
  float VtxY_;
  float VtxZ_;
  float VtxXErr_;
  float VtxYErr_;
  float VtxZErr_;
  float VtxNChi2_;

  float IanVtxX_;
  float IanVtxY_;
  float IanVtxZ_;
  float IanVtxXErr_;
  float IanVtxYErr_;
  float IanVtxZErr_;
  float IanVtxNChi2_;

  float ReFitVtxX_;
  float ReFitVtxY_;
  float ReFitVtxZ_;
  float ReFitVtxXErr_;
  float ReFitVtxYErr_;
  float ReFitVtxZErr_;
  int   ReFitVtxRho_;
  float ReFitVtxNdof_;
  float ReFitVtxNChi2_;

  float ReFitVtxGenMatchX_;
  float ReFitVtxGenMatchY_;
  float ReFitVtxGenMatchZ_;
  float ReFitVtxGenMatchXErr_;
  float ReFitVtxGenMatchYErr_;
  float ReFitVtxGenMatchZErr_;
  int   ReFitVtxGenMatchRho_;
  float ReFitVtxGenMatchNdof_;
  float ReFitVtxGenMatchNChi2_;

  std::vector< Point3D >* diTauLegsPCA_;
  std::vector< Point3D >* diTauLegsPCAM2_;
  std::vector< Point3D >* diTauLegsPCAM2_GenMatchTks_;
  std::vector< Point3D >* diTauLegsPCAOPV_;
  std::vector< Point3D >* diTauLegsPCAGen_;
  std::vector< Point3D >* diTauLegsPCAGenM2_;
  std::vector< Point3D >* diTauLegsPCABS_;
  std::vector< Point3D >* diTauLegsPCAIan_;
  
  std::vector< SymMatrix33 >* diTauLegsPCAM2Cov_;
  std::vector< SymMatrix33 >* diTauLegsPCAGenM2Cov_;

  std::vector< Point3D >* VtxPos_;
  std::vector< Point3D >* BSPos_;
  std::vector< Point3D >* ReFitVtxPos_;
  std::vector< Point3D >* ReFitVtxGenMatchPos_;
  std::vector< Point3D >* IanVtxPos_;
  std::vector< Point3D >* HiggsGenVtx_;
  std::vector< Point3D >* TausGenVtx_;
  std::vector< Point3D >* TauPSonsGenVtx_;
  std::vector< Point3D >* TauNSonsGenVtx_;

  std::vector< Vector3D >* diTauLegsLchP3AtPCA_;
  std::vector< Vector3D >* diTauLegsLchP3AtPCAOPV_;
  std::vector< Vector3D >* diTauLegsLchP3AtPCABS_;
  std::vector< Vector3D >* diTauLegsLchP3AtPCAGen_;

  std::vector< Vector3D >* diTauLegsIPAtPCA_;
  std::vector< Vector3D >* diTauLegsIPAtPCA_GenMatchTks_;
  std::vector< Vector3D >* diTauLegsIPAtPCAV2_;
  std::vector< Vector3D >* diTauLegsIPAtPCAOPV_;
  std::vector< Vector3D >* diTauLegsIPAtPCAOPVV2_;
  std::vector< Vector3D >* diTauLegsIPAtPCABS_;
  std::vector< Vector3D >* diTauLegsIPAtPCABSV2_;
  std::vector< Vector3D >* diTauLegsIPAtPCAGen_;
  std::vector< Vector3D >* diTauLegsIPAtPCAGenV2_;
  std::vector< Vector3D >* diTauLegsIPAtPCATauVtx_;

  float trackPtLeg1_;
  float trackPtErrLeg1_;
  int nMisingHitsLeg1_;
  int nHitsLeg1_;
  int nTkHitsLeg1_;
  int nPxlHitsLeg1_;
  int hasFirstPxlHitLeg1_;
  int TkAlgoLeg1_;
  float trackPtLeg2_;
  float trackPtErrLeg2_;
  int nMisingHitsLeg2_;
  int nHitsLeg2_;
  int nTkHitsLeg2_;
  int nPxlHitsLeg2_;
  int hasFirstPxlHitLeg2_;
  int TkAlgoLeg2_;

  float dxyLeg1_;
  float dxyLeg2_;
  float dxyLeg1Err_;
  float dxyLeg2Err_;
  float dzLeg1_;
  float dzLeg2_;
  float dzLeg1Err_;
  float dzLeg2Err_;

  float dxyOPVLeg1_;
  float dxyOPVLeg2_;
  float dxyOPVLeg1Err_;
  float dxyOPVLeg2Err_;
  float dzOPVLeg1_;
  float dzOPVLeg2_;
  float dzOPVLeg1Err_;
  float dzOPVLeg2Err_;

  float dxyBSLeg1_;
  float dxyBSLeg2_;
  float dxyBSLeg1Err_;
  float dxyBSLeg2Err_;
  float dzBSLeg1_;
  float dzBSLeg2_;
  float dzBSLeg1Err_;
  float dzBSLeg2Err_;

  float dxyGenLeg1_;
  float dxyGenLeg2_;
  float dxyGenLeg1Err_;
  float dxyGenLeg2Err_;
  float dzGenLeg1_;
  float dzGenLeg2_;
  float dzGenLeg1Err_;
  float dzGenLeg2Err_;

  float genVtxX_;
  float genVtxY_;
  float genVtxZ_;
};

VertexAnalyzer::VertexAnalyzer(const edm::ParameterSet& iConfig):
  TauTag_(iConfig.getParameter<edm::InputTag>("tauTag")),
  PVTag_(iConfig.getParameter<edm::InputTag>("PVTag")),
  beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpot")),
  tauPVTag_(iConfig.getParameter<edm::InputTag>("tauPVTag")),
  TrackCollectionTag_(iConfig.getParameter<edm::InputTag>("TrackCollectionTag")),
  jetTag_(iConfig.getParameter<edm::InputTag>("JetTag")),
  useBeamSpot_(iConfig.getParameter<bool>("useBeamSpot")),
  genParticlesTag_(iConfig.getParameter<edm::InputTag>("genParticles"))
{

  pfjetIDFunctor_ = PFJetIDSelectionFunctor(iConfig.getParameter<edm::ParameterSet>("PFJetId") );
}

VertexAnalyzer::~VertexAnalyzer(){
  delete diTauLegsP4_; delete genDiTauLegsP4_; delete genTausP4_; delete JetsP4_;
  delete diTauLegsPCA_; delete diTauLegsLchP4_; delete diTauLegsPCAGen_;
  delete diTauLegsPCAM2_; delete diTauLegsPCAGenM2_; delete diTauLegsPCAIan_;
  delete diTauLegsPCAM2_GenMatchTks_; delete diTauLegsIPAtPCA_GenMatchTks_;
  delete diTauLegsPCAM2Cov_; delete diTauLegsPCAGenM2Cov_;
  delete diTauLegsLchP3AtPCA_; delete diTauLegsLchP3AtPCAOPV_; delete diTauLegsLchP3AtPCABS_; 
  delete diTauLegsLchP3AtPCAGen_; delete diTauLegsPCAOPV_; delete diTauLegsPCABS_;
  delete diTauLegsIPAtPCAGen_; delete diTauLegsIPAtPCAGenV2_;
  delete diTauLegsIPAtPCA_; delete diTauLegsIPAtPCAV2_;
  delete diTauLegsIPAtPCAOPV_; delete diTauLegsIPAtPCAOPVV2_;
  delete diTauLegsIPAtPCABS_; delete diTauLegsIPAtPCABSV2_; delete diTauLegsIPAtPCATauVtx_;
  delete genVP4_; delete BSPos_;
  delete VtxPos_; delete ReFitVtxPos_; delete IanVtxPos_; delete ReFitVtxGenMatchPos_;
  delete HiggsGenVtx_; delete TausGenVtx_; delete TauPSonsGenVtx_;
  delete TauNSonsGenVtx_;

  delete genTausPid_; delete genTausCharge_; delete genTausStatus_;
  delete genTauPSonsP4_; delete genTauNSonsP4_; delete genTauPSonsPid_;
  delete genTauPSonsCharge_; delete genTauPSonsStatus_;
  delete genTauNSonsPid_; delete genTauNSonsCharge_;
  delete genTauNSonsStatus_;
}

void VertexAnalyzer::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup){

  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.luminosityBlock();
  
  // Obtain Collections
  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
  
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByLabel(TauTag_,taus);
  //const pat::TauCollection* taus = taus.product();

  edm::Handle<reco::VertexCollection > PV;
  iEvent.getByLabel(PVTag_,PV);
  const reco::VertexCollection* vertexes = PV.product();

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(beamSpotTag_,beamSpot);

  edm::Handle<reco::VertexCollection > tauPV;
  iEvent.getByLabel(tauPVTag_,tauPV);
  const reco::VertexCollection* tauvertexes = tauPV.product();

  edm::Handle<reco::TrackCollection> trackCollection;
  iEvent.getByLabel(TrackCollectionTag_,trackCollection);

  edm::Handle<reco::GenParticleCollection> genHandle;
  const reco::GenParticleCollection* genParticles = 0;
  iEvent.getByLabel(genParticlesTag_,genHandle); 

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel(jetTag_,jets);

  //Store Jets
  JetsP4_->clear();
  pat::strbitset pfjetid = pfjetIDFunctor_.getBitTemplate();
  for(pat::JetCollection::size_type ijet = 0; ijet < jets->size(); ijet++) {
    pat::Jet jet( (*jets)[ijet] );
    
    bool pfId_ = pfjetIDFunctor_(jet, pfjetid);
    if(pfId_){
      JetsP4_->push_back(jet.p4());
    }
  }

  genTausCharge_->clear();
  genTausPid_->clear();
  genTausP4_->clear();
  genTausStatus_->clear();
  genVP4_->clear();
  genTauPSonsP4_->clear();
  genTauPSonsPid_->clear();
  genTauPSonsCharge_->clear();
  genTauPSonsStatus_->clear();
  genTauNSonsP4_->clear();
  genTauNSonsPid_->clear();
  genTauNSonsCharge_->clear();
  genTauNSonsStatus_->clear();
  HiggsGenVtx_->clear();
  TausGenVtx_->clear();
  TauPSonsGenVtx_->clear();
  TauNSonsGenVtx_->clear();
  
  genVtxX_ = -99.;
  genVtxY_ = -99.;
  genVtxZ_ = -99.;
  genVPid_ = -9999;

  if(genHandle.isValid()){
    genParticles = genHandle.product();
    /*
    int index1 = -99;
    int index2 = -99;
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( fabs((*genParticles)[k].pdgId()) != 15 || (*genParticles)[k].status()!=2) continue;
      if(index1<0) 
	index1 = k;
      else 
	index2 = k;
    }
    if(index1>=0 && index2>=0){
      if((*genParticles)[index1].pt()<(*genParticles)[index2].pt()){
	int indexBkp = index1;
	index1 = index2;
	index2 = indexBkp;
      }
      genTausP4_->push_back((*genParticles)[index1].p4());
      genTausP4_->push_back((*genParticles)[index2].p4());
    }
    */
    //Find Higgs production vertex
    for(unsigned int k = 0; k < genParticles->size(); k ++){
      if( ((*genParticles)[k].pdgId() == 23 || (*genParticles)[k].pdgId() == 25 || 
	   (*genParticles)[k].pdgId() == 35  || (*genParticles)[k].pdgId() == 36) && 
	  (*genParticles)[k].status()==3)
	{
	  genVP4_->push_back( (*genParticles)[k].p4() );
	  HiggsGenVtx_->push_back( Point3D((*genParticles)[k].vx(), (*genParticles)[k].vy(), (*genParticles)[k].vz()));
	  genVtxX_ = (*genParticles)[k].vx();
	  genVtxY_ = (*genParticles)[k].vy();
	  genVtxZ_ = (*genParticles)[k].vz();
	  genVPid_ = (*genParticles)[k].pdgId();
	  
	  for(unsigned int j = 0; j< (*genParticles)[k].numberOfDaughters(); j++){
	    const reco::Candidate * Hson = (*genParticles)[k].daughter(j);
	    if( abs(Hson->pdgId()) == 15 && (Hson->status() == 3 || Hson->status() == 2) ){
	      genTausP4_->push_back(Hson->p4());
	      genTausPid_->push_back(Hson->pdgId());
	      genTausCharge_->push_back(Hson->charge());
	      genTausStatus_->push_back(Hson->status());
	      TausGenVtx_->push_back(Point3D(Hson->vx(), Hson->vy(), Hson->vz()));
	      
	      //get the tau decay products
	      if(Hson->charge() > 0){ //Tau+
		if(Hson->numberOfDaughters() > 0){
		  for(unsigned int d = 0; d<Hson->numberOfDaughters(); d++) {
		    const reco::Candidate * TauSon = Hson->daughter(d);
		    if(abs(TauSon->pdgId()) != 15  && TauSon->status() == 1){
		      genTauPSonsP4_->push_back(TauSon->p4());
		      genTauPSonsPid_->push_back(TauSon->pdgId());
		      genTauPSonsCharge_->push_back(TauSon->charge());
		      genTauPSonsStatus_->push_back(TauSon->status());
		      TauPSonsGenVtx_->push_back(Point3D(TauSon->vx(), TauSon->vy(), TauSon->vz()));
		    }
		    else{
		      for(unsigned int e = 0; e<TauSon->numberOfDaughters(); e++) {
			const reco::Candidate * TauSon2 = TauSon->daughter(e);
			if(abs(TauSon2->pdgId()) != 15 && TauSon2->status() == 1){
			  genTauPSonsP4_->push_back(TauSon2->p4());
			  genTauPSonsPid_->push_back(TauSon2->pdgId());
			  genTauPSonsCharge_->push_back(TauSon2->charge());
			  genTauPSonsStatus_->push_back(TauSon2->status());
			  TauPSonsGenVtx_->push_back(Point3D(TauSon2->vx(), TauSon2->vy(), TauSon2->vz()));
			}
			else{
			  for(unsigned int f = 0; f<TauSon2->numberOfDaughters(); f++) {
			    const reco::Candidate * TauSon3 = TauSon2->daughter(f);
			    if(abs(TauSon3->pdgId()) != 15 && TauSon3->status() == 1){
			      genTauPSonsP4_->push_back(TauSon3->p4());
			      genTauPSonsPid_->push_back(TauSon3->pdgId());
			      genTauPSonsCharge_->push_back(TauSon3->charge());
			      genTauPSonsStatus_->push_back(TauSon3->status());
			      TauPSonsGenVtx_->push_back(Point3D(TauSon3->vx(), TauSon3->vy(), TauSon3->vz()));
			    }
			    else{
			      for(unsigned int g = 0; g<TauSon3->numberOfDaughters(); g++) {
				const reco::Candidate * TauSon4 = TauSon3->daughter(g);
				if(abs(TauSon4->pdgId()) != 15 && TauSon4->status() == 1){
				  genTauPSonsP4_->push_back(TauSon4->p4());
				  genTauPSonsPid_->push_back(TauSon4->pdgId());
				  genTauPSonsCharge_->push_back(TauSon4->charge());
				  genTauPSonsStatus_->push_back(TauSon4->status());
				  TauPSonsGenVtx_->push_back(Point3D(TauSon4->vx(), TauSon4->vy(), TauSon4->vz()));
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }//end of tau+
	      if(Hson->charge() < 0){ //Tau-
                if(Hson->numberOfDaughters() > 0){
                  for(unsigned int d = 0; d<Hson->numberOfDaughters(); d++) {
                    const reco::Candidate * TauSon = Hson->daughter(d);
                    if(abs(TauSon->pdgId()) != 15 && TauSon->status() == 1){
                      genTauNSonsP4_->push_back(TauSon->p4());
                      genTauNSonsPid_->push_back(TauSon->pdgId());
                      genTauNSonsCharge_->push_back(TauSon->charge());
		      genTauNSonsStatus_->push_back(TauSon->status());
		      TauNSonsGenVtx_->push_back(Point3D(TauSon->vx(), TauSon->vy(), TauSon->vz()));
                    }
                    else{
                      for(unsigned int e = 0; e<TauSon->numberOfDaughters(); e++) {
                        const reco::Candidate * TauSon2 = TauSon->daughter(e);
                        if(abs(TauSon2->pdgId()) != 15 && TauSon2->status() == 1){
                          genTauNSonsP4_->push_back(TauSon2->p4());
                          genTauNSonsPid_->push_back(TauSon2->pdgId());
                          genTauNSonsCharge_->push_back(TauSon2->charge());
			  genTauNSonsStatus_->push_back(TauSon2->status());
			  TauNSonsGenVtx_->push_back(Point3D(TauSon2->vx(), TauSon2->vy(), TauSon2->vz()));
                        }
			else{
			  for(unsigned int f = 0; f<TauSon2->numberOfDaughters(); f++) {
			    const reco::Candidate * TauSon3 = TauSon2->daughter(f);
			    if(abs(TauSon3->pdgId()) != 15 && TauSon3->status() == 1){
			      genTauNSonsP4_->push_back(TauSon3->p4());
			      genTauNSonsPid_->push_back(TauSon3->pdgId());
			      genTauNSonsCharge_->push_back(TauSon3->charge());
			      genTauNSonsStatus_->push_back(TauSon3->status());
			      TauNSonsGenVtx_->push_back(Point3D(TauSon3->vx(), TauSon3->vy(), TauSon3->vz()));
			    }
			    else{
			      for(unsigned int g = 0; g<TauSon3->numberOfDaughters(); g++) {
				const reco::Candidate * TauSon4 = TauSon3->daughter(g);
				if(abs(TauSon4->pdgId()) != 15 && TauSon4->status() == 1){
				  genTauNSonsP4_->push_back(TauSon4->p4());
				  genTauNSonsPid_->push_back(TauSon4->pdgId());
				  genTauNSonsCharge_->push_back(TauSon4->charge());
				  genTauNSonsStatus_->push_back(TauSon4->status());
				  TauNSonsGenVtx_->push_back(Point3D(TauSon4->vx(), TauSon4->vy(), TauSon4->vz()));
				}
			      }
			    }
			  }
			}
                      }
                    }
                  }
                }
	      }//end of tau-

	    }
	  } //end of taus
	} //end of Higgs

    } 
  } //end of genParticles

  
  std::vector<DiTauInfo>sortDiTauInfos; sortDiTauInfos.clear();

  //Apply selection to select tau pair
  if(taus->size()<2) return;
  for(pat::TauCollection::size_type itau = 0; itau < taus->size(); itau++) {
    pat::Tau aTauLeg1( (*taus)[itau] );
    
    if(aTauLeg1.pt() < 10.)continue;
    if(aTauLeg1.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue;
    if(aTauLeg1.tauID("againstMuonTight")<0.5) continue;
    if(aTauLeg1.tauID("againstElectronTight")<0.5) continue;
    
    for(pat::TauCollection::size_type jtau = itau+1; jtau < taus->size(); jtau++) {
      if(jtau<=itau) continue;
      pat::Tau aTauLeg2( (*taus)[jtau] );
      
      if(aTauLeg2.pt() < 10.)continue;
      if(aTauLeg2.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")<0.5) continue;
      if(aTauLeg2.tauID("againstMuonTight")<0.5) continue;
      if(aTauLeg2.tauID("againstElectronTight")<0.5) continue;
      
      float sumpt = aTauLeg1.pt() + aTauLeg2.pt();
      float pairCharge = aTauLeg1.charge()*aTauLeg2.charge();

      DiTauInfo sortDiTauInfo; 
      sortDiTauInfo.index1_ = itau; 
      sortDiTauInfo.index2_ = jtau;
      sortDiTauInfo.sumPt_ = sumpt;
      sortDiTauInfo.diTauCharge_ = pairCharge; 
      sortDiTauInfos.push_back(sortDiTauInfo); 
    }
  }

  //sort diTaus, first OS and then according to sumPt  
  std::sort(sortDiTauInfos.begin(), sortDiTauInfos.end(), SortDiTauPairs()); 

  //Intialize variables
  VtxX_ = -99.;
  VtxY_ = -99.;
  VtxZ_ = -99.;
  VtxXErr_ = -99.;
  VtxYErr_ = -99.;
  VtxZErr_ = -99.;

  // Get Primary vertex
  reco::Vertex thePV;
  bool vtxFound = false;
  if(vertexes->size() > 0){
    for(size_t ivtx = 0; ivtx < vertexes->size(); ivtx++){
      reco::Vertex theVtx = vertexes->at(ivtx);
      if(!theVtx.isValid()) continue;
      if(theVtx.ndof() < 4)continue;
      if(theVtx.z() < -24 || theVtx.z() > 24)continue;
      if(theVtx.position().Rho() > 2.)continue;
      thePV=vertexes->at(ivtx);
      vtxFound=true;
      break;
      }
    }
  
  if(!vtxFound) return;
  
  VtxX_ = thePV.x();
  VtxY_ = thePV.y();
  VtxZ_ = thePV.z();
  VtxXErr_ = thePV.xError();
  VtxYErr_ = thePV.yError();
  VtxZErr_ = thePV.zError();
  VtxNChi2_ = thePV.normalizedChi2();

  VtxPos_->clear();
  VtxPos_->push_back(Point3D(thePV.x(), thePV.y(), thePV.z()));
  BSPos_->clear();
  BSPos_->push_back((*beamSpot).position());

  // Get Refitted Primary vertex from Ian
  IanVtxX_ = -99.;
  IanVtxY_ = -99.;
  IanVtxZ_ = -99.;
  IanVtxXErr_ = -99.;
  IanVtxYErr_ = -99.;
  IanVtxZErr_ = -99.;
  IanVtxNChi2_ = -99.;
  IanVtxPos_->clear();

  if(tauvertexes->size() > 0){
    IanVtxX_ = (*tauvertexes)[0].x();
    IanVtxY_ = (*tauvertexes)[0].y();
    IanVtxZ_ = (*tauvertexes)[0].z();
    IanVtxXErr_ = (*tauvertexes)[0].xError();
    IanVtxYErr_ = (*tauvertexes)[0].yError();
    IanVtxZErr_ = (*tauvertexes)[0].zError();
    IanVtxNChi2_ = (*tauvertexes)[0].normalizedChi2();
    IanVtxPos_->push_back(Point3D((*tauvertexes)[0].x(), (*tauvertexes)[0].y(), (*tauvertexes)[0].z()));
  }
    
  ///////////////////////////////////////////////////////////////////////////////////////////////
  //Loop over tau-pairs and re-fit the vertex
  int diTauCounter = -1;
  for(std::vector<DiTauInfo>::iterator iter = sortDiTauInfos.begin();  iter != sortDiTauInfos.end() ; iter++){ 
    diTauCounter++;

    index_   = diTauCounter;

    ReFitVtxX_ = -99.;
    ReFitVtxY_ = -99.;
    ReFitVtxZ_ = -99.;
    ReFitVtxXErr_ = -99.;
    ReFitVtxYErr_ = -99.;
    ReFitVtxZErr_ = -99.;

    ReFitVtxGenMatchX_ = -99.;
    ReFitVtxGenMatchY_ = -99.;
    ReFitVtxGenMatchZ_ = -99.;
    ReFitVtxGenMatchXErr_ = -99.;
    ReFitVtxGenMatchYErr_ = -99.;
    ReFitVtxGenMatchZErr_ = -99.;

    dxyLeg1_ = -99.;
    dxyLeg2_ = -99.;
    dxyLeg1Err_ = -99.;
    dxyLeg2Err_ = -99.;
    dzLeg1_ = -99.;
    dzLeg2_ = -99.;
    dzLeg1Err_ = -99.;
    dzLeg2Err_ = -99.;

    dxyOPVLeg1_ = -99.;
    dxyOPVLeg2_ = -99.;
    dxyOPVLeg1Err_ = -99.;
    dxyOPVLeg2Err_ = -99.;
    dzOPVLeg1_ = -99.;
    dzOPVLeg2_ = -99.;
    dzOPVLeg1Err_ = -99.;
    dzOPVLeg2Err_ = -99.;

    diTauLegsP4_->clear();
    genDiTauLegsP4_->clear();
    diTauLegsLchP4_->clear();
    diTauLegsPCA_->clear();
    diTauLegsPCAM2_->clear();
    diTauLegsPCAM2_GenMatchTks_->clear();
    diTauLegsPCAOPV_->clear();
    diTauLegsPCABS_->clear();
    diTauLegsPCAGen_->clear();
    diTauLegsPCAGenM2_->clear();
    diTauLegsPCAM2Cov_->clear();
    diTauLegsPCAGenM2Cov_->clear();
    ReFitVtxPos_->clear();
    ReFitVtxGenMatchPos_->clear();

    diTauLegsLchP3AtPCA_->clear();
    diTauLegsLchP3AtPCAOPV_->clear();
    diTauLegsLchP3AtPCABS_->clear();
    diTauLegsLchP3AtPCAGen_->clear();

    diTauLegsIPAtPCA_->clear();
    diTauLegsIPAtPCA_GenMatchTks_->clear();
    diTauLegsIPAtPCAV2_->clear();
    diTauLegsIPAtPCAOPV_->clear();
    diTauLegsIPAtPCAOPVV2_->clear();
    diTauLegsIPAtPCABS_->clear();
    diTauLegsIPAtPCABSV2_->clear();
    diTauLegsIPAtPCAGen_->clear();
    diTauLegsIPAtPCAGenV2_->clear();

    trackPtLeg1_ = -99.;
    trackPtErrLeg1_ = 999.;
    nMisingHitsLeg1_ = 99;
    nHitsLeg1_ = -99;
    nTkHitsLeg1_ = -9;
    nPxlHitsLeg1_ = -99;
    hasFirstPxlHitLeg1_ = -99;
    TkAlgoLeg1_ = -99;
    trackPtLeg2_ = -99.;
    trackPtErrLeg2_ = 999.;
    nMisingHitsLeg2_ = 99;
    nHitsLeg2_ = -99;
    nTkHitsLeg2_ = -9;
    nPxlHitsLeg2_ = -99;
    hasFirstPxlHitLeg2_ = -99;
    TkAlgoLeg2_= -99;

    pat::Tau aTauLeg1( (*taus)[iter->index1_] );
    pat::Tau aTauLeg2( (*taus)[iter->index2_] );

    tightestHPSDB3HWPLeg1_ = -1;
    if(aTauLeg1.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWPLeg1_=0; 
    if(aTauLeg1.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")>0.5) tightestHPSDB3HWPLeg1_=1; 
    if(aTauLeg1.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWPLeg1_=2;
    hpsDB3HLeg1_  = aTauLeg1.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    tightestHPSDB3HWPLeg2_ = -1;
    if(aTauLeg2.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWPLeg2_=0;
    if(aTauLeg2.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")>0.5) tightestHPSDB3HWPLeg2_=1;
    if(aTauLeg2.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWPLeg2_=2;
    hpsDB3HLeg2_  = aTauLeg2.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    diTauCharge_ = aTauLeg1.charge()+aTauLeg2.charge();
    chargeLeg1_ = aTauLeg1.charge();
    
    if((aTauLeg1.signalPFChargedHadrCands()).size()==1 && (aTauLeg1.signalPFGammaCands()).size()==0) decayModeLeg1_ = 0; 
    else if((aTauLeg1.signalPFChargedHadrCands()).size()==1 && (aTauLeg1.signalPFGammaCands()).size()>0)  decayModeLeg1_ = 1; 
    else if((aTauLeg1.signalPFChargedHadrCands()).size()==2 && (aTauLeg1.signalPFGammaCands()).size()==0)  decayModeLeg1_ = 2; 
    else if((aTauLeg1.signalPFChargedHadrCands()).size()==2 && (aTauLeg1.signalPFGammaCands()).size()>0)  decayModeLeg1_ = 3; 
    else if((aTauLeg1.signalPFChargedHadrCands()).size()==3) decayModeLeg1_ = 4; 
    else  decayModeLeg1_ = -99;

    if((aTauLeg2.signalPFChargedHadrCands()).size()==1 && (aTauLeg2.signalPFGammaCands()).size()==0) decayModeLeg2_ = 0;
    else if((aTauLeg2.signalPFChargedHadrCands()).size()==1 && (aTauLeg2.signalPFGammaCands()).size()>0)  decayModeLeg2_ = 1;
    else if((aTauLeg2.signalPFChargedHadrCands()).size()==2 && (aTauLeg2.signalPFGammaCands()).size()==0)  decayModeLeg2_ = 2;
    else if((aTauLeg2.signalPFChargedHadrCands()).size()==2 && (aTauLeg2.signalPFGammaCands()).size()>0)  decayModeLeg2_ = 3;
    else if((aTauLeg2.signalPFChargedHadrCands()).size()==3) decayModeLeg2_ = 4;
    else  decayModeLeg2_ = -99;
    
    //gendecayMode
    genDecayModeLeg1_ = -99;
    if( (aTauLeg1.genParticleById( 15,0,true)).isNonnull() || 
	(aTauLeg1.genParticleById(-15,0,true)).isNonnull() ) {
      if(aTauLeg1.genJet()){
	genDiTauLegsP4_->push_back(aTauLeg1.genJet()->p4());
	string genTauDecay = JetMCTagUtils::genTauDecayMode( *(aTauLeg1.genJet()) );
	if( genTauDecay.find("oneProng0Pi0")!=string::npos ) 
	  genDecayModeLeg1_ = 0;
	else if( genTauDecay.find("oneProng1Pi0")!=string::npos )
	  genDecayModeLeg1_ = 1;
	else if( genTauDecay.find("oneProng2Pi0")!=string::npos )
	  genDecayModeLeg1_ = 2;
	else if( genTauDecay.find("oneProngOther")!=string::npos )
	  genDecayModeLeg1_ = 3;
	else if( genTauDecay.find("threeProng0Pi0")!=string::npos )
	  genDecayModeLeg1_ = 4;
	else if( genTauDecay.find("threeProng1Pi0")!=string::npos )
	  genDecayModeLeg1_ = 5;
	else if( genTauDecay.find("threeProngOther")!=string::npos )
	  genDecayModeLeg1_ = 6;
	else if( genTauDecay.find("rare")!=string::npos )
	  genDecayModeLeg1_ = 7;
	else
	  genDecayModeLeg1_ = 99;
      }
      else genDiTauLegsP4_->push_back(math::XYZTLorentzVectorD(0,0,0,0));
    }
    else genDiTauLegsP4_->push_back(math::XYZTLorentzVectorD(0,0,0,0));

    genDecayModeLeg2_ = -99;
    if( (aTauLeg2.genParticleById( 15,0,true)).isNonnull() ||
        (aTauLeg2.genParticleById(-15,0,true)).isNonnull() ) {
      if(aTauLeg2.genJet()){
	genDiTauLegsP4_->push_back(aTauLeg2.genJet()->p4());
        string genTauDecay = JetMCTagUtils::genTauDecayMode( *(aTauLeg2.genJet()) );
        if( genTauDecay.find("oneProng0Pi0")!=string::npos )
          genDecayModeLeg2_ = 0;
        else if( genTauDecay.find("oneProng1Pi0")!=string::npos )
          genDecayModeLeg2_ = 1;
        else if( genTauDecay.find("oneProng2Pi0")!=string::npos )
          genDecayModeLeg2_ = 2;
        else if( genTauDecay.find("oneProngOther")!=string::npos )
          genDecayModeLeg2_ = 3;
        else if( genTauDecay.find("threeProng0Pi0")!=string::npos )
          genDecayModeLeg2_ = 4;
        else if( genTauDecay.find("threeProng1Pi0")!=string::npos )
          genDecayModeLeg2_ = 5;
        else if( genTauDecay.find("threeProngOther")!=string::npos )
          genDecayModeLeg2_ = 6;
        else if( genTauDecay.find("rare")!=string::npos )
          genDecayModeLeg2_ = 7;
        else
          genDecayModeLeg2_ = 99;
      }
      else genDiTauLegsP4_->push_back(math::XYZTLorentzVectorD(0,0,0,0));
    }
    else genDiTauLegsP4_->push_back(math::XYZTLorentzVectorD(0,0,0,0));
    
    diTauLegsP4_->push_back(aTauLeg1.p4());
    diTauLegsP4_->push_back(aTauLeg2.p4());

    if(aTauLeg1.leadPFChargedHadrCand().isNonnull()){
      diTauLegsLchP4_->push_back(aTauLeg1.leadPFChargedHadrCand()->p4());
      if(aTauLeg1.leadPFChargedHadrCand()->trackRef().isNonnull()){
	reco::TrackRef track1 = aTauLeg1.leadPFChargedHadrCand()->trackRef();
	trackPtLeg1_ = track1->pt();
	trackPtErrLeg1_ = track1->ptError();
	nMisingHitsLeg1_ = track1->trackerExpectedHitsInner().numberOfHits();
	nHitsLeg1_ = track1->hitPattern().numberOfValidHits();
	nTkHitsLeg1_ = track1->hitPattern().numberOfValidTrackerHits();
	nPxlHitsLeg1_ = track1->hitPattern().numberOfValidPixelHits();
	hasFirstPxlHitLeg1_ = int(track1->hitPattern().hasValidHitInFirstPixelBarrel() ||  track1->hitPattern().hasValidHitInFirstPixelEndcap());
	TkAlgoLeg1_ = track1->algo();
      }
    }
    if(aTauLeg2.leadPFChargedHadrCand().isNonnull()){
      diTauLegsLchP4_->push_back(aTauLeg2.leadPFChargedHadrCand()->p4());
      if(aTauLeg2.leadPFChargedHadrCand()->trackRef().isNonnull()){
	reco::TrackRef track2 = aTauLeg2.leadPFChargedHadrCand()->trackRef();
        trackPtLeg2_ = track2->pt();
        trackPtErrLeg2_ = track2->ptError();
        nMisingHitsLeg2_ = track2->trackerExpectedHitsInner().numberOfHits();
        nHitsLeg2_ = track2->hitPattern().numberOfValidHits();
        nTkHitsLeg2_ = track2->hitPattern().numberOfValidTrackerHits();
        nPxlHitsLeg2_ = track2->hitPattern().numberOfValidPixelHits();
	hasFirstPxlHitLeg2_ = int(track2->hitPattern().hasValidHitInFirstPixelBarrel() ||  track2->hitPattern().hasValidHitInFirstPixelEndcap());
	TkAlgoLeg2_ = track2->algo();
      } 
    }
    if(aTauLeg2.leadPFChargedHadrCand().isNonnull())diTauLegsLchP4_->push_back(aTauLeg2.leadPFChargedHadrCand()->p4());

    //Get signal tracks from Taus
    std::vector<reco::TrackBaseRef> SignalTracks;

    // Get tracks form PFTau daugthers of tauLeg1
    const std::vector<reco::PFCandidatePtr>& candsLeg1 = aTauLeg1.signalPFChargedHadrCands(); 
    for (std::vector<reco::PFCandidatePtr>::const_iterator iter = candsLeg1.begin(); iter!=candsLeg1.end(); iter++){
      if(iter->get()->trackRef().isNonnull()) SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
      else if(iter->get()->gsfTrackRef().isNonnull()){SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
    }

    // Get tracks form PFTau daugthers of tauLeg2
    const std::vector<reco::PFCandidatePtr>& candsLeg2 = aTauLeg2.signalPFChargedHadrCands();
    for (std::vector<reco::PFCandidatePtr>::const_iterator iter = candsLeg2.begin(); iter!=candsLeg2.end(); iter++){
      if(iter->get()->trackRef().isNonnull()) SignalTracks.push_back(reco::TrackBaseRef((*iter)->trackRef()));
      else if(iter->get()->gsfTrackRef().isNonnull()){SignalTracks.push_back(reco::TrackBaseRef(((*iter)->gsfTrackRef())));}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Get Non-Tau tracks 
    reco::TrackCollection nonTauTracks;
    reco::TrackCollection nonTauTracks_genMatched;
    if (trackCollection.isValid()) {
      // remove tau tracks and only tracks associated with the vertex
      unsigned int idx = 0;
      for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) {
	reco::TrackRef tmpRef(trackCollection, idx);
	reco::TrackRef tmpRefForBase=tmpRef;
	bool isSigTrk = false;
	bool fromVertex=false;
	for (unsigned int sigTrk = 0; sigTrk < SignalTracks.size(); sigTrk++) {
	  if (reco::TrackBaseRef(tmpRefForBase)==SignalTracks.at(sigTrk)){isSigTrk = true; /*std::cout<<"matches signal track "<<std::endl;*/ break;}
	}
	for(std::vector<reco::TrackBaseRef>::const_iterator vtxTrkRef=thePV.tracks_begin();vtxTrkRef<thePV.tracks_end();vtxTrkRef++){
	  if(thePV.trackWeight(*vtxTrkRef)>0 ){
	    if((*vtxTrkRef)==reco::TrackBaseRef(tmpRefForBase)){fromVertex=true; /*std::cout<<"part of vertex "<<std::endl;*/ break;}
	  }
	}
	if (!isSigTrk && fromVertex){
	  nonTauTracks.push_back(*iTrk);

	  //Check if matches to genParticle
	  for(unsigned int k = 0; k < genParticles->size(); k ++){
	    if( (*genParticles)[k].status() == 1 && (*genParticles)[k].charge() != 0 ){
	      if(deltaR((*genParticles)[k].p4(), iTrk->momentum()) < 0.3){
		nonTauTracks_genMatched.push_back(*iTrk);
		break;
	      }
	    }
	  }
	}
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////
    reco::Vertex primaryVertexReFit;
    bool RefitOK_ = false;
    if(nonTauTracks.size() >= 2){
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

      //reco::Vertex primaryVertexReFit;
      if(FitOk){
	primaryVertexReFit=transVtx;
	RefitOK_ = FitOk;
	
	ReFitVtxX_ = primaryVertexReFit.x();
	ReFitVtxY_ = primaryVertexReFit.y();
	ReFitVtxZ_ = primaryVertexReFit.z();
	ReFitVtxXErr_ = primaryVertexReFit.xError();
	ReFitVtxYErr_ = primaryVertexReFit.yError();
	ReFitVtxZErr_ = primaryVertexReFit.zError();
	ReFitVtxRho_ = primaryVertexReFit.position().Rho();
	ReFitVtxNdof_ = primaryVertexReFit.ndof();
	ReFitVtxNChi2_ = primaryVertexReFit.normalizedChi2(); //normalisedChiSquared();
	ReFitVtxPos_->push_back(Point3D(primaryVertexReFit.x(), primaryVertexReFit.y(), primaryVertexReFit.z()));
      }
    }
  
    reco::Vertex pvReFit_genMatched;
    bool Refit_genMatched_OK_ = false;
    if(nonTauTracks_genMatched.size() >= 2){
      // Refit the vertex
      TransientVertex transVtx;
      std::vector<reco::TransientTrack> transTracks;
      for (reco::TrackCollection::iterator iter=nonTauTracks_genMatched.begin(); iter!=nonTauTracks_genMatched.end(); ++iter){
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

      if(FitOk){
	pvReFit_genMatched=transVtx;
        Refit_genMatched_OK_ = FitOk;

        ReFitVtxGenMatchX_ = pvReFit_genMatched.x();
        ReFitVtxGenMatchY_ = pvReFit_genMatched.y();
        ReFitVtxGenMatchZ_ = pvReFit_genMatched.z();
        ReFitVtxGenMatchXErr_ = pvReFit_genMatched.xError();
        ReFitVtxGenMatchYErr_ = pvReFit_genMatched.yError();
        ReFitVtxGenMatchZErr_ = pvReFit_genMatched.zError();
        ReFitVtxGenMatchRho_ = pvReFit_genMatched.position().Rho();
        ReFitVtxGenMatchNdof_ = pvReFit_genMatched.ndof();
        ReFitVtxGenMatchNChi2_ = pvReFit_genMatched.normalizedChi2(); //normalisedChiSquared();
        ReFitVtxGenMatchPos_->push_back(Point3D(pvReFit_genMatched.x(), pvReFit_genMatched.y(), pvReFit_genMatched.z()));
      }
    }

    //Compute IP vectors for the signal tracks
    if(aTauLeg1.leadPFChargedHadrCand().isNonnull()){
      if(aTauLeg1.leadPFChargedHadrCand()->trackRef().isNonnull()){
	reco::TransientTrack transTrk=transTrackBuilder->build(aTauLeg1.leadPFChargedHadrCand()->trackRef());
	AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
	if(RefitOK_){
	  GlobalPoint pv(primaryVertexReFit.position().x(),primaryVertexReFit.position().y(),primaryVertexReFit.position().z());
	  TrajectoryStateClosestToPoint pca_ = transTrk.trajectoryStateClosestToPoint(pv);
	  float dxy=pca_.perigeeParameters().transverseImpactParameter();
	  float dz=pca_.perigeeParameters().longitudinalImpactParameter();
	  float dxy_err=pca_.perigeeError().transverseImpactParameterError();
	  float dz_err=pca_.perigeeError().longitudinalImpactParameterError();
	  GlobalPoint pos=pca_.position();
	  Point3D pcaLeg1=reco::Vertex::Point(pos.x(),pos.y(),pos.z());
	  dxyLeg1_ = dxy;
	  dxyLeg1Err_ = dxy_err;
	  dzLeg1_ = dz;
	  dzLeg1Err_ = dz_err;
	  diTauLegsPCA_->push_back(pcaLeg1);
	  diTauLegsLchP3AtPCA_->push_back(Vector3D(pca_.momentum().x(), pca_.momentum().y(), pca_.momentum().z()));
	  //Use a separate method to do linear extrapolation in Z.
	  //Extrapolate to closest point on transverse plane
	  //AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
	  TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(primaryVertexReFit.position()));
	  diTauLegsPCAM2_->push_back(Point3D(closestIn3DSpaceState.globalPosition().x(), closestIn3DSpaceState.globalPosition().y(), closestIn3DSpaceState.globalPosition().z()));
	  GlobalVector ip = IPTools::linearImpactParameter(closestIn3DSpaceState, pv);
	  diTauLegsIPAtPCA_->push_back(Vector3D(ip.x(), ip.y(), ip.z()));
	  GlobalVector ipV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pv);
	  diTauLegsIPAtPCAV2_->push_back(Vector3D(ipV2.x(), ipV2.y(), ipV2.z()));
	  //Get Covariance matrix
	  const SymMatrix33 PCACov_ = closestIn3DSpaceState.cartesianError().position().matrix_new();
	  diTauLegsPCAM2Cov_->push_back(PCACov_);
	  //std::cout<<" cov. 11 "<<PCACov_(0,0)<<" 22 "<<PCACov_(1,1)<<" 33 "<<PCACov_(2,2)<<std::endl;
	}
	if(Refit_genMatched_OK_){
	  GlobalPoint pv(pvReFit_genMatched.position().x(),pvReFit_genMatched.position().y(),pvReFit_genMatched.position().z());
	  TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(pvReFit_genMatched.position()));
          diTauLegsPCAM2_GenMatchTks_->push_back(Point3D(closestIn3DSpaceState.globalPosition().x(), closestIn3DSpaceState.globalPosition().y(), closestIn3DSpaceState.globalPosition().z()));
          GlobalVector ip = IPTools::linearImpactParameter(closestIn3DSpaceState, pv);
          diTauLegsIPAtPCA_GenMatchTks_->push_back(Vector3D(ip.x(), ip.y(), ip.z()));
	}

	//wrt original vertex
	GlobalPoint pvOrig(thePV.position().x(),thePV.position().y(),thePV.position().z());
	TrajectoryStateClosestToPoint pcaopv_ = transTrk.trajectoryStateClosestToPoint(pvOrig);
	float dxyOPV=pcaopv_.perigeeParameters().transverseImpactParameter();
	float dzOPV=pcaopv_.perigeeParameters().longitudinalImpactParameter();
	float dxyOPV_err=pcaopv_.perigeeError().transverseImpactParameterError();
	float dzOPV_err=pcaopv_.perigeeError().longitudinalImpactParameterError();
	GlobalPoint posOPV=pcaopv_.position();
	Point3D pcaLeg1_opv=reco::Vertex::Point(posOPV.x(),posOPV.y(),posOPV.z());
	dxyOPVLeg1_ = dxyOPV;
	dxyOPVLeg1Err_ = dxyOPV_err;
	dzOPVLeg1_ = dzOPV;
	dzOPVLeg1Err_ = dzOPV_err;
	diTauLegsPCAOPV_->push_back(pcaLeg1_opv);
	diTauLegsLchP3AtPCAOPV_->push_back(Vector3D(pcaopv_.momentum().x(), pcaopv_.momentum().y(), pcaopv_.momentum().z()));
	//Use a separate method to do linear extrapolation in Z.
	//Extrapolate to closest point on transverse plane
	TrajectoryStateOnSurface closestIn3DSpaceStateOPV = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(thePV.position()));
	GlobalVector ipOPV = IPTools::linearImpactParameter(closestIn3DSpaceStateOPV, pvOrig);
	diTauLegsIPAtPCAOPV_->push_back(Vector3D(ipOPV.x(), ipOPV.y(), ipOPV.z()));
	GlobalVector ipOPVV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pvOrig);
	diTauLegsIPAtPCAOPVV2_->push_back(Vector3D(ipOPVV2.x(), ipOPVV2.y(), ipOPVV2.z()));
	
	//wrt BeamSpot
	const GlobalPoint pvBs((*beamSpot).x0(),(*beamSpot).y0(),(*beamSpot).z0());
	TrajectoryStateClosestToPoint pcabs_ = transTrk.trajectoryStateClosestToPoint(pvBs);
	float dxyBS=pcabs_.perigeeParameters().transverseImpactParameter();
	float dzBS=pcabs_.perigeeParameters().longitudinalImpactParameter();
	float dxyBS_err=pcabs_.perigeeError().transverseImpactParameterError();
	float dzBS_err=pcabs_.perigeeError().longitudinalImpactParameterError();
	GlobalPoint posBS=pcabs_.position();
	Point3D pcaLeg1_bs=reco::Vertex::Point(posBS.x(),posBS.y(),posBS.z());
	dxyBSLeg1_ = dxyBS;
	dxyBSLeg1Err_ = dxyBS_err;
	dzBSLeg1_ = dzBS;
	dzBSLeg1Err_ = dzBS_err;
	diTauLegsPCABS_->push_back(pcaLeg1_bs);
	diTauLegsLchP3AtPCABS_->push_back(Vector3D(pcabs_.momentum().x(), pcabs_.momentum().y(), pcabs_.momentum().z()));
	//Use a separate method to do linear extrapolation in Z.
	//Extrapolate to closest point on transverse plane
	TrajectoryStateOnSurface closestIn3DSpaceStateBS = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos((*beamSpot).position()));
	GlobalVector ipBS = IPTools::linearImpactParameter(closestIn3DSpaceStateBS, pvBs);
	diTauLegsIPAtPCABS_->push_back(Vector3D(ipBS.x(), ipBS.y(), ipBS.z()));
	GlobalVector ipBSV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pvBs);
	diTauLegsIPAtPCABSV2_->push_back(Vector3D(ipBSV2.x(), ipBSV2.y(), ipBSV2.z()));

	//wrt refitted vertex from Ian
	if(tauvertexes->size() > 0){
	  //GlobalPoint pvIan((*tauvertexes)[0].position().x(),(*tauvertexes)[0].position().y(),(*tauvertexes)[0].position().z());
	  //TrajectoryStateClosestToPoint pcaIan_ = transTrk.trajectoryStateClosestToPoint(pvIan);
	  //GlobalPoint posIan=pcaIan_.position();
	  //Point3D pcaLeg1_Ian=reco::Vertex::Point(posIan.x(),posIan.y(),posIan.z());
	  TrajectoryStateOnSurface closestIn3DSpaceStateIan = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos((*tauvertexes)[0].position()));
	  diTauLegsPCAIan_->push_back(Point3D(closestIn3DSpaceStateIan.globalPosition().x(), closestIn3DSpaceStateIan.globalPosition().y(), closestIn3DSpaceStateIan.globalPosition().z()));
	}
	
	//wrt generated Higgs vertex
	if(HiggsGenVtx_->size() > 0){
	  GlobalPoint pvGen(((*HiggsGenVtx_)[0]).x(), ((*HiggsGenVtx_)[0]).y(), ((*HiggsGenVtx_)[0]).z());
	  TrajectoryStateClosestToPoint pcagen_ = transTrk.trajectoryStateClosestToPoint(pvGen);
	  GlobalPoint posGen=pcagen_.position();
	  Point3D pcaLeg1_gen=reco::Vertex::Point(posGen.x(),posGen.y(),posGen.z());
	  diTauLegsPCAGen_->push_back(pcaLeg1_gen);
	  diTauLegsLchP3AtPCAGen_->push_back(Vector3D(pcagen_.momentum().x(), pcagen_.momentum().y(), pcagen_.momentum().z()));
	  float dxyGen=pcagen_.perigeeParameters().transverseImpactParameter();
	  float dzGen=pcagen_.perigeeParameters().longitudinalImpactParameter();
	  float dxyGen_err=pcagen_.perigeeError().transverseImpactParameterError();
	  float dzGen_err=pcagen_.perigeeError().longitudinalImpactParameterError();
	  dxyGenLeg1_ = dxyGen;
	  dxyGenLeg1Err_ = dxyGen_err;
	  dzGenLeg1_ = dzGen;
	  dzGenLeg1Err_ = dzGen_err;
	  
	  //Use a separate method to do linear extrapolation in Z.
	  //Extrapolate to closest point on transverse plane
	  //AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
	  TrajectoryStateOnSurface closestIn3DSpaceStateGen = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos((*HiggsGenVtx_)[0]));
	  diTauLegsPCAGenM2_->push_back(Point3D(closestIn3DSpaceStateGen.globalPosition().x(), closestIn3DSpaceStateGen.globalPosition().y(), closestIn3DSpaceStateGen.globalPosition().z()));
	  GlobalVector ipGen = IPTools::linearImpactParameter(closestIn3DSpaceStateGen, pvGen);
	  diTauLegsIPAtPCAGen_->push_back(Vector3D(ipGen.x(), ipGen.y(), ipGen.z()));
	  
	  GlobalVector ipGenV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pvGen);
	  diTauLegsIPAtPCAGenV2_->push_back(Vector3D(ipGenV2.x(), ipGenV2.y(), ipGenV2.z()));
	  //Get Covariance matrix
	  const SymMatrix33 PCACovGen_ = closestIn3DSpaceStateGen.cartesianError().position().matrix_new();
	  diTauLegsPCAGenM2Cov_->push_back(PCACovGen_);
	  //std::cout<<" gen cov. 11 "<<PCACovGen_(0,0)<<" 22 "<<PCACovGen_(1,1)<<" 33 "<<PCACovGen_(2,2)<<std::endl;
	}

	//wrt associated tau vertex
	GlobalPoint pvTauVtx(aTauLeg1.primaryVertexPos().x(), aTauLeg1.primaryVertexPos().y(), aTauLeg1.primaryVertexPos().z());
	TrajectoryStateOnSurface closestIn3DSpaceStateTauVtx = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(aTauLeg1.primaryVertexPos()));
	GlobalVector ipTauVtx = IPTools::linearImpactParameter(closestIn3DSpaceStateTauVtx, pvTauVtx);
	diTauLegsIPAtPCATauVtx_->push_back(Vector3D(ipTauVtx.x(), ipTauVtx.y(), ipTauVtx.z()));
	
      }
    }
    
    if(aTauLeg2.leadPFChargedHadrCand().isNonnull()){
      if(aTauLeg2.leadPFChargedHadrCand()->trackRef().isNonnull()){
	reco::TransientTrack transTrk=transTrackBuilder->build(aTauLeg2.leadPFChargedHadrCand()->trackRef());
	AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
	if(RefitOK_){
	  GlobalPoint pv(primaryVertexReFit.position().x(),primaryVertexReFit.position().y(),primaryVertexReFit.position().z());
	  TrajectoryStateClosestToPoint pca_ = transTrk.trajectoryStateClosestToPoint(pv);
	  float dxy=-pca_.perigeeParameters().transverseImpactParameter();
	  float dz=pca_.perigeeParameters().longitudinalImpactParameter();
	  float dxy_err=pca_.perigeeError().transverseImpactParameterError();
	  float dz_err=pca_.perigeeError().longitudinalImpactParameterError();
	  GlobalPoint pos=pca_.position();
	  Point3D pcaLeg2=reco::Vertex::Point(pos.x(),pos.y(),pos.z());
	  dxyLeg2_ = dxy;
	  dxyLeg2Err_ = dxy_err;
	  dzLeg2_ = dz;
	  dzLeg2Err_ = dz_err;
	  diTauLegsPCA_->push_back(pcaLeg2);
	  diTauLegsLchP3AtPCA_->push_back(Vector3D(pca_.momentum().x(), pca_.momentum().y(), pca_.momentum().z()));
	  //Use a separate method to do linear extrapolation in Z.
	  //Extrapolate to closest point on transverse plane
	  //AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
	  TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(primaryVertexReFit.position()));
	  diTauLegsPCAM2_->push_back(Point3D(closestIn3DSpaceState.globalPosition().x(), closestIn3DSpaceState.globalPosition().y(), closestIn3DSpaceState.globalPosition().z()));
	  GlobalVector ip = IPTools::linearImpactParameter(closestIn3DSpaceState, pv);
	  diTauLegsIPAtPCA_->push_back(Vector3D(ip.x(), ip.y(), ip.z()));
	  GlobalVector ipV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pv);
	  diTauLegsIPAtPCAV2_->push_back(Vector3D(ipV2.x(), ipV2.y(), ipV2.z()));
	  //Get Covariance matrix
	  const SymMatrix33 PCACov_ = closestIn3DSpaceState.cartesianError().position().matrix_new();
	  diTauLegsPCAM2Cov_->push_back(PCACov_);
	  //std::cout<<" cov. 11 "<<PCACov_(0,0)<<" 22 "<<PCACov_(1,1)<<" 33 "<<PCACov_(2,2)<<std::endl;
	}
	if(Refit_genMatched_OK_){
	  GlobalPoint pv(pvReFit_genMatched.position().x(),pvReFit_genMatched.position().y(),pvReFit_genMatched.position().z());
          TrajectoryStateOnSurface closestIn3DSpaceState = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(pvReFit_genMatched.position()));
          diTauLegsPCAM2_GenMatchTks_->push_back(Point3D(closestIn3DSpaceState.globalPosition().x(), closestIn3DSpaceState.globalPosition().y(), closestIn3DSpaceState.globalPosition().z()));
          GlobalVector ip = IPTools::linearImpactParameter(closestIn3DSpaceState, pv);
          diTauLegsIPAtPCA_GenMatchTks_->push_back(Vector3D(ip.x(), ip.y(), ip.z()));
        }

	//wrt original vertex
	GlobalPoint pvOrig(thePV.position().x(),thePV.position().y(),thePV.position().z());
	TrajectoryStateClosestToPoint pcaopv_ = transTrk.trajectoryStateClosestToPoint(pvOrig);
	float dxyOPV=pcaopv_.perigeeParameters().transverseImpactParameter();
	float dzOPV=pcaopv_.perigeeParameters().longitudinalImpactParameter();
	float dxyOPV_err=pcaopv_.perigeeError().transverseImpactParameterError();
	float dzOPV_err=pcaopv_.perigeeError().longitudinalImpactParameterError();
	GlobalPoint posOPV=pcaopv_.position();
	Point3D pcaLeg2_opv=reco::Vertex::Point(posOPV.x(),posOPV.y(),posOPV.z());
	dxyOPVLeg2_ = dxyOPV;
	dxyOPVLeg2Err_ = dxyOPV_err;
	dzOPVLeg2_ = dzOPV;
	dzOPVLeg2Err_ = dzOPV_err;
	diTauLegsPCAOPV_->push_back(pcaLeg2_opv);
	diTauLegsLchP3AtPCAOPV_->push_back(Vector3D(pcaopv_.momentum().x(), pcaopv_.momentum().y(), pcaopv_.momentum().z()));
	//Use a separate method to do linear extrapolation in Z.
	//Extrapolate to closest point on transverse plane
	TrajectoryStateOnSurface closestIn3DSpaceStateOPV = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(thePV.position()));
	GlobalVector ipOPV = IPTools::linearImpactParameter(closestIn3DSpaceStateOPV, pvOrig);
	diTauLegsIPAtPCAOPV_->push_back(Vector3D(ipOPV.x(), ipOPV.y(), ipOPV.z()));
	GlobalVector ipOPVV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pvOrig);
	diTauLegsIPAtPCAOPVV2_->push_back(Vector3D(ipOPVV2.x(), ipOPVV2.y(), ipOPVV2.z()));
	
	//wrt BeamSpot
	const GlobalPoint pvBs((*beamSpot).x0(),(*beamSpot).y0(),(*beamSpot).z0());
	TrajectoryStateClosestToPoint pcabs_ = transTrk.trajectoryStateClosestToPoint(pvBs);
	float dxyBS=pcabs_.perigeeParameters().transverseImpactParameter();
	float dzBS=pcabs_.perigeeParameters().longitudinalImpactParameter();
	float dxyBS_err=pcabs_.perigeeError().transverseImpactParameterError();
	float dzBS_err=pcabs_.perigeeError().longitudinalImpactParameterError();
	GlobalPoint posBS=pcabs_.position();
	Point3D pcaLeg2_bs=reco::Vertex::Point(posBS.x(),posBS.y(),posBS.z());
	dxyBSLeg2_ = dxyBS;
	dxyBSLeg2Err_ = dxyBS_err;
	dzBSLeg2_ = dzBS;
	dzBSLeg2Err_ = dzBS_err;
	diTauLegsPCABS_->push_back(pcaLeg2_bs);
	diTauLegsLchP3AtPCABS_->push_back(Vector3D(pcabs_.momentum().x(), pcabs_.momentum().y(), pcabs_.momentum().z()));
	//Use a separate method to do linear extrapolation in Z.
	//Extrapolate to closest point on transverse plane
	TrajectoryStateOnSurface closestIn3DSpaceStateBS = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos((*beamSpot).position()));
	GlobalVector ipBS = IPTools::linearImpactParameter(closestIn3DSpaceStateBS, pvBs);
	diTauLegsIPAtPCABS_->push_back(Vector3D(ipBS.x(), ipBS.y(), ipBS.z()));
	GlobalVector ipBSV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pvBs);
	diTauLegsIPAtPCABSV2_->push_back(Vector3D(ipBSV2.x(), ipBSV2.y(), ipBSV2.z()));

	//wrt refitted vertex from Ian
	if(tauvertexes->size() > 0){
	  TrajectoryStateOnSurface closestIn3DSpaceStateIan = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos((*tauvertexes)[0].position()));
          diTauLegsPCAIan_->push_back(Point3D(closestIn3DSpaceStateIan.globalPosition().x(), closestIn3DSpaceStateIan.globalPosition().y(), closestIn3DSpaceStateIan.globalPosition().z()));
        }
	
	//wrt generated Higgs vertex
	if(HiggsGenVtx_->size() > 0){
	  GlobalPoint pvGen(((*HiggsGenVtx_)[0]).x(), ((*HiggsGenVtx_)[0]).y(), ((*HiggsGenVtx_)[0]).z());
	  TrajectoryStateClosestToPoint pcagen_ = transTrk.trajectoryStateClosestToPoint(pvGen);
	  GlobalPoint posGen=pcagen_.position();
	  Point3D pcaLeg2_gen=reco::Vertex::Point(posGen.x(),posGen.y(),posGen.z());
	  diTauLegsPCAGen_->push_back(pcaLeg2_gen);
	  diTauLegsLchP3AtPCAGen_->push_back(Vector3D(pcagen_.momentum().x(), pcagen_.momentum().y(), pcagen_.momentum().z()));
	  float dxyGen=pcagen_.perigeeParameters().transverseImpactParameter();
	  float dzGen=pcagen_.perigeeParameters().longitudinalImpactParameter();
	  float dxyGen_err=pcagen_.perigeeError().transverseImpactParameterError();
	  float dzGen_err=pcagen_.perigeeError().longitudinalImpactParameterError();
	  dxyGenLeg2_ = dxyGen;
	  dxyGenLeg2Err_ = dxyGen_err;
	  dzGenLeg2_ = dzGen;
	  dzGenLeg2Err_ = dzGen_err;
	  
	  //Use a separate method to do linear extrapolation in Z.
	  //Extrapolate to closest point on transverse plane
	  //AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
	  TrajectoryStateOnSurface closestIn3DSpaceStateGen = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos((*HiggsGenVtx_)[0]));
	  diTauLegsPCAGenM2_->push_back(Point3D(closestIn3DSpaceStateGen.globalPosition().x(), closestIn3DSpaceStateGen.globalPosition().y(), closestIn3DSpaceStateGen.globalPosition().z()));
	  GlobalVector ipGen = IPTools::linearImpactParameter(closestIn3DSpaceStateGen, pvGen);
	  diTauLegsIPAtPCAGen_->push_back(Vector3D(ipGen.x(), ipGen.y(), ipGen.z()));
	  
	  GlobalVector ipGenV2 = IPTools::linearImpactParameter(transTrk.impactPointState(), pvGen);
	  diTauLegsIPAtPCAGenV2_->push_back(Vector3D(ipGenV2.x(), ipGenV2.y(), ipGenV2.z()));
	  //Get Covariance matrix
	  const SymMatrix33 PCACovGen_ = closestIn3DSpaceStateGen.cartesianError().position().matrix_new();
	  diTauLegsPCAGenM2Cov_->push_back(PCACovGen_);
	  //std::cout<<" gen cov. 11 "<<PCACovGen_(0,0)<<" 22 "<<PCACovGen_(1,1)<<" 33 "<<PCACovGen_(2,2)<<std::endl;
	}
	
	//wrt associated tau vertex
	GlobalPoint pvTauVtx(aTauLeg2.primaryVertexPos().x(), aTauLeg2.primaryVertexPos().y(), aTauLeg2.primaryVertexPos().z());
	TrajectoryStateOnSurface closestIn3DSpaceStateTauVtx = extrapolator.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(aTauLeg2.primaryVertexPos()));
	GlobalVector ipTauVtx = IPTools::linearImpactParameter(closestIn3DSpaceStateTauVtx, pvTauVtx);
	diTauLegsIPAtPCATauVtx_->push_back(Vector3D(ipTauVtx.x(), ipTauVtx.y(), ipTauVtx.z()));
	
      }
      
    }

    tree_->Fill();
  }
}

void VertexAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","ditau tree");

  diTauLegsP4_    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genDiTauLegsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genTausP4_      = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genVP4_      = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  diTauLegsLchP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genTauPSonsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  genTauNSonsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();
  JetsP4_ = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >();

  diTauLegsLchP3AtPCA_ = new std::vector< Vector3D >();
  diTauLegsLchP3AtPCAOPV_ = new std::vector< Vector3D >();
  diTauLegsLchP3AtPCABS_ = new std::vector< Vector3D >();
  diTauLegsLchP3AtPCAGen_ = new std::vector< Vector3D >();

  diTauLegsIPAtPCA_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCA_GenMatchTks_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCAV2_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCAOPV_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCAOPVV2_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCABS_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCABSV2_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCAGen_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCAGenV2_ = new std::vector< Vector3D >();
  diTauLegsIPAtPCATauVtx_ = new std::vector< Vector3D >();

  diTauLegsPCA_ = new std::vector< Point3D >();
  diTauLegsPCAM2_ = new std::vector< Point3D >();
  diTauLegsPCAM2_GenMatchTks_ = new std::vector< Point3D >();
  diTauLegsPCAOPV_ = new std::vector< Point3D >();
  diTauLegsPCAGen_ = new std::vector< Point3D >();
  diTauLegsPCAGenM2_ = new std::vector< Point3D >();
  diTauLegsPCABS_ = new std::vector< Point3D >();
  diTauLegsPCAIan_ = new std::vector< Point3D >();
  diTauLegsPCAM2Cov_ = new std::vector< SymMatrix33 >();
  diTauLegsPCAGenM2Cov_ = new std::vector< SymMatrix33 >();
  
  VtxPos_ = new std::vector< Point3D >();
  BSPos_ = new std::vector< Point3D >();
  ReFitVtxPos_ = new std::vector< Point3D >();
  ReFitVtxGenMatchPos_ = new std::vector< Point3D >();
  IanVtxPos_ = new std::vector< Point3D >();
  HiggsGenVtx_ = new std::vector< Point3D >();
  TausGenVtx_= new std::vector< Point3D >();
  TauPSonsGenVtx_ = new std::vector< Point3D >();
  TauNSonsGenVtx_ = new std::vector< Point3D >();

  genTausPid_ = new std::vector< int >();
  genTausCharge_ = new std::vector< int >();
  genTausStatus_ = new std::vector< int >();
  genTauPSonsPid_ = new std::vector< int >();
  genTauPSonsCharge_ = new std::vector< int >();
  genTauPSonsStatus_ = new std::vector< int >();
  genTauNSonsPid_ = new std::vector< int >();
  genTauNSonsCharge_ = new std::vector< int >();
  genTauNSonsStatus_ = new std::vector< int >();
  

  tree_->Branch("diTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauLegsP4_);
  tree_->Branch("genDiTauLegsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genDiTauLegsP4_);
  tree_->Branch("genTausP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genTausP4_);
  tree_->Branch("diTauLegsLchP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&diTauLegsLchP4_);
  tree_->Branch("genVP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genVP4_);
  tree_->Branch("genTauPSonsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genTauPSonsP4_);
  tree_->Branch("genTauNSonsP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&genTauNSonsP4_);

  tree_->Branch("JetsP4", "std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&JetsP4_);

  tree_->Branch("diTauLegsLchP3AtPCA","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsLchP3AtPCA_);
  tree_->Branch("diTauLegsLchP3AtPCAOPV","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsLchP3AtPCAOPV_);
  tree_->Branch("diTauLegsLchP3AtPCABS","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsLchP3AtPCABS_);
  tree_->Branch("diTauLegsLchP3AtPCAGen","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsLchP3AtPCAGen_);

  tree_->Branch("diTauLegsIPAtPCA","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCA_);
  tree_->Branch("diTauLegsIPAtPCA_GenMatchTks","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCA_GenMatchTks_);
  tree_->Branch("diTauLegsIPAtPCAV2","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCAV2_);
  tree_->Branch("diTauLegsIPAtPCAOPV","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCAOPV_);
  tree_->Branch("diTauLegsIPAtPCAOPVV2","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCAOPVV2_);
  tree_->Branch("diTauLegsIPAtPCABS","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCABS_);
  tree_->Branch("diTauLegsIPAtPCABSV2","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCABSV2_);
  tree_->Branch("diTauLegsIPAtPCAGen","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCAGen_);
  tree_->Branch("diTauLegsIPAtPCAGenV2","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCAGenV2_);
  tree_->Branch("diTauLegsIPAtPCATauVtx","std::vector< ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >",&diTauLegsIPAtPCATauVtx_);

  tree_->Branch("genVPid", &genVPid_, "genVPid/I");
  tree_->Branch("genTausPid", "std::vector< int >", &genTausPid_);
  tree_->Branch("genTausCharge", "std::vector< int >", &genTausCharge_);
  tree_->Branch("genTausStatus", "std::vector< int >", &genTausStatus_);
  tree_->Branch("genTauPSonsPid", "std::vector< int >", &genTauPSonsPid_);
  tree_->Branch("genTauPSonsCharge", "std::vector< int >", &genTauPSonsCharge_);
  tree_->Branch("genTauPSonsStatus", "std::vector< int >", &genTauPSonsStatus_);
  tree_->Branch("genTauNSonsPid", "std::vector< int >", &genTauNSonsPid_);
  tree_->Branch("genTauNSonsCharge", "std::vector< int >", &genTauNSonsCharge_);
  tree_->Branch("genTauNSonsStatus", "std::vector< int >", &genTauNSonsStatus_);

  tree_->Branch("run",&run_,"run/l");
  tree_->Branch("event",&event_,"event/l");
  tree_->Branch("lumi",&lumi_,"lumi/l");
  tree_->Branch("genVtxZ",&genVtxZ_,"genVtxZ/F");
  tree_->Branch("genVtxX",&genVtxX_,"genVtxX/F");
  tree_->Branch("genVtxY",&genVtxY_,"genVtxY/F");
  tree_->Branch("VtxZ",&VtxZ_,"VtxZ/F");
  tree_->Branch("VtxX",&VtxX_,"VtxX/F");
  tree_->Branch("VtxY",&VtxY_,"VtxY/F");
  tree_->Branch("VtxZErr",&VtxZErr_,"VtxZErr/F");
  tree_->Branch("VtxXErr",&VtxXErr_,"VtxXErr/F");
  tree_->Branch("VtxYErr",&VtxYErr_,"VtxYErr/F");
  tree_->Branch("VtxNChi2", &VtxNChi2_, "VtxNChi2/F");
  tree_->Branch("ReFitVtxZ",&ReFitVtxZ_,"ReFitVtxZ/F");
  tree_->Branch("ReFitVtxX",&ReFitVtxX_,"ReFitVtxX/F");
  tree_->Branch("ReFitVtxY",&ReFitVtxY_,"ReFitVtxY/F");
  tree_->Branch("ReFitVtxZErr",&ReFitVtxZErr_,"ReFitVtxZErr/F");
  tree_->Branch("ReFitVtxXErr",&ReFitVtxXErr_,"ReFitVtxXErr/F");
  tree_->Branch("ReFitVtxYErr",&ReFitVtxYErr_,"ReFitVtxYErr/F");
  tree_->Branch("ReFitVtxNdof",&ReFitVtxNdof_,"ReFitVtxNdof/I");
  tree_->Branch("ReFitVtxRho",&ReFitVtxRho_,"ReFitVtxRho/F");
  tree_->Branch("ReFitVtxNChi2", &ReFitVtxNChi2_,"ReFitVtxNChi2/F");
  tree_->Branch("IanVtxZ",&IanVtxZ_,"IanVtxZ/F");
  tree_->Branch("IanVtxX",&IanVtxX_,"IanVtxX/F");
  tree_->Branch("IanVtxY",&IanVtxY_,"IanVtxY/F");
  tree_->Branch("IanVtxZErr",&IanVtxZErr_,"IanVtxZErr/F");
  tree_->Branch("IanVtxXErr",&IanVtxXErr_,"IanVtxXErr/F");
  tree_->Branch("IanVtxYErr",&IanVtxYErr_,"IanVtxYErr/F");
  tree_->Branch("IanVtxNChi2", &IanVtxNChi2_,"IanVtxNChi2/F");
  tree_->Branch("ReFitVtxGenMatchZ",&ReFitVtxGenMatchZ_,"ReFitVtxGenMatchZ/F");
  tree_->Branch("ReFitVtxGenMatchX",&ReFitVtxGenMatchX_,"ReFitVtxGenMatchX/F");
  tree_->Branch("ReFitVtxGenMatchY",&ReFitVtxGenMatchY_,"ReFitVtxGenMatchY/F");
  tree_->Branch("ReFitVtxGenMatchZErr",&ReFitVtxGenMatchZErr_,"ReFitVtxGenMatchZErr/F");
  tree_->Branch("ReFitVtxGenMatchXErr",&ReFitVtxGenMatchXErr_,"ReFitVtxGenMatchXErr/F");
  tree_->Branch("ReFitVtxGenMatchYErr",&ReFitVtxGenMatchYErr_,"ReFitVtxGenMatchYErr/F");
  tree_->Branch("ReFitVtxGenMatchNdof",&ReFitVtxGenMatchNdof_,"ReFitVtxGenMatchNdof/I");
  tree_->Branch("ReFitVtxGenMatchRho",&ReFitVtxGenMatchRho_,"ReFitVtxGenMatchRho/F");
  tree_->Branch("ReFitVtxGenMatchNChi2", &ReFitVtxGenMatchNChi2_,"ReFitVtxGenMatchNChi2/F");

  tree_->Branch("diTauLegsPCA", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCA_);
  tree_->Branch("diTauLegsPCAM2", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCAM2_);
  tree_->Branch("diTauLegsPCAM2_GenMatchTks", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCAM2_GenMatchTks_);
  tree_->Branch("diTauLegsPCAOPV", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCAOPV_);
  tree_->Branch("diTauLegsPCAGen", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCAGen_);
  tree_->Branch("diTauLegsPCAGenM2", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCAGenM2_);
  tree_->Branch("diTauLegsPCABS", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCABS_);
  tree_->Branch("diTauLegsPCAIan", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &diTauLegsPCAIan_);

  tree_->Branch("diTauLegsPCAM2Cov", "std::vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > >", &diTauLegsPCAM2Cov_);
  tree_->Branch("diTauLegsPCAGenM2Cov", "std::vector<ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > >", &diTauLegsPCAGenM2Cov_);
  //tree_->Branch("diTauLegsPCA", "std::vector< Point3D >", &diTauLegsPCA_);
  //tree_->Branch("diTauLegsPCAOPV", "std::vector< Point3D >", &diTauLegsPCAOPV_);

  tree_->Branch("VtxPos", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &VtxPos_);
  tree_->Branch("BSPos", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &BSPos_);
  tree_->Branch("ReFitVtxPos", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &ReFitVtxPos_);
  tree_->Branch("ReFitVtxGenMatchPos", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &ReFitVtxGenMatchPos_);
  tree_->Branch("IanVtxPos", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &IanVtxPos_);
  tree_->Branch("HiggsGenVtx", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &HiggsGenVtx_);
  tree_->Branch("TausGenVtx", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &TausGenVtx_);
  tree_->Branch("TauPSonsGenVtx", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &TauPSonsGenVtx_);
  tree_->Branch("TauNSonsGenVtx", "std::vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >", &TauNSonsGenVtx_);

  tree_->Branch("trackPtLeg1", &trackPtLeg1_, "trackPtLeg1/F");
  tree_->Branch("trackPtErrLeg1", &trackPtErrLeg1_, "trackPtErrLeg1/F");
  tree_->Branch("nMisingHitsLeg1", &nMisingHitsLeg1_, "nMisingHitsLeg1/I");
  tree_->Branch("nHitsLeg1", &nHitsLeg1_, "nHitsLeg1/I");
  tree_->Branch("nTkHitsLeg1", &nTkHitsLeg1_, "nTkHitsLeg1/I");
  tree_->Branch("nPxlHitsLeg1", &nPxlHitsLeg1_, "nPxlHitsLeg1/I");
  tree_->Branch("hasFirstPxlHitLeg1", &hasFirstPxlHitLeg1_, "hasFirstPxlHitLeg1/I");
  tree_->Branch("TkAlgoLeg1", &TkAlgoLeg1_, "TkAlgoLeg1/I");
  tree_->Branch("trackPtLeg2", &trackPtLeg2_, "trackPtLeg2/F");
  tree_->Branch("trackPtErrLeg2", &trackPtErrLeg2_, "trackPtErrLeg2/F");
  tree_->Branch("nMisingHitsLeg2", &nMisingHitsLeg2_, "nMisingHitsLeg2/I");
  tree_->Branch("nHitsLeg2", &nHitsLeg2_, "nHitsLeg2/I");
  tree_->Branch("nTkHitsLeg2", &nTkHitsLeg2_, "nTkHitsLeg2/I");
  tree_->Branch("nPxlHitsLeg2", &nPxlHitsLeg2_, "nPxlHitsLeg2/I");
  tree_->Branch("hasFirstPxlHitLeg2", &hasFirstPxlHitLeg2_, "hasFirstPxlHitLeg2/I");
  tree_->Branch("TkAlgoLeg2", &TkAlgoLeg2_, "TkAlgoLeg2/I");

  tree_->Branch("dxyLeg1",&dxyLeg1_, "dxyLeg1/F");
  tree_->Branch("dxyLeg2",&dxyLeg2_, "dxyLeg2/F");
  tree_->Branch("dxyLeg1Err",&dxyLeg1Err_, "dxyLeg1Err/F");
  tree_->Branch("dxyLeg2Err",&dxyLeg2Err_, "dxyLeg2Err/F");
  tree_->Branch("dzLeg1",&dzLeg1_, "dzLeg1/F");
  tree_->Branch("dzLeg2",&dzLeg2_, "dzLeg2/F");
  tree_->Branch("dzLeg1Err",&dzLeg1Err_, "dzLeg1Err/F");
  tree_->Branch("dzLeg2Err",&dzLeg2Err_, "dzLeg2Err/F");

  tree_->Branch("dxyOPVLeg1",&dxyOPVLeg1_, "dxyOPVLeg1/F");
  tree_->Branch("dxyOPVLeg2",&dxyOPVLeg2_, "dxyOPVLeg2/F");
  tree_->Branch("dxyOPVLeg1Err",&dxyOPVLeg1Err_, "dxyOPVLeg1Err/F");
  tree_->Branch("dxyOPVLeg2Err",&dxyOPVLeg2Err_, "dxyOPVLeg2Err/F");
  tree_->Branch("dzOPVLeg1",&dzOPVLeg1_, "dzOPVLeg1/F");
  tree_->Branch("dzOPVLeg2",&dzOPVLeg2_, "dzOPVLeg2/F");
  tree_->Branch("dzOPVLeg1Err",&dzOPVLeg1Err_, "dzOPVLeg1Err/F");
  tree_->Branch("dzOPVLeg2Err",&dzOPVLeg2Err_, "dzOPVLeg2Err/F");

  tree_->Branch("dxyBSLeg1",&dxyBSLeg1_, "dxyBSLeg1/F");
  tree_->Branch("dxyBSLeg2",&dxyBSLeg2_, "dxyBSLeg2/F");
  tree_->Branch("dxyBSLeg1Err",&dxyBSLeg1Err_, "dxyBSLeg1Err/F");
  tree_->Branch("dxyBSLeg2Err",&dxyBSLeg2Err_, "dxyBSLeg2Err/F");
  tree_->Branch("dzBSLeg1",&dzBSLeg1_, "dzBSLeg1/F");
  tree_->Branch("dzBSLeg2",&dzBSLeg2_, "dzBSLeg2/F");
  tree_->Branch("dzBSLeg1Err",&dzBSLeg1Err_, "dzBSLeg1Err/F");
  tree_->Branch("dzBSLeg2Err",&dzBSLeg2Err_, "dzBSLeg2Err/F");

  tree_->Branch("dxyGenLeg1",&dxyGenLeg1_, "dxyGenLeg1/F");
  tree_->Branch("dxyGenLeg2",&dxyGenLeg2_, "dxyGenLeg2/F");
  tree_->Branch("dxyGenLeg1Err",&dxyGenLeg1Err_, "dxyGenLeg1Err/F");
  tree_->Branch("dxyGenLeg2Err",&dxyGenLeg2Err_, "dxyGenLeg2Err/F");
  tree_->Branch("dzGenLeg1",&dzGenLeg1_, "dzGenLeg1/F");
  tree_->Branch("dzGenLeg2",&dzGenLeg2_, "dzGenLeg2/F");
  tree_->Branch("dzGenLeg1Err",&dzGenLeg1Err_, "dzGenLeg1Err/F");
  tree_->Branch("dzGenLeg2Err",&dzGenLeg2Err_, "dzGenLeg2Err/F");

  tree_->Branch("decayModeLeg1",&decayModeLeg1_,"decayModeLeg1/I");
  tree_->Branch("decayModeLeg2",&decayModeLeg2_,"decayModeLeg2/I");
  tree_->Branch("diTauCharge",&diTauCharge_,"diTauCharge/F");
  tree_->Branch("chargeLeg1",&chargeLeg1_,"chargeLeg1/F");
  tree_->Branch("genDecayModeLeg1",&genDecayModeLeg1_,"genDecayModeLeg1/I");
  tree_->Branch("genDecayModeLeg2",&genDecayModeLeg2_,"genDecayModeLeg2/I");

  tree_->Branch("tightestHPSDB3HWPLeg1",&tightestHPSDB3HWPLeg1_,"tightestHPSDB3HWPLeg1/I");
  tree_->Branch("tightestHPSDB3HWPLeg2",&tightestHPSDB3HWPLeg2_,"tightestHPSDB3HWPLeg2/I");
  tree_->Branch("hpsDB3HLeg1",&hpsDB3HLeg1_,"hpsDB3HLeg1/F");
  tree_->Branch("hpsDB3HLeg2",&hpsDB3HLeg2_,"hpsDB3HLeg2/F");
  tree_->Branch("index",&index_,"index/I");
  
}

void VertexAnalyzer::endJob(){}

DEFINE_FWK_MODULE(VertexAnalyzer);
