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

//typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point3D; 
//typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Vector3D;
//typedef ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> > SymMatrix33;
//ROOT::Math::SMatrixSym3D
//typedef ROOT::Math::SMatrixSym3D SymMatrix33;

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

class  TauIDCrossCheck : public edm::EDAnalyzer {
public:

  struct more {
    bool operator() (const double& lhs, const double& rhs) const
    {return lhs>rhs;}
  };

  explicit  TauIDCrossCheck(const edm::ParameterSet& iConfig);
  ~ TauIDCrossCheck();
  void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void beginJob() ;
  void endJob() ;
  
private:
  edm::InputTag TauTag_;

  TTree* tree_;

  unsigned long run_,event_,lumi_;

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *TauP4_;
  std::vector< int > decayModeFinding_;
  std::vector< int > decayModeFindingNewDM_;
  std::vector< int > decayModeFindingOldDM_;
  std::vector< int > tightestAntiECutWP_;
  std::vector< int > tightestAntiEMVA5WP_;
  std::vector< int > tightestAntiMuWP_;
  std::vector< int > tightestAntiMu2WP_;
  std::vector< int > tightestAntiMu3WP_;
  std::vector< int > tightestAntiMuMVAWP_;
  std::vector< int > tightestHPSDB3HWP_;
  std::vector< int > tightestHPSMVA3newDMwLTWP_;
  std::vector< int > tightestHPSMVA3newDMwoLTWP_;
  std::vector< int > tightestHPSMVA3oldDMwLTWP_;
  std::vector< int > tightestHPSMVA3oldDMwoLTWP_;
  std::vector< float > AntiEMVA5raw_;
  std::vector< int > AntiEMVA5category_;
  std::vector< float > AntiMuMVAraw_;
  std::vector< float > hpsDB3H_;
  std::vector< float > hpsMVA3newDMwLT_;
  std::vector< float > hpsMVA3newDMwoLT_;
  std::vector< float > hpsMVA3oldDMwLT_;
  std::vector< float > hpsMVA3oldDMwoLT_;
  
};

 TauIDCrossCheck:: TauIDCrossCheck(const edm::ParameterSet& iConfig):
  TauTag_(iConfig.getParameter<edm::InputTag>("tauTag"))
{

}

TauIDCrossCheck::~ TauIDCrossCheck(){
  delete TauP4_;
}

void  TauIDCrossCheck::analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup){

  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.luminosityBlock();
  
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByLabel(TauTag_,taus);
  //const pat::TauCollection* taus = taus.product();

  TauP4_->clear();
  decayModeFinding_.clear();
  decayModeFindingNewDM_.clear();
  decayModeFindingOldDM_.clear();
  tightestAntiECutWP_.clear();
  tightestAntiEMVA5WP_.clear();
  tightestAntiMuWP_.clear();
  tightestAntiMu2WP_.clear();
  tightestAntiMu3WP_.clear();
  tightestAntiMuMVAWP_.clear();
  tightestHPSDB3HWP_.clear();
  tightestHPSMVA3newDMwLTWP_.clear();
  tightestHPSMVA3newDMwoLTWP_.clear();
  tightestHPSMVA3oldDMwLTWP_.clear();
  tightestHPSMVA3oldDMwoLTWP_.clear();
  AntiEMVA5raw_.clear();
  AntiEMVA5category_.clear();
  AntiMuMVAraw_.clear();
  hpsDB3H_.clear();
  hpsMVA3newDMwLT_.clear();
  hpsMVA3newDMwoLT_.clear();
  hpsMVA3oldDMwLT_.clear();
  hpsMVA3oldDMwoLT_.clear();

  //collect taus
  for(pat::TauCollection::size_type itau = 0; itau < taus->size(); itau++) {
    pat::Tau aTau( (*taus)[itau] );
    
    if(aTau.pt() < 20 || TMath::Abs(aTau.eta()) > 2.3) continue;
    TauP4_->push_back(aTau.p4());
    
    //TauID
    decayModeFinding_.push_back(aTau.tauID("decayModeFinding"));
    decayModeFindingNewDM_.push_back(aTau.tauID("decayModeFindingNewDMs"));
    decayModeFindingOldDM_.push_back(aTau.tauID("decayModeFindingOldDMs"));

    //int tightestHPSDBWP = -1;
    //if(aTau.tauID("byVLooseCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP=0;
    //if(aTau.tauID("byLooseCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP=1;
    //if(aTau.tauID("byMediumCombinedIsolationDeltaBetaCorr")>0.5) tightestHPSDBWP=2;
    //if(aTau.tauID("byTightCombinedIsolationDeltaBetaCorr")>0.5)  tightestHPSDBWP=3;
    //tightestHPSDBWP_.push_back(tightestHPSDBWP);

    int tightestHPSDB3HWP = -1; 
    if(aTau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWP=0; 
    if(aTau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")>0.5) tightestHPSDB3HWP=1; 
    if(aTau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")>0.5)  tightestHPSDB3HWP=2;
    tightestHPSDB3HWP_.push_back(tightestHPSDB3HWP);
    hpsDB3H_.push_back(aTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));

    int tightestHPSMVA3newDMwLTWP = -1;
    if(aTau.tauID("byVLooseIsolationMVA3newDMwLT") >0.5) tightestHPSMVA3newDMwLTWP=0;
    if(aTau.tauID("byLooseIsolationMVA3newDMwLT")>0.5)   tightestHPSMVA3newDMwLTWP=1;
    if(aTau.tauID("byMediumIsolationMVA3newDMwLT") >0.5) tightestHPSMVA3newDMwLTWP=2;
    if(aTau.tauID("byTightIsolationMVA3newDMwLT") >0.5)  tightestHPSMVA3newDMwLTWP=3;
    if(aTau.tauID("byVTightIsolationMVA3newDMwLT") >0.5) tightestHPSMVA3newDMwLTWP=4;
    if(aTau.tauID("byVVTightIsolationMVA3newDMwLT") >0.5)tightestHPSMVA3newDMwLTWP=5;
    tightestHPSMVA3newDMwLTWP_.push_back(tightestHPSMVA3newDMwLTWP);
    hpsMVA3newDMwLT_.push_back(aTau.tauID("byIsolationMVA3newDMwLTraw"));

    int tightestHPSMVA3newDMwoLTWP = -1;
    if(aTau.tauID("byVLooseIsolationMVA3newDMwoLT") >0.5) tightestHPSMVA3newDMwoLTWP=0;
    if(aTau.tauID("byLooseIsolationMVA3newDMwoLT")>0.5)   tightestHPSMVA3newDMwoLTWP=1;
    if(aTau.tauID("byMediumIsolationMVA3newDMwoLT") >0.5) tightestHPSMVA3newDMwoLTWP=2;
    if(aTau.tauID("byTightIsolationMVA3newDMwoLT") >0.5)  tightestHPSMVA3newDMwoLTWP=3;
    if(aTau.tauID("byVTightIsolationMVA3newDMwoLT") >0.5) tightestHPSMVA3newDMwoLTWP=4;
    if(aTau.tauID("byVVTightIsolationMVA3newDMwoLT") >0.5)tightestHPSMVA3newDMwoLTWP=5;
    tightestHPSMVA3newDMwoLTWP_.push_back(tightestHPSMVA3newDMwoLTWP);
    hpsMVA3newDMwoLT_.push_back(aTau.tauID("byIsolationMVA3newDMwoLTraw"));

    int tightestHPSMVA3oldDMwLTWP = -1;
    if(aTau.tauID("byVLooseIsolationMVA3oldDMwLT") >0.5) tightestHPSMVA3oldDMwLTWP=0;
    if(aTau.tauID("byLooseIsolationMVA3oldDMwLT")>0.5)   tightestHPSMVA3oldDMwLTWP=1;
    if(aTau.tauID("byMediumIsolationMVA3oldDMwLT") >0.5) tightestHPSMVA3oldDMwLTWP=2;
    if(aTau.tauID("byTightIsolationMVA3oldDMwLT") >0.5)  tightestHPSMVA3oldDMwLTWP=3;
    if(aTau.tauID("byVTightIsolationMVA3oldDMwLT") >0.5) tightestHPSMVA3oldDMwLTWP=4;
    if(aTau.tauID("byVVTightIsolationMVA3oldDMwLT") >0.5)tightestHPSMVA3oldDMwLTWP=5;
    tightestHPSMVA3oldDMwLTWP_.push_back(tightestHPSMVA3oldDMwLTWP);
    hpsMVA3oldDMwLT_.push_back(aTau.tauID("byIsolationMVA3oldDMwLTraw"));

    int tightestHPSMVA3oldDMwoLTWP = -1;
    if(aTau.tauID("byVLooseIsolationMVA3oldDMwoLT") >0.5) tightestHPSMVA3oldDMwoLTWP=0;
    if(aTau.tauID("byLooseIsolationMVA3oldDMwoLT")>0.5)   tightestHPSMVA3oldDMwoLTWP=1;
    if(aTau.tauID("byMediumIsolationMVA3oldDMwoLT") >0.5) tightestHPSMVA3oldDMwoLTWP=2;
    if(aTau.tauID("byTightIsolationMVA3oldDMwoLT") >0.5)  tightestHPSMVA3oldDMwoLTWP=3;
    if(aTau.tauID("byVTightIsolationMVA3oldDMwoLT") >0.5) tightestHPSMVA3oldDMwoLTWP=4;
    if(aTau.tauID("byVVTightIsolationMVA3oldDMwoLT") >0.5)tightestHPSMVA3oldDMwoLTWP=5;
    tightestHPSMVA3oldDMwoLTWP_.push_back(tightestHPSMVA3oldDMwoLTWP);
    hpsMVA3oldDMwoLT_.push_back(aTau.tauID("byIsolationMVA3oldDMwoLTraw"));

    int tightestAntiMuWP = 0; 
    if( aTau.tauID("againstMuonLoose")>0.5 )tightestAntiMuWP = 1; 
    if( aTau.tauID("againstMuonMedium")>0.5 )tightestAntiMuWP = 2; 
    if( aTau.tauID("againstMuonTight")>0.5 )tightestAntiMuWP = 3; 
    tightestAntiMuWP_.push_back(tightestAntiMuWP);

    int tightestAntiMu2WP = 0;   
    if( aTau.tauID("againstMuonLoose2")>0.5 )tightestAntiMu2WP = 1;   
    if( aTau.tauID("againstMuonMedium2")>0.5 )tightestAntiMu2WP = 2;   
    if( aTau.tauID("againstMuonTight2")>0.5 )tightestAntiMu2WP = 3;   
    tightestAntiMu2WP_.push_back(tightestAntiMu2WP);

    int tightestAntiMu3WP = 0;   
    if( aTau.tauID("againstMuonLoose3")>0.5 )tightestAntiMu3WP = 1;   
    if( aTau.tauID("againstMuonTight3")>0.5 )tightestAntiMu3WP = 2;   
    tightestAntiMu3WP_.push_back(tightestAntiMu3WP);

    int tightestAntiMuMVAWP = 0;   
    if( aTau.tauID("againstMuonLooseMVA")>0.5 )tightestAntiMuMVAWP = 1;   
    if( aTau.tauID("againstMuonMediumMVA")>0.5 )tightestAntiMuMVAWP = 2;   
    if( aTau.tauID("againstMuonTightMVA")>0.5 )tightestAntiMuMVAWP = 3;   
    tightestAntiMuMVAWP_.push_back(tightestAntiMuMVAWP);
    AntiMuMVAraw_.push_back(aTau.tauID("againstMuonMVAraw"));

    int tightestAntiECutWP = 0;
    if( aTau.tauID("againstElectronLoose")>0.5 )tightestAntiECutWP = 1;
    if( aTau.tauID("againstElectronMedium")>0.5 )tightestAntiECutWP = 2;
    if( aTau.tauID("againstElectronTight")>0.5 )tightestAntiECutWP = 3;
    tightestAntiECutWP_.push_back(tightestAntiECutWP);

    int tightestAntiEMVA5WP = 0;
    if( aTau.tauID("againstElectronVLooseMVA5")>0.5) tightestAntiEMVA5WP  = 1;
    if( aTau.tauID("againstElectronLooseMVA5")>0.5)  tightestAntiEMVA5WP  = 2;
    if( aTau.tauID("againstElectronMediumMVA5")>0.5) tightestAntiEMVA5WP  = 3;
    if( aTau.tauID("againstElectronTightMVA5")>0.5)  tightestAntiEMVA5WP  = 4;
    if( aTau.tauID("againstElectronVTightMVA5")>0.5) tightestAntiEMVA5WP  = 5;
    tightestAntiEMVA5WP_.push_back(tightestAntiEMVA5WP);
    AntiEMVA5raw_.push_back(aTau.tauID("againstElectronMVA5raw"));
    AntiEMVA5category_.push_back(aTau.tauID("againstElectronMVA5category"));

    
    }
  tree_->Fill();

}

void  TauIDCrossCheck::beginJob()
{
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tau tree");

  tree_->Branch("run",&run_,"run/l");
  tree_->Branch("event",&event_,"event/l");
  tree_->Branch("lumi",&lumi_,"lumi/l");

  tree_->Branch("TauP4","std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",&TauP4_);
  tree_->Branch("decayModeFinding", "std::vector< int >", &decayModeFinding_);
  tree_->Branch("decayModeFindingNewDM", "std::vector< int >", &decayModeFindingNewDM_);
  tree_->Branch("decayModeFindingOldDM", "std::vector< int >", &decayModeFindingOldDM_);
  tree_->Branch("tightestAntiECutWP", "std::vector< int >", &tightestAntiECutWP_);
  tree_->Branch("tightestAntiEMVA5WP", "std::vector< int >", &tightestAntiEMVA5WP_);
  tree_->Branch("tightestAntiMuWP", "std::vector< int >", &tightestAntiMuWP_);
  tree_->Branch("tightestAntiMu2WP", "std::vector< int >", &tightestAntiMu2WP_);
  tree_->Branch("tightestAntiMu3WP", "std::vector< int >", &tightestAntiMu3WP_);
  tree_->Branch("tightestAntiMuMVAWP", "std::vector< int >", &tightestAntiMuMVAWP_);
  tree_->Branch("tightestHPSDB3HWP", "std::vector< int >", &tightestHPSDB3HWP_);
  tree_->Branch("tightestHPSMVA3newDMwLTWP", "std::vector< int >", &tightestHPSMVA3newDMwLTWP_);
  tree_->Branch("tightestHPSMVA3newDMwoLTWP", "std::vector< int >", &tightestHPSMVA3newDMwoLTWP_);
  tree_->Branch("tightestHPSMVA3oldDMwLTWP", "std::vector< int >", &tightestHPSMVA3oldDMwLTWP_);
  tree_->Branch("tightestHPSMVA3oldDMwoLTWP", "std::vector< int >", &tightestHPSMVA3oldDMwoLTWP_);
  tree_->Branch("AntiEMVA5category", "std::vector< int >", &AntiEMVA5category_);
  tree_->Branch("AntiEMVA5raw", "std::vector< float >", &AntiEMVA5raw_);
  tree_->Branch("AntiMuMVAraw", "std::vector< float >", &AntiMuMVAraw_);
  tree_->Branch("hpsDB3H", "std::vector< float >", &hpsDB3H_);
  tree_->Branch("hpsMVA3newDMwLT", "std::vector< float >", &hpsMVA3newDMwLT_);
  tree_->Branch("hpsMVA3newDMwoLT", "std::vector< float >", &hpsMVA3newDMwoLT_);
  tree_->Branch("hpsMVA3oldDMwLT", "std::vector< float >", &hpsMVA3oldDMwLT_);
  tree_->Branch("hpsMVA3oldDMwoLT", "std::vector< float >", &hpsMVA3oldDMwoLT_);
  
}

void  TauIDCrossCheck::endJob(){}

DEFINE_FWK_MODULE( TauIDCrossCheck);
