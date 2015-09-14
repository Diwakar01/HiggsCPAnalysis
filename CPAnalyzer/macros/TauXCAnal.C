#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TString.h"
#include "TPad.h"
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCut.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/Vector2D.h"
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <Math/SMatrixDfwd.h>

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;

using namespace std;

void TauXCAnal(TString inputFilePath_)
{

  TChain *chain = new TChain("xcAnalyzer/tree");
  chain->Add(inputFilePath_+"/*.root");

  long int run_, lumi_, event_;
  std::vector< LV > *TauP4_ = new std::vector< LV > (); 

  //std::vector< int > decayModeFinding_;
  //std::vector< int > decayModeFindingNewDM_;
  std::vector< int > *decayModeFindingOldDM_ = new std::vector< int >();
  //std::vector< int > tightestAntiECutWP_;
  //std::vector< int > tightestAntiEMVA5WP_;
  //std::vector< int > tightestAntiMuWP_;
  //std::vector< int > tightestAntiMu2WP_;
  //std::vector< int > tightestAntiMu3WP_;
  //std::vector< int > tightestAntiMuMVAWP_;
  //std::vector< int > tightestHPSDB3HWP_;
  //std::vector< int > tightestHPSMVA3newDMwLTWP_;
  //std::vector< int > tightestHPSMVA3newDMwoLTWP_;
  std::vector< int > *tightestHPSMVA3oldDMwLTWP_ = new std::vector< int >();
  std::vector< int > *tightestHPSMVA3oldDMwoLTWP_ = new std::vector< int >();
  //std::vector< float > AntiEMVA5raw_;
  //std::vector< int > AntiEMVA5category_;
  //std::vector< float > AntiMuMVAraw_;
  std::vector< float > *hpsDB3H_ = new std::vector< float >();
  //std::vector< float > hpsMVA3newDMwLT_;
  //std::vector< float > hpsMVA3newDMwoLT_;
  //std::vector< float > hpsMVA3oldDMwLT_;
  //std::vector< float > hpsMVA3oldDMwoLT_;

  //chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("run", 1);
  chain->SetBranchStatus("lumi", 1);
  chain->SetBranchStatus("event", 1);
  chain->SetBranchStatus("TauP4", 1);
  chain->SetBranchStatus("decayModeFindingOldDM", 1);
  chain->SetBranchStatus("tightestHPSMVA3oldDMwLTWP", 1);
  chain->SetBranchStatus("tightestHPSMVA3oldDMwoLTWP", 1);

  chain->SetBranchAddress("run", &run_);
  chain->SetBranchAddress("lumi", &lumi_);
  chain->SetBranchAddress("event", &event_);
  chain->SetBranchAddress("TauP4", &TauP4_);
  chain->SetBranchAddress("decayModeFindingOldDM", &decayModeFindingOldDM_);
  chain->SetBranchAddress("tightestHPSMVA3oldDMwLTWP", &tightestHPSMVA3oldDMwLTWP_);
  chain->SetBranchAddress("tightestHPSMVA3oldDMwoLTWP", &tightestHPSMVA3oldDMwoLTWP_);
  chain->SetBranchAddress("hpsDB3H", &hpsDB3H_);

  int nEntries    = chain->GetEntries() ;
  cout<< " Number of events "<<nEntries<<endl;
  for (int n = 0; n < nEntries ; n++) {
    chain->GetEntry(n);

    if(event_ != 640031 
       && event_ != 642622
       && event_ != 643786
       && event_ != 644270
       && event_ != 644462
       && event_ != 644653
       && event_ != 646058
       && event_ != 647533
       && event_ != 648784
       && event_ != 649196
       && event_ != 649436
       && event_ != 649514
       && event_ != 651131
       && event_ != 653385
       && event_ != 653426
       && event_ != 655176
       && event_ != 655188
       && event_ != 655390
       && event_ != 657857
       && event_ != 657936
       && event_ != 660152
       && event_ != 660605
       && event_ != 661614
       && event_ != 5728
       && event_ != 666375
       && event_ != 669421
       && event_ != 669809
       && event_ != 671833
       && event_ != 673308
       && event_ != 673536
       && event_ != 674170
       && event_ != 674440)
      continue;
    
    cout<<"run  "<<run_<<" lumi "<<lumi_<<" event "<<event_<<endl;
    
    for(size_t i = 0; i < TauP4_->size(); i++){
      if((*TauP4_)[i].pt() > 45 && TMath::Abs((*TauP4_)[i].eta()) < 2.1) {
	if((*decayModeFindingOldDM_)[i] > 0.5){
	  double pt = (*TauP4_)[i].pt();
	  double eta = (*TauP4_)[i].eta(); 
	  std::cout<<"Tau  pt "<< pt <<" eta "<< eta <<" tightestHPSMVA3oldDMwLTWP "<<(*tightestHPSMVA3oldDMwLTWP_)[i]<<" tightestHPSMVA3oldDMwoLTWP "<<(*tightestHPSMVA3oldDMwoLTWP_)[i]<<std::endl;
	  std::cout<<"hps DB3H "<<(*hpsDB3H_)[i]<<std::endl;
	} 
      }
    }



  }

  delete TauP4_; delete decayModeFindingOldDM_; delete tightestHPSMVA3oldDMwLTWP_;
  delete tightestHPSMVA3oldDMwoLTWP_;
}

void TauXCAnalAll()
{
  TauXCAnal("/nfs/dust/cms/user/anayak/CMS/Ntuple_HttAnalysis/XCEventsSUSYGGH130V2/");
  //TauXCAnal("/nfs/dust/cms/user/anayak/CMS/CMSSW_5_3_13_patch3/src/MyRootMaker/MyRootMaker/test/");
}
