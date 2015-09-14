//--------------------------------------------------------------------------------------------------
// VertexAnalysis
//
// Macro to study Vertex parameters for Higgs CP-studies
// Authors: A. Nayak
//--------------------------------------------------------------------------------------------------
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
#include "TRandom3.h"

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
typedef ROOT::Math::XYZPointD Point3D;
typedef ROOT::Math::XYZVectorD PV;
typedef ROOT::Math::XYVectorD Vector2D;
typedef ROOT::Math::SMatrixSym3D SymMatrix33;
typedef ROOT::Math::SVector<double,3> AlgebraicVector3;

int verbosity = 0;
using namespace std;

TRandom3 randx(123456789);
TRandom3 randy(234567897);
TRandom3 randz(345678919);
TRandom3 randipx(456789123);
TRandom3 randipy(567891235);
TRandom3 randinx(678912357);
TRandom3 randiny(789123459);
TRandom3 randipz(891234567);
TRandom3 randinz(912345679);
TRandom3 randpcapx(987654329);
TRandom3 randpcapy(987654327);
TRandom3 randpcapz(987654325);
TRandom3 randpcanx(987654323);
TRandom3 randpcany(987654319);
TRandom3 randpcanz(987654317);

TH1F * DrawOverflow(TH1F *h)
{
  // This function paint the histogram h with an extra bin for overflows
  UInt_t nx    = h->GetNbinsX()+1;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++)
    xbins[i]=h->GetBinLowEdge(i+1);
  xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
  char *tempName= new char[strlen(h->GetName())+10];
  sprintf(tempName,"%swtOverFlow",h->GetName());
  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
  // Reset the axis labels
  htmp->SetXTitle(h->GetXaxis()->GetTitle());
  htmp->SetYTitle(h->GetYaxis()->GetTitle());
  // Fill the new hitogram including the extra bin for overflows
  for (UInt_t i=1; i<=nx; i++){
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
    htmp->SetBinError(i, h->GetBinError(i));
  }

  // Fill the underflows
  htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
  htmp->SetBinError(0, h->GetBinError(0));
  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());

  // FillStyle and color
  htmp->SetLineWidth(h->GetLineWidth());
  htmp->SetLineColor(h->GetLineColor());
  htmp->SetFillStyle(h->GetFillStyle());
  htmp->SetFillColor(h->GetFillColor());

  return htmp;
}

std::map<std::string, TH1F*>  createHistograms(){

  std::map<std::string, TH1F*> histos_;

  histos_["DeltaPVX"] = new TH1F("DeltaPVX", "", 1000, -0.1, 0.1); 
  histos_["DeltaPVY"] = new TH1F("DeltaPVY", "", 1000, -0.1, 0.1); 
  histos_["DeltaPVZ"] = new TH1F("DeltaPVZ", "", 1000, -0.1, 0.1); 

  histos_["DeltaPVwrtGenX"] = new TH1F("DeltaPVwrtGenX", "", 1000, -0.1, 0.1); 
  histos_["DeltaPVwrtGenY"] = new TH1F("DeltaPVwrtGenY", "", 1000, -0.1, 0.1); 
  histos_["DeltaPVwrtGenZ"] = new TH1F("DeltaPVwrtGenZ", "", 1000, -0.1, 0.1); 

  histos_["DeltaRfPVwrtGenX"] = new TH1F("DeltaRfPVwrtGenX", "", 1000, -0.1, 0.1); 
  histos_["DeltaRfPVwrtGenY"] = new TH1F("DeltaRfPVwrtGenY", "", 1000, -0.1, 0.1); 
  histos_["DeltaRfPVwrtGenZ"] = new TH1F("DeltaRfPVwrtGenZ", "", 1000, -0.1, 0.1); 

  histos_["DeltaIanPVX"] = new TH1F("DeltaIanPVX", "", 1000, -0.1, 0.1);
  histos_["DeltaIanPVY"] = new TH1F("DeltaIanPVY", "", 1000, -0.1, 0.1);
  histos_["DeltaIanPVZ"] = new TH1F("DeltaIanPVZ", "", 1000, -0.1, 0.1);

  histos_["DeltaIanPVwrtGenX"] = new TH1F("DeltaIanPVwrtGenX", "", 1000, -0.1, 0.1);
  histos_["DeltaIanPVwrtGenY"] = new TH1F("DeltaIanPVwrtGenY", "", 1000, -0.1, 0.1);
  histos_["DeltaIanPVwrtGenZ"] = new TH1F("DeltaIanPVwrtGenZ", "", 1000, -0.1, 0.1);

  histos_["ResPVwrtGenX"] = new TH1F("ResPVwrtGenX", "", 100, -10.0, 10.0); 
  histos_["ResPVwrtGenY"] = new TH1F("ResPVwrtGenY", "", 100, -10.0, 10.0); 
  histos_["ResPVwrtGenZ"] = new TH1F("ResPVwrtGenZ", "", 100, -10.0, 10.0); 

  histos_["ResRfPVwrtGenX"] = new TH1F("ResRfPVwrtGenX", "", 100, -10.0, 10.0); 
  histos_["ResRfPVwrtGenY"] = new TH1F("ResRfPVwrtGenY", "", 100, -10.0, 10.0); 
  histos_["ResRfPVwrtGenZ"] = new TH1F("ResRfPVwrtGenZ", "", 100, -10.0, 10.0); 

  histos_["VtxXRes"] = new TH1F("VtxXRes",  "", 1000, -0.1, 0.1); 
  histos_["VtxYRes"] = new TH1F("VtxYRes",  "", 1000, -0.1, 0.1); 
  histos_["VtxZRes"] = new TH1F("VtxZRes",  "", 1000, -0.1, 0.1); 

  histos_["ReFitVtxXRes"] = new TH1F("ReFitVtxXRes",  "", 1000, -0.1, 0.1); 
  histos_["ReFitVtxYRes"] = new TH1F("ReFitVtxYRes",  "", 1000, -0.1, 0.1); 
  histos_["ReFitVtxZRes"] = new TH1F("ReFitVtxZRes",  "", 1000, -0.1, 0.1); 
  
  histos_["VtxXErr"] = new TH1F("VtxXErr",  "", 1000, -0.05, 0.05); 
  histos_["VtxYErr"] = new TH1F("VtxYErr",  "", 1000, -0.05, 0.05); 
  histos_["VtxZErr"] = new TH1F("VtxZErr",  "", 1000, -0.05, 0.05); 

  histos_["ReFitVtxXErr"] = new TH1F("ReFitVtxXErr",  "", 1000, -0.05, 0.05); 
  histos_["ReFitVtxYErr"] = new TH1F("ReFitVtxYErr",  "", 1000, -0.05, 0.05); 
  histos_["ReFitVtxZErr"] = new TH1F("ReFitVtxZErr",  "", 1000, -0.05, 0.05); 

  histos_["VtxNChi2"] = new TH1F("VtxNChi2", "", 1000, 0., 100.);
  histos_["ReFitVtxNChi2"] =new TH1F("ReFitVtxNChi2", "", 1000, 0., 100.);
  histos_["IanVtxNChi2"] =new TH1F("IanVtxNChi2", "", 1000, 0., 100.);

  histos_["DeltaX_gen_reco_pca_genvtx"] = new TH1F("DeltaX_gen_reco_pca_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaY_gen_reco_pca_genvtx"] = new TH1F("DeltaY_gen_reco_pca_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaZ_gen_reco_pca_genvtx"] = new TH1F("DeltaZ_gen_reco_pca_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaX_gen_reco_pcaLE1_genvtx"] = new TH1F("DeltaX_gen_reco_pcaLE1_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaY_gen_reco_pcaLE1_genvtx"] = new TH1F("DeltaY_gen_reco_pcaLE1_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaZ_gen_reco_pcaLE1_genvtx"] = new TH1F("DeltaZ_gen_reco_pcaLE1_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaX_gen_reco_pcaLE2_genvtx"] = new TH1F("DeltaX_gen_reco_pcaLE2_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaY_gen_reco_pcaLE2_genvtx"] = new TH1F("DeltaY_gen_reco_pcaLE2_genvtx", "", 1000, -0.1, 0.1); 
  histos_["DeltaZ_gen_reco_pcaLE2_genvtx"] = new TH1F("DeltaZ_gen_reco_pcaLE2_genvtx", "", 1000, -0.1, 0.1); 

  histos_["CPPhiStar"] = new TH1F("CPPhiStar",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab"] = new TH1F("CPPhiLab",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv"] = new TH1F("CPPhiStar_opv",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_opv"] = new TH1F("CPPhiLab_opv",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_bs"] = new TH1F("CPPhiStar_bs",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_bs"] = new TH1F("CPPhiLab_bs",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx"] = new TH1F("CPPhiStar_genvtx",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genvtx"] = new TH1F("CPPhiLab_genvtx",  "", 18, 0.0, 3.14159); 
  
  histos_["CPPhiStar_NP"] = new TH1F("CPPhiStar_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_NP"] = new TH1F("CPPhiLab_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_NP"] = new TH1F("CPPhiStar_opv_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_opv_NP"] = new TH1F("CPPhiLab_opv_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_bs_NP"] = new TH1F("CPPhiStar_bs_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_bs_NP"] = new TH1F("CPPhiLab_bs_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_NP"] = new TH1F("CPPhiStar_genvtx_NP",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genvtx_NP"] = new TH1F("CPPhiLab_genvtx_NP",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_LE1"] = new TH1F("CPPhiStar_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_LE1"] = new TH1F("CPPhiLab_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_LE2"] = new TH1F("CPPhiStar_LE2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_LE2"] = new TH1F("CPPhiLab_LE2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_LE1"] = new TH1F("CPPhiStar_opv_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_opv_LE1"] = new TH1F("CPPhiLab_opv_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_LE2"] = new TH1F("CPPhiStar_opv_LE2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_opv_LE2"] = new TH1F("CPPhiLab_opv_LE2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_LE1"] = new TH1F("CPPhiStar_genvtx_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genvtx_LE1"] = new TH1F("CPPhiLab_genvtx_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_LE2"] = new TH1F("CPPhiStar_genvtx_LE2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genvtx_LE2"] = new TH1F("CPPhiLab_genvtx_LE2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_tauvtx_LE1"] = new TH1F("CPPhiStar_tauvtx_LE1",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_tauvtx_LE1"] = new TH1F("CPPhiLab_tauvtx_LE1",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_M2"] = new TH1F("CPPhiStar_M2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_M2"] = new TH1F("CPPhiLab_M2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_Ian"] = new TH1F("CPPhiStar_Ian",  "", 18, 0.0, 3.14159);

  histos_["CPPhiStar_CP"] = new TH1F("CPPhiStar_CP",  "", 35, 0.0, 7.0); 
  histos_["CPPhiStar_CP_LE1"] = new TH1F("CPPhiStar_CP_LE1",  "", 35, 0.0, 7.0); 
  histos_["CPPhiStar_opv_CP"] = new TH1F("CPPhiStar_opv_CP",  "", 35, 0.0, 7.0); 
  histos_["CPPhiStar_opv_CP_LE1"] = new TH1F("CPPhiStar_opv_CP_LE1",  "", 35, 0.0, 7.0); 
  histos_["CPPhiStar_genvtx_CP"] = new TH1F("CPPhiStar_genvtx_CP",  "", 35, 0.0, 7.0); 
  histos_["CPPhiStar_genvtx_CP_LE1"] = new TH1F("CPPhiStar_genvtx_CP_LE1",  "", 35, 0.0, 7.0); 
  histos_["CPPhiStar_tauvtx_CP_LE1"] = new TH1F("CPPhiStar_tauvtx_CP_LE1",  "", 35, 0.0, 7.0); 

  histos_["CPPhiStar_P4AtPCA"] = new TH1F("CPPhiStar_P4AtPCA",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_P4AtPCA"] = new TH1F("CPPhiStar_opv_P4AtPCA",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_bs_P4AtPCA"] = new TH1F("CPPhiStar_bs_P4AtPCA",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_P4AtPCA"] = new TH1F("CPPhiStar_genvtx_P4AtPCA",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_genMatch"] = new TH1F("CPPhiStar_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genMatch"] = new TH1F("CPPhiLab_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_genMatch"] = new TH1F("CPPhiStar_opv_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_bs_genMatch"] = new TH1F("CPPhiStar_bs_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_genMatch"] = new TH1F("CPPhiStar_genvtx_genMatch",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_NP_genMatch"] = new TH1F("CPPhiStar_NP_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_NP_genMatch"] = new TH1F("CPPhiLab_NP_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_NP_genMatch"] = new TH1F("CPPhiStar_opv_NP_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_bs_NP_genMatch"] = new TH1F("CPPhiStar_bs_NP_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_NP_genMatch"] = new TH1F("CPPhiStar_genvtx_NP_genMatch",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_LE1_genMatch"] = new TH1F("CPPhiStar_LE1_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_LE2_genMatch"] = new TH1F("CPPhiStar_LE2_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_LE1_genMatch"] = new TH1F("CPPhiStar_opv_LE1_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_LE2_genMatch"] = new TH1F("CPPhiStar_opv_LE2_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_LE1_genMatch"] = new TH1F("CPPhiStar_genvtx_LE1_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_LE2_genMatch"] = new TH1F("CPPhiStar_genvtx_LE2_genMatch",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_P4AtPCA_genMatch"] = new TH1F("CPPhiStar_P4AtPCA_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_P4AtPCA_genMatch"] = new TH1F("CPPhiStar_opv_P4AtPCA_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_bs_P4AtPCA_genMatch"] = new TH1F("CPPhiStar_bs_P4AtPCA_genMatch",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_P4AtPCA_genMatch"] = new TH1F("CPPhiStar_genvtx_P4AtPCA_genMatch",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_ipcut"] = new TH1F("CPPhiStar_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_ipcut"] = new TH1F("CPPhiLab_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_ipcut"] = new TH1F("CPPhiStar_opv_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_opv_ipcut"] = new TH1F("CPPhiLab_opv_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_ipcut"] = new TH1F("CPPhiStar_genvtx_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genvtx_ipcut"] = new TH1F("CPPhiLab_genvtx_ipcut",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_ipcut_TkQCut"] = new TH1F("CPPhiStar_ipcut_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_ipcut_TkQCut"] = new TH1F("CPPhiLab_ipcut_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_ipcut_TkQCutV2"] = new TH1F("CPPhiStar_ipcut_TkQCutV2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_ipcut_TkQCutV3"] = new TH1F("CPPhiStar_ipcut_TkQCutV3",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_ipcut30"] = new TH1F("CPPhiStar_ipcut30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_ipcut30"] = new TH1F("CPPhiLab_ipcut30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_ipcut30_TkQCut"] = new TH1F("CPPhiStar_ipcut30_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_ipcut30_TkQCut"] = new TH1F("CPPhiLab_ipcut30_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_ipcut30_TkQCutV2"] = new TH1F("CPPhiStar_ipcut30_TkQCutV2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ipcut_TkQCut"] = new TH1F("CPPhiStar_genvtx_ipcut_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_genvtx_ipcut_TkQCut"] = new TH1F("CPPhiLab_genvtx_ipcut_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ipcut_TkQCutV2"] = new TH1F("CPPhiStar_genvtx_ipcut_TkQCutV2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ipcut_TkQCutV3"] = new TH1F("CPPhiStar_genvtx_ipcut_TkQCutV3",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ipcut30"] = new TH1F("CPPhiStar_genvtx_ipcut30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_genvtx_ipcut30"] = new TH1F("CPPhiLab_genvtx_ipcut30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ipcut30_TkQCut"] = new TH1F("CPPhiStar_genvtx_ipcut30_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiLab_genvtx_ipcut30_TkQCut"] = new TH1F("CPPhiLab_genvtx_ipcut30_TkQCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ipcut30_TkQCutV2"] = new TH1F("CPPhiStar_genvtx_ipcut30_TkQCutV2",  "", 18, 0.0, 3.14159);

  histos_["CPPhiStar_M2_gen1p0pi0"] = new TH1F("CPPhiStar_M2_gen1p0pi0",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_LE1_gen1p0pi0"] = new TH1F("CPPhiStar_LE1_gen1p0pi0",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_gen1p0pi0"] = new TH1F("CPPhiStar_genvtx_gen1p0pi0",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_M2_gentau"] = new TH1F("CPPhiStar_M2_gentau",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_LE1_gentau"] = new TH1F("CPPhiStar_LE1_gentau",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_gentau"] = new TH1F("CPPhiStar_genvtx_gentau",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_tau30"] = new TH1F("CPPhiStar_tau30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_tau30"] = new TH1F("CPPhiStar_genvtx_tau30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_tau40"] = new TH1F("CPPhiStar_tau40",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_tau40"] = new TH1F("CPPhiStar_genvtx_tau40",  "", 18, 0.0, 3.14159);
  //histos_["CPPhiStar_LE1_tau30"] = new TH1F("CPPhiStar_LE1_tau30",  "", 18, 0.0, 3.14159);
  //histos_["CPPhiStar_LE1_tau40"] = new TH1F("CPPhiStar_LE1_tau40",  "", 18, 0.0, 3.14159);
  //histos_["CPPhiStar_LE1_tau50"] = new TH1F("CPPhiStar_LE1_tau50",  "", 18, 0.0, 3.14159);

  histos_["CPPhiStar_ResCut"] = new TH1F("CPPhiStar_ResCut",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_genvtx_ResCut"] = new TH1F("CPPhiStar_genvtx_ResCut",  "", 18, 0.0, 3.14159);

  histos_["CPPhiStar_NP_ipcut"] = new TH1F("CPPhiStar_NP_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_NP_ipcut"] = new TH1F("CPPhiLab_NP_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_opv_NP_ipcut"] = new TH1F("CPPhiStar_opv_NP_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_opv_NP_ipcut"] = new TH1F("CPPhiLab_opv_NP_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_genvtx_NP_ipcut"] = new TH1F("CPPhiStar_genvtx_NP_ipcut",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiLab_genvtx_NP_ipcut"] = new TH1F("CPPhiLab_genvtx_NP_ipcut",  "", 18, 0.0, 3.14159); 

  histos_["CPPhiStar_gen"] = new TH1F("CPPhiStar_gen",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_gen_CP"] = new TH1F("CPPhiStar_gen_CP",  "", 35, 0.0, 7.0); 
  histos_["CPPhiLab_gen"] = new TH1F("CPPhiLab_gen",  "", 18, 0.0, 3.14159); 
  histos_["CPPhi_gen"] = new TH1F("CPPhi_gen",  "", 18, 0.0, 3.14159); 
  histos_["CPPhi_gen_v2"] = new TH1F("CPPhi_gen_v2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_v2"] = new TH1F("CPPhiStar_gen_v2",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiPrime_gen"] = new TH1F("CPPhiPrime_gen",  "", 18, 0.0, 3.14159); 
  histos_["CPPsiPrime_gen"] = new TH1F("CPPsiPrime_gen",  "", 18, 0.0, 3.14159); 
  histos_["CPPhiStar_gen_vs"] = new TH1F("CPPhiStar_gen_vs",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vs_ip30"] = new TH1F("CPPhiStar_gen_vs_ip30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vs_ip50"] = new TH1F("CPPhiStar_gen_vs_ip50",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vst"] = new TH1F("CPPhiStar_gen_vst",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vst_ip30"] = new TH1F("CPPhiStar_gen_vst_ip30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vst_ip50"] = new TH1F("CPPhiStar_gen_vst_ip50",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vips"] = new TH1F("CPPhiStar_gen_vips",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vips_ip30"] = new TH1F("CPPhiStar_gen_vips_ip30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vips_ip50"] = new TH1F("CPPhiStar_gen_vips_ip50",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2"] = new TH1F("CPPhiStar_gen_vipsv2",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_ip30"] = new TH1F("CPPhiStar_gen_vipsv2_ip30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_ip50"] = new TH1F("CPPhiStar_gen_vipsv2_ip50",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_tau20"] = new TH1F("CPPhiStar_gen_vipsv2_tau20",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_tau30"] = new TH1F("CPPhiStar_gen_vipsv2_tau30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_tau40"] = new TH1F("CPPhiStar_gen_vipsv2_tau40",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_tau50"] = new TH1F("CPPhiStar_gen_vipsv2_tau50",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_ip30_tau40"] = new TH1F("CPPhiStar_gen_vipsv2_ip30_tau40",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_vipsv2_ip50_tau40"] = new TH1F("CPPhiStar_gen_vipsv2_ip50_tau40",  "", 18, 0.0, 3.14159);

  histos_["CPPhiStar_gen_pcas"] = new TH1F("CPPhiStar_gen_pcas",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_pcas_ip30"] = new TH1F("CPPhiStar_gen_pcas_ip30",  "", 18, 0.0, 3.14159);
  histos_["CPPhiStar_gen_pcas_ip50"] = new TH1F("CPPhiStar_gen_pcas_ip50",  "", 18, 0.0, 3.14159);

  histos_["xcheck_angle1"] = new TH1F("xcheck_angle1", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle2"] = new TH1F("xcheck_angle2", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3"] = new TH1F("xcheck_angle3", "", 40, 0.0, 4.0); 

  histos_["xcheck_angle3_reco"] = new TH1F("xcheck_angle3_reco", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_opv"] = new TH1F("xcheck_angle3_reco_opv", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_bs"] = new TH1F("xcheck_angle3_reco_bs", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genvtx"] = new TH1F("xcheck_angle3_reco_genvtx", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_pca"] = new TH1F("xcheck_angle3_reco_pca", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_opv_pca"] = new TH1F("xcheck_angle3_reco_opv_pca", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_bs_pca"] = new TH1F("xcheck_angle3_reco_bs_pca", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genvtx_pca"] = new TH1F("xcheck_angle3_reco_genvtx_pca", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_2d"] = new TH1F("xcheck_angle3_reco_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_opv_2d"] = new TH1F("xcheck_angle3_reco_opv_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_bs_2d"] = new TH1F("xcheck_angle3_reco_bs_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genvtx_2d"] = new TH1F("xcheck_angle3_reco_genvtx_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_pca_2d"] = new TH1F("xcheck_angle3_reco_pca_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_opv_pca_2d"] = new TH1F("xcheck_angle3_reco_opv_pca_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_bs_pca_2d"] = new TH1F("xcheck_angle3_reco_bs_pca_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genvtx_pca_2d"] = new TH1F("xcheck_angle3_reco_genvtx_pca_2d", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genP4"] = new TH1F("xcheck_angle3_reco_genP4", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_opv_genP4"] = new TH1F("xcheck_angle3_reco_opv_genP4", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_bs_genP4"] = new TH1F("xcheck_angle3_reco_bs_genP4", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genvtx_genP4"] = new TH1F("xcheck_angle3_reco_genvtx_genP4", "", 40, 0.0, 4.0); 
  histos_["xcheck_angle3_reco_genIP"] = new TH1F("xcheck_angle3_reco_genIP", "", 400, 0.0, 4.0);

  histos_["angle_dp_refitVtx"] = new TH1F("angle_dp_refitVtx", "", 200, 0.0, 4.0);
  histos_["angle_dp_origVtx"] = new TH1F("angle_dp_origVtx", "", 200, 0.0, 4.0);
  histos_["angle_dp_genVtx"] = new TH1F("angle_dp_genVtx", "", 200, 0.0, 4.0);
  histos_["angle_dp_refitVtx_genP4"] = new TH1F("angle_dp_refitVtx_genP4", "", 200, 0.0, 4.0);
  histos_["angle_dp_origVtx_genP4"] = new TH1F("angle_dp_origVtx_genP4", "", 200, 0.0, 4.0);
  histos_["angle_dp_genVtx_genP4"] = new TH1F("angle_dp_genVtx_genP4", "", 200, 0.0, 4.0);
  histos_["angle_dp_genIP"] = new TH1F("angle_dp_genIP", "", 400, 0.0, 4.0);
  histos_["angle_dp_refitVtx_NP"] = new TH1F("angle_dp_refitVtx_NP", "", 200, 0.0, 4.0);
  histos_["angle_dp_genVtx_NP"] = new TH1F("angle_dp_genVtx_NP", "", 200, 0.0, 4.0);
  histos_["delphi_dp_genVtx"] = new TH1F("delphi_dp_genVtx", "", 200, 0.0, 4.0);

  histos_["angle_dp_refitVtx_LE1"] = new TH1F("angle_dp_refitVtx_LE1", "", 400, -4.0, 4.0);
  histos_["angle_dp_origVtx_LE1"] = new TH1F("angle_dp_origVtx_LE1", "", 400, -4.0, 4.0);
  histos_["angle_dp_genVtx_LE1"] = new TH1F("angle_dp_genVtx_LE1", "", 400, -4.0, 4.0);

  histos_["DeltaPt_Gen_Reco"] = new TH1F("DeltaPt_Gen_Reco", "", 100, 0.0, 10.0);
  histos_["DeltaP_Gen_Reco"] = new TH1F("DeltaP_Gen_Reco", "", 100, 0.0, 10.0);
  histos_["DeltaPtOverPt_Gen_Reco"] = new TH1F("DeltaPtOverPt_Gen_Reco", "", 1000, 0.0, 1.0);
  histos_["DeltaPOverP_Gen_Reco"] = new TH1F("DeltaPOverP_Gen_Reco", "", 1000, 0.0, 1.0);

  histos_["angle_pion_momentum"] = new TH1F("angle_pion_momentum", "", 1000, 0.0, 0.1);

  histos_["DeltaTkLength_gen_reco_genvtx_le1_l"] = new TH1F("DeltaTkLength_gen_reco_genvtx_le1_l", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_genvtx_le1_t"] = new TH1F("DeltaTkLength_gen_reco_genvtx_le1_t", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_genvtx_le2_l"] = new TH1F("DeltaTkLength_gen_reco_genvtx_le2_l", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_genvtx_le2_t"] = new TH1F("DeltaTkLength_gen_reco_genvtx_le2_t", "", 2000, -0.1, 0.1);
  histos_["AngleTkLength_gen_reco_genvtx_le1_l"] = new TH1F("AngleTkLength_gen_reco_genvtx_le1_l", "", 400, 0.0, 4.0);
  histos_["AngleTkLength_gen_reco_genvtx_le1_t"] = new TH1F("AngleTkLength_gen_reco_genvtx_le1_t", "", 400, 0.0, 4.0);
  histos_["AngleTkLength_gen_reco_genvtx_le2_l"] = new TH1F("AngleTkLength_gen_reco_genvtx_le2_l", "", 400, 0.0, 4.0);
  histos_["AngleTkLength_gen_reco_genvtx_le2_t"] = new TH1F("AngleTkLength_gen_reco_genvtx_le2_t", "", 400, 0.0, 4.0);
  histos_["TkLength_gen_reco_genvtx_le1_tdp"] = new TH1F("TkLength_gen_reco_genvtx_le1_tdp", "", 2000, -0.1, 0.1);
  
  histos_["DeltaTkLength_gen_reco_le1_l"] = new TH1F("DeltaTkLength_gen_reco_le1_l", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_le1_t"] = new TH1F("DeltaTkLength_gen_reco_le1_t", "", 2000, -0.1, 0.1);
  histos_["AngleTkLength_gen_reco_le1_l"] = new TH1F("AngleTkLength_gen_reco_le1_l", "", 400, 0.0, 4.0);
  histos_["AngleTkLength_gen_reco_le1_t"] = new TH1F("AngleTkLength_gen_reco_le1_t", "", 400, 0.0, 4.0);
  histos_["TkLength_gen_reco_le1_tdp"] = new TH1F("TkLength_gen_reco_le1_tdp", "", 2000, -0.1, 0.1);
  
  histos_["DeltaTkLength_gen_reco_m2_l"] = new TH1F("DeltaTkLength_gen_reco_m2_l", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_m2_t"] = new TH1F("DeltaTkLength_gen_reco_m2_t", "", 2000, -0.1, 0.1);
  histos_["AngleTkLength_gen_reco_m2_l"] = new TH1F("AngleTkLength_gen_reco_m2_l", "", 400, 0.0, 4.0);
  histos_["AngleTkLength_gen_reco_m2_t"] = new TH1F("AngleTkLength_gen_reco_m2_t", "", 400, 0.0, 4.0);
  histos_["TkLength_gen_reco_m2_tdp"] = new TH1F("TkLength_gen_reco_m2_tdp", "", 2000, -0.1, 0.1);

  histos_["DeltaTkLength_gen_reco_ip30_l_tkQCut"] = new TH1F("DeltaTkLength_gen_reco_ip30_l_tkQCut", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_ip30_t_tkQCut"] = new TH1F("DeltaTkLength_gen_reco_ip30_t_tkQCut", "", 2000, -0.1, 0.1);
  histos_["TkLength_gen_reco_ip30_tdp_tkQCut"] = new TH1F("TkLength_gen_reco_ip30_tdp_tkQCut", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_ip50_l_tkQCut"] = new TH1F("DeltaTkLength_gen_reco_ip50_l_tkQCut", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_ip50_t_tkQCut"] = new TH1F("DeltaTkLength_gen_reco_ip50_t_tkQCut", "", 2000, -0.1, 0.1);
  histos_["TkLength_gen_reco_ip50_tdp_tkQCut"] = new TH1F("TkLength_gen_reco_ip50_tdp_tkQCut", "", 2000, -0.1, 0.1);

  histos_["DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut"] = new TH1F("DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut", "", 2000, -0.1, 0.1);
  histos_["TkLength_gen_reco_genvtx_ip30_tdp_tkQCut"] = new TH1F("TkLength_gen_reco_genvtx_ip30_tdp_tkQCut", "", 2000, -0.1, 0.1);
  histos_["DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut"] = new TH1F("DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut", "", 2000, -0.1, 0.1);
  histos_["TkLength_gen_reco_genvtx_ip50_tdp_tkQCut"] = new TH1F("TkLength_gen_reco_genvtx_ip50_tdp_tkQCut", "", 2000, -0.1, 0.1);

  histos_["ip_error_m2_l"] = new TH1F("ip_error_m2_l", "", 2000, -0.1, 0.1);
  histos_["ip_error_m2_t"] = new TH1F("ip_error_m2_t", "", 2000, -0.1, 0.1);
  histos_["ip_error_m2_tdp"] = new TH1F("ip_error_m2_tdp", "", 2000, -0.1, 0.1);

  for(std::map<std::string, TH1F*>::const_iterator it = histos_.begin(); it != histos_.end(); it++){
    (it->second)->Sumw2();
  }
  return histos_;
}
  
void VertexAnalysis(TString inputFilePath_, TString outFileName_, int HPid_)
{

  TChain *chain = new TChain("vertexAnalyzer/tree");
  chain->Add(inputFilePath_+"/*.root");

  TFile *outFile_ = TFile::Open(outFileName_+".root", "RECREATE");
  outFile_->SetCompressionLevel( 9 );

  //std::map<std::string, TH1F*> histos_;
  std::map<std::string, TH1F*> histos_ = createHistograms();
  
  unsigned long run_,event_,lumi_;
  //int index_;

  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *diTauLegsP4_ 
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > (); 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *diTauLegsLchP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCA_ 
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCAM2_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCAOPV_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCABS_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCAGen_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCAGenM2_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCAIan_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTausP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genVP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTauPSonsP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTauNSonsP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  
  std::vector< SymMatrix33 >* diTauLegsPCAM2Cov_ 
    = new std::vector< SymMatrix33 > ();
  std::vector< SymMatrix33 >* diTauLegsPCAGenM2Cov_ 
    = new std::vector< SymMatrix33 > ();

  std::vector< Point3D >* VtxPos_ = new std::vector< Point3D >();
  std::vector< Point3D >* BSPos_ = new std::vector< Point3D >();
  std::vector< Point3D >* ReFitVtxPos_ = new std::vector< Point3D >();
  std::vector< Point3D >* IanVtxPos_ = new std::vector< Point3D >();
  std::vector< Point3D >* HiggsGenVtx_ = new std::vector< Point3D >();
  std::vector< Point3D >* TausGenVtx_ = new std::vector< Point3D >();
  std::vector< Point3D >* TauPSonsGenVtx_ = new std::vector< Point3D >();
  std::vector< Point3D >* TauNSonsGenVtx_ = new std::vector< Point3D >();

  std::vector< PV >* diTauLegsLchP3AtPCA_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsLchP3AtPCAOPV_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsLchP3AtPCABS_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsLchP3AtPCAGen_ = new std::vector< PV >();

  std::vector< PV >* diTauLegsIPAtPCA_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsIPAtPCAV2_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsIPAtPCAOPV_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsIPAtPCAOPVV2_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsIPAtPCAGen_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsIPAtPCAGenV2_ = new std::vector< PV >();
  std::vector< PV >* diTauLegsIPAtPCATauVtx_ = new std::vector< PV >();

  std::vector< int > *genTausPid_ = new std::vector< int >();
  std::vector< int > *genTausCharge_ = new std::vector< int >();
  std::vector< int > *genTausStatus_ = new std::vector< int >();
  std::vector< int > *genTauPSonsPid_ = new std::vector< int >();
  std::vector< int > *genTauPSonsCharge_ = new std::vector< int >();
  std::vector< int > *genTauPSonsStatus_ = new std::vector< int >();
  std::vector< int > *genTauNSonsPid_ = new std::vector< int >();
  std::vector< int > *genTauNSonsCharge_ = new std::vector< int >();
  std::vector< int > *genTauNSonsStatus_ = new std::vector< int >();

  int genVPid_;
  int tightestHPSDB3HWPLeg1_, tightestHPSDB3HWPLeg2_;
  float hpsDB3HLeg1_, hpsDB3HLeg2_, diTauCharge_, chargeLeg1_;
  int decayModeLeg1_, decayModeLeg2_;
  int genDecayModeLeg1_, genDecayModeLeg2_;

  float VtxX_, VtxY_, VtxZ_, VtxXErr_, VtxYErr_, VtxZErr_, VtxNChi2_;
  float genVtxX_, genVtxY_, genVtxZ_;
  
  float ReFitVtxX_, ReFitVtxY_, ReFitVtxZ_, ReFitVtxXErr_, ReFitVtxYErr_, ReFitVtxZErr_;
  int   ReFitVtxRho_;
  float ReFitVtxNdof_, ReFitVtxNChi2_;
  float IanVtxX_, IanVtxY_, IanVtxZ_, IanVtxXErr_, IanVtxYErr_, IanVtxZErr_, IanVtxNChi2_;

  float trackPtLeg1_, trackPtErrLeg1_;
  int nMisingHitsLeg1_, nHitsLeg1_, nTkHitsLeg1_;
  int nPxlHitsLeg1_, hasFirstPxlHitLeg1_;
  float trackPtLeg2_, trackPtErrLeg2_;
  int nMisingHitsLeg2_, nHitsLeg2_, nTkHitsLeg2_;
  int nPxlHitsLeg2_, hasFirstPxlHitLeg2_;

  float dxyLeg1_, dxyLeg2_, dxyLeg1Err_, dxyLeg2Err_;
  float dzLeg1_, dzLeg2_, dzLeg1Err_, dzLeg2Err_;
  float dxyOPVLeg1_, dxyOPVLeg2_, dxyOPVLeg1Err_, dxyOPVLeg2Err_;
  float dzOPVLeg1_, dzOPVLeg2_, dzOPVLeg1Err_, dzOPVLeg2Err_;

  float dxyBSLeg1_, dxyBSLeg2_, dxyBSLeg1Err_, dxyBSLeg2Err_;
  float dzBSLeg1_, dzBSLeg2_, dzBSLeg1Err_, dzBSLeg2Err_;

  float dxyGenLeg1_, dxyGenLeg2_, dxyGenLeg1Err_, dxyGenLeg2Err_;
  float dzGenLeg1_, dzGenLeg2_, dzGenLeg1Err_, dzGenLeg2Err_;

  //Set branches
  chain->SetBranchStatus("run", 1);
  chain->SetBranchStatus("event", 1);
  chain->SetBranchStatus("lumi", 1);
  chain->SetBranchStatus("diTauLegsP4", 1);
  chain->SetBranchStatus("diTauLegsLchP4", 1);
  chain->SetBranchStatus("diTauLegsPCA", 1);
  chain->SetBranchStatus("diTauLegsPCAOPV", 1);
  chain->SetBranchStatus("diTauLegsPCABS", 1);
  chain->SetBranchStatus("diTauLegsPCAGen", 1);
  chain->SetBranchStatus("genTausP4", 1);
  chain->SetBranchStatus("genVP4", 1);
  chain->SetBranchStatus("genVPid", 1);
  chain->SetBranchStatus("genTausPid", 1);
  chain->SetBranchStatus("genTausCharge", 1);
  chain->SetBranchStatus("genTausStatus", 1);
  chain->SetBranchStatus("genTauPSonsP4", 1);
  chain->SetBranchStatus("genTauNSonsP4", 1);
  chain->SetBranchStatus("genTauPSonsPid", 1);
  chain->SetBranchStatus("genTauPSonsCharge", 1);
  chain->SetBranchStatus("genTauPSonsStatus", 1);
  chain->SetBranchStatus("genTauNSonsPid", 1);
  chain->SetBranchStatus("genTauNSonsCharge", 1);
  chain->SetBranchStatus("genTauNSonsStatus", 1);
  chain->SetBranchStatus("VtxPos", 1);
  chain->SetBranchStatus("BSPos", 1);
  chain->SetBranchStatus("ReFitVtxPos", 1);
  chain->SetBranchStatus("HiggsGenVtx", 1);
  chain->SetBranchStatus("TausGenVtx", 1);
  chain->SetBranchStatus("TauPSonsGenVtx", 1);
  chain->SetBranchStatus("TauNSonsGenVtx", 1);
  chain->SetBranchStatus("tightestHPSDB3HWPLeg1", 1);
  chain->SetBranchStatus("tightestHPSDB3HWPLeg2", 1);
  chain->SetBranchStatus("hpsDB3HLeg1", 1);
  chain->SetBranchStatus("hpsDB3HLeg2", 1);
  chain->SetBranchStatus("diTauCharge", 1);
  chain->SetBranchStatus("chargeLeg1", 1);
  chain->SetBranchStatus("decayModeLeg1", 1);
  chain->SetBranchStatus("decayModeLeg2", 1);
  chain->SetBranchStatus("genVtxX", 1);
  chain->SetBranchStatus("genVtxY", 1);
  chain->SetBranchStatus("genVtxZ", 1);
  chain->SetBranchStatus("VtxX", 1);
  chain->SetBranchStatus("VtxY", 1);
  chain->SetBranchStatus("VtxZ", 1);
  chain->SetBranchStatus("VtxXErr", 1);
  chain->SetBranchStatus("VtxYErr", 1);
  chain->SetBranchStatus("VtxZErr", 1);
  chain->SetBranchStatus("VtxNChi2", 1);
  chain->SetBranchStatus("ReFitVtxX", 1);
  chain->SetBranchStatus("ReFitVtxY", 1);
  chain->SetBranchStatus("ReFitVtxZ", 1);
  chain->SetBranchStatus("ReFitVtxXErr", 1);
  chain->SetBranchStatus("ReFitVtxYErr", 1);
  chain->SetBranchStatus("ReFitVtxZErr", 1);
  chain->SetBranchStatus("ReFitVtxRho", 1);
  chain->SetBranchStatus("ReFitVtxNdof", 1);
  chain->SetBranchStatus("ReFitVtxNChi2", 1);
  chain->SetBranchStatus("IanVtxX", 1);
  chain->SetBranchStatus("IanVtxY", 1);
  chain->SetBranchStatus("IanVtxZ", 1);
  chain->SetBranchStatus("IanVtxXErr", 1);
  chain->SetBranchStatus("IanVtxYErr", 1);
  chain->SetBranchStatus("IanVtxZErr", 1);
  chain->SetBranchStatus("IanVtxNChi2", 1);
  chain->SetBranchStatus("dxyLeg1", 1);
  chain->SetBranchStatus("dxyLeg2", 1);
  chain->SetBranchStatus("dxyLeg1Err", 1);
  chain->SetBranchStatus("dxyLeg2Err", 1);
  chain->SetBranchStatus("dzLeg1", 1);
  chain->SetBranchStatus("dzLeg2", 1);
  chain->SetBranchStatus("dzLeg1Err", 1);
  chain->SetBranchStatus("dzLeg2Err", 1);
  chain->SetBranchStatus("dxyOPVLeg1", 1);
  chain->SetBranchStatus("dxyOPVLeg2", 1);
  chain->SetBranchStatus("dxyOPVLeg1Err", 1);
  chain->SetBranchStatus("dxyOPVLeg2Err", 1);
  chain->SetBranchStatus("dzOPVLeg1", 1);
  chain->SetBranchStatus("dzOPVLeg2", 1);
  chain->SetBranchStatus("dzOPVLeg1Err", 1);
  chain->SetBranchStatus("dzOPVLeg2Err", 1);
  chain->SetBranchStatus("dxyBSLeg1", 1);
  chain->SetBranchStatus("dxyBSLeg2", 1);
  chain->SetBranchStatus("dxyBSLeg1Err", 1);
  chain->SetBranchStatus("dxyBSLeg2Err", 1);
  chain->SetBranchStatus("dzBSLeg1", 1);
  chain->SetBranchStatus("dzBSLeg2", 1);
  chain->SetBranchStatus("dzBSLeg1Err", 1);
  chain->SetBranchStatus("dzBSLeg2Err", 1);
  chain->SetBranchStatus("dxyGenLeg1", 1);
  chain->SetBranchStatus("dxyGenLeg2", 1);
  chain->SetBranchStatus("dxyGenLeg1Err", 1);
  chain->SetBranchStatus("dxyGenLeg2Err", 1);
  chain->SetBranchStatus("dzGenLeg1", 1);
  chain->SetBranchStatus("dzGenLeg2", 1);
  chain->SetBranchStatus("dzGenLeg1Err", 1);
  chain->SetBranchStatus("dzGenLeg2Err", 1);

  chain->SetBranchAddress("run", &run_);
  chain->SetBranchAddress("event", &event_);
  chain->SetBranchAddress("lumi", &lumi_);
  chain->SetBranchAddress("diTauLegsP4", &diTauLegsP4_);
  chain->SetBranchAddress("diTauLegsLchP4", &diTauLegsLchP4_);
  chain->SetBranchAddress("diTauLegsPCA", &diTauLegsPCA_);
  chain->SetBranchAddress("diTauLegsPCAM2", &diTauLegsPCAM2_);
  chain->SetBranchAddress("diTauLegsPCAOPV", &diTauLegsPCAOPV_);
  chain->SetBranchAddress("diTauLegsPCABS", &diTauLegsPCABS_);
  chain->SetBranchAddress("diTauLegsPCAGen", &diTauLegsPCAGen_);
  chain->SetBranchAddress("diTauLegsPCAGenM2", &diTauLegsPCAGenM2_);
  chain->SetBranchAddress("diTauLegsPCAIan", &diTauLegsPCAIan_);
  chain->SetBranchAddress("diTauLegsPCAM2Cov", &diTauLegsPCAM2Cov_);
  chain->SetBranchAddress("diTauLegsPCAGenM2Cov", &diTauLegsPCAGenM2Cov_);
  chain->SetBranchAddress("diTauLegsLchP3AtPCA", &diTauLegsLchP3AtPCA_);
  chain->SetBranchAddress("diTauLegsLchP3AtPCAOPV", &diTauLegsLchP3AtPCAOPV_);
  chain->SetBranchAddress("diTauLegsLchP3AtPCABS", &diTauLegsLchP3AtPCABS_);
  chain->SetBranchAddress("diTauLegsLchP3AtPCAGen", &diTauLegsLchP3AtPCAGen_);
  chain->SetBranchAddress("diTauLegsIPAtPCA", &diTauLegsIPAtPCA_);
  chain->SetBranchAddress("diTauLegsIPAtPCAV2", &diTauLegsIPAtPCAV2_);
  chain->SetBranchAddress("diTauLegsIPAtPCAOPV", &diTauLegsIPAtPCAOPV_);
  chain->SetBranchAddress("diTauLegsIPAtPCAOPVV2", &diTauLegsIPAtPCAOPVV2_);
  chain->SetBranchAddress("diTauLegsIPAtPCAGen", &diTauLegsIPAtPCAGen_);
  chain->SetBranchAddress("diTauLegsIPAtPCAGenV2", &diTauLegsIPAtPCAGenV2_);
  chain->SetBranchAddress("diTauLegsIPAtPCATauVtx", &diTauLegsIPAtPCATauVtx_);
  chain->SetBranchAddress("genTausP4", &genTausP4_);
  chain->SetBranchAddress("genVP4", &genVP4_);
  chain->SetBranchAddress("genVPid", &genVPid_);
  chain->SetBranchAddress("genTausPid", &genTausPid_);
  chain->SetBranchAddress("genTausCharge", &genTausCharge_);
  chain->SetBranchAddress("genTausStatus", &genTausStatus_);
  chain->SetBranchAddress("genTauPSonsP4", &genTauPSonsP4_);
  chain->SetBranchAddress("genTauNSonsP4", &genTauNSonsP4_);
  chain->SetBranchAddress("genTauPSonsPid", &genTauPSonsPid_);
  chain->SetBranchAddress("genTauPSonsCharge", &genTauPSonsCharge_);
  chain->SetBranchAddress("genTauPSonsStatus", &genTauPSonsStatus_);
  chain->SetBranchAddress("genTauNSonsPid", &genTauNSonsPid_);
  chain->SetBranchAddress("genTauNSonsCharge", &genTauNSonsCharge_);
  chain->SetBranchAddress("genTauNSonsStatus", &genTauNSonsStatus_);
  chain->SetBranchAddress("VtxPos", &VtxPos_);
  chain->SetBranchAddress("BSPos", &BSPos_);
  chain->SetBranchAddress("ReFitVtxPos", &ReFitVtxPos_);
  chain->SetBranchAddress("IanVtxPos", &IanVtxPos_);
  chain->SetBranchAddress("HiggsGenVtx", &HiggsGenVtx_);
  chain->SetBranchAddress("TausGenVtx", &TausGenVtx_);
  chain->SetBranchAddress("TauPSonsGenVtx", &TauPSonsGenVtx_);
  chain->SetBranchAddress("TauNSonsGenVtx", &TauNSonsGenVtx_);
  chain->SetBranchAddress("tightestHPSDB3HWPLeg1", &tightestHPSDB3HWPLeg1_);
  chain->SetBranchAddress("tightestHPSDB3HWPLeg2", &tightestHPSDB3HWPLeg2_);
  chain->SetBranchAddress("hpsDB3HLeg1", &hpsDB3HLeg1_);
  chain->SetBranchAddress("hpsDB3HLeg2", &hpsDB3HLeg2_);
  chain->SetBranchAddress("diTauCharge", &diTauCharge_);
  chain->SetBranchAddress("chargeLeg1", &chargeLeg1_);
  chain->SetBranchAddress("decayModeLeg1", &decayModeLeg1_);
  chain->SetBranchAddress("decayModeLeg2", &decayModeLeg2_);
  chain->SetBranchAddress("genDecayModeLeg1", &genDecayModeLeg1_);
  chain->SetBranchAddress("genDecayModeLeg2", &genDecayModeLeg2_);
  chain->SetBranchAddress("genVtxX", &genVtxX_);
  chain->SetBranchAddress("genVtxY", &genVtxY_);
  chain->SetBranchAddress("genVtxZ", &genVtxZ_);
  chain->SetBranchAddress("VtxX", &VtxX_);
  chain->SetBranchAddress("VtxY", &VtxY_);
  chain->SetBranchAddress("VtxZ", &VtxZ_);
  chain->SetBranchAddress("VtxXErr", &VtxXErr_);
  chain->SetBranchAddress("VtxYErr", &VtxYErr_);
  chain->SetBranchAddress("VtxZErr", &VtxZErr_);
  chain->SetBranchAddress("VtxNChi2", &VtxNChi2_);
  chain->SetBranchAddress("ReFitVtxX", &ReFitVtxX_);
  chain->SetBranchAddress("ReFitVtxY", &ReFitVtxY_);
  chain->SetBranchAddress("ReFitVtxZ", &ReFitVtxZ_);
  chain->SetBranchAddress("ReFitVtxXErr", &ReFitVtxXErr_);
  chain->SetBranchAddress("ReFitVtxYErr", &ReFitVtxYErr_);
  chain->SetBranchAddress("ReFitVtxZErr", &ReFitVtxZErr_);
  chain->SetBranchAddress("ReFitVtxRho", &ReFitVtxRho_);
  chain->SetBranchAddress("ReFitVtxNdof", &ReFitVtxNdof_);
  chain->SetBranchAddress("ReFitVtxNChi2", &ReFitVtxNChi2_);
  chain->SetBranchAddress("IanVtxX", &IanVtxX_);
  chain->SetBranchAddress("IanVtxY", &IanVtxY_);
  chain->SetBranchAddress("IanVtxZ", &IanVtxZ_);
  chain->SetBranchAddress("IanVtxXErr", &IanVtxXErr_);
  chain->SetBranchAddress("IanVtxYErr", &IanVtxYErr_);
  chain->SetBranchAddress("IanVtxZErr", &IanVtxZErr_);
  chain->SetBranchAddress("IanVtxNChi2", &IanVtxNChi2_);
  chain->SetBranchAddress("dxyLeg1", &dxyLeg1_);
  chain->SetBranchAddress("dxyLeg2", &dxyLeg2_);
  chain->SetBranchAddress("dxyLeg1Err", &dxyLeg1Err_);
  chain->SetBranchAddress("dxyLeg2Err", &dxyLeg2Err_);
  chain->SetBranchAddress("dzLeg1", &dzLeg1_);
  chain->SetBranchAddress("dzLeg2", &dzLeg2_);
  chain->SetBranchAddress("dzLeg1Err", &dzLeg1Err_);
  chain->SetBranchAddress("dzLeg2Err", &dzLeg2Err_);
  chain->SetBranchAddress("dxyOPVLeg1", &dxyOPVLeg1_);
  chain->SetBranchAddress("dxyOPVLeg2", &dxyOPVLeg2_);
  chain->SetBranchAddress("dxyOPVLeg1Err", &dxyOPVLeg1Err_);
  chain->SetBranchAddress("dxyOPVLeg2Err", &dxyOPVLeg2Err_);
  chain->SetBranchAddress("dzOPVLeg1", &dzOPVLeg1_);
  chain->SetBranchAddress("dzOPVLeg2", &dzOPVLeg2_);
  chain->SetBranchAddress("dzOPVLeg1Err", &dzOPVLeg1Err_);
  chain->SetBranchAddress("dzOPVLeg2Err", &dzOPVLeg2Err_);
  chain->SetBranchAddress("dxyBSLeg1", &dxyBSLeg1_);
  chain->SetBranchAddress("dxyBSLeg2", &dxyBSLeg2_);
  chain->SetBranchAddress("dxyBSLeg1Err", &dxyBSLeg1Err_);
  chain->SetBranchAddress("dxyBSLeg2Err", &dxyBSLeg2Err_);
  chain->SetBranchAddress("dzBSLeg1", &dzBSLeg1_);
  chain->SetBranchAddress("dzBSLeg2", &dzBSLeg2_);
  chain->SetBranchAddress("dzBSLeg1Err", &dzBSLeg1Err_);
  chain->SetBranchAddress("dzBSLeg2Err", &dzBSLeg2Err_);
  chain->SetBranchAddress("dxyGenLeg1", &dxyGenLeg1_);
  chain->SetBranchAddress("dxyGenLeg2", &dxyGenLeg2_);
  chain->SetBranchAddress("dxyGenLeg1Err", &dxyGenLeg1Err_);
  chain->SetBranchAddress("dxyGenLeg2Err", &dxyGenLeg2Err_);
  chain->SetBranchAddress("dzGenLeg1", &dzGenLeg1_);
  chain->SetBranchAddress("dzGenLeg2", &dzGenLeg2_);
  chain->SetBranchAddress("dzGenLeg1Err", &dzGenLeg1Err_);
  chain->SetBranchAddress("dzGenLeg2Err", &dzGenLeg2Err_);
  chain->SetBranchAddress("trackPtLeg1", &trackPtLeg1_);
  chain->SetBranchAddress("trackPtErrLeg1", &trackPtErrLeg1_);
  chain->SetBranchAddress("nMisingHitsLeg1", &nMisingHitsLeg1_);
  chain->SetBranchAddress("nHitsLeg1", &nHitsLeg1_);
  chain->SetBranchAddress("nTkHitsLeg1", &nTkHitsLeg1_);
  chain->SetBranchAddress("nPxlHitsLeg1", &nPxlHitsLeg1_);
  chain->SetBranchAddress("hasFirstPxlHitLeg1", &hasFirstPxlHitLeg1_);
  chain->SetBranchAddress("trackPtLeg2", &trackPtLeg2_);
  chain->SetBranchAddress("trackPtErrLeg2", &trackPtErrLeg2_);
  chain->SetBranchAddress("nMisingHitsLeg2", &nMisingHitsLeg2_);
  chain->SetBranchAddress("nHitsLeg2", &nHitsLeg2_);
  chain->SetBranchAddress("nTkHitsLeg2", &nTkHitsLeg2_);
  chain->SetBranchAddress("nPxlHitsLeg2", &nPxlHitsLeg2_);
  chain->SetBranchAddress("hasFirstPxlHitLeg2", &hasFirstPxlHitLeg2_);
  
  ///////////////////////
  // LOOP OVER ENTRIES //
  ///////////////////////

  //Sync check with Christian using VBF events without tau - polarization
  vector<int>EventsCheck; EventsCheck.clear();
  EventsCheck.push_back(77941); EventsCheck.push_back(154153);
  EventsCheck.push_back(11531); EventsCheck.push_back(109713);
  EventsCheck.push_back(186146); EventsCheck.push_back(186204);
  EventsCheck.push_back(126861); EventsCheck.push_back(73063);
  EventsCheck.push_back(81845); EventsCheck.push_back(154618);
  EventsCheck.push_back(9244);
 
  int nEntries    = chain->GetEntries() ;
  unsigned int lastRun_, lastLumi_, lastEvent_=0;
  bool foundDiTauPair_ = false;
  for (int n = 0; n < nEntries; n++) {
    
    chain->GetEntry(n);
    if(n%1000==0) std::cout << n <<"/"<<nEntries<< std::endl;

    //apply kinematic cuts
    double ptL1     = (*diTauLegsP4_)[0].Pt();
    double ptL2     = (*diTauLegsP4_)[1].Pt();
    double etaL1    = (*diTauLegsP4_)[0].Eta();
    double etaL2    = (*diTauLegsP4_)[1].Eta();
    if(ptL1 < 20 || ptL2 < 20 || TMath::Abs(etaL1) > 2.1 || TMath::Abs(etaL2) > 2.1)
      continue;
    
    //if(std::find(EventsCheck.begin(), EventsCheck.end(), event_)==EventsCheck.end()) continue;
    if(verbosity)std::cout<<std::endl;
    if(verbosity)std::cout<<"run : "<<run_<<", lumi : "<<lumi_<<", event : "<<event_<<std::endl;  
    //compute decay plane in generator level
    if(event_ != lastEvent_){ //take one entry per event
	if(genVPid_ == HPid_){
	  //rest frame of Higgs
	  ROOT::Math::Boost boost_to_rf_h((*genVP4_)[0].BoostToCM());
	  
	  LV gentauP_lab(0, 0, 0, 0); LV gentauN_lab(0, 0, 0, 0);
	  Point3D gentauPVtx_(0, 0, 0); Point3D gentauNVtx_(0, 0, 0);
	  if(genTausP4_->size() == 2 && genTausCharge_->size() == 2){
	    if((*genTausCharge_)[0] > 0 && (*genTausCharge_)[1] < 0){
	      gentauP_lab = (*genTausP4_)[0];
	      gentauN_lab = (*genTausP4_)[1];
	      gentauPVtx_ = (*TausGenVtx_)[0];
	      gentauNVtx_ = (*TausGenVtx_)[1];
	    }
	    else if((*genTausCharge_)[0] < 0 && (*genTausCharge_)[1] > 0){
	      gentauP_lab = (*genTausP4_)[1];
              gentauN_lab = (*genTausP4_)[0];
	      gentauPVtx_ = (*TausGenVtx_)[1];
              gentauNVtx_ = (*TausGenVtx_)[0];
            } 
	  }
	  
	  LV gentauP_rf = gentauP_lab.pt()>0 ? boost_to_rf_h(gentauP_lab) : LV(0, 0, 0, 0);
	  LV gentauN_rf = gentauN_lab.pt()>0 ? boost_to_rf_h(gentauN_lab) : LV(0, 0, 0, 0);

	  LV genpionP_lab(0, 0, 0, 0); LV genpionN_lab(0, 0, 0, 0);
	  Point3D genpionPVtx_(0, 0, 0); Point3D genpionNVtx_(0, 0, 0);
	  if(genTauPSonsP4_->size()>0 && genTauPSonsP4_->size()<3 && genTauPSonsPid_->size() > 0){ //take only one prong
	    for(size_t i = 0; i<genTauPSonsP4_->size(); i++){
	      if(abs((*genTauPSonsPid_)[i]) == 211){
		genpionP_lab = (*genTauPSonsP4_)[i];
		genpionPVtx_ = (*TauPSonsGenVtx_)[i];
	      }
	    }
	  }
	  if(genTauNSonsP4_->size()>0 && genTauNSonsP4_->size()<3 && genTauNSonsPid_->size() > 0){ //take only one prong
            for(size_t i = 0; i<genTauNSonsP4_->size(); i++){
              if(abs((*genTauNSonsPid_)[i]) == 211){
		genpionN_lab = (*genTauNSonsP4_)[i];
		genpionNVtx_ = (*TauNSonsGenVtx_)[i];
	      }
            }
          }

	  ROOT::Math::Boost boost_to_rf_gentauP(gentauP_lab.BoostToCM());
	  ROOT::Math::Boost boost_to_rf_gentauN(gentauN_lab.BoostToCM());

	  //LV genpionP_rf = genpionP_lab.pt()>0 ? boost_to_rf_h(genpionP_lab) : LV(0, 0, 0, 0);
          //LV genpionN_rf = genpionN_lab.pt()>0 ? boost_to_rf_h(genpionN_lab) : LV(0, 0, 0, 0);
	  LV genpionP_rf = genpionP_lab.pt()>0 ? boost_to_rf_gentauP(genpionP_lab) : LV(0, 0, 0, 0);
          LV genpionN_rf = genpionN_lab.pt()>0 ? boost_to_rf_gentauN(genpionN_lab) : LV(0, 0, 0, 0);

	  Point3D GenHVtx_ =  (*HiggsGenVtx_)[0];
	  if(verbosity)std::cout<<"Higgs Vertex ("<<GenHVtx_.x()<<","<<GenHVtx_.y()<<","<<GenHVtx_.z()<<")"<<std::endl;
	  if(verbosity){
	    std::cout<<"In Lab frame :"<<std::endl;
	    std::cout<<" Gen Tau+ Pt = "<<gentauP_lab.pt()<<" eta = "<<gentauP_lab.eta()<<" phi = "<<gentauP_lab.phi()<<" mass = "<<gentauP_lab.mass()<<std::endl;
	    std::cout<<" Gen Tau- Pt = "<<gentauN_lab.pt()<<" eta = "<<gentauN_lab.eta()<<" phi = "<<gentauN_lab.phi()<<" mass = "<<gentauN_lab.mass()<<std::endl;
	    std::cout<<" Gen Pion+ Pt = "<<genpionP_lab.pt()<<" eta = "<<genpionP_lab.eta()<<" phi = "<<genpionP_lab.phi()<<" mass = "<<genpionP_lab.mass()<<std::endl;
	    std::cout<<" Gen Pion- Pt = "<<genpionN_lab.pt()<<" eta = "<<genpionN_lab.eta()<<" phi = "<<genpionN_lab.phi()<<" mass = "<<genpionN_lab.mass()<<std::endl;
	  }
	  if(verbosity){
	    std::cout<<"In Higgs rest frame :"<<std::endl;
	    std::cout<<" Gen Tau+ Pt = "<<gentauP_rf.pt()<<" eta = "<<gentauP_rf.eta()<<" phi = "<<gentauP_rf.phi()<<" mass = "<<gentauP_rf.mass()<<std::endl;
	    std::cout<<" Gen Tau- Pt = "<<gentauN_rf.pt()<<" eta = "<<gentauN_rf.eta()<<" phi = "<<gentauN_rf.phi()<<" mass = "<<gentauN_rf.mass()<<std::endl;
	    std::cout<<" Gen Pion+ Pt = "<<genpionP_rf.pt()<<" eta = "<<genpionP_rf.eta()<<" phi = "<<genpionP_rf.phi()<<" mass = "<<genpionP_rf.mass()<<std::endl;
	    std::cout<<" Gen Pion- Pt = "<<genpionN_rf.pt()<<" eta = "<<genpionN_rf.eta()<<" phi = "<<genpionN_rf.phi()<<" mass = "<<genpionN_rf.mass()<<std::endl;
          }
	  if(gentauP_rf.pt()>0 && gentauN_rf.pt()>0 && genpionP_rf.pt()>0 && genpionN_rf.pt()>0){
	    //Get only the position component of the vector
	    PV gentauP_rf_3d = gentauP_rf.Vect();
	    PV gentauN_rf_3d = gentauN_rf.Vect();
	    PV genpionP_rf_3d = genpionP_rf.Vect();
            PV genpionN_rf_3d = genpionN_rf.Vect();

	    //Get the normal to the decay plane
	    PV dplane_normal_1 = gentauP_rf_3d.Cross(genpionP_rf_3d);
	    PV dplane_normal_2 = genpionN_rf_3d.Cross(gentauN_rf_3d); //gentauN_rf_3d.Cross(genpionN_rf_3d);
	    PV dplane_normal_1_u = dplane_normal_1/sqrt(dplane_normal_1.mag2());
	    PV dplane_normal_2_u = dplane_normal_2/sqrt(dplane_normal_2.mag2());

	    //Get CP angle in ZMF
	    double phi_cp_gen = TMath::ACos(dplane_normal_1_u.Dot(dplane_normal_2_u));
	    if(verbosity)std::cout<<"phi* in Higgs rest-frame: method1 = "<<phi_cp_gen<<std::endl;
	    
	    histos_["CPPhi_gen"]->Fill(phi_cp_gen);
	    histos_["CPPhi_gen_v2"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(dplane_normal_1_u, dplane_normal_2_u));
	  }

	  //compute phi_star in pi-pi+ rf.
	  //Point3D GenHVtx_ =  (*HiggsGenVtx_)[0];
	  if(genpionP_lab.pt() > 0 && genpionN_lab.pt() > 0){
	    LV gen_pion_pair_lab = genpionP_lab + genpionN_lab;
	    ROOT::Math::Boost boost_to_rf_gen(gen_pion_pair_lab.BoostToCM());

	    //if(verbosity)std::cout<<"Higgs Vertex ("<<GenHVtx_.x()<<","<<GenHVtx_.y()<<","<<GenHVtx_.z()<<")"<<std::endl;
	    //if(verbosity)std::cout<<"lab frame : pionP pt = "<<genpionP_lab.pt()<<" pionN pt "<<genpionN_lab.pt()<<std::endl;

	    //boost to rf of two pion system
	    LV genpionP_rf_2p = boost_to_rf_gen(genpionP_lab);
	    LV genpionN_rf_2p = boost_to_rf_gen(genpionN_lab);
	    if(verbosity)std::cout<<"p1+p2 rest frame : pionP pt = "<<genpionP_rf_2p.pt()<<" pionN pt "<<genpionN_rf_2p.pt()<<std::endl;

	    //get unit momentum vectors
	    PV genpionP_lab_3d = genpionP_lab.Vect();
	    PV genpionN_lab_3d = genpionN_lab.Vect();
	    PV genpionP_lab_3d_u = genpionP_lab_3d/sqrt(genpionP_lab_3d.mag2());
	    PV genpionN_lab_3d_u = genpionN_lab_3d/sqrt(genpionN_lab_3d.mag2());

	    if(verbosity)std::cout<<"genTau+ Vertex ("<<gentauPVtx_.x()<<","<<gentauPVtx_.y()<<","<<gentauPVtx_.z()<<")"<<std::endl;
            if(verbosity)std::cout<<"genTau- Vertex ("<<gentauNVtx_.x()<<","<<gentauNVtx_.y()<<","<<gentauNVtx_.z()<<")"<<std::endl;

	    if(verbosity)std::cout<<"genPion+ Vertex ("<<genpionPVtx_.x()<<","<<genpionPVtx_.y()<<","<<genpionPVtx_.z()<<")"<<std::endl;
	    if(verbosity)std::cout<<"genPion- Vertex ("<<genpionNVtx_.x()<<","<<genpionNVtx_.y()<<","<<genpionNVtx_.z()<<")"<<std::endl;
	    
	    //Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
	    double tP = genpionP_lab_3d_u.x()*(GenHVtx_.x() - genpionPVtx_.x()) + genpionP_lab_3d_u.y()*(GenHVtx_.y() - genpionPVtx_.y()) + genpionP_lab_3d_u.z()*(GenHVtx_.z() - genpionPVtx_.z());
	    double tN = genpionN_lab_3d_u.x()*(GenHVtx_.x() - genpionNVtx_.x()) + genpionN_lab_3d_u.y()*(GenHVtx_.y() - genpionNVtx_.y()) + genpionN_lab_3d_u.z()*(GenHVtx_.z() - genpionNVtx_.z());
	      
	    Point3D pcaP(genpionPVtx_.x() + genpionP_lab_3d_u.x()*tP, genpionPVtx_.y() + genpionP_lab_3d_u.y()*tP, genpionPVtx_.z() + genpionP_lab_3d_u.z()*tP);
	    Point3D pcaN(genpionNVtx_.x() + genpionN_lab_3d_u.x()*tN, genpionNVtx_.y() + genpionN_lab_3d_u.y()*tN, genpionNVtx_.z() + genpionN_lab_3d_u.z()*tN);

	    //if(verbosity)std::cout<<"PCA Pion+ ("<<pcaP.x()<<","<<pcaP.y()<<","<<pcaP.z()<<")"<<std::endl;
	    //if(verbosity)std::cout<<"PCA Pion- ("<<pcaN.x()<<","<<pcaN.y()<<","<<pcaN.z()<<")"<<std::endl;

	    //Get the normalized IP vectors
	    PV ipvP_lab_gen_3d(pcaP.x() - GenHVtx_.x(), pcaP.y() - GenHVtx_.y(), pcaP.z() - GenHVtx_.z());
	    PV ipvN_lab_gen_3d(pcaN.x() - GenHVtx_.x(), pcaN.y() - GenHVtx_.y(), pcaN.z() - GenHVtx_.z());

	    if(verbosity)std::cout<<"IPV lab Pion+ ("<<ipvP_lab_gen_3d.x()<<","<<ipvP_lab_gen_3d.y()<<","<<ipvP_lab_gen_3d.z()<<")"<<std::endl;
            if(verbosity)std::cout<<"IPV lab Pion- ("<<ipvN_lab_gen_3d.x()<<","<<ipvN_lab_gen_3d.y()<<","<<ipvN_lab_gen_3d.z()<<")"<<std::endl;

	    LV ipvP_lab_gen(ipvP_lab_gen_3d.x()/sqrt(ipvP_lab_gen_3d.mag2()), ipvP_lab_gen_3d.y()/sqrt(ipvP_lab_gen_3d.mag2()), ipvP_lab_gen_3d.z()/sqrt(ipvP_lab_gen_3d.mag2()), 0);
	    LV ipvN_lab_gen(ipvN_lab_gen_3d.x()/sqrt(ipvN_lab_gen_3d.mag2()), ipvN_lab_gen_3d.y()/sqrt(ipvN_lab_gen_3d.mag2()), ipvN_lab_gen_3d.z()/sqrt(ipvN_lab_gen_3d.mag2()), 0);
	    
	    //boost to ZMF
	    LV ipvP_rf_gen = boost_to_rf_gen(ipvP_lab_gen);
	    LV ipvN_rf_gen = boost_to_rf_gen(ipvN_lab_gen);

	    //Get only the position component of the vector
	    PV ipvP_rf_gen_3d = ipvP_rf_gen.Vect();
	    PV ipvN_rf_gen_3d = ipvN_rf_gen.Vect();

	    //if(verbosity)std::cout<<"IPV Pion+ ("<<ipvP_rf_gen_3d.x()<<","<<ipvP_rf_gen_3d.y()<<","<<ipvP_rf_gen_3d.z()<<")"<<std::endl;
	    //if(verbosity)std::cout<<"IPV Pion- ("<<ipvN_rf_gen_3d.x()<<","<<ipvN_rf_gen_3d.y()<<","<<ipvN_rf_gen_3d.z()<<")"<<std::endl;

	    PV genpionP_rf_2p_3d = genpionP_rf_2p.Vect();
	    PV genpionN_rf_2p_3d = genpionN_rf_2p.Vect();

	    //Get the unit (normalized) vector along pion momentum
	    PV genpionP_rf_2p_3d_u = genpionP_rf_2p_3d/TMath::Sqrt(genpionP_rf_2p_3d.mag2());
	    PV genpionN_rf_2p_3d_u = genpionN_rf_2p_3d/TMath::Sqrt(genpionN_rf_2p_3d.mag2());

	    //Get the longitudinal component of IP vector parallel to pion momenta
	    PV ipvP_rf_gen_3d_l = ipvP_rf_gen_3d.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
	    PV ipvN_rf_gen_3d_l = ipvN_rf_gen_3d.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;

	    //if(verbosity)std::cout<<"LIPV Pion+ ("<<ipvP_rf_gen_3d_l.x()<<","<<ipvP_rf_gen_3d_l.y()<<","<<ipvP_rf_gen_3d_l.z()<<")"<<std::endl;
	    //if(verbosity)std::cout<<"LIPV Pion- ("<<ipvN_rf_gen_3d_l.x()<<","<<ipvN_rf_gen_3d_l.y()<<","<<ipvN_rf_gen_3d_l.z()<<")"<<std::endl;

	    //Get IP vector normal to pion momenta
	    PV ipvP_rf_gen_3d_t = ipvP_rf_gen_3d - ipvP_rf_gen_3d_l;
	    PV ipvN_rf_gen_3d_t = ipvN_rf_gen_3d - ipvN_rf_gen_3d_l;

	    //if(verbosity)std::cout<<"TIPV Pion+ ("<<ipvP_rf_gen_3d_t.x()<<","<<ipvP_rf_gen_3d_t.y()<<","<<ipvP_rf_gen_3d_t.z()<<")"<<std::endl;
	    //if(verbosity)std::cout<<"TIPV Pion- ("<<ipvN_rf_gen_3d_t.x()<<","<<ipvN_rf_gen_3d_t.y()<<","<<ipvN_rf_gen_3d_t.z()<<")"<<std::endl;

	    //Get normalized normal IP vector
	    PV ipvP_rf_gen_3d_t_u = ipvP_rf_gen_3d_t/TMath::Sqrt(ipvP_rf_gen_3d_t.mag2());
	    PV ipvN_rf_gen_3d_t_u = ipvN_rf_gen_3d_t/TMath::Sqrt(ipvN_rf_gen_3d_t.mag2());

	    //Get CP angle in ZMF
	    double phi_star_gen = TMath::ACos(ipvP_rf_gen_3d_t_u.Dot(ipvN_rf_gen_3d_t_u));
	    histos_["CPPhiStar_gen"]->Fill(phi_star_gen);
	    //Get CP angle in lab frame
	    double phi_lab_gen = TMath::ACos(ipvP_lab_gen_3d.Dot(ipvN_lab_gen_3d)/(sqrt(ipvP_lab_gen_3d.mag2())*sqrt(ipvN_lab_gen_3d.mag2())));
	    histos_["CPPhiLab_gen"]->Fill(phi_lab_gen);
	    double OstarCP_gen = genpionN_rf_2p_3d_u.Dot(ipvP_rf_gen_3d_t_u.Cross(ipvN_rf_gen_3d_t_u));
	    double phi_star_gen_cp = (OstarCP_gen >= 0) ? phi_star_gen : (TMath::TwoPi() - phi_star_gen);
	    histos_["CPPhiStar_gen_CP"]->Fill(phi_star_gen_cp);

	    //smear Vertex, 10 micron in X/Y, 30 micron in Z.
	    double spv_x = GenHVtx_.x() + randx.Gaus(0., 0.001); // randx.Uniform(-0.001, 0.001);
	    double spv_y = GenHVtx_.y() + randy.Gaus(0., 0.001); // randy.Uniform(-0.001, 0.001);
	    double spv_z = GenHVtx_.z() + randz.Gaus(0., 0.003); // randz.Uniform(-0.003, 0.003);
	    //cout<<" randx.Gaus(0., 0.001) "<<randx.Gaus(0., 0.001)<<" randy.Gaus(0., 0.001) "<<randy.Gaus(0., 0.001)<<" randz.Gaus(0., 0.003) "<<randz.Gaus(0., 0.003)<<endl;
	    
	    double tP_s = genpionP_lab_3d_u.x()*(spv_x - genpionPVtx_.x()) + genpionP_lab_3d_u.y()*(spv_y - genpionPVtx_.y()) + genpionP_lab_3d_u.z()*(spv_z - genpionPVtx_.z());
            double tN_s = genpionN_lab_3d_u.x()*(spv_x - genpionNVtx_.x()) + genpionN_lab_3d_u.y()*(spv_y - genpionNVtx_.y()) + genpionN_lab_3d_u.z()*(spv_z - genpionNVtx_.z());

            Point3D pcaP_s(genpionPVtx_.x() + genpionP_lab_3d_u.x()*tP_s, genpionPVtx_.y() + genpionP_lab_3d_u.y()*tP_s, genpionPVtx_.z() + genpionP_lab_3d_u.z()*tP_s);
            Point3D pcaN_s(genpionNVtx_.x() + genpionN_lab_3d_u.x()*tN_s, genpionNVtx_.y() + genpionN_lab_3d_u.y()*tN_s, genpionNVtx_.z() + genpionN_lab_3d_u.z()*tN_s);

	    PV ipvP_lab_gen_3d_s(pcaP_s.x() - spv_x, pcaP_s.y() - spv_y, pcaP_s.z() - spv_z);
            PV ipvN_lab_gen_3d_s(pcaN_s.x() - spv_x, pcaN_s.y() - spv_y, pcaN_s.z() - spv_z);

	    LV ipvP_lab_gen_s(ipvP_lab_gen_3d_s.Unit().x(), ipvP_lab_gen_3d_s.Unit().y(), ipvP_lab_gen_3d_s.Unit().z(), 0);
            LV ipvN_lab_gen_s(ipvN_lab_gen_3d_s.Unit().x(), ipvN_lab_gen_3d_s.Unit().y(), ipvN_lab_gen_3d_s.Unit().z(), 0);

            //boost to ZMF
            LV ipvP_rf_gen_s = boost_to_rf_gen(ipvP_lab_gen_s);
            LV ipvN_rf_gen_s = boost_to_rf_gen(ipvN_lab_gen_s);

            //Get only the position component of the vector
            PV ipvP_rf_gen_3d_s = ipvP_rf_gen_s.Vect();
            PV ipvN_rf_gen_3d_s = ipvN_rf_gen_s.Vect();

            //Get IP vector normal to pion momenta
            PV ipvP_rf_gen_3d_t_s = ipvP_rf_gen_3d_s - ipvP_rf_gen_3d_s.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
            PV ipvN_rf_gen_3d_t_s = ipvN_rf_gen_3d_s - ipvN_rf_gen_3d_s.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;

	    //Get normalized normal IP vector
            PV ipvP_rf_gen_3d_t_u_s = ipvP_rf_gen_3d_t_s.Unit();
            PV ipvN_rf_gen_3d_t_u_s = ipvN_rf_gen_3d_t_s.Unit();

            //Get CP angle in ZMF
            histos_["CPPhiStar_gen_vs"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_s.Dot(ipvN_rf_gen_3d_t_u_s)));
	    if(sqrt(ipvP_lab_gen_3d_s.mag2()) > 0.003 && sqrt(ipvN_lab_gen_3d_s.mag2()) > 0.003)
	      histos_["CPPhiStar_gen_vs_ip30"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_s.Dot(ipvN_rf_gen_3d_t_u_s)));
	    if(sqrt(ipvP_lab_gen_3d_s.mag2()) > 0.005 && sqrt(ipvN_lab_gen_3d_s.mag2()) > 0.005)
              histos_["CPPhiStar_gen_vs_ip50"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_s.Dot(ipvN_rf_gen_3d_t_u_s)));

	    //Use only transverse smearing
	    double tP_st = genpionP_lab_3d_u.x()*(spv_x - genpionPVtx_.x()) + genpionP_lab_3d_u.y()*(spv_y - genpionPVtx_.y()) + genpionP_lab_3d_u.z()*(GenHVtx_.z() - genpionPVtx_.z());
            double tN_st = genpionN_lab_3d_u.x()*(spv_x - genpionNVtx_.x()) + genpionN_lab_3d_u.y()*(spv_y - genpionNVtx_.y()) + genpionN_lab_3d_u.z()*(GenHVtx_.z() - genpionNVtx_.z());

            Point3D pcaP_st(genpionPVtx_.x() + genpionP_lab_3d_u.x()*tP_st, genpionPVtx_.y() + genpionP_lab_3d_u.y()*tP_st, genpionPVtx_.z() + genpionP_lab_3d_u.z()*tP_st);
            Point3D pcaN_st(genpionNVtx_.x() + genpionN_lab_3d_u.x()*tN_st, genpionNVtx_.y() + genpionN_lab_3d_u.y()*tN_st, genpionNVtx_.z() + genpionN_lab_3d_u.z()*tN_st);

            PV ipvP_lab_gen_3d_st(pcaP_st.x() - spv_x, pcaP_st.y() - spv_y, pcaP_st.z() - GenHVtx_.z());
            PV ipvN_lab_gen_3d_st(pcaN_st.x() - spv_x, pcaN_st.y() - spv_y, pcaN_st.z() - GenHVtx_.z());

            LV ipvP_lab_gen_st(ipvP_lab_gen_3d_st.Unit().x(), ipvP_lab_gen_3d_st.Unit().y(), ipvP_lab_gen_3d_st.Unit().z(), 0);
            LV ipvN_lab_gen_st(ipvN_lab_gen_3d_st.Unit().x(), ipvN_lab_gen_3d_st.Unit().y(), ipvN_lab_gen_3d_st.Unit().z(), 0);

            //boost to ZMF
            LV ipvP_rf_gen_st = boost_to_rf_gen(ipvP_lab_gen_st);
            LV ipvN_rf_gen_st = boost_to_rf_gen(ipvN_lab_gen_st);

            //Get only the position component of the vector
            PV ipvP_rf_gen_3d_st = ipvP_rf_gen_st.Vect();
            PV ipvN_rf_gen_3d_st = ipvN_rf_gen_st.Vect();

	    //Get IP vector normal to pion momenta
            PV ipvP_rf_gen_3d_t_st = ipvP_rf_gen_3d_st - ipvP_rf_gen_3d_st.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
            PV ipvN_rf_gen_3d_t_st = ipvN_rf_gen_3d_st - ipvN_rf_gen_3d_st.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;

            //Get normalized normal IP vector
            PV ipvP_rf_gen_3d_t_u_st = ipvP_rf_gen_3d_t_st.Unit();
            PV ipvN_rf_gen_3d_t_u_st = ipvN_rf_gen_3d_t_st.Unit();

            //Get CP angle in ZMF
            histos_["CPPhiStar_gen_vst"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_st.Dot(ipvN_rf_gen_3d_t_u_st)));
            if(sqrt(ipvP_lab_gen_3d_st.mag2()) > 0.003 && sqrt(ipvN_lab_gen_3d_st.mag2()) > 0.003)
              histos_["CPPhiStar_gen_vst_ip30"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_st.Dot(ipvN_rf_gen_3d_t_u_st)));
            if(sqrt(ipvP_lab_gen_3d_st.mag2()) > 0.005 && sqrt(ipvN_lab_gen_3d_st.mag2()) > 0.005)
              histos_["CPPhiStar_gen_vst_ip50"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_st.Dot(ipvN_rf_gen_3d_t_u_st)));

	    //smear Tr. IP by 20 micron
	    double spcaP_x = (pcaP_s.x()+randipx.Gaus(0., 0.002));
	    double spcaP_y = (pcaP_s.y()+randipy.Gaus(0., 0.002));
	    double spcaN_x = (pcaN_s.x()+randinx.Gaus(0., 0.002));
	    double spcaN_y = (pcaN_s.y()+randiny.Gaus(0., 0.002));

	    PV ipvP_lab_gen_3d_vips(spcaP_x - spv_x, spcaP_y - spv_y, pcaP_s.z() - spv_z);
            PV ipvN_lab_gen_3d_vips(spcaN_x - spv_x, spcaN_y - spv_y, pcaN_s.z() - spv_z);
	    
	    LV ipvP_lab_gen_vips(ipvP_lab_gen_3d_vips.Unit().x(), ipvP_lab_gen_3d_vips.Unit().y(), ipvP_lab_gen_3d_vips.Unit().z(), 0);
            LV ipvN_lab_gen_vips(ipvN_lab_gen_3d_vips.Unit().x(), ipvN_lab_gen_3d_vips.Unit().y(), ipvN_lab_gen_3d_vips.Unit().z(), 0);

            //boost to ZMF
            LV ipvP_rf_gen_vips = boost_to_rf_gen(ipvP_lab_gen_vips);
            LV ipvN_rf_gen_vips = boost_to_rf_gen(ipvN_lab_gen_vips);

            //Get only the position component of the vector
            PV ipvP_rf_gen_3d_vips = ipvP_rf_gen_vips.Vect();
            PV ipvN_rf_gen_3d_vips = ipvN_rf_gen_vips.Vect();

            //Get IP vector normal to pion momenta
            PV ipvP_rf_gen_3d_t_vips = ipvP_rf_gen_3d_vips - ipvP_rf_gen_3d_vips.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
            PV ipvN_rf_gen_3d_t_vips = ipvN_rf_gen_3d_vips - ipvN_rf_gen_3d_vips.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;

            //Get normalized normal IP vector
            PV ipvP_rf_gen_3d_t_u_vips = ipvP_rf_gen_3d_t_vips.Unit();
            PV ipvN_rf_gen_3d_t_u_vips = ipvN_rf_gen_3d_t_vips.Unit();

            //Get CP angle in ZMF
            histos_["CPPhiStar_gen_vips"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vips.Dot(ipvN_rf_gen_3d_t_u_vips)));
	    if(sqrt(ipvP_lab_gen_3d_vips.mag2()) > 0.003 && sqrt(ipvN_lab_gen_3d_vips.mag2()) > 0.003)
	      histos_["CPPhiStar_gen_vips_ip30"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vips.Dot(ipvN_rf_gen_3d_t_u_vips)));
	    if(sqrt(ipvP_lab_gen_3d_vips.mag2()) > 0.005 && sqrt(ipvN_lab_gen_3d_vips.mag2()) > 0.005)
              histos_["CPPhiStar_gen_vips_ip50"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vips.Dot(ipvN_rf_gen_3d_t_u_vips)));

	    //smear also IPz
	    double spcaP_z = (pcaP_s.z()+randipz.Gaus(0., 0.002));
	    double spcaN_z = (pcaN_s.z()+randinz.Gaus(0., 0.002));

	    PV ipvP_lab_gen_3d_vipsv2(spcaP_x - spv_x, spcaP_y - spv_y, spcaP_z - spv_z);
            PV ipvN_lab_gen_3d_vipsv2(spcaN_x - spv_x, spcaN_y - spv_y, spcaN_z - spv_z);

            LV ipvP_lab_gen_vipsv2(ipvP_lab_gen_3d_vipsv2.Unit().x(), ipvP_lab_gen_3d_vipsv2.Unit().y(), ipvP_lab_gen_3d_vipsv2.Unit().z(), 0);
            LV ipvN_lab_gen_vipsv2(ipvN_lab_gen_3d_vipsv2.Unit().x(), ipvN_lab_gen_3d_vipsv2.Unit().y(), ipvN_lab_gen_3d_vipsv2.Unit().z(), 0);

            //boost to ZMF
            LV ipvP_rf_gen_vipsv2 = boost_to_rf_gen(ipvP_lab_gen_vipsv2);
            LV ipvN_rf_gen_vipsv2 = boost_to_rf_gen(ipvN_lab_gen_vipsv2);

            //Get only the position component of the vector
            PV ipvP_rf_gen_3d_vipsv2 = ipvP_rf_gen_vipsv2.Vect();
            PV ipvN_rf_gen_3d_vipsv2 = ipvN_rf_gen_vipsv2.Vect();

            //Get IP vector normal to pion momenta
            PV ipvP_rf_gen_3d_t_vipsv2 = ipvP_rf_gen_3d_vipsv2 - ipvP_rf_gen_3d_vipsv2.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
            PV ipvN_rf_gen_3d_t_vipsv2 = ipvN_rf_gen_3d_vipsv2 - ipvN_rf_gen_3d_vipsv2.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;

            //Get normalized normal IP vector
            PV ipvP_rf_gen_3d_t_u_vipsv2 = ipvP_rf_gen_3d_t_vipsv2.Unit();
            PV ipvN_rf_gen_3d_t_u_vipsv2 = ipvN_rf_gen_3d_t_vipsv2.Unit();

            //Get CP angle in ZMF
            histos_["CPPhiStar_gen_vipsv2"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
            if(sqrt(ipvP_lab_gen_3d_vipsv2.mag2()) > 0.003 && sqrt(ipvN_lab_gen_3d_vipsv2.mag2()) > 0.003)
              histos_["CPPhiStar_gen_vipsv2_ip30"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(sqrt(ipvP_lab_gen_3d_vipsv2.mag2()) > 0.005 && sqrt(ipvN_lab_gen_3d_vipsv2.mag2()) > 0.005)
              histos_["CPPhiStar_gen_vipsv2_ip50"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(genpionP_lab.pt() > 20 && genpionN_lab.pt() > 20)
	      histos_["CPPhiStar_gen_vipsv2_tau20"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(genpionP_lab.pt() > 30 && genpionN_lab.pt() > 30)
	      histos_["CPPhiStar_gen_vipsv2_tau30"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(genpionP_lab.pt() > 40 && genpionN_lab.pt() > 40)
	      histos_["CPPhiStar_gen_vipsv2_tau40"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(genpionP_lab.pt() > 50 && genpionN_lab.pt() > 50)
	      histos_["CPPhiStar_gen_vipsv2_tau50"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(genpionP_lab.pt() > 40 && genpionN_lab.pt() > 40 && sqrt(ipvP_lab_gen_3d_vipsv2.mag2()) > 0.003 && sqrt(ipvN_lab_gen_3d_vipsv2.mag2()) > 0.003)
	      histos_["CPPhiStar_gen_vipsv2_ip30_tau40"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));
	    if(genpionP_lab.pt() > 40 && genpionN_lab.pt() > 40 && sqrt(ipvP_lab_gen_3d_vipsv2.mag2()) > 0.005 && sqrt(ipvN_lab_gen_3d_vipsv2.mag2()) > 0.005)
	      histos_["CPPhiStar_gen_vipsv2_ip50_tau40"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_vipsv2.Dot(ipvN_rf_gen_3d_t_u_vipsv2)));

	    //Smear only PCA, not vertex
	    double s_pcaP_x = (pcaP.x()+randpcapx.Gaus(0., 0.002));
            double s_pcaP_y = (pcaP.y()+randpcapy.Gaus(0., 0.002));
            double s_pcaN_x = (pcaN.x()+randpcanx.Gaus(0., 0.002));
            double s_pcaN_y = (pcaN.y()+randpcany.Gaus(0., 0.002));
	    double s_pcaP_z = (pcaP.z()+randpcapz.Gaus(0., 0.002));
            double s_pcaN_z = (pcaN.z()+randpcanz.Gaus(0., 0.002));

            PV ipvP_lab_gen_3d_ips(s_pcaP_x - GenHVtx_.x(), s_pcaP_y - GenHVtx_.y(), s_pcaP_z - GenHVtx_.z());
            PV ipvN_lab_gen_3d_ips(s_pcaN_x - GenHVtx_.x(), s_pcaN_y - GenHVtx_.y(), s_pcaN_z - GenHVtx_.z());

            LV ipvP_lab_gen_ips(ipvP_lab_gen_3d_ips.Unit().x(), ipvP_lab_gen_3d_ips.Unit().y(), ipvP_lab_gen_3d_ips.Unit().z(), 0);
            LV ipvN_lab_gen_ips(ipvN_lab_gen_3d_ips.Unit().x(), ipvN_lab_gen_3d_ips.Unit().y(), ipvN_lab_gen_3d_ips.Unit().z(), 0);

            //boost to ZMF
            LV ipvP_rf_gen_ips = boost_to_rf_gen(ipvP_lab_gen_ips);
            LV ipvN_rf_gen_ips = boost_to_rf_gen(ipvN_lab_gen_ips);

            //Get only the position component of the vector
            PV ipvP_rf_gen_3d_ips = ipvP_rf_gen_ips.Vect();
            PV ipvN_rf_gen_3d_ips = ipvN_rf_gen_ips.Vect();

            //Get IP vector normal to pion momenta
            PV ipvP_rf_gen_3d_t_ips = ipvP_rf_gen_3d_ips - ipvP_rf_gen_3d_ips.Dot(genpionP_rf_2p_3d_u)*genpionP_rf_2p_3d_u;
            PV ipvN_rf_gen_3d_t_ips = ipvN_rf_gen_3d_ips - ipvN_rf_gen_3d_ips.Dot(genpionN_rf_2p_3d_u)*genpionN_rf_2p_3d_u;

            //Get normalized normal IP vector
            PV ipvP_rf_gen_3d_t_u_ips = ipvP_rf_gen_3d_t_ips.Unit();
            PV ipvN_rf_gen_3d_t_u_ips = ipvN_rf_gen_3d_t_ips.Unit();

            //Get CP angle in ZMF
            histos_["CPPhiStar_gen_pcas"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_ips.Dot(ipvN_rf_gen_3d_t_u_ips)));
            if(sqrt(ipvP_lab_gen_3d_ips.mag2()) > 0.003 && sqrt(ipvN_lab_gen_3d_ips.mag2()) > 0.003)
              histos_["CPPhiStar_gen_pcas_ip30"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_ips.Dot(ipvN_rf_gen_3d_t_u_ips)));
	    if(sqrt(ipvP_lab_gen_3d_ips.mag2()) > 0.005 && sqrt(ipvN_lab_gen_3d_ips.mag2()) > 0.005)
              histos_["CPPhiStar_gen_pcas_ip50"]->Fill(TMath::ACos(ipvP_rf_gen_3d_t_u_ips.Dot(ipvN_rf_gen_3d_t_u_ips)));
	    
	    //Compute the angle between pi+-pi- momentum in respective tau+/- rf.  check arxiv:1108.0670
	    //It should also be same as the angle bewteen the decay planes. 
	    //ROOT::Math::Boost boost_to_rf_gentauP(gentauP_lab.BoostToCM());
	    //ROOT::Math::Boost boost_to_rf_gentauN(gentauN_lab.BoostToCM());
	    
	    //boost to rf of respective taus
            LV genpionP_rf_tau = boost_to_rf_gentauP(genpionP_lab);
            LV genpionN_rf_tau = boost_to_rf_gentauN(genpionN_lab);

	    PV genpionP_rf_tau_3d = genpionP_rf_tau.Vect();
	    PV genpionN_rf_tau_3d = genpionN_rf_tau.Vect();

	    PV genpionP_rf_tau_3d_u = genpionP_rf_tau_3d/sqrt(genpionP_rf_tau_3d.mag2());
	    PV genpionN_rf_tau_3d_u = genpionN_rf_tau_3d/sqrt(genpionN_rf_tau_3d.mag2());

	    double phi_prime_gen = TMath::ACos(genpionP_rf_tau_3d_u.Dot(genpionN_rf_tau_3d_u));
	    
	    PV gentauN_rf_3d = gentauN_rf.Vect();
	    PV gentauN_rf_3d_u = gentauN_rf_3d/sqrt(gentauN_rf_3d.mag2());
	    double psi_prime_gen = TMath::ACos(gentauN_rf_3d_u.Dot(genpionP_rf_tau_3d_u.Cross(genpionN_rf_tau_3d_u)));

	    histos_["CPPhiPrime_gen"]->Fill(phi_prime_gen);
	    histos_["CPPsiPrime_gen"]->Fill(psi_prime_gen);
	    

	    //Check a method used by Pooja, using the tau momenta normal to pion momenta
	    PV gentauP_lab_3d = gentauP_lab.Vect();
	    PV gentauN_lab_3d = gentauN_lab.Vect();
	    PV gentauP_lab_3d_n = gentauP_lab_3d - (gentauP_lab_3d.Dot(genpionP_lab_3d_u))*genpionP_lab_3d_u;
	    PV gentauN_lab_3d_n = gentauN_lab_3d - (gentauN_lab_3d.Dot(genpionN_lab_3d_u))*genpionN_lab_3d_u;
	    PV gentauP_lab_3d_n_u = gentauP_lab_3d_n/sqrt(gentauP_lab_3d_n.mag2());
            PV gentauN_lab_3d_n_u = gentauN_lab_3d_n/sqrt(gentauN_lab_3d_n.mag2());
	    LV gentauP_lab_n_u = LV(gentauP_lab_3d_n_u.x(), gentauP_lab_3d_n_u.y(), gentauP_lab_3d_n_u.z(), 0);
	    LV gentauN_lab_n_u = LV(gentauN_lab_3d_n_u.x(), gentauN_lab_3d_n_u.y(), gentauN_lab_3d_n_u.z(), 0);
	    //boost to pi+pi- RF
	    LV gentauP_rf_2p_n = boost_to_rf_gen(gentauP_lab_n_u);
	    LV gentauN_rf_2p_n = boost_to_rf_gen(gentauN_lab_n_u);
	    PV gentauP_rf_2p_n_3d = gentauP_rf_2p_n.Vect();
	    PV gentauN_rf_2p_n_3d = gentauN_rf_2p_n.Vect(); 
	    PV gentauP_rf_2p_n_3d_t = gentauP_rf_2p_n_3d - (gentauP_rf_2p_n_3d.Dot(genpionP_rf_2p_3d_u))*genpionP_rf_2p_3d_u;
	    PV gentauN_rf_2p_n_3d_t = gentauN_rf_2p_n_3d - (gentauN_rf_2p_n_3d.Dot(genpionN_rf_2p_3d_u))*genpionN_rf_2p_3d_u;
	    PV gentauP_rf_2p_n_3d_t_u = gentauP_rf_2p_n_3d_t/sqrt(gentauP_rf_2p_n_3d_t.mag2());
	    PV gentauN_rf_2p_n_3d_t_u = gentauN_rf_2p_n_3d_t/sqrt(gentauN_rf_2p_n_3d_t.mag2());

	    double phi_star_gen_v2 = TMath::ACos(gentauP_rf_2p_n_3d_t_u.Dot(gentauN_rf_2p_n_3d_t_u));
	    histos_["CPPhiStar_gen_v2"]->Fill(phi_star_gen_v2);
	    
	    //cross check of IP vs Momenta
	    //check if tau vertex direction is parallel to tau momenta
	    PV tauPDecayVect(genpionPVtx_.x() - GenHVtx_.x(), genpionPVtx_.y() - GenHVtx_.y(), genpionPVtx_.z() - GenHVtx_.z());
	    PV tauNDecayVect(genpionNVtx_.x() - GenHVtx_.x(), genpionNVtx_.y() - GenHVtx_.y(), genpionNVtx_.z() - GenHVtx_.z());
	    histos_["xcheck_angle1"]->Fill(TMath::ACos(tauPDecayVect.Dot(gentauP_lab_3d)/(sqrt(tauPDecayVect.mag2())*sqrt(gentauP_lab_3d.mag2()))));
	    histos_["xcheck_angle1"]->Fill(TMath::ACos(tauNDecayVect.Dot(gentauN_lab_3d)/(sqrt(tauNDecayVect.mag2())*sqrt(gentauN_lab_3d.mag2()))));

	    if(verbosity)std::cout<<"tau+ Dir. Vect ("<<tauPDecayVect.Unit().x()<<", "<<tauPDecayVect.Unit().y()<<", "<<tauPDecayVect.Unit().z()<<")"<<std::endl;
	    if(verbosity)std::cout<<"tau+ Mom. Vect ("<<gentauP_lab_3d.Unit().x()<<", "<<gentauP_lab_3d.Unit().y()<<", "<<gentauP_lab_3d.Unit().z()<<")"<<std::endl;
	    if(verbosity)std::cout<<"tau- Dir. Vect ("<<tauNDecayVect.Unit().x()<<", "<<tauNDecayVect.Unit().y()<<", "<<tauNDecayVect.Unit().z()<<")"<<std::endl;
            if(verbosity)std::cout<<"tau- Mom. Vect ("<<gentauN_lab_3d.Unit().x()<<", "<<gentauN_lab_3d.Unit().y()<<", "<<gentauN_lab_3d.Unit().z()<<")"<<std::endl;
	    
	    //check if IP of pion is parallel to normal componennt of tau vector
	    histos_["xcheck_angle2"]->Fill(TMath::ACos(ipvP_lab_gen_3d.Dot(gentauP_lab_3d_n_u)/sqrt(ipvP_lab_gen_3d.mag2())));
	    histos_["xcheck_angle2"]->Fill(TMath::ACos(ipvN_lab_gen_3d.Dot(gentauN_lab_3d_n_u)/sqrt(ipvN_lab_gen_3d.mag2())));
	    
	    //check if IP of pion is parallel to normal componennt of tau vector, in pi-pi+ RF
	    histos_["xcheck_angle3"]->Fill(TMath::ACos(ipvP_rf_gen_3d.Dot(gentauP_rf_2p_n_3d)/(sqrt(ipvP_rf_gen_3d.mag2())*sqrt(gentauP_rf_2p_n_3d.mag2()))));
            histos_["xcheck_angle3"]->Fill(TMath::ACos(ipvN_rf_gen_3d.Dot(gentauN_rf_2p_n_3d)/(sqrt(ipvN_rf_gen_3d.mag2())*sqrt(gentauN_rf_2p_n_3d.mag2()))));

	    if(verbosity)std::cout<<"phi* in pi+ pi- rest-frame: method3 = "<<phi_star_gen<<std::endl;
	  }//end of gen pi-pi+ system
	}
    }//end of Generator studies

    //Take only one pair per events
    //tau pairs are ordered in 1st OS, the sum pT
    if(foundDiTauPair_ && run_ == lastRun_ &&
       lumi_ == lastLumi_ && event_ == lastEvent_)
      continue;

    foundDiTauPair_ = false;
    lastRun_ = run_;
    lastLumi_ = lumi_; 
    lastEvent_ = event_;
    
    bool passDecayMode = false;
    if(decayModeLeg1_ == 0 && decayModeLeg2_ == 0) passDecayMode = true;

    //Apply Medium DeltaBeta3Hits isolation
    bool passIsolation = false;
    if(tightestHPSDB3HWPLeg1_ >= 0 && tightestHPSDB3HWPLeg2_ >= 0)
      passIsolation = true;

    //apply kinematic cuts
    //double ptL1     = (*diTauLegsP4_)[0].Pt();
    //double ptL2     = (*diTauLegsP4_)[1].Pt();
    //double etaL1    = (*diTauLegsP4_)[0].Eta();
    //double etaL2    = (*diTauLegsP4_)[1].Eta();
    bool passKineCuts = false;
    if(ptL1 > 20 && ptL2 > 20 && TMath::Abs(etaL1) < 2.1 && TMath::Abs(etaL2) < 2.1)
      passKineCuts = true;

    bool passLeadChHadCut = false;
    if((*diTauLegsLchP4_)[0].Pt() > 0.0 && (*diTauLegsLchP4_)[1].Pt() > 0.0)
      passLeadChHadCut = true;

    //Apply a dZ and dxy cut wrt original vertex
    //Information not available in the current ntuple

    bool passDiTauCharge = false;
    if(diTauCharge_ == 0) passDiTauCharge = true;

    //use this pair only
    if(passDecayMode && passIsolation)
      foundDiTauPair_ = true;

    if(!passDecayMode) continue;
    if(!passIsolation) continue;
    if(!passKineCuts) continue;
    if(!passLeadChHadCut) continue;
    if(!passDiTauCharge) continue;

    //Quality cuts on track
    bool passQualityLeg1_ = false, passQualityLeg2_ = false;
    if(genDecayModeLeg1_ == 0 && trackPtLeg1_ > 20 && trackPtErrLeg1_/trackPtLeg1_ < 0.10 &&
       nHitsLeg1_ >= 8 && nPxlHitsLeg1_ >= 2 && nMisingHitsLeg1_ == 0 && TMath::Abs(etaL1) < 1.5)
      passQualityLeg1_ = true;
    if(genDecayModeLeg2_ == 0 && trackPtLeg2_ >20 && trackPtErrLeg2_/trackPtLeg2_ < 0.10 &&
       nHitsLeg2_ >= 8 && nPxlHitsLeg2_ >= 2 && nMisingHitsLeg2_ == 0 && TMath::Abs(etaL2) < 1.5)
      passQualityLeg2_ = true;
    
    bool passQualityLeg1V2_ = false, passQualityLeg2V2_ = false;
    if(genDecayModeLeg1_ == 0 && trackPtLeg1_ > 20 && trackPtErrLeg1_/trackPtLeg1_ < 0.10 &&
       nHitsLeg1_ >= 10 && nPxlHitsLeg1_ >= 3 && nMisingHitsLeg1_ == 0 && TMath::Abs(etaL1) < 1.5)
      passQualityLeg1V2_ = true;
    if(genDecayModeLeg2_ == 0 && trackPtLeg2_ > 20 && trackPtErrLeg2_/trackPtLeg2_ < 0.10 &&
       nHitsLeg2_ >= 10 && nPxlHitsLeg2_ >= 3 && nMisingHitsLeg2_ == 0 && TMath::Abs(etaL2) < 1.5)
      passQualityLeg2V2_ = true;

    bool passQualityLeg1V3_ = false, passQualityLeg2V3_ = false;
    if(trackPtLeg1_ > 20 && trackPtErrLeg1_/trackPtLeg1_ < 0.10 &&
       nHitsLeg1_ >= 8 && nPxlHitsLeg1_ >= 2 && nMisingHitsLeg1_ == 0 && TMath::Abs(etaL1) < 1.5)
      passQualityLeg1V3_ = true;
    if(trackPtLeg2_ >20 && trackPtErrLeg2_/trackPtLeg2_ < 0.10 &&
       nHitsLeg2_ >= 8 && nPxlHitsLeg2_ >= 2 && nMisingHitsLeg2_ == 0 && TMath::Abs(etaL2) < 1.5)
      passQualityLeg2V3_ = true;

    //Fill Histograms

    //Vertex Resolution
    histos_["VtxXRes"]->Fill(VtxXErr_/VtxX_);
    histos_["VtxYRes"]->Fill(VtxYErr_/VtxY_);
    histos_["VtxZRes"]->Fill(VtxZErr_/VtxZ_);

    histos_["ReFitVtxXRes"]->Fill(ReFitVtxXErr_/ReFitVtxX_);
    histos_["ReFitVtxYRes"]->Fill(ReFitVtxYErr_/ReFitVtxY_);
    histos_["ReFitVtxZRes"]->Fill(ReFitVtxZErr_/ReFitVtxZ_);
    
    histos_["VtxXErr"]->Fill(VtxXErr_);
    histos_["VtxYErr"]->Fill(VtxYErr_);
    histos_["VtxZErr"]->Fill(VtxZErr_);

    histos_["ReFitVtxXErr"]->Fill(ReFitVtxXErr_);
    histos_["ReFitVtxYErr"]->Fill(ReFitVtxYErr_);
    histos_["ReFitVtxZErr"]->Fill(ReFitVtxZErr_);
    
    histos_["DeltaPVX"]->Fill(VtxX_ - ReFitVtxX_);
    histos_["DeltaPVY"]->Fill(VtxY_ - ReFitVtxY_);
    histos_["DeltaPVZ"]->Fill(VtxZ_ - ReFitVtxZ_);

    histos_["DeltaPVwrtGenX"]->Fill(genVtxX_ - VtxX_);
    histos_["DeltaPVwrtGenY"]->Fill(genVtxY_ - VtxY_);
    histos_["DeltaPVwrtGenZ"]->Fill(genVtxZ_ - VtxZ_);

    histos_["DeltaRfPVwrtGenX"]->Fill(genVtxX_ - ReFitVtxX_);
    histos_["DeltaRfPVwrtGenY"]->Fill(genVtxY_ - ReFitVtxY_);
    histos_["DeltaRfPVwrtGenZ"]->Fill(genVtxZ_ - ReFitVtxZ_);

    histos_["ResPVwrtGenX"]->Fill((genVtxX_ - VtxX_)/VtxXErr_);
    histos_["ResPVwrtGenY"]->Fill((genVtxY_ - VtxY_)/VtxYErr_);
    histos_["ResPVwrtGenZ"]->Fill((genVtxZ_ - VtxZ_)/VtxZErr_);

    histos_["ResRfPVwrtGenX"]->Fill((genVtxX_ - ReFitVtxX_)/ReFitVtxXErr_);
    histos_["ResRfPVwrtGenY"]->Fill((genVtxY_ - ReFitVtxY_)/ReFitVtxYErr_);
    histos_["ResRfPVwrtGenZ"]->Fill((genVtxZ_ - ReFitVtxZ_)/ReFitVtxZErr_);

    histos_["DeltaIanPVX"]->Fill(VtxX_ - IanVtxX_);
    histos_["DeltaIanPVY"]->Fill(VtxY_ - IanVtxY_);
    histos_["DeltaIanPVZ"]->Fill(VtxZ_ - IanVtxZ_);

    histos_["DeltaIanPVwrtGenX"]->Fill(genVtxX_ - IanVtxX_);
    histos_["DeltaIanPVwrtGenY"]->Fill(genVtxY_ - IanVtxY_);
    histos_["DeltaIanPVwrtGenZ"]->Fill(genVtxZ_ - IanVtxZ_);

    histos_["VtxNChi2"]->Fill(VtxNChi2_);
    histos_["ReFitVtxNChi2"]->Fill(ReFitVtxNChi2_);
    histos_["IanVtxNChi2"]->Fill(IanVtxNChi2_);

    //if(verbosity)std::cout<<" diTau pT"<<(*diTauLegsP4_)[0].Pt()+(*diTauLegsP4_)[1].Pt()<<std::endl;
    
    // Compute the CP angle at Detector level, using re-fitted vetex///////////////////////////////////////////////
    //Make Pi+Pi- Zero-momentum-frame
    LV pion1_lab = ((*diTauLegsLchP4_)[0]);
    LV pion2_lab = ((*diTauLegsLchP4_)[1]);
    LV pion_pair_lab = pion1_lab + pion2_lab;
    
    ROOT::Math::Boost boost_to_rf(pion_pair_lab.BoostToCM());
    //ROOT::Math::Boost boost_to_lab(boost_to_rf.Inverse());

    if(verbosity)std::cout<<" pion1 pt "<<pion1_lab.Pt()<<" pion2 pt "<<pion2_lab.Pt()<<std::endl;
    LV pion1_rf = boost_to_rf(pion1_lab);
    LV pion2_rf = boost_to_rf(pion2_lab);

    //Define the normalized impact parameter vectors
    PV ipv1_lab_3d((*diTauLegsPCA_)[0].x() - ReFitVtxX_, (*diTauLegsPCA_)[0].y() - ReFitVtxY_, (*diTauLegsPCA_)[0].z() - ReFitVtxZ_);
    PV ipv2_lab_3d((*diTauLegsPCA_)[1].x() - ReFitVtxX_, (*diTauLegsPCA_)[1].y() - ReFitVtxY_, (*diTauLegsPCA_)[1].z() - ReFitVtxZ_);

    if(verbosity)std::cout<<"Reco Refit Vertex ("<<ReFitVtxX_<<","<<ReFitVtxY_<<","<<ReFitVtxZ_<<")"<<std::endl;
    if(verbosity)std::cout<<"PCA wrt Refit Vertex pion 1 ("<<(*diTauLegsPCA_)[0].x()<<","<<(*diTauLegsPCA_)[0].y()<<","<<(*diTauLegsPCA_)[0].z()<<")"<<std::endl;
    
    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco"]->Fill(TMath::ACos(ipv1_lab_3d.Unit().Dot(pion1_lab.Vect().Unit())));
    histos_["xcheck_angle3_reco"]->Fill(TMath::ACos(ipv2_lab_3d.Unit().Dot(pion2_lab.Vect().Unit())));
    Vector2D ipv1_lab_xy(ipv1_lab_3d.x(), ipv1_lab_3d.y());
    Vector2D ipv2_lab_xy(ipv2_lab_3d.x(), ipv2_lab_3d.y());
    Vector2D pion1_lab_xy(pion1_lab.px(), pion1_lab.py());
    Vector2D pion2_lab_xy(pion2_lab.px(), pion2_lab.py());
    histos_["xcheck_angle3_reco_2d"]->Fill(TMath::ACos(ipv1_lab_xy.Unit().Dot(pion1_lab_xy.Unit())));
    histos_["xcheck_angle3_reco_2d"]->Fill(TMath::ACos(ipv2_lab_xy.Unit().Dot(pion2_lab_xy.Unit())));

    LV ipv1_lab(ipv1_lab_3d.Unit().x(), ipv1_lab_3d.Unit().y(), ipv1_lab_3d.Unit().z(), 0); 
    LV ipv2_lab(ipv2_lab_3d.Unit().x(), ipv2_lab_3d.Unit().y(), ipv2_lab_3d.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_rf = boost_to_rf(ipv1_lab);
    LV ipv2_rf = boost_to_rf(ipv2_lab);
    
    //Get only the position component of the vector
    PV ipv1_rf_3d = ipv1_rf.Vect();
    PV ipv2_rf_3d = ipv2_rf.Vect();

    PV pion1_rf_3d = pion1_rf.Vect();
    PV pion2_rf_3d = pion2_rf.Vect();

    //Get the unit (normalized) vector along pion momentum 
    PV pion1_rf_3d_u = pion1_rf_3d/TMath::Sqrt(pion1_rf_3d.mag2());
    PV pion2_rf_3d_u = pion2_rf_3d/TMath::Sqrt(pion2_rf_3d.mag2());

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_rf_3d_l = ipv1_rf_3d.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_l = ipv2_rf_3d.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_t = ipv1_rf_3d - ipv1_rf_3d_l;
    PV ipv2_rf_3d_t = ipv2_rf_3d - ipv2_rf_3d_l;

    //Get normalized normal IP vector
    PV ipv1_rf_3d_t_u = ipv1_rf_3d_t/TMath::Sqrt(ipv1_rf_3d_t.mag2());
    PV ipv2_rf_3d_t_u = ipv2_rf_3d_t/TMath::Sqrt(ipv2_rf_3d_t.mag2());

    //Get CP angle in ZMF
    double phi_star = TMath::ACos(ipv1_rf_3d_t_u.Dot(ipv2_rf_3d_t_u));
    histos_["CPPhiStar"]->Fill(phi_star);
    //Get CP angle in lab frame
    double phi_lab = TMath::ACos(ipv1_lab_3d.Unit().Dot(ipv2_lab_3d.Unit()));
    histos_["CPPhiLab"]->Fill(phi_lab);
    //if(sqrt(ipv1_lab_3d.mag2()) > 0.005 && sqrt(ipv2_lab_3d.mag2()) > 0.005){
    //  histos_["CPPhiStar_ipcut"]->Fill(phi_star);
    //  histos_["CPPhiLab_ipcut"]->Fill(phi_lab);
    //}
    double OstarCP_reco = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_rf_3d_t_u.Cross(ipv2_rf_3d_t_u))
      : pion2_rf_3d_u.Dot(ipv2_rf_3d_t_u.Cross(ipv1_rf_3d_t_u));
    double phi_star_cp = (OstarCP_reco >= 0) ? phi_star : (TMath::TwoPi() - phi_star);
    histos_["CPPhiStar_CP"]->Fill(phi_star_cp);

    //replace IP vector by the component of IP normal to pion momenta
    PV ipv1_lab_3d_NP = ipv1_lab_3d - ipv1_lab_3d.Dot(pion1_lab.Vect().Unit())*pion1_lab.Vect().Unit();
    PV ipv2_lab_3d_NP = ipv2_lab_3d - ipv2_lab_3d.Dot(pion2_lab.Vect().Unit())*pion2_lab.Vect().Unit();

    LV ipv1_lab_NP(ipv1_lab_3d_NP.Unit().x(), ipv1_lab_3d_NP.Unit().y(), ipv1_lab_3d_NP.Unit().z(), 0);
    LV ipv2_lab_NP(ipv2_lab_3d_NP.Unit().x(), ipv2_lab_3d_NP.Unit().y(), ipv2_lab_3d_NP.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_rf_NP = boost_to_rf(ipv1_lab_NP);
    LV ipv2_rf_NP = boost_to_rf(ipv2_lab_NP);

    //Get only the position component of the vector
    PV ipv1_rf_3d_NP = ipv1_rf_NP.Vect();
    PV ipv2_rf_3d_NP = ipv2_rf_NP.Vect();

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_NP_t = ipv1_rf_3d_NP - ipv1_rf_3d_NP.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_NP_t = ipv2_rf_3d_NP - ipv2_rf_3d_NP.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get normalized normal IP vector
    PV ipv1_rf_3d_NP_t_u = ipv1_rf_3d_NP_t.Unit();
    PV ipv2_rf_3d_NP_t_u = ipv2_rf_3d_NP_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_NP"]->Fill(TMath::ACos(ipv1_rf_3d_NP_t_u.Dot(ipv2_rf_3d_NP_t_u)));
    //Get CP angle in lab frame
    histos_["CPPhiLab_NP"]->Fill(TMath::ACos(ipv1_lab_3d_NP.Unit().Dot(ipv2_lab_3d_NP.Unit())));
    if(sqrt(ipv1_lab_3d.mag2()) > 0.005 && sqrt(ipv2_lab_3d.mag2()) > 0.005){
      histos_["CPPhiStar_NP_ipcut"]->Fill(TMath::ACos(ipv1_rf_3d_NP_t_u.Dot(ipv2_rf_3d_NP_t_u)));
      histos_["CPPhiLab_NP_ipcut"]->Fill(TMath::ACos(ipv1_lab_3d_NP.Unit().Dot(ipv2_lab_3d_NP.Unit())));
    }

    ///////Use momentum of charged hadron at PCA ////////////////////
    double m_pi = 0.13957018;
    LV pion1_lab_pca((*diTauLegsLchP3AtPCA_)[0].x(), (*diTauLegsLchP3AtPCA_)[0].y(), (*diTauLegsLchP3AtPCA_)[0].z(), TMath::Sqrt((*diTauLegsLchP3AtPCA_)[0].mag2() + m_pi*m_pi) );
    LV pion2_lab_pca((*diTauLegsLchP3AtPCA_)[1].x(), (*diTauLegsLchP3AtPCA_)[1].y(), (*diTauLegsLchP3AtPCA_)[1].z(), TMath::Sqrt((*diTauLegsLchP3AtPCA_)[1].mag2() + m_pi*m_pi) );
    LV pion_pair_lab_pca = pion1_lab_pca + pion2_lab_pca;
    if(verbosity)std::cout<<" At PCA, pion1 pt "<<pion1_lab_pca.Pt()<<" pion2 pt "<<pion2_lab_pca.Pt()<<std::endl;
    ROOT::Math::Boost boost_to_rf_pca(pion_pair_lab_pca.BoostToCM());

    LV pion1_rf_pca = boost_to_rf_pca(pion1_lab_pca);
    LV pion2_rf_pca = boost_to_rf_pca(pion2_lab_pca);

    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_pca"]->Fill(TMath::ACos(ipv1_lab_3d.Unit().Dot(pion1_lab_pca.Vect().Unit())));
    histos_["xcheck_angle3_reco_pca"]->Fill(TMath::ACos(ipv2_lab_3d.Unit().Dot(pion2_lab_pca.Vect().Unit())));
    Vector2D pion1_lab_pca_xy(pion1_lab_pca.px(), pion1_lab_pca.py());
    Vector2D pion2_lab_pca_xy(pion2_lab_pca.px(), pion2_lab_pca.py());
    histos_["xcheck_angle3_reco_pca_2d"]->Fill(TMath::ACos(ipv1_lab_xy.Unit().Dot(pion1_lab_pca_xy.Unit())));
    histos_["xcheck_angle3_reco_pca_2d"]->Fill(TMath::ACos(ipv2_lab_xy.Unit().Dot(pion2_lab_pca_xy.Unit())));

    //boost to ZMF
    LV ipv1_rf_pca = boost_to_rf_pca(ipv1_lab);
    LV ipv2_rf_pca = boost_to_rf_pca(ipv2_lab);

    //Get only the position component of the vector
    PV ipv1_rf_pca_3d = ipv1_rf_pca.Vect();
    PV ipv2_rf_pca_3d = ipv2_rf_pca.Vect();

    PV pion1_rf_pca_3d = pion1_rf_pca.Vect();
    PV pion2_rf_pca_3d = pion2_rf_pca.Vect();

    //Get the unit (normalized) vector along pion momentum
    PV pion1_rf_pca_3d_u = pion1_rf_pca_3d.Unit();
    PV pion2_rf_pca_3d_u = pion2_rf_pca_3d.Unit();

    //Get IP vector normal to pion momenta
    PV ipv1_rf_pca_3d_t = ipv1_rf_3d - ipv1_rf_pca_3d.Dot(pion1_rf_pca_3d_u)*pion1_rf_pca_3d_u;
    PV ipv2_rf_pca_3d_t = ipv2_rf_3d - ipv2_rf_pca_3d.Dot(pion2_rf_pca_3d_u)*pion2_rf_pca_3d_u;

    //Get normalized normal IP vector
    PV ipv1_rf_pca_3d_t_u = ipv1_rf_pca_3d_t.Unit();
    PV ipv2_rf_pca_3d_t_u = ipv2_rf_pca_3d_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_P4AtPCA"]->Fill(TMath::ACos(ipv1_rf_pca_3d_t_u.Dot(ipv2_rf_pca_3d_t_u)));

    //Use PCA with Analytical extrapolator
    //Define the normalized impact parameter vectors
    PV ipv1_lab_3d_m2((*diTauLegsPCAM2_)[0].x() - ReFitVtxX_, (*diTauLegsPCAM2_)[0].y() - ReFitVtxY_, (*diTauLegsPCAM2_)[0].z() - ReFitVtxZ_);
    PV ipv2_lab_3d_m2((*diTauLegsPCAM2_)[1].x() - ReFitVtxX_, (*diTauLegsPCAM2_)[1].y() - ReFitVtxY_, (*diTauLegsPCAM2_)[1].z() - ReFitVtxZ_);

    LV ipv1_lab_m2(ipv1_lab_3d_m2.Unit().x(), ipv1_lab_3d_m2.Unit().y(), ipv1_lab_3d_m2.Unit().z(), 0);
    LV ipv2_lab_m2(ipv2_lab_3d_m2.Unit().x(), ipv2_lab_3d_m2.Unit().y(), ipv2_lab_3d_m2.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_rf_m2 = boost_to_rf(ipv1_lab_m2);
    LV ipv2_rf_m2 = boost_to_rf(ipv2_lab_m2);

    //Get only the position component of the vector
    PV ipv1_rf_3d_m2 = ipv1_rf_m2.Vect();
    PV ipv2_rf_3d_m2 = ipv2_rf_m2.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_rf_3d_m2_l = ipv1_rf_3d_m2.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_m2_l = ipv2_rf_3d_m2.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_m2_t = ipv1_rf_3d_m2 - ipv1_rf_3d_m2_l;
    PV ipv2_rf_3d_m2_t = ipv2_rf_3d_m2 - ipv2_rf_3d_m2_l;

    //Get normalized normal IP vector
    PV ipv1_rf_3d_m2_t_u = ipv1_rf_3d_m2_t.unit();
    PV ipv2_rf_3d_m2_t_u = ipv2_rf_3d_m2_t.unit();
    histos_["CPPhiStar_M2"]->Fill(TMath::ACos(ipv1_rf_3d_m2_t_u.Dot(ipv2_rf_3d_m2_t_u)));
    histos_["CPPhiLab_M2"]->Fill(TMath::ACos(ipv1_lab_3d_m2.Unit().Dot(ipv2_lab_3d_m2.Unit())));
    if(genDecayModeLeg1_ == 0 && genDecayModeLeg2_ == 0)
      histos_["CPPhiStar_M2_gen1p0pi0"]->Fill(TMath::ACos(ipv1_rf_3d_m2_t_u.Dot(ipv2_rf_3d_m2_t_u)));
    if(genDecayModeLeg1_ >= 0 && genDecayModeLeg2_ >= 0)
      histos_["CPPhiStar_M2_gentau"]->Fill(TMath::ACos(ipv1_rf_3d_m2_t_u.Dot(ipv2_rf_3d_m2_t_u)));

    //Compute error on PCA using the cov. matrix
    SymMatrix33 covMatrix1 = (*diTauLegsPCAM2Cov_)[0];
    SymMatrix33 covMatrix2 = (*diTauLegsPCAM2Cov_)[1];
    AlgebraicVector3 dir_pion1, dir_pion2;
    dir_pion1[0] = pion1_lab.Vect().Unit().x();
    dir_pion1[1] = pion1_lab.Vect().Unit().y();
    dir_pion1[2] = pion1_lab.Vect().Unit().z();
    dir_pion2[0] = pion2_lab.Vect().Unit().x();
    dir_pion2[1] = pion2_lab.Vect().Unit().y();
    dir_pion2[2] = pion2_lab.Vect().Unit().z();

    AlgebraicVector3 dir_ip1, dir_ip2;
    dir_ip1[0] = ipv1_lab_3d_m2.Unit().x();
    dir_ip1[1] = ipv1_lab_3d_m2.Unit().y();
    dir_ip1[2] = ipv1_lab_3d_m2.Unit().z();
    dir_ip2[0] = ipv2_lab_3d_m2.Unit().x();
    dir_ip2[1] = ipv2_lab_3d_m2.Unit().y();
    dir_ip2[2] = ipv2_lab_3d_m2.Unit().z();
    
    PV normal_decayPlane_pion1 = pion1_lab.Vect().Cross(ipv1_lab_3d_m2);
    PV normal_decayPlane_pion2 = pion2_lab.Vect().Cross(ipv2_lab_3d_m2);
    AlgebraicVector3 dir_dp1, dir_dp2;
    dir_dp1[0] = normal_decayPlane_pion1.Unit().x();
    dir_dp1[1] = normal_decayPlane_pion1.Unit().y();
    dir_dp1[2] = normal_decayPlane_pion1.Unit().z();
    dir_dp2[0] = normal_decayPlane_pion2.Unit().x();
    dir_dp2[1] = normal_decayPlane_pion2.Unit().y();
    dir_dp2[2] = normal_decayPlane_pion2.Unit().z();

    double err_p1_l = ROOT::Math::Similarity(dir_pion1, covMatrix1);
    double err_p2_l = ROOT::Math::Similarity(dir_pion2, covMatrix2);
    double err_p1_t = ROOT::Math::Similarity(dir_ip1, covMatrix1);
    double err_p2_t = ROOT::Math::Similarity(dir_ip2, covMatrix2);
    double err_p1_tdp = ROOT::Math::Similarity(dir_dp1, covMatrix1);
    double err_p2_tdp = ROOT::Math::Similarity(dir_dp2, covMatrix2);
    
    histos_["ip_error_m2_l"]->Fill(err_p1_l);
    histos_["ip_error_m2_l"]->Fill(err_p2_l);
    histos_["ip_error_m2_t"]->Fill(err_p1_t);
    histos_["ip_error_m2_t"]->Fill(err_p2_t);
    histos_["ip_error_m2_tdp"]->Fill(err_p1_tdp);
    histos_["ip_error_m2_tdp"]->Fill(err_p2_tdp);
    
    /////////////Check using PCA with 3d extrapolation method////////////////////////////
    PV ipv1_lab_3d_le1 =  (*diTauLegsIPAtPCA_)[0];
    PV ipv2_lab_3d_le1 =  (*diTauLegsIPAtPCA_)[1];

    LV ipv1_lab_le1(ipv1_lab_3d_le1.unit().x(), ipv1_lab_3d_le1.unit().y(), ipv1_lab_3d_le1.unit().z(), 0);
    LV ipv2_lab_le1(ipv2_lab_3d_le1.unit().x(), ipv2_lab_3d_le1.unit().y(), ipv2_lab_3d_le1.unit().z(), 0);

    //boost to ZMF
    LV ipv1_rf_le1 = boost_to_rf(ipv1_lab_le1);
    LV ipv2_rf_le1 = boost_to_rf(ipv2_lab_le1);

    //Get only the position component of the vector
    PV ipv1_rf_3d_le1 = ipv1_rf_le1.Vect();
    PV ipv2_rf_3d_le1 = ipv2_rf_le1.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_rf_3d_le1_l = ipv1_rf_3d_le1.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_le1_l = ipv2_rf_3d_le1.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_le1_t = ipv1_rf_3d_le1 - ipv1_rf_3d_le1_l;
    PV ipv2_rf_3d_le1_t = ipv2_rf_3d_le1 - ipv2_rf_3d_le1_l;

       //Get normalized normal IP vector
    PV ipv1_rf_3d_le1_t_u = ipv1_rf_3d_le1_t.unit();
    PV ipv2_rf_3d_le1_t_u = ipv2_rf_3d_le1_t.unit();
    histos_["CPPhiStar_LE1"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
    histos_["CPPhiLab_LE1"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v1 "<<TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u))<<std::endl;
    double OstarCP_reco_le1 = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_rf_3d_le1_t_u.Cross(ipv2_rf_3d_le1_t_u))
      : pion2_rf_3d_u.Dot(ipv2_rf_3d_le1_t_u.Cross(ipv1_rf_3d_le1_t_u));
    double phi_star_cp_le1 = (OstarCP_reco_le1 >= 0) ? TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)) 
      : (TMath::TwoPi() - TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
    histos_["CPPhiStar_CP_LE1"]->Fill(phi_star_cp_le1);

    if(genDecayModeLeg1_ == 0 && genDecayModeLeg2_ == 0)
      histos_["CPPhiStar_LE1_gen1p0pi0"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
    if(genDecayModeLeg1_ >= 0 && genDecayModeLeg2_ >= 0)
      histos_["CPPhiStar_LE1_gentau"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));

    if(sqrt(ipv1_lab_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_3d_le1.mag2()) > 0.005){
      histos_["CPPhiStar_ipcut"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
      histos_["CPPhiLab_ipcut"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
      if(passQualityLeg1_ && passQualityLeg2_){
	histos_["CPPhiStar_ipcut_TkQCut"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
	histos_["CPPhiLab_ipcut_TkQCut"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
      }
      if(passQualityLeg1V2_ && passQualityLeg2V2_){
	histos_["CPPhiStar_ipcut_TkQCutV2"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
        //histos_["CPPhiLab_ipcut_TkQCutV2"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
      }
      if(passQualityLeg1V3_ && passQualityLeg2V3_)
	histos_["CPPhiStar_ipcut_TkQCutV3"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
    }
    if(sqrt(ipv1_lab_3d_le1.mag2()) > 0.003 && sqrt(ipv2_lab_3d_le1.mag2()) > 0.003){
      histos_["CPPhiStar_ipcut30"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
      histos_["CPPhiLab_ipcut30"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
      if(passQualityLeg1_ && passQualityLeg2_){
	histos_["CPPhiStar_ipcut30_TkQCut"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
        histos_["CPPhiLab_ipcut30_TkQCut"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
      } 
      if(passQualityLeg1V2_ && passQualityLeg2V2_){
        histos_["CPPhiStar_ipcut30_TkQCutV2"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
        //histos_["CPPhiLab_ipcut30_TkQCutV2"]->Fill(TMath::ACos(ipv1_lab_3d_le1.Unit().Dot(ipv2_lab_3d_le1.Unit())));
      }
    }
    if(ptL1 > 40 && ptL2 > 40)
      histos_["CPPhiStar_tau40"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
    if(ptL1 > 30 && ptL2 > 30)
      histos_["CPPhiStar_tau30"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));

    PV ipv1_lab_3d_le2 =  (*diTauLegsIPAtPCAV2_)[0];
    PV ipv2_lab_3d_le2 =  (*diTauLegsIPAtPCAV2_)[1];

    LV ipv1_lab_le2(ipv1_lab_3d_le2.unit().x(), ipv1_lab_3d_le2.unit().y(), ipv1_lab_3d_le2.unit().z(), 0);
    LV ipv2_lab_le2(ipv2_lab_3d_le2.unit().x(), ipv2_lab_3d_le2.unit().y(), ipv2_lab_3d_le2.unit().z(), 0);

    //boost to ZMF
    LV ipv1_rf_le2 = boost_to_rf(ipv1_lab_le2);
    LV ipv2_rf_le2 = boost_to_rf(ipv2_lab_le2);

    //Get only the position component of the vector
    PV ipv1_rf_3d_le2 = ipv1_rf_le2.Vect();
    PV ipv2_rf_3d_le2 = ipv2_rf_le2.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_rf_3d_le2_l = ipv1_rf_3d_le2.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_le2_l = ipv2_rf_3d_le2.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_le2_t = ipv1_rf_3d_le2 - ipv1_rf_3d_le2_l;
    PV ipv2_rf_3d_le2_t = ipv2_rf_3d_le2 - ipv2_rf_3d_le2_l;
    
    //Get normalized normal IP vector
    PV ipv1_rf_3d_le2_t_u = ipv1_rf_3d_le2_t.unit();
    PV ipv2_rf_3d_le2_t_u = ipv2_rf_3d_le2_t.unit();
    histos_["CPPhiStar_LE2"]->Fill(TMath::ACos(ipv1_rf_3d_le2_t_u.Dot(ipv2_rf_3d_le2_t_u)));
    histos_["CPPhiLab_LE2"]->Fill(TMath::ACos(ipv1_lab_3d_le2.Unit().Dot(ipv2_lab_3d_le2.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v2 "<<TMath::ACos(ipv1_rf_3d_le2_t_u.Dot(ipv2_rf_3d_le2_t_u))<<std::endl;
    
    ///////Use Vertex with Ian's method ///////////
    //Define the normalized impact parameter vectors
    PV ipv1_lab_3d_ian((*diTauLegsPCAIan_)[0].x() - IanVtxX_, (*diTauLegsPCAIan_)[0].y() - IanVtxY_, (*diTauLegsPCAIan_)[0].z() - IanVtxZ_);
    PV ipv2_lab_3d_ian((*diTauLegsPCAIan_)[1].x() - IanVtxX_, (*diTauLegsPCAIan_)[1].y() - IanVtxY_, (*diTauLegsPCAIan_)[1].z() - IanVtxZ_);

    LV ipv1_lab_ian(ipv1_lab_3d_ian.Unit().x(), ipv1_lab_3d_ian.Unit().y(), ipv1_lab_3d_ian.Unit().z(), 0);
    LV ipv2_lab_ian(ipv2_lab_3d_ian.Unit().x(), ipv2_lab_3d_ian.Unit().y(), ipv2_lab_3d_ian.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_rf_ian = boost_to_rf(ipv1_lab_ian);
    LV ipv2_rf_ian = boost_to_rf(ipv2_lab_ian);

    //Get only the position component of the vector
    PV ipv1_rf_3d_ian = ipv1_rf_ian.Vect();
    PV ipv2_rf_3d_ian = ipv2_rf_ian.Vect();

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_ian_t = ipv1_rf_3d_ian - ipv1_rf_3d_ian.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_ian_t = ipv2_rf_3d_ian - ipv2_rf_3d_ian.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get normalized normal IP vector
    PV ipv1_rf_3d_ian_t_u = ipv1_rf_3d_ian_t.unit();
    PV ipv2_rf_3d_ian_t_u = ipv2_rf_3d_ian_t.unit();
    histos_["CPPhiStar_Ian"]->Fill(TMath::ACos(ipv1_rf_3d_ian_t_u.Dot(ipv2_rf_3d_ian_t_u)));

    /////// Compute the CP angle at Detector level, using Original vertex///////////////////////////////////////////////
    //Define the normalized impact parameter vectors
    PV ipv1_lab_opv_3d((*diTauLegsPCAOPV_)[0].x() - VtxX_, (*diTauLegsPCAOPV_)[0].y() - VtxY_, (*diTauLegsPCAOPV_)[0].z() - VtxZ_);
    PV ipv2_lab_opv_3d((*diTauLegsPCAOPV_)[1].x() - VtxX_, (*diTauLegsPCAOPV_)[1].y() - VtxY_, (*diTauLegsPCAOPV_)[1].z() - VtxZ_);

    LV ipv1_lab_opv(ipv1_lab_opv_3d.Unit().x(), ipv1_lab_opv_3d.Unit().y(), ipv1_lab_opv_3d.Unit().z(), 0);
    LV ipv2_lab_opv(ipv2_lab_opv_3d.Unit().x(), ipv2_lab_opv_3d.Unit().y(), ipv2_lab_opv_3d.Unit().z(), 0);

    if(verbosity)std::cout<<"Reco P. Vertex ("<<VtxX_<<","<<VtxY_<<","<<VtxZ_<<")"<<std::endl;
    if(verbosity)std::cout<<"PCA wrt P. Vertex pion 1 ("<<(*diTauLegsPCAOPV_)[0].x()<<","<<(*diTauLegsPCAOPV_)[0].y()<<","<<(*diTauLegsPCAOPV_)[0].z()<<")"<<std::endl;
    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_opv"]->Fill(TMath::ACos(ipv1_lab_opv_3d.Unit().Dot(pion1_lab.Vect().Unit())));
    histos_["xcheck_angle3_reco_opv"]->Fill(TMath::ACos(ipv2_lab_opv_3d.Unit().Dot(pion2_lab.Vect().Unit())));
    Vector2D ipv1_lab_opv_xy(ipv1_lab_opv_3d.x(), ipv1_lab_opv_3d.y());
    Vector2D ipv2_lab_opv_xy(ipv2_lab_opv_3d.x(), ipv2_lab_opv_3d.y());
    histos_["xcheck_angle3_reco_opv_2d"]->Fill(TMath::ACos(ipv1_lab_opv_xy.Unit().Dot(pion1_lab_xy.Unit())));
    histos_["xcheck_angle3_reco_opv_2d"]->Fill(TMath::ACos(ipv2_lab_opv_xy.Unit().Dot(pion2_lab_xy.Unit())));

    //boost to ZMF
    LV ipv1_opv_rf = boost_to_rf(ipv1_lab_opv);
    LV ipv2_opv_rf = boost_to_rf(ipv2_lab_opv);

    //Get only the position component of the vector
    PV ipv1_opv_rf_3d = ipv1_opv_rf.Vect();
    PV ipv2_opv_rf_3d = ipv2_opv_rf.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_opv_rf_3d_l = ipv1_opv_rf_3d.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_opv_rf_3d_l = ipv2_opv_rf_3d.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_opv_rf_3d_t = ipv1_opv_rf_3d - ipv1_opv_rf_3d_l;
    PV ipv2_opv_rf_3d_t = ipv2_opv_rf_3d - ipv2_opv_rf_3d_l;

    //Get normalized normal IP vector
    PV ipv1_opv_rf_3d_t_u = ipv1_opv_rf_3d_t/TMath::Sqrt(ipv1_opv_rf_3d_t.mag2());
    PV ipv2_opv_rf_3d_t_u = ipv2_opv_rf_3d_t/TMath::Sqrt(ipv2_opv_rf_3d_t.mag2());

    //Get CP angle in ZMF
    double phi_star_opv = TMath::ACos(ipv1_opv_rf_3d_t_u.Dot(ipv2_opv_rf_3d_t_u));
    histos_["CPPhiStar_opv"]->Fill(phi_star_opv);
    histos_["CPPhiLab_opv"]->Fill(TMath::ACos(ipv1_lab_opv_3d.Unit().Dot(ipv2_lab_opv_3d.Unit())));
    //if(sqrt(ipv1_lab_opv_3d.mag2()) > 0.005 && sqrt(ipv2_lab_opv_3d.mag2()) > 0.005){
    //  histos_["CPPhiStar_opv_ipcut"]->Fill(phi_star_opv);
    //  histos_["CPPhiLab_opv_ipcut"]->Fill(TMath::ACos(ipv1_lab_opv_3d.Unit().Dot(ipv2_lab_opv_3d.Unit())));
    //}
    double OstarCP_opv = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_opv_rf_3d_t_u.Cross(ipv2_opv_rf_3d_t_u))
      : pion2_rf_3d_u.Dot(ipv2_opv_rf_3d_t_u.Cross(ipv1_opv_rf_3d_t_u));
    double phi_star_opv_cp = (OstarCP_opv >= 0) ? phi_star_opv : (TMath::TwoPi() - phi_star_opv);
    histos_["CPPhiStar_opv_CP"]->Fill(phi_star_opv_cp);

    //replace IP vector by the component of IP normal to pion momenta
    PV ipv1_lab_opv_3d_NP = ipv1_lab_opv_3d - ipv1_lab_opv_3d.Dot(pion1_lab.Vect().Unit())*pion1_lab.Vect().Unit();
    PV ipv2_lab_opv_3d_NP = ipv2_lab_opv_3d - ipv2_lab_opv_3d.Dot(pion2_lab.Vect().Unit())*pion2_lab.Vect().Unit();

    LV ipv1_lab_opv_NP(ipv1_lab_opv_3d_NP.Unit().x(), ipv1_lab_opv_3d_NP.Unit().y(), ipv1_lab_opv_3d_NP.Unit().z(), 0);
    LV ipv2_lab_opv_NP(ipv2_lab_opv_3d_NP.Unit().x(), ipv2_lab_opv_3d_NP.Unit().y(), ipv2_lab_opv_3d_NP.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_opv_rf_NP = boost_to_rf(ipv1_lab_opv_NP);
    LV ipv2_opv_rf_NP = boost_to_rf(ipv2_lab_opv_NP);

    //Get only the position component of the vector
    PV ipv1_opv_rf_3d_NP = ipv1_opv_rf_NP.Vect();
    PV ipv2_opv_rf_3d_NP = ipv2_opv_rf_NP.Vect();

    //Get IP vector normal to pion momenta
    PV ipv1_opv_rf_3d_NP_t = ipv1_opv_rf_3d_NP - ipv1_opv_rf_3d_NP.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_opv_rf_3d_NP_t = ipv2_opv_rf_3d_NP - ipv2_opv_rf_3d_NP.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get normalized normal IP vector
    PV ipv1_opv_rf_3d_NP_t_u = ipv1_opv_rf_3d_NP_t.Unit();
    PV ipv2_opv_rf_3d_NP_t_u = ipv2_opv_rf_3d_NP_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_opv_NP"]->Fill(TMath::ACos(ipv1_opv_rf_3d_NP_t_u.Dot(ipv2_opv_rf_3d_NP_t_u)));
    //Get CP angle in lab frame
    histos_["CPPhiLab_opv_NP"]->Fill(TMath::ACos(ipv1_lab_opv_3d_NP.Unit().Dot(ipv2_lab_opv_3d_NP.Unit())));
    if(sqrt(ipv1_lab_opv_3d.mag2()) > 0.005 && sqrt(ipv2_lab_opv_3d.mag2()) > 0.005){
      histos_["CPPhiStar_opv_NP_ipcut"]->Fill(TMath::ACos(ipv1_opv_rf_3d_NP_t_u.Dot(ipv2_opv_rf_3d_NP_t_u)));
      histos_["CPPhiLab_opv_NP_ipcut"]->Fill(TMath::ACos(ipv1_lab_opv_3d_NP.Unit().Dot(ipv2_lab_opv_3d_NP.Unit())));
    }

    ///////Use momentum of charged hadron at PCA ////////////////////
    LV pion1_lab_opv_pca((*diTauLegsLchP3AtPCAOPV_)[0].x(), (*diTauLegsLchP3AtPCAOPV_)[0].y(), (*diTauLegsLchP3AtPCAOPV_)[0].z(), TMath::Sqrt((*diTauLegsLchP3AtPCAOPV_)[0].mag2() + m_pi*m_pi) );
    LV pion2_lab_opv_pca((*diTauLegsLchP3AtPCAOPV_)[1].x(), (*diTauLegsLchP3AtPCAOPV_)[1].y(), (*diTauLegsLchP3AtPCAOPV_)[1].z(), TMath::Sqrt((*diTauLegsLchP3AtPCAOPV_)[1].mag2() + m_pi*m_pi) );
    LV pion_pair_lab_opv_pca = pion1_lab_opv_pca + pion2_lab_opv_pca;
    if(verbosity)std::cout<<" At PCA OPV, pion1 pt "<<pion1_lab_opv_pca.Pt()<<" pion2 pt "<<pion2_lab_opv_pca.Pt()<<std::endl;
    ROOT::Math::Boost boost_to_rf_opv_pca(pion_pair_lab_opv_pca.BoostToCM());

    LV pion1_rf_opv_pca = boost_to_rf_opv_pca(pion1_lab_opv_pca);
    LV pion2_rf_opv_pca = boost_to_rf_opv_pca(pion2_lab_opv_pca);

    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_opv_pca"]->Fill(TMath::ACos(ipv1_lab_opv_3d.Unit().Dot(pion1_lab_opv_pca.Vect().Unit())));
    histos_["xcheck_angle3_reco_opv_pca"]->Fill(TMath::ACos(ipv2_lab_opv_3d.Unit().Dot(pion2_lab_opv_pca.Vect().Unit())));
    Vector2D pion1_lab_opv_pca_xy(pion1_lab_opv_pca.px(), pion1_lab_opv_pca.py());
    Vector2D pion2_lab_opv_pca_xy(pion2_lab_opv_pca.px(), pion2_lab_opv_pca.py());
    histos_["xcheck_angle3_reco_opv_pca_2d"]->Fill(TMath::ACos(ipv1_lab_opv_xy.Unit().Dot(pion1_lab_opv_pca_xy.Unit())));
    histos_["xcheck_angle3_reco_opv_pca_2d"]->Fill(TMath::ACos(ipv2_lab_opv_xy.Unit().Dot(pion2_lab_opv_pca_xy.Unit())));

    //boost to ZMF
    LV ipv1_opv_rf_pca = boost_to_rf_opv_pca(ipv1_lab_opv);
    LV ipv2_opv_rf_pca = boost_to_rf_opv_pca(ipv2_lab_opv);

    //Get only the position component of the vector
    PV ipv1_opv_rf_pca_3d = ipv1_opv_rf_pca.Vect();
    PV ipv2_opv_rf_pca_3d = ipv2_opv_rf_pca.Vect();

    PV pion1_rf_opv_pca_3d = pion1_rf_opv_pca.Vect();
    PV pion2_rf_opv_pca_3d = pion2_rf_opv_pca.Vect();

    //Get the unit (normalized) vector along pion momentum
    PV pion1_rf_opv_pca_3d_u = pion1_rf_opv_pca_3d.Unit();
    PV pion2_rf_opv_pca_3d_u = pion2_rf_opv_pca_3d.Unit();

    //Get IP vector normal to pion momenta
    PV ipv1_opv_rf_pca_3d_t = ipv1_rf_3d - ipv1_opv_rf_pca_3d.Dot(pion1_rf_opv_pca_3d_u)*pion1_rf_opv_pca_3d_u;
    PV ipv2_opv_rf_pca_3d_t = ipv2_rf_3d - ipv2_opv_rf_pca_3d.Dot(pion2_rf_opv_pca_3d_u)*pion2_rf_opv_pca_3d_u;

    //Get normalized normal IP vector
    PV ipv1_opv_rf_pca_3d_t_u = ipv1_opv_rf_pca_3d_t.Unit();
    PV ipv2_opv_rf_pca_3d_t_u = ipv2_opv_rf_pca_3d_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_opv_P4AtPCA"]->Fill(TMath::ACos(ipv1_opv_rf_pca_3d_t_u.Dot(ipv2_opv_rf_pca_3d_t_u)));

    /////////////Check using PCA with 3d extrapolation method////////////////////////////
    PV ipv1_lab_opv_3d_le1 =  (*diTauLegsIPAtPCAOPV_)[0];
    PV ipv2_lab_opv_3d_le1 =  (*diTauLegsIPAtPCAOPV_)[1];

    LV ipv1_lab_opv_le1(ipv1_lab_opv_3d_le1.unit().x(), ipv1_lab_opv_3d_le1.unit().y(), ipv1_lab_opv_3d_le1.unit().z(), 0);
    LV ipv2_lab_opv_le1(ipv2_lab_opv_3d_le1.unit().x(), ipv2_lab_opv_3d_le1.unit().y(), ipv2_lab_opv_3d_le1.unit().z(), 0);

    //boost to ZMF
    LV ipv1_opv_rf_le1 = boost_to_rf(ipv1_lab_opv_le1);
    LV ipv2_opv_rf_le1 = boost_to_rf(ipv2_lab_opv_le1);

    //Get only the position component of the vector
    PV ipv1_opv_rf_3d_le1 = ipv1_opv_rf_le1.Vect();
    PV ipv2_opv_rf_3d_le1 = ipv2_opv_rf_le1.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_opv_rf_3d_le1_l = ipv1_opv_rf_3d_le1.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_opv_rf_3d_le1_l = ipv2_opv_rf_3d_le1.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_opv_rf_3d_le1_t = ipv1_opv_rf_3d_le1 - ipv1_opv_rf_3d_le1_l;
    PV ipv2_opv_rf_3d_le1_t = ipv2_opv_rf_3d_le1 - ipv2_opv_rf_3d_le1_l;

    //Get normalized normal IP vector
    PV ipv1_opv_rf_3d_le1_t_u = ipv1_opv_rf_3d_le1_t.unit();
    PV ipv2_opv_rf_3d_le1_t_u = ipv2_opv_rf_3d_le1_t.unit();
    histos_["CPPhiStar_opv_LE1"]->Fill(TMath::ACos(ipv1_opv_rf_3d_le1_t_u.Dot(ipv2_opv_rf_3d_le1_t_u)));
    histos_["CPPhiLab_opv_LE1"]->Fill(TMath::ACos(ipv1_lab_opv_3d_le1.Unit().Dot(ipv2_lab_opv_3d_le1.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v1 (OPV) "<<TMath::ACos(ipv1_opv_rf_3d_le1_t_u.Dot(ipv2_opv_rf_3d_le1_t_u))<<std::endl;
    double OstarCP_opv_le1 = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_opv_rf_3d_le1_t_u.Cross(ipv2_opv_rf_3d_le1_t_u))
      : pion2_rf_3d_u.Dot(ipv2_opv_rf_3d_le1_t_u.Cross(ipv1_opv_rf_3d_le1_t_u));
    double phi_star_opv_cp_le1 = (OstarCP_opv_le1 >= 0) ? TMath::ACos(ipv1_opv_rf_3d_le1_t_u.Dot(ipv2_opv_rf_3d_le1_t_u)) 
      : (TMath::TwoPi() - TMath::ACos(ipv1_opv_rf_3d_le1_t_u.Dot(ipv2_opv_rf_3d_le1_t_u)));
    histos_["CPPhiStar_opv_CP_LE1"]->Fill(phi_star_opv_cp_le1);

    if(sqrt(ipv1_lab_opv_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_opv_3d_le1.mag2()) > 0.005){
      histos_["CPPhiStar_opv_ipcut"]->Fill(TMath::ACos(ipv1_opv_rf_3d_le1_t_u.Dot(ipv2_opv_rf_3d_le1_t_u)));
      histos_["CPPhiLab_opv_ipcut"]->Fill(TMath::ACos(ipv1_lab_opv_3d_le1.Unit().Dot(ipv2_lab_opv_3d_le1.Unit())));
    }

    PV ipv1_lab_opv_3d_le2 =  (*diTauLegsIPAtPCAOPVV2_)[0];
    PV ipv2_lab_opv_3d_le2 =  (*diTauLegsIPAtPCAOPVV2_)[1];

    LV ipv1_lab_opv_le2(ipv1_lab_opv_3d_le2.unit().x(), ipv1_lab_opv_3d_le2.unit().y(), ipv1_lab_opv_3d_le2.unit().z(), 0);
    LV ipv2_lab_opv_le2(ipv2_lab_opv_3d_le2.unit().x(), ipv2_lab_opv_3d_le2.unit().y(), ipv2_lab_opv_3d_le2.unit().z(), 0);

    //boost to ZMF
    LV ipv1_opv_rf_le2 = boost_to_rf(ipv1_lab_opv_le2);
    LV ipv2_opv_rf_le2 = boost_to_rf(ipv2_lab_opv_le2);

    //Get only the position component of the vector
    PV ipv1_opv_rf_3d_le2 = ipv1_opv_rf_le2.Vect();
    PV ipv2_opv_rf_3d_le2 = ipv2_opv_rf_le2.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_opv_rf_3d_le2_l = ipv1_opv_rf_3d_le2.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_opv_rf_3d_le2_l = ipv2_opv_rf_3d_le2.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_opv_rf_3d_le2_t = ipv1_opv_rf_3d_le2 - ipv1_opv_rf_3d_le2_l;
    PV ipv2_opv_rf_3d_le2_t = ipv2_opv_rf_3d_le2 - ipv2_opv_rf_3d_le2_l;

    //Get normalized normal IP vector
    PV ipv1_opv_rf_3d_le2_t_u = ipv1_opv_rf_3d_le2_t.unit();
    PV ipv2_opv_rf_3d_le2_t_u = ipv2_opv_rf_3d_le2_t.unit();
    histos_["CPPhiStar_opv_LE2"]->Fill(TMath::ACos(ipv1_opv_rf_3d_le2_t_u.Dot(ipv2_opv_rf_3d_le2_t_u)));
    histos_["CPPhiLab_opv_LE2"]->Fill(TMath::ACos(ipv1_lab_opv_3d_le2.Unit().Dot(ipv2_lab_opv_3d_le2.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v2 (OPV) "<<TMath::ACos(ipv1_opv_rf_3d_le2_t_u.Dot(ipv2_opv_rf_3d_le2_t_u))<<std::endl;

    ////////////////// Compute the CP angle at Detector level, using BeamSpot///////////////////////////////////////////////
    //Define the normalized impact parameter vectors
    PV ipv1_lab_bs_3d((*diTauLegsPCABS_)[0].x() - (*BSPos_)[0].x(), (*diTauLegsPCABS_)[0].y() - (*BSPos_)[0].y(), (*diTauLegsPCABS_)[0].z() - (*BSPos_)[0].z());
    PV ipv2_lab_bs_3d((*diTauLegsPCABS_)[1].x() - (*BSPos_)[0].x(), (*diTauLegsPCABS_)[1].y() - (*BSPos_)[0].y(), (*diTauLegsPCABS_)[1].z() - (*BSPos_)[0].z());

    LV ipv1_lab_bs(ipv1_lab_bs_3d.Unit().x(), ipv1_lab_bs_3d.Unit().y(), ipv1_lab_bs_3d.Unit().z(), 0);
    LV ipv2_lab_bs(ipv2_lab_bs_3d.Unit().x(), ipv2_lab_bs_3d.Unit().y(), ipv2_lab_bs_3d.Unit().z(), 0);

    if(verbosity)std::cout<<"BeamSpot ("<<((*BSPos_)[0]).x()<<","<<((*BSPos_)[0]).y()<<","<<((*BSPos_)[0]).z()<<")"<<std::endl;
    if(verbosity)std::cout<<"PCA wrt BeamSpot pion 1 ("<<(*diTauLegsPCABS_)[0].x()<<","<<(*diTauLegsPCABS_)[0].y()<<","<<(*diTauLegsPCABS_)[0].z()<<")"<<std::endl;
    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_bs"]->Fill(TMath::ACos(ipv1_lab_bs_3d.Unit().Dot(pion1_lab.Vect().Unit())));
    histos_["xcheck_angle3_reco_bs"]->Fill(TMath::ACos(ipv2_lab_bs_3d.Unit().Dot(pion2_lab.Vect().Unit())));
    Vector2D ipv1_lab_bs_xy(ipv1_lab_bs_3d.x(), ipv1_lab_bs_3d.y());
    Vector2D ipv2_lab_bs_xy(ipv2_lab_bs_3d.x(), ipv2_lab_bs_3d.y());
    histos_["xcheck_angle3_reco_bs_2d"]->Fill(TMath::ACos(ipv1_lab_bs_xy.Unit().Dot(pion1_lab_xy.Unit())));
    histos_["xcheck_angle3_reco_bs_2d"]->Fill(TMath::ACos(ipv2_lab_bs_xy.Unit().Dot(pion2_lab_xy.Unit())));

    //boost to ZMF
    LV ipv1_bs_rf = boost_to_rf(ipv1_lab_bs);
    LV ipv2_bs_rf = boost_to_rf(ipv2_lab_bs);

    //Get only the position component of the vector
    PV ipv1_bs_rf_3d = ipv1_bs_rf.Vect();
    PV ipv2_bs_rf_3d = ipv2_bs_rf.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_bs_rf_3d_l = ipv1_bs_rf_3d.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_bs_rf_3d_l = ipv2_bs_rf_3d.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_bs_rf_3d_t = ipv1_bs_rf_3d - ipv1_bs_rf_3d_l;
    PV ipv2_bs_rf_3d_t = ipv2_bs_rf_3d - ipv2_bs_rf_3d_l;

    //Get normalized normal IP vector
    PV ipv1_bs_rf_3d_t_u = ipv1_bs_rf_3d_t/TMath::Sqrt(ipv1_bs_rf_3d_t.mag2());
    PV ipv2_bs_rf_3d_t_u = ipv2_bs_rf_3d_t/TMath::Sqrt(ipv2_bs_rf_3d_t.mag2());

    //Get CP angle in ZMF
    double phi_star_bs = TMath::ACos(ipv1_bs_rf_3d_t_u.Dot(ipv2_bs_rf_3d_t_u));
    histos_["CPPhiStar_bs"]->Fill(phi_star_bs);
    histos_["CPPhiLab_bs"]->Fill(TMath::ACos(ipv1_lab_bs_3d.Unit().Dot(ipv2_lab_bs_3d.Unit())));
    
    //replace IP vector by the component of IP normal to pion momenta
    PV ipv1_lab_bs_3d_NP = ipv1_lab_bs_3d - ipv1_lab_bs_3d.Dot(pion1_lab.Vect().Unit())*pion1_lab.Vect().Unit();
    PV ipv2_lab_bs_3d_NP = ipv2_lab_bs_3d - ipv2_lab_bs_3d.Dot(pion2_lab.Vect().Unit())*pion2_lab.Vect().Unit();

    LV ipv1_lab_bs_NP(ipv1_lab_bs_3d_NP.Unit().x(), ipv1_lab_bs_3d_NP.Unit().y(), ipv1_lab_bs_3d_NP.Unit().z(), 0);
    LV ipv2_lab_bs_NP(ipv2_lab_bs_3d_NP.Unit().x(), ipv2_lab_bs_3d_NP.Unit().y(), ipv2_lab_bs_3d_NP.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_bs_rf_NP = boost_to_rf(ipv1_lab_bs_NP);
    LV ipv2_bs_rf_NP = boost_to_rf(ipv2_lab_bs_NP);

    //Get only the position component of the vector
    PV ipv1_bs_rf_3d_NP = ipv1_bs_rf_NP.Vect();
    PV ipv2_bs_rf_3d_NP = ipv2_bs_rf_NP.Vect();

    //Get IP vector normal to pion momenta
    PV ipv1_bs_rf_3d_NP_t = ipv1_bs_rf_3d_NP - ipv1_bs_rf_3d_NP.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_bs_rf_3d_NP_t = ipv2_bs_rf_3d_NP - ipv2_bs_rf_3d_NP.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get normalized normal IP vector
    PV ipv1_bs_rf_3d_NP_t_u = ipv1_bs_rf_3d_NP_t.Unit();
    PV ipv2_bs_rf_3d_NP_t_u = ipv2_bs_rf_3d_NP_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_bs_NP"]->Fill(TMath::ACos(ipv1_bs_rf_3d_NP_t_u.Dot(ipv2_bs_rf_3d_NP_t_u)));
    //Get CP angle in lab frame
    histos_["CPPhiLab_bs_NP"]->Fill(TMath::ACos(ipv1_lab_bs_3d_NP.Unit().Dot(ipv2_lab_bs_3d_NP.Unit())));

    ///////Use momentum of charged hadron at PCA ////////////////////
    LV pion1_lab_bs_pca((*diTauLegsLchP3AtPCABS_)[0].x(), (*diTauLegsLchP3AtPCABS_)[0].y(), (*diTauLegsLchP3AtPCABS_)[0].z(), TMath::Sqrt((*diTauLegsLchP3AtPCABS_)[0].mag2() + m_pi*m_pi) );
    LV pion2_lab_bs_pca((*diTauLegsLchP3AtPCABS_)[1].x(), (*diTauLegsLchP3AtPCABS_)[1].y(), (*diTauLegsLchP3AtPCABS_)[1].z(), TMath::Sqrt((*diTauLegsLchP3AtPCABS_)[1].mag2() + m_pi*m_pi) );
    LV pion_pair_lab_bs_pca = pion1_lab_bs_pca + pion2_lab_bs_pca;
    if(verbosity)std::cout<<" At PCA BS, pion1 pt "<<pion1_lab_bs_pca.Pt()<<" pion2 pt "<<pion2_lab_bs_pca.Pt()<<std::endl;
    ROOT::Math::Boost boost_to_rf_bs_pca(pion_pair_lab_bs_pca.BoostToCM());

    LV pion1_rf_bs_pca = boost_to_rf_bs_pca(pion1_lab_bs_pca);
    LV pion2_rf_bs_pca = boost_to_rf_bs_pca(pion2_lab_bs_pca);

    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_bs_pca"]->Fill(TMath::ACos(ipv1_lab_bs_3d.Unit().Dot(pion1_lab_bs_pca.Vect().Unit())));
    histos_["xcheck_angle3_reco_bs_pca"]->Fill(TMath::ACos(ipv2_lab_bs_3d.Unit().Dot(pion2_lab_bs_pca.Vect().Unit())));
    Vector2D pion1_lab_bs_pca_xy(pion1_lab_bs_pca.px(), pion1_lab_bs_pca.py());
    Vector2D pion2_lab_bs_pca_xy(pion2_lab_bs_pca.px(), pion2_lab_bs_pca.py());
    histos_["xcheck_angle3_reco_bs_pca_2d"]->Fill(TMath::ACos(ipv1_lab_bs_xy.Unit().Dot(pion1_lab_bs_pca_xy.Unit())));
    histos_["xcheck_angle3_reco_bs_pca_2d"]->Fill(TMath::ACos(ipv2_lab_bs_xy.Unit().Dot(pion2_lab_bs_pca_xy.Unit())));

    //boost to ZMF
    LV ipv1_bs_rf_pca = boost_to_rf_bs_pca(ipv1_lab_bs);
    LV ipv2_bs_rf_pca = boost_to_rf_bs_pca(ipv2_lab_bs);

    //Get only the position component of the vector
    PV ipv1_bs_rf_pca_3d = ipv1_bs_rf_pca.Vect();
    PV ipv2_bs_rf_pca_3d = ipv2_bs_rf_pca.Vect();

    PV pion1_rf_bs_pca_3d = pion1_rf_bs_pca.Vect();
    PV pion2_rf_bs_pca_3d = pion2_rf_bs_pca.Vect();

    //Get the unit (normalized) vector along pion momentum
    PV pion1_rf_bs_pca_3d_u = pion1_rf_bs_pca_3d.Unit();
    PV pion2_rf_bs_pca_3d_u = pion2_rf_bs_pca_3d.Unit();

    //Get IP vector normal to pion momenta
    PV ipv1_bs_rf_pca_3d_t = ipv1_rf_3d - ipv1_bs_rf_pca_3d.Dot(pion1_rf_bs_pca_3d_u)*pion1_rf_bs_pca_3d_u;
    PV ipv2_bs_rf_pca_3d_t = ipv2_rf_3d - ipv2_bs_rf_pca_3d.Dot(pion2_rf_bs_pca_3d_u)*pion2_rf_bs_pca_3d_u;
    
    //Get normalized normal IP vector
    PV ipv1_bs_rf_pca_3d_t_u = ipv1_bs_rf_pca_3d_t.Unit();
    PV ipv2_bs_rf_pca_3d_t_u = ipv2_bs_rf_pca_3d_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_bs_P4AtPCA"]->Fill(TMath::ACos(ipv1_bs_rf_pca_3d_t_u.Dot(ipv2_bs_rf_pca_3d_t_u)));

    ///////////////////// Compute the CP angle at Detector level, using Gen Higgs Vertex///////////////////////////////////////////////
    //Define the normalized impact parameter vectors
    PV ipv1_lab_genvtx_3d((*diTauLegsPCAGen_)[0].x() - ((*HiggsGenVtx_)[0]).x(), (*diTauLegsPCAGen_)[0].y() - ((*HiggsGenVtx_)[0]).y(), (*diTauLegsPCAGen_)[0].z() - ((*HiggsGenVtx_)[0]).z());
    PV ipv2_lab_genvtx_3d((*diTauLegsPCAGen_)[1].x() - ((*HiggsGenVtx_)[0]).x(), (*diTauLegsPCAGen_)[1].y() - ((*HiggsGenVtx_)[0]).y(), (*diTauLegsPCAGen_)[1].z() - ((*HiggsGenVtx_)[0]).z());

    LV ipv1_lab_genvtx(ipv1_lab_genvtx_3d.Unit().x(), ipv1_lab_genvtx_3d.Unit().y(), ipv1_lab_genvtx_3d.Unit().z(), 0);
    LV ipv2_lab_genvtx(ipv2_lab_genvtx_3d.Unit().x(), ipv2_lab_genvtx_3d.Unit().y(), ipv2_lab_genvtx_3d.Unit().z(), 0);

    if(verbosity)std::cout<<"Gen Higgs Vertex ("<<((*HiggsGenVtx_)[0]).x()<<","<<((*HiggsGenVtx_)[0]).y()<<","<<((*HiggsGenVtx_)[0]).z()<<")"<<std::endl;
    if(verbosity)std::cout<<"PCA wrt Gen Vertex pion 1 ("<<(*diTauLegsPCAGen_)[0].x()<<","<<(*diTauLegsPCAGen_)[0].y()<<","<<(*diTauLegsPCAGen_)[0].z()<<")"<<std::endl;
    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_genvtx"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d.Unit().Dot(pion1_lab.Vect().Unit())));
    histos_["xcheck_angle3_reco_genvtx"]->Fill(TMath::ACos(ipv2_lab_genvtx_3d.Unit().Dot(pion2_lab.Vect().Unit())));
    Vector2D ipv1_lab_genvtx_xy(ipv1_lab_genvtx_3d.x(), ipv1_lab_genvtx_3d.y());
    Vector2D ipv2_lab_genvtx_xy(ipv2_lab_genvtx_3d.x(), ipv2_lab_genvtx_3d.y());
    histos_["xcheck_angle3_reco_genvtx_2d"]->Fill(TMath::ACos(ipv1_lab_genvtx_xy.Unit().Dot(pion1_lab_xy.Unit())));
    histos_["xcheck_angle3_reco_genvtx_2d"]->Fill(TMath::ACos(ipv2_lab_genvtx_xy.Unit().Dot(pion2_lab_xy.Unit())));

    //boost to ZMF
    LV ipv1_genvtx_rf = boost_to_rf(ipv1_lab_genvtx);
    LV ipv2_genvtx_rf = boost_to_rf(ipv2_lab_genvtx);

    //Get only the position component of the vector
    PV ipv1_genvtx_rf_3d = ipv1_genvtx_rf.Vect();
    PV ipv2_genvtx_rf_3d = ipv2_genvtx_rf.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_genvtx_rf_3d_l = ipv1_genvtx_rf_3d.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_genvtx_rf_3d_l = ipv2_genvtx_rf_3d.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_genvtx_rf_3d_t = ipv1_genvtx_rf_3d - ipv1_genvtx_rf_3d_l;
    PV ipv2_genvtx_rf_3d_t = ipv2_genvtx_rf_3d - ipv2_genvtx_rf_3d_l;

    //Get normalized normal IP vector
    PV ipv1_genvtx_rf_3d_t_u = ipv1_genvtx_rf_3d_t/TMath::Sqrt(ipv1_genvtx_rf_3d_t.mag2());
    PV ipv2_genvtx_rf_3d_t_u = ipv2_genvtx_rf_3d_t/TMath::Sqrt(ipv2_genvtx_rf_3d_t.mag2());

    //Get CP angle in ZMF
    double phi_star_genvtx = TMath::ACos(ipv1_genvtx_rf_3d_t_u.Dot(ipv2_genvtx_rf_3d_t_u));
    histos_["CPPhiStar_genvtx"]->Fill(phi_star_genvtx);
    histos_["CPPhiLab_genvtx"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d.Unit().Dot(ipv2_lab_genvtx_3d.Unit())));
    //if(sqrt(ipv1_lab_genvtx_3d.mag2()) > 0.005 && sqrt(ipv2_lab_genvtx_3d.mag2()) > 0.005){
    //  histos_["CPPhiStar_genvtx_ipcut"]->Fill(phi_star_genvtx);
    //  histos_["CPPhiLab_genvtx_ipcut"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d.Unit().Dot(ipv2_lab_genvtx_3d.Unit())));
    //}
    if(verbosity)std::cout<<"CP phi* using IP at GenVtx "<<phi_star_genvtx<<std::endl;
    double OstarCP_genvtx = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_genvtx_rf_3d_t_u.Cross(ipv2_genvtx_rf_3d_t_u))
      : pion2_rf_3d_u.Dot(ipv2_genvtx_rf_3d_t_u.Cross(ipv1_genvtx_rf_3d_t_u));
    double phi_star_genvtx_cp = (OstarCP_genvtx >= 0) ? phi_star_genvtx : (TMath::TwoPi() - phi_star_genvtx);
    histos_["CPPhiStar_genvtx_CP"]->Fill(phi_star_genvtx_cp);

    //replace IP vector by the component of IP normal to pion momenta
    PV ipv1_lab_genvtx_3d_NP = ipv1_lab_genvtx_3d - ipv1_lab_genvtx_3d.Dot(pion1_lab.Vect().Unit())*pion1_lab.Vect().Unit();
    PV ipv2_lab_genvtx_3d_NP = ipv2_lab_genvtx_3d - ipv2_lab_genvtx_3d.Dot(pion2_lab.Vect().Unit())*pion2_lab.Vect().Unit();

    LV ipv1_lab_genvtx_NP(ipv1_lab_genvtx_3d_NP.Unit().x(), ipv1_lab_genvtx_3d_NP.Unit().y(), ipv1_lab_genvtx_3d_NP.Unit().z(), 0);
    LV ipv2_lab_genvtx_NP(ipv2_lab_genvtx_3d_NP.Unit().x(), ipv2_lab_genvtx_3d_NP.Unit().y(), ipv2_lab_genvtx_3d_NP.Unit().z(), 0);

    //boost to ZMF
    LV ipv1_genvtx_rf_NP = boost_to_rf(ipv1_lab_genvtx_NP);
    LV ipv2_genvtx_rf_NP = boost_to_rf(ipv2_lab_genvtx_NP);

    //Get only the position component of the vector
    PV ipv1_genvtx_rf_3d_NP = ipv1_genvtx_rf_NP.Vect();
    PV ipv2_genvtx_rf_3d_NP = ipv2_genvtx_rf_NP.Vect();

    //Get IP vector normal to pion momenta
    PV ipv1_genvtx_rf_3d_NP_t = ipv1_genvtx_rf_3d_NP - ipv1_genvtx_rf_3d_NP.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_genvtx_rf_3d_NP_t = ipv2_genvtx_rf_3d_NP - ipv2_genvtx_rf_3d_NP.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get normalized normal IP vector
    PV ipv1_genvtx_rf_3d_NP_t_u = ipv1_genvtx_rf_3d_NP_t.Unit();
    PV ipv2_genvtx_rf_3d_NP_t_u = ipv2_genvtx_rf_3d_NP_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_genvtx_NP"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_NP_t_u.Dot(ipv2_genvtx_rf_3d_NP_t_u)));
    //Get CP angle in lab frame
    histos_["CPPhiLab_genvtx_NP"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_NP.Unit().Dot(ipv2_lab_genvtx_3d_NP.Unit())));
    if(sqrt(ipv1_lab_genvtx_3d.mag2()) > 0.005 && sqrt(ipv2_lab_genvtx_3d.mag2()) > 0.005){
      histos_["CPPhiStar_genvtx_NP_ipcut"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_NP_t_u.Dot(ipv2_genvtx_rf_3d_NP_t_u)));
      histos_["CPPhiLab_genvtx_NP_ipcut"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_NP.Unit().Dot(ipv2_lab_genvtx_3d_NP.Unit())));
    }
    if(verbosity)std::cout<<"CP phi* using normal of IP "<<TMath::ACos(ipv1_genvtx_rf_3d_NP_t_u.Dot(ipv2_genvtx_rf_3d_NP_t_u))<<std::endl;

    ///////Use momentum of charged hadron at PCA ////////////////////
    LV pion1_lab_genvtx_pca((*diTauLegsLchP3AtPCAGen_)[0].x(), (*diTauLegsLchP3AtPCAGen_)[0].y(), (*diTauLegsLchP3AtPCAGen_)[0].z(), TMath::Sqrt((*diTauLegsLchP3AtPCAGen_)[0].mag2() + m_pi*m_pi) );
    LV pion2_lab_genvtx_pca((*diTauLegsLchP3AtPCAGen_)[1].x(), (*diTauLegsLchP3AtPCAGen_)[1].y(), (*diTauLegsLchP3AtPCAGen_)[1].z(), TMath::Sqrt((*diTauLegsLchP3AtPCAGen_)[1].mag2() + m_pi*m_pi) );
    LV pion_pair_lab_genvtx_pca = pion1_lab_genvtx_pca + pion2_lab_genvtx_pca;
    if(verbosity)std::cout<<" At PCA GevVtx, pion1 pt "<<pion1_lab_genvtx_pca.Pt()<<" pion2 pt "<<pion2_lab_genvtx_pca.Pt()<<std::endl;
    ROOT::Math::Boost boost_to_rf_genvtx_pca(pion_pair_lab_genvtx_pca.BoostToCM());

    LV pion1_rf_genvtx_pca = boost_to_rf_genvtx_pca(pion1_lab_genvtx_pca);
    LV pion2_rf_genvtx_pca = boost_to_rf_genvtx_pca(pion2_lab_genvtx_pca);

    //Check the angle between pion momenta and IP vector
    histos_["xcheck_angle3_reco_genvtx_pca"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d.Unit().Dot(pion1_lab_genvtx_pca.Vect().Unit())));
    histos_["xcheck_angle3_reco_genvtx_pca"]->Fill(TMath::ACos(ipv2_lab_genvtx_3d.Unit().Dot(pion2_lab_genvtx_pca.Vect().Unit())));
    Vector2D pion1_lab_genvtx_pca_xy(pion1_lab_genvtx_pca.px(), pion1_lab_genvtx_pca.py());
    Vector2D pion2_lab_genvtx_pca_xy(pion2_lab_genvtx_pca.px(), pion2_lab_genvtx_pca.py());
    histos_["xcheck_angle3_reco_genvtx_pca_2d"]->Fill(TMath::ACos(ipv1_lab_genvtx_xy.Unit().Dot(pion1_lab_genvtx_pca_xy.Unit())));
    histos_["xcheck_angle3_reco_genvtx_pca_2d"]->Fill(TMath::ACos(ipv2_lab_genvtx_xy.Unit().Dot(pion2_lab_genvtx_pca_xy.Unit())));

    if(verbosity)std::cout<<"angle IP and Pion 1 "<<TMath::ACos(((*diTauLegsIPAtPCAGen_)[0]).unit().Dot(pion1_lab_genvtx_pca.Vect().Unit()))<<std::endl;
    if(verbosity)std::cout<<"angle IP and Pion 2 "<<TMath::ACos(((*diTauLegsIPAtPCAGen_)[1]).unit().Dot(pion2_lab_genvtx_pca.Vect().Unit()))<<std::endl;
    if(verbosity)std::cout<<"angle IP and Pion 1 V2 "<<TMath::ACos(((*diTauLegsIPAtPCAGenV2_)[0]).unit().Dot(pion1_lab_genvtx_pca.Vect().Unit()))<<std::endl;
    if(verbosity)std::cout<<"angle IP and Pion 2 V2 "<<TMath::ACos(((*diTauLegsIPAtPCAGenV2_)[1]).unit().Dot(pion2_lab_genvtx_pca.Vect().Unit()))<<std::endl;

    //boost to ZMF
    LV ipv1_genvtx_rf_pca = boost_to_rf_genvtx_pca(ipv1_lab_genvtx);
    LV ipv2_genvtx_rf_pca = boost_to_rf_genvtx_pca(ipv2_lab_genvtx);

    //Get only the position component of the vector
    PV ipv1_genvtx_rf_pca_3d = ipv1_genvtx_rf_pca.Vect();
    PV ipv2_genvtx_rf_pca_3d = ipv2_genvtx_rf_pca.Vect();

    PV pion1_rf_genvtx_pca_3d = pion1_rf_genvtx_pca.Vect();
    PV pion2_rf_genvtx_pca_3d = pion2_rf_genvtx_pca.Vect();

    //Get the unit (normalized) vector along pion momentum
    PV pion1_rf_genvtx_pca_3d_u = pion1_rf_genvtx_pca_3d.Unit();
    PV pion2_rf_genvtx_pca_3d_u = pion2_rf_genvtx_pca_3d.Unit();

    //Get IP vector normal to pion momenta
    PV ipv1_genvtx_rf_pca_3d_t = ipv1_rf_3d - ipv1_genvtx_rf_pca_3d.Dot(pion1_rf_genvtx_pca_3d_u)*pion1_rf_genvtx_pca_3d_u;
    PV ipv2_genvtx_rf_pca_3d_t = ipv2_rf_3d - ipv2_genvtx_rf_pca_3d.Dot(pion2_rf_genvtx_pca_3d_u)*pion2_rf_genvtx_pca_3d_u;

    //Get normalized normal IP vector
    PV ipv1_genvtx_rf_pca_3d_t_u = ipv1_genvtx_rf_pca_3d_t.Unit();
    PV ipv2_genvtx_rf_pca_3d_t_u = ipv2_genvtx_rf_pca_3d_t.Unit();

    //Get CP angle in ZMF
    histos_["CPPhiStar_genvtx_P4AtPCA"]->Fill(TMath::ACos(ipv1_genvtx_rf_pca_3d_t_u.Dot(ipv2_genvtx_rf_pca_3d_t_u)));

    /////////////Check using PCA with 3d extrapolation method////////////////////////////
    PV ipv1_lab_genvtx_3d_le1 =  (*diTauLegsIPAtPCAGen_)[0];
    PV ipv2_lab_genvtx_3d_le1 =  (*diTauLegsIPAtPCAGen_)[1];
    
    LV ipv1_lab_genvtx_le1(ipv1_lab_genvtx_3d_le1.unit().x(), ipv1_lab_genvtx_3d_le1.unit().y(), ipv1_lab_genvtx_3d_le1.unit().z(), 0);
    LV ipv2_lab_genvtx_le1(ipv2_lab_genvtx_3d_le1.unit().x(), ipv2_lab_genvtx_3d_le1.unit().y(), ipv2_lab_genvtx_3d_le1.unit().z(), 0);

    //boost to ZMF
    LV ipv1_genvtx_rf_le1 = boost_to_rf(ipv1_lab_genvtx_le1);
    LV ipv2_genvtx_rf_le1 = boost_to_rf(ipv2_lab_genvtx_le1);

    //Get only the position component of the vector
    PV ipv1_genvtx_rf_3d_le1 = ipv1_genvtx_rf_le1.Vect();
    PV ipv2_genvtx_rf_3d_le1 = ipv2_genvtx_rf_le1.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_genvtx_rf_3d_le1_l = ipv1_genvtx_rf_3d_le1.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_genvtx_rf_3d_le1_l = ipv2_genvtx_rf_3d_le1.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_genvtx_rf_3d_le1_t = ipv1_genvtx_rf_3d_le1 - ipv1_genvtx_rf_3d_le1_l;
    PV ipv2_genvtx_rf_3d_le1_t = ipv2_genvtx_rf_3d_le1 - ipv2_genvtx_rf_3d_le1_l;

    //Get normalized normal IP vector
    PV ipv1_genvtx_rf_3d_le1_t_u = ipv1_genvtx_rf_3d_le1_t.unit();
    PV ipv2_genvtx_rf_3d_le1_t_u = ipv2_genvtx_rf_3d_le1_t.unit();
    histos_["CPPhiStar_genvtx_LE1"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
    histos_["CPPhiLab_genvtx_LE1"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v1 "<<TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u))<<std::endl;
    double OstarCP_genvtx_le1 = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_genvtx_rf_3d_le1_t_u.Cross(ipv2_genvtx_rf_3d_le1_t_u))
      : pion2_rf_3d_u.Dot(ipv2_genvtx_rf_3d_le1_t_u.Cross(ipv1_genvtx_rf_3d_le1_t_u));
    double phi_star_genvtx_cp_le1 = (OstarCP_genvtx_le1 >= 0) ? TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u))
      : (TMath::TwoPi() - TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
    histos_["CPPhiStar_genvtx_CP_LE1"]->Fill(phi_star_genvtx_cp_le1);
    if(sqrt(ipv1_lab_genvtx_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_genvtx_3d_le1.mag2()) > 0.005){
      histos_["CPPhiStar_genvtx_ipcut"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
      histos_["CPPhiLab_genvtx_ipcut"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
      if(passQualityLeg1_ && passQualityLeg2_){
	histos_["CPPhiStar_genvtx_ipcut_TkQCut"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
        histos_["CPPhiLab_genvtx_ipcut_TkQCut"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
      }
      if(passQualityLeg1V2_ && passQualityLeg2V2_){
        histos_["CPPhiStar_genvtx_ipcut_TkQCutV2"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
        //histos_["CPPhiLab_genvtx_ipcut_TkQCutV2"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
      }
      if(passQualityLeg1V3_ && passQualityLeg2V3_)
	histos_["CPPhiStar_genvtx_ipcut_TkQCutV3"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
    }
    if(sqrt(ipv1_lab_genvtx_3d_le1.mag2()) > 0.003 && sqrt(ipv2_lab_genvtx_3d_le1.mag2()) > 0.003){
      histos_["CPPhiStar_genvtx_ipcut30"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
      histos_["CPPhiLab_genvtx_ipcut30"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
      if(passQualityLeg1_ && passQualityLeg2_){
        histos_["CPPhiStar_genvtx_ipcut30_TkQCut"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
        histos_["CPPhiLab_genvtx_ipcut30_TkQCut"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
      }
      if(passQualityLeg1V2_ && passQualityLeg2V2_){
        histos_["CPPhiStar_genvtx_ipcut30_TkQCutV2"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
        //histos_["CPPhiLab_genvtx_ipcut30_TkQCutV2"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le1.Unit().Dot(ipv2_lab_genvtx_3d_le1.Unit())));
      }
    }
    if(ptL1 > 40 && ptL2 > 40)
      histos_["CPPhiStar_genvtx_tau40"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
    if(ptL1 > 30 && ptL2 > 30)
      histos_["CPPhiStar_genvtx_tau30"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
    if(genDecayModeLeg1_ == 0 && genDecayModeLeg2_ == 0)
      histos_["CPPhiStar_genvtx_gen1p0pi0"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
    if(genDecayModeLeg1_ >= 0 && genDecayModeLeg2_ >= 0)
      histos_["CPPhiStar_genvtx_gentau"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));

    PV ipv1_lab_genvtx_3d_le2 =  (*diTauLegsIPAtPCAGenV2_)[0];
    PV ipv2_lab_genvtx_3d_le2 =  (*diTauLegsIPAtPCAGenV2_)[1];

    LV ipv1_lab_genvtx_le2(ipv1_lab_genvtx_3d_le2.unit().x(), ipv1_lab_genvtx_3d_le2.unit().y(), ipv1_lab_genvtx_3d_le2.unit().z(), 0);
    LV ipv2_lab_genvtx_le2(ipv2_lab_genvtx_3d_le2.unit().x(), ipv2_lab_genvtx_3d_le2.unit().y(), ipv2_lab_genvtx_3d_le2.unit().z(), 0);

    //boost to ZMF
    LV ipv1_genvtx_rf_le2 = boost_to_rf(ipv1_lab_genvtx_le2);
    LV ipv2_genvtx_rf_le2 = boost_to_rf(ipv2_lab_genvtx_le2);

    //Get only the position component of the vector
    PV ipv1_genvtx_rf_3d_le2 = ipv1_genvtx_rf_le2.Vect();
    PV ipv2_genvtx_rf_3d_le2 = ipv2_genvtx_rf_le2.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_genvtx_rf_3d_le2_l = ipv1_genvtx_rf_3d_le2.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_genvtx_rf_3d_le2_l = ipv2_genvtx_rf_3d_le2.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_genvtx_rf_3d_le2_t = ipv1_genvtx_rf_3d_le2 - ipv1_genvtx_rf_3d_le2_l;
    PV ipv2_genvtx_rf_3d_le2_t = ipv2_genvtx_rf_3d_le2 - ipv2_genvtx_rf_3d_le2_l;

    //Get normalized normal IP vector
    PV ipv1_genvtx_rf_3d_le2_t_u = ipv1_genvtx_rf_3d_le2_t.unit();
    PV ipv2_genvtx_rf_3d_le2_t_u = ipv2_genvtx_rf_3d_le2_t.unit();
    histos_["CPPhiStar_genvtx_LE2"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le2_t_u.Dot(ipv2_genvtx_rf_3d_le2_t_u)));
    histos_["CPPhiLab_genvtx_LE2"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d_le2.Unit().Dot(ipv2_lab_genvtx_3d_le2.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v2 "<<TMath::ACos(ipv1_genvtx_rf_3d_le2_t_u.Dot(ipv2_genvtx_rf_3d_le2_t_u))<<std::endl;

    //Vertex associated to tau , using Ian's method
    PV ipv1_lab_tauvtx_3d_le1 =  (*diTauLegsIPAtPCATauVtx_)[0];
    PV ipv2_lab_tauvtx_3d_le1 =  (*diTauLegsIPAtPCATauVtx_)[1];

    LV ipv1_lab_tauvtx_le1(ipv1_lab_tauvtx_3d_le1.unit().x(), ipv1_lab_tauvtx_3d_le1.unit().y(), ipv1_lab_tauvtx_3d_le1.unit().z(), 0);
    LV ipv2_lab_tauvtx_le1(ipv2_lab_tauvtx_3d_le1.unit().x(), ipv2_lab_tauvtx_3d_le1.unit().y(), ipv2_lab_tauvtx_3d_le1.unit().z(), 0);

    //boost to ZMF
    LV ipv1_tauvtx_rf_le1 = boost_to_rf(ipv1_lab_tauvtx_le1);
    LV ipv2_tauvtx_rf_le1 = boost_to_rf(ipv2_lab_tauvtx_le1);

    //Get only the position component of the vector
    PV ipv1_tauvtx_rf_3d_le1 = ipv1_tauvtx_rf_le1.Vect();
    PV ipv2_tauvtx_rf_3d_le1 = ipv2_tauvtx_rf_le1.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_tauvtx_rf_3d_le1_l = ipv1_tauvtx_rf_3d_le1.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_tauvtx_rf_3d_le1_l = ipv2_tauvtx_rf_3d_le1.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;
    
    //Get IP vector normal to pion momenta
    PV ipv1_tauvtx_rf_3d_le1_t = ipv1_tauvtx_rf_3d_le1 - ipv1_tauvtx_rf_3d_le1_l;
    PV ipv2_tauvtx_rf_3d_le1_t = ipv2_tauvtx_rf_3d_le1 - ipv2_tauvtx_rf_3d_le1_l;

    //Get normalized normal IP vector
    PV ipv1_tauvtx_rf_3d_le1_t_u = ipv1_tauvtx_rf_3d_le1_t.unit();
    PV ipv2_tauvtx_rf_3d_le1_t_u = ipv2_tauvtx_rf_3d_le1_t.unit();
    histos_["CPPhiStar_tauvtx_LE1"]->Fill(TMath::ACos(ipv1_tauvtx_rf_3d_le1_t_u.Dot(ipv2_tauvtx_rf_3d_le1_t_u)));
    histos_["CPPhiLab_tauvtx_LE1"]->Fill(TMath::ACos(ipv1_lab_tauvtx_3d_le1.Unit().Dot(ipv2_lab_tauvtx_3d_le1.Unit())));
    if(verbosity)std::cout<<"CP phi* using new IP v1 "<<TMath::ACos(ipv1_tauvtx_rf_3d_le1_t_u.Dot(ipv2_tauvtx_rf_3d_le1_t_u))<<std::endl;
    double OstarCP_tauvtx_le1 = (chargeLeg1_ < 0) ? pion1_rf_3d_u.Dot(ipv1_tauvtx_rf_3d_le1_t_u.Cross(ipv2_tauvtx_rf_3d_le1_t_u))
      : pion2_rf_3d_u.Dot(ipv2_tauvtx_rf_3d_le1_t_u.Cross(ipv1_tauvtx_rf_3d_le1_t_u));
    double phi_star_tauvtx_cp_le1 = (OstarCP_tauvtx_le1 >= 0) ? TMath::ACos(ipv1_tauvtx_rf_3d_le1_t_u.Dot(ipv2_tauvtx_rf_3d_le1_t_u))
      : (TMath::TwoPi() - TMath::ACos(ipv1_tauvtx_rf_3d_le1_t_u.Dot(ipv2_tauvtx_rf_3d_le1_t_u)));
    histos_["CPPhiStar_tauvtx_CP_LE1"]->Fill(phi_star_tauvtx_cp_le1);
    

    //Compare to genLevel info
    if(genVPid_ == HPid_){

      LV gentauP_lab(0, 0, 0, 0); LV gentauN_lab(0, 0, 0, 0);
      if(genTausP4_->size() == 2 && genTausCharge_->size() == 2){
	if((*genTausCharge_)[0] > 0 && (*genTausCharge_)[1] < 0){
	  gentauP_lab = (*genTausP4_)[0];
	  gentauN_lab = (*genTausP4_)[1];
	}
	else if((*genTausCharge_)[0] < 0 && (*genTausCharge_)[1] > 0){
	  gentauP_lab = (*genTausP4_)[1];
	  gentauN_lab = (*genTausP4_)[0];
	}
      }

      LV genpionP_lab(0, 0, 0, 0); LV genpionN_lab(0, 0, 0, 0);
      Point3D genpionPVtx_(0, 0, 0); Point3D genpionNVtx_(0, 0, 0);
      if(genTauPSonsP4_->size()>0 && genTauPSonsP4_->size()<3 && genTauPSonsPid_->size() > 0){ //take only one prong
	for(size_t i = 0; i<genTauPSonsP4_->size(); i++){
	  if(abs((*genTauPSonsPid_)[i]) == 211){
	    genpionP_lab = (*genTauPSonsP4_)[i];
	    genpionPVtx_ = (*TauPSonsGenVtx_)[i];
	  }
	}
      }
      if(genTauNSonsP4_->size()>0 && genTauNSonsP4_->size()<3 && genTauNSonsPid_->size() > 0){ //take only one prong
	for(size_t i = 0; i<genTauNSonsP4_->size(); i++){
	  if(abs((*genTauNSonsPid_)[i]) == 211){
	    genpionN_lab = (*genTauNSonsP4_)[i];
	    genpionNVtx_ = (*TauNSonsGenVtx_)[i];
	  }
	}
      }

      bool matchPion1 = false; bool matchPion2 = false;
      if(genpionP_lab.pt() > 0 && genpionN_lab.pt() > 0 && gentauP_lab.pt() > 0 && gentauN_lab.pt() > 0){
	if(ROOT::Math::VectorUtil::DeltaR(genpionP_lab, pion1_lab) < 0.3 || ROOT::Math::VectorUtil::DeltaR(genpionP_lab, pion2_lab) < 0.3 )matchPion1 = true;
	if(ROOT::Math::VectorUtil::DeltaR(genpionN_lab, pion1_lab) < 0.3 || ROOT::Math::VectorUtil::DeltaR(genpionN_lab, pion2_lab) < 0.3 )matchPion2 = true;
	if(matchPion1 && matchPion2){
	  //Fill the reco distributions again
	  histos_["CPPhiStar_genMatch"]->Fill(phi_star);
	  histos_["CPPhiLab_genMatch"]->Fill(phi_lab);
	  histos_["CPPhiStar_opv_genMatch"]->Fill(phi_star_opv);
	  histos_["CPPhiStar_bs_genMatch"]->Fill(phi_star_bs);
	  histos_["CPPhiStar_genvtx_genMatch"]->Fill(phi_star_genvtx);
	  histos_["CPPhiStar_NP_genMatch"]->Fill(TMath::ACos(ipv1_rf_3d_NP_t_u.Dot(ipv2_rf_3d_NP_t_u)));
	  histos_["CPPhiLab_NP_genMatch"]->Fill(TMath::ACos(ipv1_lab_3d_NP.Unit().Dot(ipv2_lab_3d_NP.Unit())));
	  histos_["CPPhiStar_opv_NP_genMatch"]->Fill(TMath::ACos(ipv1_opv_rf_3d_NP_t_u.Dot(ipv2_opv_rf_3d_NP_t_u)));
	  histos_["CPPhiStar_bs_NP_genMatch"]->Fill(TMath::ACos(ipv1_bs_rf_3d_NP_t_u.Dot(ipv2_bs_rf_3d_NP_t_u)));
	  histos_["CPPhiStar_genvtx_NP_genMatch"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_NP_t_u.Dot(ipv2_genvtx_rf_3d_NP_t_u)));
	  histos_["CPPhiStar_P4AtPCA_genMatch"]->Fill(TMath::ACos(ipv1_rf_pca_3d_t_u.Dot(ipv2_rf_pca_3d_t_u)));
	  histos_["CPPhiStar_opv_P4AtPCA_genMatch"]->Fill(TMath::ACos(ipv1_opv_rf_pca_3d_t_u.Dot(ipv2_opv_rf_pca_3d_t_u)));
	  histos_["CPPhiStar_bs_P4AtPCA_genMatch"]->Fill(TMath::ACos(ipv1_bs_rf_pca_3d_t_u.Dot(ipv2_bs_rf_pca_3d_t_u)));
	  histos_["CPPhiStar_genvtx_P4AtPCA_genMatch"]->Fill(TMath::ACos(ipv1_genvtx_rf_pca_3d_t_u.Dot(ipv2_genvtx_rf_pca_3d_t_u)));

	  histos_["CPPhiStar_LE1_genMatch"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
          histos_["CPPhiStar_LE2_genMatch"]->Fill(TMath::ACos(ipv1_rf_3d_le2_t_u.Dot(ipv2_rf_3d_le2_t_u)));
	  histos_["CPPhiStar_opv_LE1_genMatch"]->Fill(TMath::ACos(ipv1_opv_rf_3d_le1_t_u.Dot(ipv2_opv_rf_3d_le1_t_u)));
          histos_["CPPhiStar_opv_LE2_genMatch"]->Fill(TMath::ACos(ipv1_opv_rf_3d_le2_t_u.Dot(ipv2_opv_rf_3d_le2_t_u)));
	  histos_["CPPhiStar_genvtx_LE1_genMatch"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));
	  histos_["CPPhiStar_genvtx_LE2_genMatch"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le2_t_u.Dot(ipv2_genvtx_rf_3d_le2_t_u)));

	  //if(verbosity)std::cout<<"pT of Lch in lab "<<pion1_lab.pt()<<" "<<pion2_lab.pt()<<std::endl;
	  ////if(verbosity)std::cout<<"pT of Lch in rf  "<<pion1_rf.pt()<<" "<<pion2_rf.pt()<<std::endl;
	  //if(verbosity)std::cout<<"pT of gen pion in lab "<<genpionP_lab.pt()<<" "<<genpionN_lab.pt()<<std::endl;
	  //if(verbosity)std::cout<<"Gen Higgs Vertex ("<<((*HiggsGenVtx_)[0]).x()<<","<<((*HiggsGenVtx_)[0]).y()<<","<<((*HiggsGenVtx_)[0]).z()<<")"<<std::endl;
	  //if(verbosity)std::cout<<"Reco P. Vertex ("<<VtxX_<<","<<VtxY_<<","<<VtxZ_<<")"<<std::endl;
	  //if(verbosity)std::cout<<"Reco Refit Vertex ("<<ReFitVtxX_<<","<<ReFitVtxY_<<","<<ReFitVtxZ_<<")"<<std::endl;

	  //get angle between generated and reconstructed decay plane 
	  PV genpionP_lab_3d_u = genpionP_lab.Vect().Unit();
	  PV genpionN_lab_3d_u = genpionN_lab.Vect().Unit();

	  Point3D GenHVtx_ = (*HiggsGenVtx_)[0];
	  //Get PCAs, by extrapolating the line represented by pion vertex and pion momentum
	  double tP = genpionP_lab_3d_u.x()*(GenHVtx_.x() - genpionPVtx_.x()) + genpionP_lab_3d_u.y()*(GenHVtx_.y() - genpionPVtx_.y()) + genpionP_lab_3d_u.z()*(GenHVtx_.z() - genpionPVtx_.z());
	  double tN = genpionN_lab_3d_u.x()*(GenHVtx_.x() - genpionNVtx_.x()) + genpionN_lab_3d_u.y()*(GenHVtx_.y() - genpionNVtx_.y()) + genpionN_lab_3d_u.z()*(GenHVtx_.z() - genpionNVtx_.z());
              
	  Point3D pcaP(genpionPVtx_.x() + genpionP_lab_3d_u.x()*tP, genpionPVtx_.y() + genpionP_lab_3d_u.y()*tP, genpionPVtx_.z() + genpionP_lab_3d_u.z()*tP);
	  Point3D pcaN(genpionNVtx_.x() + genpionN_lab_3d_u.x()*tN, genpionNVtx_.y() + genpionN_lab_3d_u.y()*tN, genpionNVtx_.z() + genpionN_lab_3d_u.z()*tN);

	  //Get the IP vectors
	  PV ipvP_lab_gen_3d(pcaP.x() - GenHVtx_.x(), pcaP.y() - GenHVtx_.y(), pcaP.z() - GenHVtx_.z());
	  PV ipvN_lab_gen_3d(pcaN.x() - GenHVtx_.x(), pcaN.y() - GenHVtx_.y(), pcaN.z() - GenHVtx_.z());
	  
	  //get normal to decay planes in generator
	  PV normal_decayPlane_P = genpionP_lab_3d_u.Cross(ipvP_lab_gen_3d);
	  PV normal_decayPlane_N = genpionN_lab_3d_u.Cross(ipvN_lab_gen_3d);

	  //get normal to decay planes in detector
	  PV normal_decayPlane_p1 = pion1_lab.Vect().Unit().Cross(ipv1_lab_3d);
	  PV normal_decayPlane_p2 = pion2_lab.Vect().Unit().Cross(ipv2_lab_3d);

	  PV normal_decayPlane_p1_opv = pion1_lab.Vect().Unit().Cross(ipv1_lab_opv_3d);
          PV normal_decayPlane_p2_opv = pion2_lab.Vect().Unit().Cross(ipv2_lab_opv_3d);

	  PV normal_decayPlane_p1_genvtx = pion1_lab.Vect().Unit().Cross(ipv1_lab_genvtx_3d);
          PV normal_decayPlane_p2_genvtx = pion2_lab.Vect().Unit().Cross(ipv2_lab_genvtx_3d);

	  PV normal_decayPlane_p1_le1 = pion1_lab.Vect().Unit().Cross(ipv1_lab_3d_le1);
          PV normal_decayPlane_p2_le1 = pion2_lab.Vect().Unit().Cross(ipv2_lab_3d_le1);

          PV normal_decayPlane_p1_opv_le1 = pion1_lab.Vect().Unit().Cross(ipv1_lab_opv_3d_le1);
          PV normal_decayPlane_p2_opv_le1 = pion2_lab.Vect().Unit().Cross(ipv2_lab_opv_3d_le1);

          PV normal_decayPlane_p1_genvtx_le1 = pion1_lab.Vect().Unit().Cross(ipv1_lab_genvtx_3d_le1);
          PV normal_decayPlane_p2_genvtx_le1 = pion2_lab.Vect().Unit().Cross(ipv2_lab_genvtx_3d_le1);

	  //PCA - genTauDecayVertex projected parallel to gen. tau lepton momentum
	  //PCA - genTauDecayVertex projected perpendicular to gen. tau lepton momentum
	  //At Gen level
	  PV pcaTogenTauDecayVtx_gen_P(pcaP.x() - genpionPVtx_.x(), pcaP.y() - genpionPVtx_.y(), pcaP.z() - genpionPVtx_.z());
	  PV pcaTogenTauDecayVtx_gen_N(pcaN.x() - genpionNVtx_.x(), pcaN.y() - genpionNVtx_.y(), pcaN.z() - genpionNVtx_.z());
	  PV pcaTogenTauDecayVtx_gen_P_l = pcaTogenTauDecayVtx_gen_P.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
	  PV pcaTogenTauDecayVtx_gen_P_t = pcaTogenTauDecayVtx_gen_P - pcaTogenTauDecayVtx_gen_P_l;
	  PV pcaTogenTauDecayVtx_gen_N_l = pcaTogenTauDecayVtx_gen_N.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
	  PV pcaTogenTauDecayVtx_gen_N_t = pcaTogenTauDecayVtx_gen_N - pcaTogenTauDecayVtx_gen_N_l;

	  if(ROOT::Math::VectorUtil::DeltaR(genpionP_lab, pion1_lab) < 0.3 && ROOT::Math::VectorUtil::DeltaR(genpionN_lab, pion2_lab) < 0.3 ){
	    

	    histos_["DeltaX_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[0].x() - pcaP.x());
	    histos_["DeltaX_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[1].x() - pcaN.x());

	    histos_["DeltaY_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[0].y() - pcaP.y());
            histos_["DeltaY_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[1].y() - pcaN.y());

	    histos_["DeltaZ_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[0].z() - pcaP.z());
            histos_["DeltaZ_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[1].z() - pcaN.z());

	    histos_["DeltaX_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.x()+ipv1_lab_genvtx_3d_le1.x()) - pcaP.x());
	    histos_["DeltaX_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.x()+ipv2_lab_genvtx_3d_le1.x()) - pcaN.x());
	    histos_["DeltaY_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.y()+ipv1_lab_genvtx_3d_le1.y()) - pcaP.y());
            histos_["DeltaY_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.y()+ipv2_lab_genvtx_3d_le1.y()) - pcaN.y());
	    histos_["DeltaZ_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.z()+ipv1_lab_genvtx_3d_le1.z()) - pcaP.z());
            histos_["DeltaZ_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.z()+ipv2_lab_genvtx_3d_le1.z()) - pcaN.z());
	    
	    histos_["DeltaX_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.x()+ipv1_lab_genvtx_3d_le2.x()) - pcaP.x());
            histos_["DeltaX_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.x()+ipv2_lab_genvtx_3d_le2.x()) - pcaN.x());
            histos_["DeltaY_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.y()+ipv1_lab_genvtx_3d_le2.y()) - pcaP.y());
            histos_["DeltaY_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.y()+ipv2_lab_genvtx_3d_le2.y()) - pcaN.y());
            histos_["DeltaZ_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.z()+ipv1_lab_genvtx_3d_le2.z()) - pcaP.z());
            histos_["DeltaZ_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.z()+ipv2_lab_genvtx_3d_le2.z()) - pcaN.z());

	    histos_["angle_dp_refitVtx"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1.Unit())));
	    histos_["angle_dp_refitVtx"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2.Unit())));

	    histos_["angle_dp_origVtx"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_opv.Unit())));
            histos_["angle_dp_origVtx"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_opv.Unit())));

	    histos_["angle_dp_genVtx"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_genvtx.Unit())));
            histos_["angle_dp_genVtx"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_genvtx.Unit())));

	    histos_["delphi_dp_genVtx"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(normal_decayPlane_P.Unit(), normal_decayPlane_p1_genvtx.Unit()));
	    histos_["delphi_dp_genVtx"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(normal_decayPlane_N.Unit(), normal_decayPlane_p2_genvtx.Unit()));

	    histos_["angle_dp_refitVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_le1.Unit())));
            histos_["angle_dp_refitVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_le1.Unit())));
            histos_["angle_dp_origVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_opv_le1.Unit())));
            histos_["angle_dp_origVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_opv_le1.Unit())));
            histos_["angle_dp_genVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_genvtx_le1.Unit())));
            histos_["angle_dp_genVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_genvtx_le1.Unit())));

	    //get normal to decay planes in detector, replacing detector Pi p4 by gen Pi p4
	    PV normal_decayPlane_p1_genP4 = genpionP_lab_3d_u.Cross(ipv1_lab_3d);
	    PV normal_decayPlane_p2_genP4 = genpionN_lab_3d_u.Cross(ipv2_lab_3d);

	    PV normal_decayPlane_p1_opv_genP4 = genpionP_lab_3d_u.Cross(ipv1_lab_opv_3d);
	    PV normal_decayPlane_p2_opv_genP4 = genpionN_lab_3d_u.Cross(ipv2_lab_opv_3d);

	    PV normal_decayPlane_p1_genvtx_genP4 = genpionP_lab_3d_u.Cross(ipv1_lab_genvtx_3d);
	    PV normal_decayPlane_p2_genvtx_genP4 = genpionN_lab_3d_u.Cross(ipv2_lab_genvtx_3d);

	    histos_["angle_dp_refitVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_genP4.Unit())));
            histos_["angle_dp_refitVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_genP4.Unit())));

            histos_["angle_dp_origVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_opv_genP4.Unit())));
            histos_["angle_dp_origVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_opv_genP4.Unit())));

            histos_["angle_dp_genVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_genvtx_genP4.Unit())));
            histos_["angle_dp_genVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_genvtx_genP4.Unit())));

	    //get normal to decay planes in detector, replacing detector Pi IP by gen Pi IP
	    PV normal_decayPlane_p1_genIP = pion1_lab.Vect().Unit().Cross(ipvP_lab_gen_3d);
	    PV normal_decayPlane_p2_genIP = pion2_lab.Vect().Unit().Cross(ipvN_lab_gen_3d);

	    histos_["angle_dp_genIP"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_genIP.Unit())));
            histos_["angle_dp_genIP"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_genIP.Unit())));

	    //get normal to decay planes in detector, replacing detector Pi IP by normal component of IP wrt pion momenta
	    PV normal_decayPlane_p1_NP = pion1_lab.Vect().Unit().Cross(ipv1_lab_3d_NP.Unit());
            PV normal_decayPlane_p2_NP = pion2_lab.Vect().Unit().Cross(ipv2_lab_3d_NP.Unit());

	    histos_["angle_dp_refitVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_NP.Unit())));
	    histos_["angle_dp_refitVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_NP.Unit())));
	    
	    PV normal_decayPlane_p1_genvtx_NP = pion1_lab.Vect().Unit().Cross(ipv1_lab_genvtx_3d_NP.Unit());
            PV normal_decayPlane_p2_genvtx_NP = pion2_lab.Vect().Unit().Cross(ipv2_lab_genvtx_3d_NP.Unit());

	    histos_["angle_dp_genVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p1_genvtx_NP.Unit())));
            histos_["angle_dp_genVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p2_genvtx_NP.Unit())));

	    //Get pT resolution of pion
	    histos_["DeltaPt_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.pt() - pion1_lab.pt()));
	    histos_["DeltaPt_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.pt() - pion2_lab.pt()));

	    histos_["DeltaP_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.P() - pion1_lab.P()));
            histos_["DeltaP_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.P() - pion2_lab.P()));
	    
	    histos_["DeltaPtOverPt_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.pt() - pion1_lab.pt())/pion1_lab.pt());
            histos_["DeltaPtOverPt_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.pt() - pion2_lab.pt())/pion2_lab.pt());

            histos_["DeltaPOverP_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.P() - pion1_lab.P())/pion1_lab.P());
            histos_["DeltaPOverP_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.P() - pion2_lab.P())/pion2_lab.P());

	    //Get angle between IP and pion momentum, use gen pion momenta
	    histos_["xcheck_angle3_reco_genP4"]->Fill(TMath::ACos(ipv1_lab_3d.Unit().Dot(genpionP_lab_3d_u)));
	    histos_["xcheck_angle3_reco_genP4"]->Fill(TMath::ACos(ipv2_lab_3d.Unit().Dot(genpionN_lab_3d_u)));

	    histos_["xcheck_angle3_reco_opv_genP4"]->Fill(TMath::ACos(ipv1_lab_opv_3d.Unit().Dot(genpionP_lab_3d_u)));
	    histos_["xcheck_angle3_reco_opv_genP4"]->Fill(TMath::ACos(ipv2_lab_opv_3d.Unit().Dot(genpionN_lab_3d_u)));

	    histos_["xcheck_angle3_reco_bs_genP4"]->Fill(TMath::ACos(ipv1_lab_bs_3d.Unit().Dot(genpionP_lab_3d_u)));
            histos_["xcheck_angle3_reco_bs_genP4"]->Fill(TMath::ACos(ipv2_lab_bs_3d.Unit().Dot(genpionN_lab_3d_u)));

	    histos_["xcheck_angle3_reco_genvtx_genP4"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d.Unit().Dot(genpionP_lab_3d_u)));
            histos_["xcheck_angle3_reco_genvtx_genP4"]->Fill(TMath::ACos(ipv2_lab_genvtx_3d.Unit().Dot(genpionN_lab_3d_u)));

	    //Get angle between IP and pion momentum, use gen pion IP
	    histos_["xcheck_angle3_reco_genIP"]->Fill(TMath::ACos(ipvP_lab_gen_3d.Unit().Dot(pion1_lab.Vect().Unit())));
            histos_["xcheck_angle3_reco_genIP"]->Fill(TMath::ACos(ipvN_lab_gen_3d.Unit().Dot(pion2_lab.Vect().Unit())));
	    
	    //Angle bwteeen gen and reco pion
	    histos_["angle_pion_momentum"]->Fill(TMath::ACos(genpionP_lab_3d_u.Dot(pion1_lab.Vect().Unit())));
	    histos_["angle_pion_momentum"]->Fill(TMath::ACos(genpionN_lab_3d_u.Dot(pion2_lab.Vect().Unit())));

	    //PCA - genTauDecayVertex projected parallel to gen. tau lepton momentum
	    //PCA - genTauDecayVertex projected perpendicular to gen. tau lepton momentum
	    //At Reco level with linear extrapolation method 1
	    PV pcaTogenTauDecayVtx_gen_1_le1((GenHVtx_.x()+ipv1_lab_genvtx_3d_le1.x()) - genpionPVtx_.x(), (GenHVtx_.y()+ipv1_lab_genvtx_3d_le1.y()) - genpionPVtx_.y(), (GenHVtx_.z()+ipv1_lab_genvtx_3d_le1.z()) - genpionPVtx_.z());
	    PV pcaTogenTauDecayVtx_gen_2_le1((GenHVtx_.x()+ipv2_lab_genvtx_3d_le1.x()) - genpionNVtx_.x(), (GenHVtx_.y()+ipv2_lab_genvtx_3d_le1.y()) - genpionNVtx_.y(), (GenHVtx_.z()+ipv2_lab_genvtx_3d_le1.z()) - genpionNVtx_.z());
	    PV pcaTogenTauDecayVtx_gen_1_le1_l = pcaTogenTauDecayVtx_gen_1_le1.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
	    PV pcaTogenTauDecayVtx_gen_1_le1_t = pcaTogenTauDecayVtx_gen_1_le1 - pcaTogenTauDecayVtx_gen_1_le1_l;
	    PV pcaTogenTauDecayVtx_gen_2_le1_l = pcaTogenTauDecayVtx_gen_2_le1.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_2_le1_t = pcaTogenTauDecayVtx_gen_2_le1 - pcaTogenTauDecayVtx_gen_2_le1_l;
	    //std::cout<<" mag pcaTogenTauDecayVtx_gen_1_le1_l "<<pcaTogenTauDecayVtx_gen_1_le1_l.mag2()<<std::endl;
	    //std::cout<<" mag pcaTogenTauDecayVtx_gen_P_l "<<pcaTogenTauDecayVtx_gen_P_l.mag2()<<std::endl;
	    //At Reco level with linear extrapolation method 2
	    PV pcaTogenTauDecayVtx_gen_1_le2((GenHVtx_.x()+ipv1_lab_genvtx_3d_le2.x()) - genpionPVtx_.x(), (GenHVtx_.y()+ipv1_lab_genvtx_3d_le2.y()) - genpionPVtx_.y(), (GenHVtx_.z()+ipv1_lab_genvtx_3d_le2.z()) - genpionPVtx_.z());
            PV pcaTogenTauDecayVtx_gen_2_le2((GenHVtx_.x()+ipv2_lab_genvtx_3d_le2.x()) - genpionNVtx_.x(), (GenHVtx_.y()+ipv2_lab_genvtx_3d_le2.y()) - genpionNVtx_.y(), (GenHVtx_.z()+ipv2_lab_genvtx_3d_le2.z()) - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_gen_1_le2_l = pcaTogenTauDecayVtx_gen_1_le2.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_1_le2_t = pcaTogenTauDecayVtx_gen_1_le2 - pcaTogenTauDecayVtx_gen_1_le2_l;
            PV pcaTogenTauDecayVtx_gen_2_le2_l = pcaTogenTauDecayVtx_gen_2_le2.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_2_le2_t = pcaTogenTauDecayVtx_gen_2_le2 - pcaTogenTauDecayVtx_gen_2_le2_l;

	    histos_["DeltaTkLength_gen_reco_genvtx_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
	    histos_["DeltaTkLength_gen_reco_genvtx_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
	    histos_["DeltaTkLength_gen_reco_genvtx_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
	    
	    histos_["DeltaTkLength_gen_reco_genvtx_le2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));

	    histos_["AngleTkLength_gen_reco_genvtx_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));
																		  
            histos_["AngleTkLength_gen_reco_genvtx_le2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le2_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le2_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le2_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le2_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));

	    //IP TkLength normal to gen decay plane
            double  pcaTogenTauDecayVtx_gen_1_le1_tdp = pcaTogenTauDecayVtx_gen_1_le1.Dot(normal_decayPlane_P.Unit());
            double  pcaTogenTauDecayVtx_gen_2_le1_tdp = pcaTogenTauDecayVtx_gen_2_le1.Dot(normal_decayPlane_N.Unit());
            histos_["TkLength_gen_reco_genvtx_le1_tdp"]->Fill(pcaTogenTauDecayVtx_gen_1_le1_tdp);
            histos_["TkLength_gen_reco_genvtx_le1_tdp"]->Fill(pcaTogenTauDecayVtx_gen_2_le1_tdp);

	    if(fabs(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2())) < 0.002 &&
               fabs(sqrt(pcaTogenTauDecayVtx_gen_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2())) < 0.002 &&
               fabs(pcaTogenTauDecayVtx_gen_1_le1_tdp) < 0.002 && fabs(pcaTogenTauDecayVtx_gen_2_le1_tdp) < 0.002)
	      histos_["CPPhiStar_genvtx_ResCut"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));

	    //with track quality cuts
            if(passQualityLeg1_ && passQualityLeg2_){
              if(sqrt(ipv1_lab_genvtx_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_genvtx_3d_le1.mag2()) > 0.005){
		histos_["DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
		histos_["TkLength_gen_reco_genvtx_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_1_le1_tdp);
		histos_["TkLength_gen_reco_genvtx_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_2_le1_tdp);
	      }
	      if(sqrt(ipv1_lab_genvtx_3d_le1.mag2()) > 0.003 && sqrt(ipv2_lab_genvtx_3d_le1.mag2()) > 0.003){
		histos_["DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
                histos_["TkLength_gen_reco_genvtx_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_1_le1_tdp);
		histos_["TkLength_gen_reco_genvtx_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_2_le1_tdp);
	      }
	    }

	    //With Re-fitted Reco Vertex,
	    PV pcaTogenTauDecayVtx_1_le1((ReFitVtxX_+ipv1_lab_3d_le1.x()) - genpionPVtx_.x(), (ReFitVtxY_+ipv1_lab_3d_le1.y()) - genpionPVtx_.y(), (ReFitVtxZ_+ipv1_lab_3d_le1.z()) - genpionPVtx_.z());
            PV pcaTogenTauDecayVtx_2_le1((ReFitVtxX_+ipv2_lab_3d_le1.x()) - genpionNVtx_.x(), (ReFitVtxY_+ipv2_lab_3d_le1.y()) - genpionNVtx_.y(), (ReFitVtxZ_+ipv2_lab_3d_le1.z()) - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_1_le1_l = pcaTogenTauDecayVtx_1_le1.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_1_le1_t = pcaTogenTauDecayVtx_1_le1 - pcaTogenTauDecayVtx_1_le1_l;
            PV pcaTogenTauDecayVtx_2_le1_l = pcaTogenTauDecayVtx_2_le1.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_2_le1_t = pcaTogenTauDecayVtx_2_le1 - pcaTogenTauDecayVtx_2_le1_l;

	    histos_["DeltaTkLength_gen_reco_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
            histos_["DeltaTkLength_gen_reco_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));

	    histos_["AngleTkLength_gen_reco_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));
            histos_["AngleTkLength_gen_reco_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));

	    //IP TkLength normal to gen decay plane
	    double  pcaTogenTauDecayVtx_1_le1_tdp = pcaTogenTauDecayVtx_1_le1.Dot(normal_decayPlane_P.Unit());
	    double  pcaTogenTauDecayVtx_2_le1_tdp = pcaTogenTauDecayVtx_2_le1.Dot(normal_decayPlane_N.Unit());
	    histos_["TkLength_gen_reco_le1_tdp"]->Fill(pcaTogenTauDecayVtx_1_le1_tdp);
	    histos_["TkLength_gen_reco_le1_tdp"]->Fill(pcaTogenTauDecayVtx_2_le1_tdp);

	    if(fabs(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2())) < 0.002 &&
	       fabs(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2())) < 0.002 &&
	       fabs(pcaTogenTauDecayVtx_1_le1_tdp) < 0.002 && fabs(pcaTogenTauDecayVtx_2_le1_tdp) < 0.002)
	      histos_["CPPhiStar_ResCut"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));
	    
	    //with track quality cuts
	    if(passQualityLeg1_ && passQualityLeg2_){
	      if(sqrt(ipv1_lab_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_3d_le1.mag2()) > 0.005){
		histos_["DeltaTkLength_gen_reco_ip50_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
		histos_["DeltaTkLength_gen_reco_ip50_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
		histos_["DeltaTkLength_gen_reco_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["DeltaTkLength_gen_reco_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
		histos_["TkLength_gen_reco_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_1_le1_tdp);
		histos_["TkLength_gen_reco_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_2_le1_tdp);
	      }
	      if(sqrt(ipv1_lab_3d_le1.mag2()) > 0.003 && sqrt(ipv2_lab_3d_le1.mag2()) > 0.003){
		histos_["DeltaTkLength_gen_reco_ip30_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
		histos_["DeltaTkLength_gen_reco_ip30_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
		histos_["DeltaTkLength_gen_reco_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["DeltaTkLength_gen_reco_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
		histos_["TkLength_gen_reco_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_1_le1_tdp);
		histos_["TkLength_gen_reco_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_2_le1_tdp);
	      }
	    }
	    
	    //With PCA without LE
	    PV pcaTogenTauDecayVtx_1_m2((*diTauLegsPCAM2_)[0].x() - genpionPVtx_.x(), (*diTauLegsPCAM2_)[0].y() - genpionPVtx_.y(), (*diTauLegsPCAM2_)[0].z() - genpionPVtx_.z());
            PV pcaTogenTauDecayVtx_2_m2((*diTauLegsPCAM2_)[1].x() - genpionNVtx_.x(), (*diTauLegsPCAM2_)[1].y() - genpionNVtx_.y(), (*diTauLegsPCAM2_)[1].z() - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_1_m2_l = pcaTogenTauDecayVtx_1_m2.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_1_m2_t = pcaTogenTauDecayVtx_1_m2 - pcaTogenTauDecayVtx_1_m2_l;
            PV pcaTogenTauDecayVtx_2_m2_l = pcaTogenTauDecayVtx_2_m2.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_2_m2_t = pcaTogenTauDecayVtx_2_m2 - pcaTogenTauDecayVtx_2_m2_l;

            histos_["DeltaTkLength_gen_reco_m2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_1_m2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_m2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_2_m2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_m2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_1_m2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
            histos_["DeltaTkLength_gen_reco_m2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_2_m2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));

            histos_["AngleTkLength_gen_reco_m2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_m2_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_m2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_m2_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_m2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_m2_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));
            histos_["AngleTkLength_gen_reco_m2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_m2_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));

            //IP TkLength normal to gen decay plane
            double  pcaTogenTauDecayVtx_1_m2_tdp = pcaTogenTauDecayVtx_1_m2.Dot(normal_decayPlane_P.Unit());
            double  pcaTogenTauDecayVtx_2_m2_tdp = pcaTogenTauDecayVtx_2_m2.Dot(normal_decayPlane_N.Unit());
            histos_["TkLength_gen_reco_m2_tdp"]->Fill(pcaTogenTauDecayVtx_1_m2_tdp);
            histos_["TkLength_gen_reco_m2_tdp"]->Fill(pcaTogenTauDecayVtx_2_m2_tdp);


	  }
	  else if(ROOT::Math::VectorUtil::DeltaR(genpionN_lab, pion1_lab) < 0.3 && ROOT::Math::VectorUtil::DeltaR(genpionP_lab, pion2_lab) < 0.3 ){

	    histos_["DeltaX_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[0].x() - pcaN.x());
            histos_["DeltaX_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[1].x() - pcaP.x());

            histos_["DeltaY_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[0].y() - pcaN.y());
            histos_["DeltaY_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[1].y() - pcaP.y());

            histos_["DeltaZ_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[0].z() - pcaN.z());
            histos_["DeltaZ_gen_reco_pca_genvtx"]->Fill((*diTauLegsPCAGen_)[1].z() - pcaP.z());

	    histos_["DeltaX_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.x()+ipv1_lab_genvtx_3d_le1.x()) - pcaN.x());
            histos_["DeltaX_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.x()+ipv2_lab_genvtx_3d_le1.x()) - pcaP.x());
            histos_["DeltaY_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.y()+ipv1_lab_genvtx_3d_le1.y()) - pcaN.y());
            histos_["DeltaY_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.y()+ipv2_lab_genvtx_3d_le1.y()) - pcaP.y());
            histos_["DeltaZ_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.z()+ipv1_lab_genvtx_3d_le1.z()) - pcaN.z());
            histos_["DeltaZ_gen_reco_pcaLE1_genvtx"]->Fill((GenHVtx_.z()+ipv2_lab_genvtx_3d_le1.z()) - pcaP.z());

            histos_["DeltaX_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.x()+ipv1_lab_genvtx_3d_le2.x()) - pcaN.x());
            histos_["DeltaX_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.x()+ipv2_lab_genvtx_3d_le2.x()) - pcaP.x());
            histos_["DeltaY_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.y()+ipv1_lab_genvtx_3d_le2.y()) - pcaN.y());
            histos_["DeltaY_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.y()+ipv2_lab_genvtx_3d_le2.y()) - pcaP.y());
            histos_["DeltaZ_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.z()+ipv1_lab_genvtx_3d_le2.z()) - pcaN.z());
            histos_["DeltaZ_gen_reco_pcaLE2_genvtx"]->Fill((GenHVtx_.z()+ipv2_lab_genvtx_3d_le2.z()) - pcaP.z());

	    histos_["angle_dp_refitVtx"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2.Unit())));
            histos_["angle_dp_refitVtx"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1.Unit())));

            histos_["angle_dp_origVtx"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_opv.Unit())));
            histos_["angle_dp_origVtx"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_opv.Unit())));

            histos_["angle_dp_genVtx"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_genvtx.Unit())));
            histos_["angle_dp_genVtx"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_genvtx.Unit())));

	    histos_["angle_dp_refitVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_le1.Unit())));
            histos_["angle_dp_refitVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_le1.Unit())));
            histos_["angle_dp_origVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_opv_le1.Unit())));
            histos_["angle_dp_origVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_opv_le1.Unit())));
            histos_["angle_dp_genVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_genvtx_le1.Unit())));
            histos_["angle_dp_genVtx_LE1"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_genvtx_le1.Unit())));

	    //get normal to decay planes in detector, replacing detector Pi p4 by gen Pi p4
            PV normal_decayPlane_p1_genP4 = genpionN_lab_3d_u.Cross(ipv1_lab_3d);
            PV normal_decayPlane_p2_genP4 = genpionP_lab_3d_u.Cross(ipv2_lab_3d);

            PV normal_decayPlane_p1_opv_genP4 = genpionN_lab_3d_u.Cross(ipv1_lab_opv_3d);
            PV normal_decayPlane_p2_opv_genP4 = genpionP_lab_3d_u.Cross(ipv2_lab_opv_3d);

            PV normal_decayPlane_p1_genvtx_genP4 = genpionN_lab_3d_u.Cross(ipv1_lab_genvtx_3d);
            PV normal_decayPlane_p2_genvtx_genP4 = genpionP_lab_3d_u.Cross(ipv2_lab_genvtx_3d);

            histos_["angle_dp_refitVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_genP4.Unit())));
            histos_["angle_dp_refitVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_genP4.Unit())));

            histos_["angle_dp_origVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_opv_genP4.Unit())));
            histos_["angle_dp_origVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_opv_genP4.Unit())));

            histos_["angle_dp_genVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_genvtx_genP4.Unit())));
            histos_["angle_dp_genVtx_genP4"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_genvtx_genP4.Unit())));

	    histos_["delphi_dp_genVtx"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(normal_decayPlane_P.Unit(), normal_decayPlane_p2_genvtx.Unit()));
            histos_["delphi_dp_genVtx"]->Fill(ROOT::Math::VectorUtil::DeltaPhi(normal_decayPlane_N.Unit(), normal_decayPlane_p1_genvtx.Unit()));

            //get normal to decay planes in detector, replacing detector Pi IP by gen Pi IP
            PV normal_decayPlane_p1_genIP = pion1_lab.Vect().Unit().Cross(ipvN_lab_gen_3d);
            PV normal_decayPlane_p2_genIP = pion2_lab.Vect().Unit().Cross(ipvP_lab_gen_3d);

            histos_["angle_dp_genIP"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_genIP.Unit())));
            histos_["angle_dp_genIP"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_genIP.Unit())));

	    //get normal to decay planes in detector, replacing detector Pi IP by normal component of IP wrt pion momenta
            PV normal_decayPlane_p1_NP = pion1_lab.Vect().Unit().Cross(ipv1_lab_3d_NP.Unit());
            PV normal_decayPlane_p2_NP = pion2_lab.Vect().Unit().Cross(ipv2_lab_3d_NP.Unit());

            histos_["angle_dp_refitVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_NP.Unit())));
            histos_["angle_dp_refitVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_NP.Unit())));

            PV normal_decayPlane_p1_genvtx_NP = pion1_lab.Vect().Unit().Cross(ipv1_lab_genvtx_3d_NP.Unit());
            PV normal_decayPlane_p2_genvtx_NP = pion2_lab.Vect().Unit().Cross(ipv2_lab_genvtx_3d_NP.Unit());

            histos_["angle_dp_genVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_P.Unit().Dot(normal_decayPlane_p2_genvtx_NP.Unit())));
            histos_["angle_dp_genVtx_NP"]->Fill(TMath::ACos(normal_decayPlane_N.Unit().Dot(normal_decayPlane_p1_genvtx_NP.Unit())));

            //Get pT resolution of pion
	    histos_["DeltaPt_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.pt() - pion2_lab.pt()));
            histos_["DeltaPt_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.pt() - pion1_lab.pt()));

            histos_["DeltaP_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.P() - pion2_lab.P()));
            histos_["DeltaP_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.P() - pion1_lab.P()));

	    histos_["DeltaPtOverPt_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.pt() - pion2_lab.pt())/pion2_lab.pt());
            histos_["DeltaPtOverPt_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.pt() - pion1_lab.pt())/pion1_lab.pt());

            histos_["DeltaPOverP_Gen_Reco"]->Fill(TMath::Abs(genpionP_lab.P() - pion2_lab.P())/pion2_lab.P());
            histos_["DeltaPOverP_Gen_Reco"]->Fill(TMath::Abs(genpionN_lab.P() - pion1_lab.P())/pion1_lab.P());

	    histos_["xcheck_angle3_reco_genP4"]->Fill(TMath::ACos(ipv1_lab_3d.Unit().Dot(genpionN_lab_3d_u)));
            histos_["xcheck_angle3_reco_genP4"]->Fill(TMath::ACos(ipv2_lab_3d.Unit().Dot(genpionP_lab_3d_u)));

            histos_["xcheck_angle3_reco_opv_genP4"]->Fill(TMath::ACos(ipv1_lab_opv_3d.Unit().Dot(genpionN_lab_3d_u)));
            histos_["xcheck_angle3_reco_opv_genP4"]->Fill(TMath::ACos(ipv2_lab_opv_3d.Unit().Dot(genpionP_lab_3d_u)));

            histos_["xcheck_angle3_reco_bs_genP4"]->Fill(TMath::ACos(ipv1_lab_bs_3d.Unit().Dot(genpionN_lab_3d_u)));
            histos_["xcheck_angle3_reco_bs_genP4"]->Fill(TMath::ACos(ipv2_lab_bs_3d.Unit().Dot(genpionP_lab_3d_u)));

            histos_["xcheck_angle3_reco_genvtx_genP4"]->Fill(TMath::ACos(ipv1_lab_genvtx_3d.Unit().Dot(genpionN_lab_3d_u)));
            histos_["xcheck_angle3_reco_genvtx_genP4"]->Fill(TMath::ACos(ipv2_lab_genvtx_3d.Unit().Dot(genpionP_lab_3d_u)));

	    //Get angle between IP and pion momentum, use gen pion IP
            histos_["xcheck_angle3_reco_genIP"]->Fill(TMath::ACos(ipvP_lab_gen_3d.Unit().Dot(pion2_lab.Vect().Unit())));
            histos_["xcheck_angle3_reco_genIP"]->Fill(TMath::ACos(ipvN_lab_gen_3d.Unit().Dot(pion1_lab.Vect().Unit())));

	    //Angle bwteeen gen and reco pion
            histos_["angle_pion_momentum"]->Fill(TMath::ACos(genpionP_lab_3d_u.Dot(pion2_lab.Vect().Unit())));
            histos_["angle_pion_momentum"]->Fill(TMath::ACos(genpionN_lab_3d_u.Dot(pion1_lab.Vect().Unit())));

	    //PCA - genTauDecayVertex projected parallel to gen. tau lepton momentum
            //PCA - genTauDecayVertex projected perpendicular to gen. tau lepton momentum
            //At Reco level with linear extrapolation method 1
            PV pcaTogenTauDecayVtx_gen_1_le1((GenHVtx_.x()+ipv1_lab_genvtx_3d_le1.x()) - genpionNVtx_.x(), (GenHVtx_.y()+ipv1_lab_genvtx_3d_le1.y()) - genpionNVtx_.y(), (GenHVtx_.z()+ipv1_lab_genvtx_3d_le1.z()) - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_gen_2_le1((GenHVtx_.x()+ipv2_lab_genvtx_3d_le1.x()) - genpionPVtx_.x(), (GenHVtx_.y()+ipv2_lab_genvtx_3d_le1.y()) - genpionPVtx_.y(), (GenHVtx_.z()+ipv2_lab_genvtx_3d_le1.z()) - genpionPVtx_.z());
            PV pcaTogenTauDecayVtx_gen_1_le1_l = pcaTogenTauDecayVtx_gen_1_le1.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_1_le1_t = pcaTogenTauDecayVtx_gen_1_le1 - pcaTogenTauDecayVtx_gen_1_le1_l;
            PV pcaTogenTauDecayVtx_gen_2_le1_l = pcaTogenTauDecayVtx_gen_2_le1.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_2_le1_t = pcaTogenTauDecayVtx_gen_2_le1 - pcaTogenTauDecayVtx_gen_2_le1_l;
            //At Reco level with linear extrapolation method 2
            PV pcaTogenTauDecayVtx_gen_1_le2((GenHVtx_.x()+ipv1_lab_genvtx_3d_le2.x()) - genpionNVtx_.x(), (GenHVtx_.y()+ipv1_lab_genvtx_3d_le2.y()) - genpionNVtx_.y(), (GenHVtx_.z()+ipv1_lab_genvtx_3d_le2.z()) - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_gen_2_le2((GenHVtx_.x()+ipv2_lab_genvtx_3d_le2.x()) - genpionPVtx_.x(), (GenHVtx_.y()+ipv2_lab_genvtx_3d_le2.y()) - genpionPVtx_.y(), (GenHVtx_.z()+ipv2_lab_genvtx_3d_le2.z()) - genpionPVtx_.z());
            PV pcaTogenTauDecayVtx_gen_1_le2_l = pcaTogenTauDecayVtx_gen_1_le2.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_1_le2_t = pcaTogenTauDecayVtx_gen_1_le2 - pcaTogenTauDecayVtx_gen_1_le2_l;
            PV pcaTogenTauDecayVtx_gen_2_le2_l = pcaTogenTauDecayVtx_gen_2_le2.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_gen_2_le2_t = pcaTogenTauDecayVtx_gen_2_le2 - pcaTogenTauDecayVtx_gen_2_le2_l;

            histos_["DeltaTkLength_gen_reco_genvtx_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));

            histos_["DeltaTkLength_gen_reco_genvtx_le2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
            histos_["DeltaTkLength_gen_reco_genvtx_le2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));

            histos_["AngleTkLength_gen_reco_genvtx_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
	    histos_["AngleTkLength_gen_reco_genvtx_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));

            histos_["AngleTkLength_gen_reco_genvtx_le2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le2_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le2_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_1_le2_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));
            histos_["AngleTkLength_gen_reco_genvtx_le2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_gen_2_le2_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));

	    //IP TkLength normal to gen decay plane
            double  pcaTogenTauDecayVtx_gen_1_le1_tdp = pcaTogenTauDecayVtx_gen_1_le1.Dot(normal_decayPlane_N.Unit());
            double  pcaTogenTauDecayVtx_gen_2_le1_tdp = pcaTogenTauDecayVtx_gen_2_le1.Dot(normal_decayPlane_P.Unit());
            histos_["TkLength_gen_reco_genvtx_le1_tdp"]->Fill(pcaTogenTauDecayVtx_gen_1_le1_tdp);
            histos_["TkLength_gen_reco_genvtx_le1_tdp"]->Fill(pcaTogenTauDecayVtx_gen_2_le1_tdp);

	    if(fabs(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2())) < 0.002 &&
               fabs(sqrt(pcaTogenTauDecayVtx_gen_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2())) < 0.002 &&
               fabs(pcaTogenTauDecayVtx_gen_1_le1_tdp) < 0.002 && fabs(pcaTogenTauDecayVtx_gen_2_le1_tdp) < 0.002)
	      histos_["CPPhiStar_genvtx_ResCut"]->Fill(TMath::ACos(ipv1_genvtx_rf_3d_le1_t_u.Dot(ipv2_genvtx_rf_3d_le1_t_u)));

	    //with track quality cuts
            if(passQualityLeg1_ && passQualityLeg2_){
              if(sqrt(ipv1_lab_genvtx_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_genvtx_3d_le1.mag2()) > 0.005){
                histos_["DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
                histos_["DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["TkLength_gen_reco_genvtx_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_1_le1_tdp);
                histos_["TkLength_gen_reco_genvtx_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_2_le1_tdp);
	      }
	      if(sqrt(ipv1_lab_genvtx_3d_le1.mag2()) > 0.003 && sqrt(ipv2_lab_genvtx_3d_le1.mag2()) > 0.003){
		histos_["DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
		histos_["DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_gen_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["TkLength_gen_reco_genvtx_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_1_le1_tdp);
                histos_["TkLength_gen_reco_genvtx_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_gen_2_le1_tdp);
              }
            }

	    //With Re-fitted Reco Vertex, ReFitVtxX_
	    PV pcaTogenTauDecayVtx_1_le1((ReFitVtxX_+ipv1_lab_3d_le1.x()) - genpionNVtx_.x(), (ReFitVtxY_+ipv1_lab_3d_le1.y()) - genpionNVtx_.y(), (ReFitVtxZ_+ipv1_lab_3d_le1.z()) - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_2_le1((ReFitVtxX_+ipv2_lab_3d_le1.x()) - genpionPVtx_.x(), (ReFitVtxY_+ipv2_lab_3d_le1.y()) - genpionPVtx_.y(), (ReFitVtxZ_+ipv2_lab_3d_le1.z()) - genpionPVtx_.z());
            PV pcaTogenTauDecayVtx_1_le1_l = pcaTogenTauDecayVtx_1_le1.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_1_le1_t = pcaTogenTauDecayVtx_1_le1 - pcaTogenTauDecayVtx_1_le1_l;
            PV pcaTogenTauDecayVtx_2_le1_l = pcaTogenTauDecayVtx_2_le1.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_2_le1_t = pcaTogenTauDecayVtx_2_le1 - pcaTogenTauDecayVtx_2_le1_l;
	    
	    histos_["DeltaTkLength_gen_reco_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_le1_l"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
            histos_["DeltaTkLength_gen_reco_le1_t"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));

	    histos_["AngleTkLength_gen_reco_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_le1_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_le1_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));
            histos_["AngleTkLength_gen_reco_le1_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_le1_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));

	    //IP TkLength normal to gen decay plane
            double  pcaTogenTauDecayVtx_1_le1_tdp = pcaTogenTauDecayVtx_1_le1.Dot(normal_decayPlane_N.Unit());
            double  pcaTogenTauDecayVtx_2_le1_tdp = pcaTogenTauDecayVtx_2_le1.Dot(normal_decayPlane_P.Unit());
            histos_["TkLength_gen_reco_le1_tdp"]->Fill(pcaTogenTauDecayVtx_1_le1_tdp);
            histos_["TkLength_gen_reco_le1_tdp"]->Fill(pcaTogenTauDecayVtx_2_le1_tdp);

	    if(fabs(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2())) < 0.002 &&
	       fabs(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2())) < 0.002 &&
	       fabs(pcaTogenTauDecayVtx_1_le1_tdp) < 0.002 && fabs(pcaTogenTauDecayVtx_2_le1_tdp) < 0.002)
              histos_["CPPhiStar_ResCut"]->Fill(TMath::ACos(ipv1_rf_3d_le1_t_u.Dot(ipv2_rf_3d_le1_t_u)));

	    //with track quality cuts
            if(passQualityLeg1_ && passQualityLeg2_){
              if(sqrt(ipv1_lab_3d_le1.mag2()) > 0.005 && sqrt(ipv2_lab_3d_le1.mag2()) > 0.005){
                histos_["DeltaTkLength_gen_reco_ip50_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
                histos_["DeltaTkLength_gen_reco_ip50_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
                histos_["DeltaTkLength_gen_reco_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
                histos_["DeltaTkLength_gen_reco_ip50_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
		histos_["TkLength_gen_reco_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_1_le1_tdp);
                histos_["TkLength_gen_reco_ip50_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_2_le1_tdp);
              }
              if(sqrt(ipv1_lab_3d_le1.mag2()) > 0.003 && sqrt(ipv2_lab_3d_le1.mag2()) > 0.003){
                histos_["DeltaTkLength_gen_reco_ip30_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
                histos_["DeltaTkLength_gen_reco_ip30_l_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
                histos_["DeltaTkLength_gen_reco_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_1_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
                histos_["DeltaTkLength_gen_reco_ip30_t_tkQCut"]->Fill(sqrt(pcaTogenTauDecayVtx_2_le1_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));
                histos_["TkLength_gen_reco_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_1_le1_tdp);
                histos_["TkLength_gen_reco_ip30_tdp_tkQCut"]->Fill(pcaTogenTauDecayVtx_2_le1_tdp);
              }
            }

	    //With PCA without LE
            PV pcaTogenTauDecayVtx_1_m2((*diTauLegsPCAM2_)[0].x() - genpionNVtx_.x(), (*diTauLegsPCAM2_)[0].y() - genpionNVtx_.y(), (*diTauLegsPCAM2_)[0].z() - genpionNVtx_.z());
            PV pcaTogenTauDecayVtx_2_m2((*diTauLegsPCAM2_)[1].x() - genpionPVtx_.x(), (*diTauLegsPCAM2_)[1].y() - genpionPVtx_.y(), (*diTauLegsPCAM2_)[1].z() - genpionPVtx_.z());
	    PV pcaTogenTauDecayVtx_1_m2_l = pcaTogenTauDecayVtx_1_m2.Dot(gentauN_lab.Vect().Unit())*gentauN_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_1_m2_t = pcaTogenTauDecayVtx_1_m2 - pcaTogenTauDecayVtx_1_m2_l;
	    PV pcaTogenTauDecayVtx_2_m2_l = pcaTogenTauDecayVtx_2_m2.Dot(gentauP_lab.Vect().Unit())*gentauP_lab.Vect().Unit();
            PV pcaTogenTauDecayVtx_2_m2_t = pcaTogenTauDecayVtx_2_m2 - pcaTogenTauDecayVtx_2_m2_l;

            histos_["DeltaTkLength_gen_reco_m2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_1_m2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_l.mag2()));
            histos_["DeltaTkLength_gen_reco_m2_l"]->Fill(sqrt(pcaTogenTauDecayVtx_2_m2_l.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_l.mag2()));
            histos_["DeltaTkLength_gen_reco_m2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_1_m2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_N_t.mag2()));
            histos_["DeltaTkLength_gen_reco_m2_t"]->Fill(sqrt(pcaTogenTauDecayVtx_2_m2_t.mag2()) - sqrt(pcaTogenTauDecayVtx_gen_P_t.mag2()));

            histos_["AngleTkLength_gen_reco_m2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_m2_l.unit().Dot(pcaTogenTauDecayVtx_gen_N_l.unit())));
            histos_["AngleTkLength_gen_reco_m2_l"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_m2_l.unit().Dot(pcaTogenTauDecayVtx_gen_P_l.unit())));
            histos_["AngleTkLength_gen_reco_m2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_1_m2_t.unit().Dot(pcaTogenTauDecayVtx_gen_N_t.unit())));
            histos_["AngleTkLength_gen_reco_m2_t"]->Fill(TMath::ACos(pcaTogenTauDecayVtx_2_m2_t.unit().Dot(pcaTogenTauDecayVtx_gen_P_t.unit())));

            //IP TkLength normal to gen decay plane
            double  pcaTogenTauDecayVtx_1_m2_tdp = pcaTogenTauDecayVtx_1_m2.Dot(normal_decayPlane_N.Unit());
            double  pcaTogenTauDecayVtx_2_m2_tdp = pcaTogenTauDecayVtx_2_m2.Dot(normal_decayPlane_P.Unit());
            histos_["TkLength_gen_reco_m2_tdp"]->Fill(pcaTogenTauDecayVtx_1_m2_tdp);
            histos_["TkLength_gen_reco_m2_tdp"]->Fill(pcaTogenTauDecayVtx_2_m2_tdp);

	  }
	  
	}
      }
    } //end gen match


  }

  outFile_->Write();
  outFile_->Close();
  
  delete diTauLegsPCA_; delete diTauLegsPCAOPV_; delete diTauLegsLchP4_; delete diTauLegsP4_;
  delete diTauLegsPCABS_; delete diTauLegsPCAGen_;
  delete genTausP4_; delete genVP4_; delete genTauPSonsP4_; delete genTauNSonsP4_;

  delete VtxPos_; delete BSPos_; delete ReFitVtxPos_; delete HiggsGenVtx_; delete IanVtxPos_;
  delete TausGenVtx_; delete TauPSonsGenVtx_; delete TauNSonsGenVtx_;

  delete diTauLegsLchP3AtPCA_; delete diTauLegsLchP3AtPCAOPV_; 
  delete diTauLegsLchP3AtPCABS_; delete diTauLegsLchP3AtPCAGen_;
  delete diTauLegsIPAtPCA_; delete diTauLegsIPAtPCAV2_;
  delete diTauLegsIPAtPCAOPV_; delete diTauLegsIPAtPCAOPVV2_;
  delete diTauLegsIPAtPCAGen_; delete diTauLegsIPAtPCAGenV2_;
  delete diTauLegsIPAtPCATauVtx_;

  delete genTausPid_; delete genTausCharge_; delete genTausStatus_;
  delete genTauPSonsPid_; delete genTauPSonsCharge_; delete genTauPSonsStatus_;
  delete genTauNSonsPid_; delete genTauNSonsCharge_; delete genTauNSonsStatus_;
}

void VertexAnalysisAll()
{
  
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-125_MC_v13_vtxWithBS/", "GluGluToHToTauTau_M-125_MC_v13_vtxWithBS_anal", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-120_MC_v17_vtxWithBS/", "GluGluToHToTauTau_M-120_MC_v17_vtxWithBS_anal", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/VBF_HToTauTau_M-125_MC_v14_vtxWithBS/", "VBF_HToTauTau_M-125_MC_v14_vtxWithBS_anal", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-120_MC_v17_vtxWithBS/", "SUSYGluGluToHToTauTau_M-120_MC_v17_vtxWithBS_anal", 36);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/DYJetsToLL_M-50_MC_v17_vtxWithBS/", "DYJetsToLL_M-50_MC_v17_vtxWithBS_anal", 23);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS/", "GluGluToHToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS_anal", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/VBF_HToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS/", "VBF_HToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS_anal", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/DYJetsToLL_M-50_tauPolarOff_MC_v13_vtxWithBS/", "DYJetsToLL_M-50_tauPolarOff_MC_v13_vtxWithBS_anal", 23);

  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/DYJetsToLL_M-50_MC_v18/", "CPWthoutPi0s/DYJetsToLL_M-50_MC_v18", 23);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-110_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-110_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-115_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-115_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-120_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-120_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-125_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-125_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-130_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-130_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-135_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-135_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-140_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-140_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-145_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-145_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-150_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-150_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-155_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-155_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-160_MC_v18/", "CPWthoutPi0s/GluGluToHToTauTau_M-160_MC_v18", 25);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-100_MC_v18/", "CPWthoutPi0s/SUSYGluGluToHToTauTau_M-100_MC_v18", 36);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-110_MC_v18/", "CPWthoutPi0s/SUSYGluGluToHToTauTau_M-110_MC_v18", 36);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-120_MC_v18/", "CPWthoutPi0s/SUSYGluGluToHToTauTau_M-120_MC_v18", 36);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-130_MC_v18/", "CPWthoutPi0s/SUSYGluGluToHToTauTau_M-130_MC_v18", 36);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-140_MC_v18/", "CPWthoutPi0s/SUSYGluGluToHToTauTau_M-140_MC_v18", 36);
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-160_MC_v18/", "CPWthoutPi0s/SUSYGluGluToHToTauTau_M-160_MC_v18", 36);


  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-120_MC_v19/", "GluGluToHToTauTau_M-120_MC_v19", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-120_MC_v19/", "SUSYGluGluToHToTauTau_M-120_MC_v19", 36);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/VBF0P_HToTauTau_M-125_JHUGenV4_MC_v18/", "VBF0P_HToTauTau_M-125_JHUGenV4_MC_v18", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/VBF0M_HToTauTau_M-125_JHUGenV4_MC_v18/", "VBF0M_HToTauTau_M-125_JHUGenV4_MC_v18", 25);

  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GGF120_8TeV_v25_noBS/", "CPWthoutPi0s/GluGluToHToTauTau_M-120_MC_v25_noBS", 25);
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GGF125_8TeV_v25_noBS/", "CPWthoutPi0s/GluGluToHToTauTau_M-125_MC_v25_noBS", 25);
  
} 
