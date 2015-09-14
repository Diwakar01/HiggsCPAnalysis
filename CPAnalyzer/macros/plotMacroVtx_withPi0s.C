void plotMacro(TString refHist,  TString refitHist, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_, TString inFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");
  
  //TFile *f1 = TFile::Open("GluGluToHToTauTau_M-125_MC_v4_anal.root", "READONLY");
  //TFile *f1 = TFile::Open("VBF_HToTauTau_M-125_MC_v4_anal.root", "READONLY");
  TFile *f1 = TFile::Open(inFileName_, "READONLY");
  
  TH1F* hPVorig = (TH1F*)f1->Get(refHist);
  hPVorig->GetXaxis()->SetRangeUser(xmin_, xmax_);
  hPVorig->SetMaximum(hPVorig->GetMaximum()*1.2);
  //hPVorig->GetYaxis()->SetRangeUser(ymin_, ymax_);
  hPVorig->SetLineColor(kRed);
  hPVorig->GetXaxis()->SetTitle(xtitle_);
  hPVorig->GetYaxis()->SetTitle(ytitle_);

  TH1F* hPVrefit = (TH1F*)f1->Get(refitHist);
  hPVrefit->SetLineColor(kBlue);
  
  hPVorig->Draw();
  hPVrefit->Draw("same");

  TLegend* leg = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hPVorig,"Original PV");
  leg->AddEntry(hPVrefit,"Re-fitted PV");
  
  leg->Draw();
  
  c1->SaveAs("plots/vertex/"+outFileName_+".png");
}

void plotMacroAll()
{
  /*
  plotMacro("VtxXRes",  "ReFitVtxXRes", 0.0, 0.04, 0., 1.0, "Resolution X-position", "Arbitrary Units", "GGFH_VtxResolution_X");
  plotMacro("VtxYRes",  "ReFitVtxYRes", 0.0, 0.04, 0., 1.0, "Resolution Y-position", "Arbitrary Units", "GGFH_VtxResolution_Y");
  plotMacro("VtxZRes",  "ReFitVtxZRes", -0.02, 0.02, 0., 1.0, "Resolution Z-position", "Arbitrary Units", "GGFH_VtxResolution_Z");
  */
  /*
  plotMacro("VtxXRes",  "ReFitVtxXRes", 0.0, 0.04, 0., 1.0, "Resolution X-position", "Arbitrary Units", "VBFH_VtxResolution_X");
  plotMacro("VtxYRes",  "ReFitVtxYRes", 0.0, 0.04, 0., 1.0, "Resolution Y-position", "Arbitrary Units", "VBFH_VtxResolution_Y");
  plotMacro("VtxZRes",  "ReFitVtxZRes", -0.02, 0.02, 0., 1.0, "Resolution Z-position", "Arbitrary Units", "VBFH_VtxResolution_Z");
  */
  /*
  //plotMacro("VtxXErr",  "ReFitVtxXErr", 0.0, 0.02, 0., 1.0, "Error X-position", "Arbitrary Units", "VBFH_VtxError_X", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  //plotMacro("VtxYErr",  "ReFitVtxYErr", 0.0, 0.02, 0., 1.0, "Error Y-position", "Arbitrary Units", "VBFH_VtxError_Y", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  //plotMacro("VtxZErr",  "ReFitVtxZErr", 0.0, 0.02, 0., 1.0, "Error Z-position", "Arbitrary Units", "VBFH_VtxError_Z", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("VtxXErr",  "ReFitVtxXErr", 0.0, 0.005, 0., 1.0, "Error X-position", "Arbitrary Units", "GGFH_VtxError_X", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  plotMacro("VtxYErr",  "ReFitVtxYErr", 0.0, 0.005, 0., 1.0, "Error Y-position", "Arbitrary Units", "GGFH_VtxError_Y", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  plotMacro("VtxZErr",  "ReFitVtxZErr", 0.0, 0.02, 0., 1.0, "Error Z-position", "Arbitrary Units", "GGFH_VtxError_Z", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");


  //plotMacro("DeltaPVwrtGenX",  "DeltaRfPVwrtGenX", -0.01, 0.01, 0., 1.0, "Gen X-position - Reco X-position", "Arbitrary Units", "VBFH_PosDiffWrtGen_X", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  //plotMacro("DeltaPVwrtGenY",  "DeltaRfPVwrtGenY", -0.01, 0.01, 0., 1.0, "Gen Y-position - Reco Y-position", "Arbitrary Units", "VBFH_PosDiffWrtGen_Y", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  //plotMacro("DeltaPVwrtGenZ",  "DeltaRfPVwrtGenZ", -0.01, 0.01, 0., 1.0, "Gen Z-position - Reco Z-position", "Arbitrary Units", "VBFH_PosDiffWrtGen_Z", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("DeltaPVwrtGenX",  "DeltaRfPVwrtGenX", -0.01, 0.01, 0., 1.0, "Gen X-position - Reco X-position", "Arbitrary Units", "GGFH_PosDiffWrtGen_X", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  plotMacro("DeltaPVwrtGenY",  "DeltaRfPVwrtGenY", -0.01, 0.01, 0., 1.0, "Gen Y-position - Reco Y-position", "Arbitrary Units", "GGFH_PosDiffWrtGen_Y", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  plotMacro("DeltaPVwrtGenZ",  "DeltaRfPVwrtGenZ", -0.01, 0.01, 0., 1.0, "Gen Z-position - Reco Z-position", "Arbitrary Units", "GGFH_PosDiffWrtGen_Z", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");

  //plotMacro("ResPVwrtGenX",  "ResRfPVwrtGenX", -5.0, 5.0, 0., 1.0, "(Gen X - Reco X) / Reco X Error", "Arbitrary Units", "VBFH_PosResWrtGen_X", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  //plotMacro("ResPVwrtGenY",  "ResRfPVwrtGenY", -5.0, 5.0, 0., 1.0, "(Gen Y - Reco Y) / Reco Y Error", "Arbitrary Units", "VBFH_PosResWrtGen_Y", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  //plotMacro("ResPVwrtGenZ",  "ResRfPVwrtGenZ", -5.0, 5.0, 0., 1.0, "(Gen Z - Reco Z) / Reco Z Error", "Arbitrary Units", "VBFH_PosResWrtGen_Z", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("ResPVwrtGenX",  "ResRfPVwrtGenX", -5.0, 5.0, 0., 1.0, "(Gen X - Reco X) / Reco X Error", "Arbitrary Units", "GGFH_PosResWrtGen_X", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  plotMacro("ResPVwrtGenY",  "ResRfPVwrtGenY", -5.0, 5.0, 0., 1.0, "(Gen Y - Reco Y) / Reco Y Error", "Arbitrary Units", "GGFH_PosResWrtGen_Y", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  plotMacro("ResPVwrtGenZ",  "ResRfPVwrtGenZ", -5.0, 5.0, 0., 1.0, "(Gen Z - Reco Z) / Reco Z Error", "Arbitrary Units", "GGFH_PosResWrtGen_Z", "GluGluToHToTauTau_M-120_MC_v11_vtxWithBS_anal.root");
  */

}

void plotMacroSample(TString histname_, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile *f1 = TFile::Open("VBF_HToTauTau_M-125_MC_v4_noBS_anal.root", "READONLY");
  TH1F* hPVvbf = (TH1F*)f1->Get(histname_);
  hPVvbf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  //hPVvbf->GetYaxis()->SetRangeUser(ymin_, ymax_);
  hPVvbf->SetLineColor(kRed);
  hPVvbf->GetXaxis()->SetTitle(xtitle_);
  hPVvbf->GetYaxis()->SetTitle(ytitle_);

  TFile *f2 = TFile::Open("GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root", "READONLY");
  TH1F* hPVggf = (TH1F*)f2->Get(histname_);
  hPVggf->SetLineColor(kBlue);

  hPVvbf->DrawNormalized();
  hPVggf->DrawNormalized("same");

  TLegend* leg = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hPVvbf,"VBF H#rightarrow#tau#tau");
  leg->AddEntry(hPVggf,"GGF H#rightarrow#tau#tau");

  leg->Draw();

  c1->SaveAs("plots/"+outFileName_+".png");
}

void plotMacroSampleAll()
{
  plotMacroSample("DeltaPVX", -0.02, 0.02, 0., 1.0, "Difference X-position", "Arbitrary Units", "Position_Difference_X");
  plotMacroSample("DeltaPVY", -0.02, 0.02, 0., 1.0, "Difference Y-position", "Arbitrary Units", "Position_Difference_Y");
  plotMacroSample("DeltaPVZ", -0.02, 0.02, 0., 1.0, "Difference Z-position", "Arbitrary Units", "Position_Difference_Z");
}
  
void plotMacroPhi(TString histname_, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TPad *pad1DEta = new TPad("pad1DEta", "",0.05,0.22,0.96,0.97);
  pad1DEta->Draw();
  pad1DEta->cd();
  pad1DEta->Range(-25,-0.1375,225,1.2375);
  pad1DEta->SetFillColor(0);
  pad1DEta->SetBorderMode(0);
  pad1DEta->SetBorderSize(2);
  pad1DEta->SetFrameBorderMode(0);
  pad1DEta->SetFrameBorderMode(0);
  //pad1DEta->SetLogy();
  //pad1DEta->SetGridy();
  //pad1DEta->SetGridx();

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile *f1 = TFile::Open("GluGluToHToTauTau_M-All_MC_v21.root", "READONLY");
  //TFile *f1 = TFile::Open("VBF_HToTauTau_M-125_MC_v13_vtxWithBS_anal_pt10.root", "READONLY");
  TH1F* hggf = (TH1F*)f1->Get(histname_);
  hggf->Scale(1./hggf->Integral());
  hggf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  hggf->GetYaxis()->SetRangeUser(ymin_, ymax_);
  //hggf->SetMaximum(hggf->GetMaximum()*2.0);
  hggf->SetLineColor(kRed);
  hggf->SetMarkerColor(kRed);
  hggf->SetMarkerSize(0.8);
  hggf->SetMarkerStyle(20);
  hggf->GetXaxis()->SetTitle(xtitle_);
  hggf->GetYaxis()->SetTitle(ytitle_);
  hggf->GetYaxis()->SetTitleOffset(1.4);

  TFile *f2 = TFile::Open("SUSYGluGluToHToTauTau_M-All_MC_v21.root", "READONLY");
  //TFile *f2 = TFile::Open("VBF_HToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS_anal_pt10.root", "READONLY");
  TH1F* hSusyggf = (TH1F*)f2->Get(histname_);
  hSusyggf->SetLineColor(kBlue);
  hSusyggf->SetMarkerColor(kBlue);
  hSusyggf->SetMarkerSize(0.8);
  hSusyggf->SetMarkerStyle(24);
  hSusyggf->Scale(1./hSusyggf->Integral());

  TFile *f3 = TFile::Open("DYJetsToLL_M-50_MC_v21.root", "READONLY");
  TH1F* hdy = (TH1F*)f3->Get(histname_);
  hdy->SetLineColor(kGreen);
  hdy->SetMarkerColor(kGreen);
  hdy->SetMarkerSize(0.8);
  hdy->SetMarkerStyle(22);
  hdy->Scale(1./hdy->Integral());

  hggf->DrawNormalized();
  hSusyggf->DrawNormalized("same");
  hdy->DrawNormalized("same");

  TLegend* leg = new TLegend(0.6,0.7,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hggf,"CP Even GGF H#rightarrow#tau#tau");
  leg->AddEntry(hSusyggf,"CP Odd GGF H#rightarrow#tau#tau");
  ////leg->AddEntry(hggf,"VBF H#rightarrow#tau#tau, With Pol.");
  ////leg->AddEntry(hSusyggf,"VBF H#rightarrow#tau#tau, Without Pol.");
  leg->AddEntry(hdy,"Z#rightarrow#tau#tau");

  leg->Draw();

  pad1DEta->Modified();
  c1->cd();

  TPad* pad2DEta = new TPad("pad2DEta", "",0.05,0.02,0.96,0.22);
  pad2DEta->Draw();
  pad2DEta->cd();
  pad2DEta->Range(-25,-0.2687085,225,0.2687085);
  pad2DEta->SetFillColor(0);
  pad2DEta->SetBorderMode(0);
  pad2DEta->SetBorderSize(2);
  pad2DEta->SetFrameBorderMode(0);
  pad2DEta->SetFrameBorderMode(0);
  pad2DEta->SetGridy();

  TH1F* h_sum_Higgs = hggf->Clone();
  h_sum_Higgs->Reset();
  h_sum_Higgs->Add(hggf, hSusyggf);
  TH1F* h_diff_Higgs = hggf->Clone();
  h_diff_Higgs->Add(hSusyggf, -1.0);
  TH1F* h_diff_dy = hggf->Clone();
  h_diff_dy->Add(hdy, -1.0);

  TH1F* h_Ratio_Higgs = h_sum_Higgs->Clone();
  h_Ratio_Higgs->Reset();
  h_Ratio_Higgs->Divide(h_diff_Higgs, h_sum_Higgs);
  TH1F* h_Ratio_dy = h_sum_Higgs->Clone();
  h_Ratio_dy->Reset();
  h_Ratio_dy->Divide(h_diff_dy, h_sum_Higgs);
  h_Ratio_Higgs->SetLineColor(kBlue);
  h_Ratio_Higgs->SetMarkerColor(kBlue);
  h_Ratio_Higgs->SetMarkerSize(0.8);
  h_Ratio_Higgs->SetMarkerStyle(24);
  h_Ratio_dy->SetLineColor(kGreen);
  h_Ratio_dy->SetMarkerColor(kGreen);
  h_Ratio_dy->SetMarkerSize(0.8);
  h_Ratio_dy->SetMarkerStyle(22);
  h_Ratio_Higgs->GetXaxis()->SetTitle(xtitle_);
  h_Ratio_Higgs->GetYaxis()->SetTitle("Ratio");
  h_Ratio_Higgs->GetYaxis()->SetTitleOffset(0.3);
  h_Ratio_Higgs->SetTitleSize(0.08, "X");
  h_Ratio_Higgs->SetTitleSize(0.15, "Y");
  h_Ratio_Higgs->SetLabelSize(0.12, "X");
  h_Ratio_Higgs->SetLabelSize(0.12, "Y");
  h_Ratio_Higgs->Draw();
  h_Ratio_dy->Draw("same");

  pad2DEta->Modified();
  c1->cd();
  c1->Modified();
  c1->cd();

  c1->SaveAs("plots/"+outFileName_+".png");
  c1->SaveAs("plots/"+outFileName_+".pdf");
}

void plotMacroPhiAll()
{

  //plotMacroPhi("CPPhiStar", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star");
  ////plotMacroPhi("CPPhiLab", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab");
  ////plotMacroPhi("CPPhiStar_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genPiMatch");
  ////plotMacroPhi("CPPhiLab_genMatch", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_genPiMatch");
  ////plotMacroPhi("CPPhiStar_NP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_NP");
  ////plotMacroPhi("CPPhiLab_NP", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_NP");
  ////plotMacroPhi("CPPhiStar_NP_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_NP_genPiMatch");
  ////plotMacroPhi("CPPhiLab_NP_genMatch", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_NP_genPiMatch");
  //plotMacroPhi("CPPhiStar_ipcut", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut50");
  ////plotMacroPhi("CPPhiLab_ipcut", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_ipcut50");
  ////plotMacroPhi("CPPhiStar_NP_ipcut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_NP_ipcut");
  ////plotMacroPhi("CPPhiLab_NP_ipcut", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_NP_ipcut");
  //plotMacroPhi("CPPhiStar_LE1", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_le1");
  ////plotMacroPhi("CPPhiLab_LE1", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_le1");
  //plotMacroPhi("CPPhiStar_ipcut30", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut30");
  //plotMacroPhi("CPPhiStar_ipcut_TkQCut", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut50_tkQuality");
  //plotMacroPhi("CPPhiStar_ipcut30_TkQCut", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut30_tkQuality");
  //plotMacroPhi("CPPhiStar_ipcut_TkQCutV2", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut50_tkQualityTight");
  //plotMacroPhi("CPPhiStar_ipcut30_TkQCutV2", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut30_tkQualityTight");
  //plotMacroPhi("CPPhiStar_ipcut_TkQCutV3", 0., 3.2, 0., 0.18, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ipcut50_tkQualityLoose");
  ////plotMacroPhi("CPPhiStar_ResCut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_ResCut");
  //plotMacroPhi("CPPhiStar_M2", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_m2");
  //plotMacroPhi("CPPhiStar_M2_gen1p0pi0", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_m2_gen1p0pi0");
  //plotMacroPhi("CPPhiStar_LE1_gen1p0pi0", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen1p0pi0");
  //plotMacroPhi("CPPhiStar_tau40", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_tau40");
  ////
  ////plotMacroPhi("CPPhiStar_opv", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv");
  ////plotMacroPhi("CPPhiLab_opv", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_opv");
  ////plotMacroPhi("CPPhiStar_opv_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_genPiMatch");
  ////plotMacroPhi("CPPhiStar_opv_NP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_NP");
  ////plotMacroPhi("CPPhiLab_opv_NP", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_opv_NP");
  ////plotMacroPhi("CPPhiStar_opv_NP_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_NP_genPiMatch");
  //plotMacroPhi("CPPhiStar_opv_ipcut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_ipcut");
  ////plotMacroPhi("CPPhiLab_opv_ipcut", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_opv_ipcut");
  ////plotMacroPhi("CPPhiStar_opv_NP_ipcut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_NP_ipcut");
  ////plotMacroPhi("CPPhiLab_opv_NP_ipcut", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_opv_NP_ipcut");
  //plotMacroPhi("CPPhiStar_opv_LE1", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_le1");
  ////plotMacroPhi("CPPhiLab_opv_LE1", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_opv_le1");
  //
  ////
  ////plotMacroPhi("CPPhiStar_bs", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_bs");
  ////plotMacroPhi("CPPhiLab_bs", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_bs");
  ////plotMacroPhi("CPPhiStar_bs_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_bs_genPiMatch");
  ////plotMacroPhi("CPPhiStar_bs_NP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_bs_NP");
  ////plotMacroPhi("CPPhiLab_bs_NP", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_bs_NP");
  ////plotMacroPhi("CPPhiStar_bs_NP_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_bs_NP_genPiMatch");
  ////
  ////plotMacroPhi("CPPhiStar_genvtx", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx");
  ////plotMacroPhi("CPPhiLab_genvtx", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_genvtx");
  ////plotMacroPhi("CPPhiStar_genvtx_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_genPiMatch");
  ////plotMacroPhi("CPPhiStar_genvtx_NP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_NP");
  ////plotMacroPhi("CPPhiLab_genvtx_NP", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_genvtx_NP");
  ////plotMacroPhi("CPPhiStar_genvtx_NP_genMatch", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_NP_genPiMatch");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut");
  ////plotMacroPhi("CPPhiLab_genvtx_ipcut", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_genvtx_ipcut");
  ////plotMacroPhi("CPPhiStar_genvtx_NP_ipcut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_NP_ipcut");
  ////plotMacroPhi("CPPhiLab_genvtx_NP_ipcut", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_genvtx_NP_ipcut");
  //
  //plotMacroPhi("CPPhiStar_genvtx_LE1", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_le1");
  ////plotMacroPhi("CPPhiLab_genvtx_LE1", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_genvtx_le1");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut30", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut30");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut_TkQCut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut50_tkQuality");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut30_TkQCut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut30_tkQuality");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut30_TkQCutV2", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut30_tkQualityTight");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut_TkQCutV2", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut50_tkQualityTight");
  //plotMacroPhi("CPPhiStar_genvtx_ipcut_TkQCutV3", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ipcut50_tkQualityLoose");
  //plotMacroPhi("CPPhiStar_genvtx_tau40", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_tau40");
  //plotMacroPhi("CPPhiStar_genvtx_gen1p0pi0", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_gen1p0pi0");
  //
  ////plotMacroPhi("CPPhiStar_genvtx_ResCut", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_ResCut");
  //
  //plotMacroPhi("CPPhiStar_tauvtx_LE1", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_tauvtx_le1");
  ////plotMacroPhi("CPPhiLab_tauvtx_LE1", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_tauvtx_le1");
  //
  ////
  ////plotMacroPhi("CPPhiStar_CP", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_cp");
  ////plotMacroPhi("CPPhiStar_opv_CP", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_opv_cp");
  //////plotMacroPhi("CPPhiStar_genvtx_CP", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_cp");
  ////plotMacroPhi("CPPhiStar_CP_LE1", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_cp_le1");
  ////plotMacroPhi("CPPhiStar_opv_CP_LE1", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_opv_cp_le1");
  //////plotMacroPhi("CPPhiStar_genvtx_CP_LE1", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_genvtx_cp_le1");
  ////plotMacroPhi("CPPhiStar_tauvtx_CP_LE1", 0., 6.6, 0., 0.15, "#phi^{*}_{CP} (Radian)", "Arbitrary Units", "cp_phi_star_tauvtx_cp_le1");
  
  
  //plotMacroPhi("xcheck_angle3_reco", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion");
  //plotMacroPhi("xcheck_angle3_reco_opv", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion_opv");
  //plotMacroPhi("xcheck_angle3_reco_bs", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion_bs");
  //plotMacroPhi("xcheck_angle3_reco_genvtx", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion_genvtx");
  //plotMacroPhi("xcheck_angle3_reco_genP4", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion_genP4");
  //plotMacroPhi("xcheck_angle3_reco_genvtx_genP4", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion_genvtx_genP4");
  //plotMacroPhi("xcheck_angle3_reco_genIP", 1.0, 2.1, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_ip_and_pion_genIP");
  //
  //plotMacroPhi("angle_dp_refitVtx", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_refitVtx");
  //plotMacroPhi("angle_dp_origVtx", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_origVtx");
  //plotMacroPhi("angle_dp_genVtx", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_genVtx");
  //plotMacroPhi("angle_dp_refitVtx_genP4", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_refitVtx_genP4");
  //plotMacroPhi("angle_dp_genVtx_genP4", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_genVtx_genP4");
  //plotMacroPhi("angle_dp_genIP", 0., 0.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_genIP");
  //plotMacroPhi("angle_dp_refitVtx_NP", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_refitVtx_NP");
  //plotMacroPhi("angle_dp_genVtx_NP", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_genVtx_NP");
  //
  //plotMacroPhi("angle_dp_refitVtx_LE1", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_refitVtx_le1");
  //plotMacroPhi("angle_dp_origVtx_LE1", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_origVtx_le1");
  //plotMacroPhi("angle_dp_genVtx_LE1", 0., 3.2, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_dp_genVtx_le1");
  
  //plotMacroPhi("DeltaPt_Gen_Reco", 0., 10.0, 0., 0.15, "DeltaP_{T} (GeV)", "Arbitrary Units", "DeltaPt_Gen_Reco");
  //plotMacroPhi("DeltaP_Gen_Reco", 0., 10.0, 0., 0.15, "DeltaP (GeV)", "Arbitrary Units", "DeltaP_Gen_Reco");
  //plotMacroPhi("DeltaPtOverPt_Gen_Reco", 0., 0.1, 0., 0.15, "DeltaP_{T}/P_{T}", "Arbitrary Units", "DeltaPtOverPt_Gen_Reco");
  //plotMacroPhi("DeltaPOverP_Gen_Reco", 0., 0.1, 0., 0.15, "DeltaP/P", "Arbitrary Units", "DeltaPOverP_Gen_Reco");
  //
  //plotMacroPhi("angle_pion_momentum", 0., 0.01, 0., 0.15, "#theta (Radian)", "Arbitrary Units", "angle_pion_momentum");
  
  //plotMacroPhi("CPPhiStar_gen", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen");
  
  
  //plotMacroPhi("CPPhi_gen", 0., 3.2, 0., 0.15, "#phi^{decayPlane} (Radian)", "Arbitrary Units", "cp_phi_gen");
  //plotMacroPhi("CPPhiStar_gen", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen");
  //////plotMacroPhi("CPPhiPrime_gen", 0., 3.2, 0., 0.15, "#phi^{prime} (Radian)", "Arbitrary Units", "cp_phi_prime_gen");
  //////plotMacroPhi("CPPhiLab_gen", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_gen");
  //////plotMacroPhi("CPPhiStar_gen_v2", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_v2");
  //plotMacroPhi("CPPhiStar_gen_vs", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vs");
  //plotMacroPhi("CPPhiStar_gen_vs_ip30", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vs_ip30");
  //plotMacroPhi("CPPhiStar_gen_vs_ip50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vs_ip50");
  //plotMacroPhi("CPPhiStar_gen_vips", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ips");
  //plotMacroPhi("CPPhiStar_gen_vips_ip30", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ips_ip30");
  //plotMacroPhi("CPPhiStar_gen_vips_ip50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ips_ip50");
  //plotMacroPhi("CPPhiStar_gen_vipsv2", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2");
  //plotMacroPhi("CPPhiStar_gen_vipsv2_ip30", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2_ip30");
  //plotMacroPhi("CPPhiStar_gen_vipsv2_ip50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2_ip50");
  //plotMacroPhi("CPPhiStar_gen_pcas", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_pcas");
  //plotMacroPhi("CPPhiStar_gen_pcas_ip30", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_pcas_ip30");
  //plotMacroPhi("CPPhiStar_gen_pcas_ip50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_pcas_ip50");
  //plotMacroPhi("CPPhiStar_gen_vst", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vst");
  //plotMacroPhi("CPPhiStar_gen_vst_ip30", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vst_ip30");
  //plotMacroPhi("CPPhiStar_gen_vst_ip50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vst_ip50");
  //plotMacroPhi("CPPhiStar_gen_ycutP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ycutP");
  //plotMacroPhi("CPPhiStar_gen_ycutN", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ycutN");
  //plotMacroPhi("CPPhiStar_gen_vs_ycutP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vs_ycutP");
  //plotMacroPhi("CPPhiStar_gen_vs_ycutN", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_vs_ycutN");
  //plotMacroPhi("CPPhiStar_gen_vipsv2_ycutP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2_ycutP");
  //plotMacroPhi("CPPhiStar_gen_vipsv2_ycutN", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2_ycutN");
  //plotMacroPhi("CPPhiStar_gen_pcas_ycutP", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_pcas_ycutP");
  //plotMacroPhi("CPPhiStar_gen_pcas_ycutN", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_pcas_ycutN");
  plotMacroPhi("CPPhiStar_gen_ycutP_visPt50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ycutP_visPt50");
  plotMacroPhi("CPPhiStar_gen_ycutN_visPt50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ycutN_visPt50");
  plotMacroPhi("CPPhiStar_gen_vipsv2_ycutP_visPt50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2_ycutP_visPt50");
  plotMacroPhi("CPPhiStar_gen_vipsv2_ycutN_visPt50", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_ipsv2_ycutN_visPt50");

  /*
  //plotMacroPhi("DeltaX_gen_reco_pca_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA X - Reco PCA X", "Arbitrary Units", "DeltaX_gen_reco_pca_genvtx");
  //plotMacroPhi("DeltaY_gen_reco_pca_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA Y - Reco PCA Y", "Arbitrary Units", "DeltaY_gen_reco_pca_genvtx");
  //plotMacroPhi("DeltaZ_gen_reco_pca_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA Z - Reco PCA Z", "Arbitrary Units", "DeltaZ_gen_reco_pca_genvtx");

  plotMacroPhi("DeltaX_gen_reco_pcaLE1_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA X - Reco PCA X", "Arbitrary Units", "DeltaX_gen_reco_pcaLE1_genvtx");
  plotMacroPhi("DeltaY_gen_reco_pcaLE1_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA Y - Reco PCA Y", "Arbitrary Units", "DeltaY_gen_reco_pcaLE1_genvtx");
  plotMacroPhi("DeltaZ_gen_reco_pcaLE1_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA Z - Reco PCA Z", "Arbitrary Units", "DeltaZ_gen_reco_pcaLE1_genvtx");

  plotMacroPhi("DeltaX_gen_reco_pcaLE2_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA X - Reco PCA X", "Arbitrary Units", "DeltaX_gen_reco_pcaLE2_genvtx");
  plotMacroPhi("DeltaY_gen_reco_pcaLE2_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA Y - Reco PCA Y", "Arbitrary Units", "DeltaY_gen_reco_pcaLE2_genvtx");
  plotMacroPhi("DeltaZ_gen_reco_pcaLE2_genvtx", -0.01, 0.01, 0., 0.15, "Gen PCA Z - Reco PCA Z", "Arbitrary Units", "DeltaZ_gen_reco_pcaLE2_genvtx");
  */
  /*
  plotMacroPhi("DeltaTkLength_gen_reco_genvtx_le1_l", -0.002, 0.002, 0., 0.15, "Delta(PCA - genTauDecayVertex) (parallel to gen tau)", "Arbitrary Units", "DeltaTkLength_gen_reco_genvtx_le1_l");
  plotMacroPhi("DeltaTkLength_gen_reco_genvtx_le1_t", -0.01, 0.01, 0., 0.15, "Delta(PCA - genTauDecayVertex) (perpendicular to gen tau)", "Arbitrary Units", "DeltaTkLength_gen_reco_genvtx_le1_t");
  plotMacroPhi("DeltaTkLength_gen_reco_genvtx_le2_l", -0.002, 0.002, 0., 0.15, "Delta(PCA - genTauDecayVertex) (parallel to gen tau)", "Arbitrary Units", "DeltaTkLength_gen_reco_genvtx_le2_l");
  plotMacroPhi("DeltaTkLength_gen_reco_genvtx_le2_t", -0.01, 0.01, 0., 0.15, "Delta(PCA - genTauDecayVertex) (perpendicular to gen tau)", "Arbitrary Units", "DeltaTkLength_gen_reco_genvtx_le2_t");
  */
  /*
  plotMacroPhi("AngleTkLength_gen_reco_genvtx_le1_l", 0., 0.05, 0., 0.15, "Angle(PCA - genTauDecayVertex) (parallel to gen tau)", "Arbitrary Units", "AngleTkLength_gen_reco_genvtx_le1_l");
  plotMacroPhi("AngleTkLength_gen_reco_genvtx_le1_t", 0., 4.0, 0., 0.15, "Angle(PCA - genTauDecayVertex) (perpendicular to gen tau)", "Arbitrary Units", "AngleTkLength_gen_reco_genvtx_le1_t");
  plotMacroPhi("AngleTkLength_gen_reco_genvtx_le2_l", 0., 0.05, 0., 0.15, "Angle(PCA - genTauDecayVertex) (parallel to gen tau)", "Arbitrary Units", "AngleTkLength_gen_reco_genvtx_le2_l");
  plotMacroPhi("AngleTkLength_gen_reco_genvtx_le2_t", 0., 4.0, 0., 0.15, "Angle(PCA - genTauDecayVertex) (perpendicular to gen tau)", "Arbitrary Units", "AngleTkLength_gen_reco_genvtx_le2_t");
  */
  /*
  plotMacroPhi("DeltaTkLength_gen_reco_le1_l", -0.01, 0.01, 0., 0.15, "Delta(PCA - genTauDecayVertex) (parallel to gen tau)", "Arbitrary Units", "DeltaTkLength_gen_reco_le1_l");
  plotMacroPhi("DeltaTkLength_gen_reco_le1_t", -0.01, 0.01, 0., 0.15, "Delta(PCA - genTauDecayVertex) (perpendicular to gen tau)", "Arbitrary Units", "DeltaTkLength_gen_reco_le1_t");
  plotMacroPhi("AngleTkLength_gen_reco_le1_l", 0., 0.05, 0., 0.15, "Angle(PCA - genTauDecayVertex) (parallel to gen tau)", "Arbitrary Units", "AngleTkLength_gen_reco_le1_l");
  plotMacroPhi("AngleTkLength_gen_reco_le1_t", 0., 4.0, 0., 0.15, "Angle(PCA - genTauDecayVertex) (perpendicular to gen tau)", "Arbitrary Units", "AngleTkLength_gen_reco_le1_t");
  */
  /*
  plotMacroPhi("TkLength_gen_reco_genvtx_le1_tdp", -0.01, 0.01, 0., 0.15, "(PCA - genTauDecayVertex) (perpendicular to gen tau decayPlane)", "Arbitrary Units", "IPLength_gen_reco_genvtx_le1_tdp");
  plotMacroPhi("TkLength_gen_reco_le1_tdp", -0.01, 0.01, 0., 0.15, "(PCA - genTauDecayVertex) (perpendicular to gen tau decayPlane)", "Arbitrary Units", "IPLength_gen_reco_le1_tdp");
  */
}

void plotMacroPhi_PolOff(TString histname_, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  //TFile *f1 = TFile::Open("DYJetsToLL_M-50_MC_v13_vtxWithBS_anal.root", "READONLY");
  //TFile *f1 = TFile::Open("GluGluToHToTauTau_M-125_MC_v13_vtxWithBS_anal.root", "READONLY");
  TFile *f1 = TFile::Open("VBF_HToTauTau_M-125_MC_v13_vtxWithBS_anal.root", "READONLY");
  TH1F* hggf = (TH1F*)f1->Get(histname_);
  hggf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  //hggf->GetYaxis()->SetRangeUser(ymin_, ymax_);
  hggf->SetMaximum(hggf->GetMaximum()*2.0);
  hggf->SetLineColor(kRed);
  hggf->GetXaxis()->SetTitle(xtitle_);
  hggf->GetYaxis()->SetTitle(ytitle_);

  //TFile *f2 = TFile::Open("DYJetsToLL_M-50_tauPolarOff_MC_v13_vtxWithBS_anal.root", "READONLY");
  //TFile *f2 = TFile::Open("GluGluToHToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS_anal.root", "READONLY");
  TFile *f2 = TFile::Open("VBF_HToTauTau_M-125_tauPolarOff_MC_v13_vtxWithBS_anal.root", "READONLY");
  TH1F* hggf_noPol = (TH1F*)f2->Get(histname_);
  hggf_noPol->SetLineColor(kBlue);

  hggf->DrawNormalized();
  hggf_noPol->DrawNormalized("same");
  
  TLegend* leg = new TLegend(0.6,0.7,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  //leg->AddEntry(hggf,"GGF H#rightarrow#tau#tau, With Tau Pol.");
  //leg->AddEntry(hggf_noPol,"GGF H#rightarrow#tau#tau, Without Tau Pol.");
  leg->AddEntry(hggf,"VBF H#rightarrow#tau#tau, With Pol.");
  leg->AddEntry(hggf_noPol,"VBF H#rightarrow#tau#tau, Without Pol.");
  //leg->AddEntry(hggf,"DY#rightarrow#tau#tau, With Pol.");
  //leg->AddEntry(hggf_noPol,"DY#rightarrow#tau#tau, Without Pol.");

  leg->Draw();

  c1->SaveAs("plots/NoSpin/"+outFileName_+".png");
}

void plotMacroPhiAll_PolOff()
{
  plotMacroPhi_PolOff("CPPhi_gen", 0., 3.2, 0., 0.15, "#phi^{decayPlane} (Radian)", "Arbitrary Units", "cp_phi_gen_NoTauPol_VBF");
  plotMacroPhi_PolOff("CPPhiStar_gen", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_NoTauPol_VBF");
  //plotMacroPhi_PolOff("CPPhiPrime_gen", 0., 3.2, 0., 0.15, "#phi^{prime} (Radian)", "Arbitrary Units", "cp_phi_prime_gen_NoTauPol");
  plotMacroPhi_PolOff("CPPhiLab_gen", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_gen_NoTauPol_VBF");
  //plotMacroPhi_PolOff("CPPhiStar_gen_v2", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_gen_NoTauPol_v2");
}

void plotMacroPhiXCheck(TString hist1_, TString hist2_, TString hist3_, TString hist4_, float xmin_, float xmax_,TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile *f1 = TFile::Open("GluGluToHToTauTau_M-120_MC_v13_vtxWithBS_anal.root", "READONLY");
  TH1F* hggf1 = (TH1F*)f1->Get(hist1_);
  hggf1->GetXaxis()->SetRangeUser(xmin_, xmax_);
  hggf1->SetMaximum(hggf1->GetMaximum()*1.6);
  hggf1->SetLineColor(kRed);
  hggf1->GetXaxis()->SetTitle(xtitle_);
  hggf1->GetYaxis()->SetTitle(ytitle_);

  TH1F* hggf2 = (TH1F*)f1->Get(hist2_);
  hggf2->SetLineColor(kBlue);
  TH1F* hggf3 = (TH1F*)f1->Get(hist3_);
  hggf3->SetLineColor(kGreen);
  TH1F* hggf4 = (TH1F*)f1->Get(hist4_);
  hggf4->SetLineColor(kViolet);

  hggf1->DrawNormalized();
  hggf2->DrawNormalized("same");
  hggf3->DrawNormalized("same");
  hggf4->DrawNormalized("same");

  TLegend* leg = new TLegend(0.3,0.7,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hggf1,"SM GGF H#rightarrow#tau#tau, IP in XY");
  leg->AddEntry(hggf2,"SM GGF H#rightarrow#tau#tau, IP in XY with normal to pi");
  leg->AddEntry(hggf3,"SM GGF H#rightarrow#tau#tau, IP linearization method1");
  leg->AddEntry(hggf4,"SM GGF H#rightarrow#tau#tau, IP linearization method2");

  leg->Draw();

  c1->SaveAs("plots/"+outFileName_+".png");
}

void plotMacroPhiXCheckAll(){

  plotMacroPhiXCheck("CPPhiStar_genvtx", "CPPhiStar_genvtx_NP", "CPPhiStar_genvtx_LE1", "CPPhiStar_genvtx_LE2", 0, 3.2, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_linearization_ggf");
  plotMacroPhiXCheck("CPPhiLab_genvtx", "CPPhiLab_genvtx_NP", "CPPhiLab_genvtx_LE1", "CPPhiLab_genvtx_LE2", 0, 3.2,"#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_linearization_ggf");
  plotMacroPhiXCheck("CPPhiStar_genvtx_genMatch", "CPPhiStar_genvtx_NP_genMatch", "CPPhiStar_genvtx_LE1_genMatch", "CPPhiStar_genvtx_LE2_genMatch", 0, 3.2,"#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_linearization_ggf_genMatch");
}


