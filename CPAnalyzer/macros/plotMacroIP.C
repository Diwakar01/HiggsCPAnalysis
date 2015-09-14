void plotMacroIP(TString histname_, float xmin_, float xmax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy();

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

  TFile *f1 = TFile::Open("GluGluToHToTauTau_M-All_MC_v18.root", "READONLY");
  //TFile *f1 = TFile::Open("CPWthoutPi0s/GluGluToHToTauTau_M-All_MC_v25_noBS.root", "READONLY");
  TH1F* hggf = (TH1F*)f1->Get(histname_);
  hggf->Scale(1./hggf->Integral());
  hggf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  hggf->SetMaximum(hggf->GetMaximum()*100.0);
  //hggf->GetYaxis()->SetRangeUser(0.0001, 10.);
  hggf->SetLineColor(1);
  hggf->SetMarkerSize(0.5);
  hggf->SetMarkerColor(1);
  hggf->SetMarkerStyle(20);
  hggf->GetXaxis()->SetTitle(xtitle_);
  hggf->GetYaxis()->SetTitle(ytitle_);
  hggf->GetYaxis()->SetTitleOffset(1.3);
  hggf->GetXaxis()->SetTitleOffset(1.25);

  hggf->Draw();
  //TF1 *func = new TF1("func", "gaus", -0.015, 0.015); //for DeltaRfPVwrtGenX,Y
  //TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.03, 0.03); //for DeltaRfPVwrtGenZ, DeltaIanPVwrtGenZ, DeltaPVwrtGenZ
  //TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.025, 0.025); //DeltaTkLength_gen_reco_le1_l
  TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.025, 0.025);  //DeltaTkLength_gen_reco_le1_t, TkLength_gen_reco_le1_tdp
  //TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.02, 0.02); //DeltaTkLength_gen_reco_ip50_l_tkQCut
  //TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.02, 0.02); //DeltaTkLength_gen_reco_ip30_l_tkQCut
  //TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.0015, 0.0015); //DeltaTkLength_gen_reco_genvtx_le1_l

  TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.02, 0.02); //for DeltaRfPVwrtGenX,Y woBS

  //hggf->Fit("gaus", "Q", "", -0.015, 0.015);
  //TF1 *fcn = hggf->GetFunction("gaus");
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.005); //for DeltaRfPVwrtGenZ
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.005); //DeltaIanPVwrtGenZ, DeltaPVwrtGenZ
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.004); //DeltaTkLength_gen_reco_le1_l, DeltaTkLength_gen_reco_ip30_l_tkQCut
  func->SetParameters(0.1, 0.001, 0.001, 0.001, 0.001, 0.002); //DeltaTkLength_gen_reco_le1_t, TkLength_gen_reco_le1_tdp
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.0045); //DeltaTkLength_gen_reco_ip30_l_tkQCut
  //func->SetParameters(0.1, 0.001, 0.001, 0.001, 0.002, 0.005); //DeltaTkLength_gen_reco_genvtx_le1_t
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.003); //TkLength_gen_reco_genvtx_ip30_tdp_tkQCut
  //func->SetParameters(0.1, 0.001, 0.0005, 0.0005, 0.0005, 0.0005); ////DeltaTkLength_gen_reco_genvtx_le1_l

  //func->SetParameters(0.1, 0.001, 0.001, 0.001, 0.001, 0.004); //for DeltaRfPVwrtGenX,Y woBS
  //func->SetParameters(10000000, 0.001, 0.001, 0.001, 0.001, 0.004);

  hggf->Fit("func", "R");
  TF1 *fcn = hggf->GetFunction("func");

  TLegend* leg = new TLegend(0.5,0.7,0.8,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->AddEntry(hggf, Form("MC: RMS = %1.2e", hggf->GetRMS()), "p");
  leg->AddEntry(fcn, Form("Gaussian fit: #sigma_{1} = %1.2e", fcn->GetParameter(2)), "l");
  leg->AddEntry(fcn, Form("Gaussian fit: #sigma_{2} = %1.2e", fcn->GetParameter(5)), "l");
  leg->Draw();

  //std::cout<<"Integral "<<fcn->Integral(0, 0.2)<<endl;
  //c1->SaveAs("plots/"+outFileName_+".png");
  //c1->SaveAs("plots/"+outFileName_+".pdf");
}

plotMacroIPAll()
{
  
  //plotMacroIP("DeltaTkLength_gen_reco_le1_l", -0.025, 0.025, "#delta PCA #parallel gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_le1_l");
  plotMacroIP("DeltaTkLength_gen_reco_le1_t", -0.025, 0.025, "#delta PCA #\perp gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_le1_t");
  //plotMacroIP("TkLength_gen_reco_le1_tdp", -0.025, 0.025, "#delta PCA #\perp gen. tau decay plane [cm]", "Events", "IPLength_gen_reco_le1_tdp");

  //plotMacroIP("DeltaTkLength_gen_reco_ip50_l_tkQCut", -0.025, 0.025, "#delta PCA #parallel gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_ip50_l_tkQCut");
  //plotMacroIP("DeltaTkLength_gen_reco_ip50_t_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_ip50_t_tkQCut");
  //plotMacroIP("TkLength_gen_reco_ip50_tdp_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau decay plane [cm]", "Events", "IPLength_gen_reco_ip50_tdp_tkQCut");
  
  //plotMacroIP("DeltaTkLength_gen_reco_ip30_l_tkQCut", -0.025, 0.025, "#delta PCA #parallel gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_ip30_l_tkQCut");
  //plotMacroIP("DeltaTkLength_gen_reco_ip30_t_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_ip30_t_tkQCut");
  //plotMacroIP("TkLength_gen_reco_ip30_tdp_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau decay plane [cm]", "Events", "IPLength_gen_reco_ip30_tdp_tkQCut");

  //plotMacroIP("DeltaTkLength_gen_reco_genvtx_le1_l", -0.0015, 0.0015, "#delta PCA #parallel gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_genvtx_le1_l");
  //plotMacroIP("DeltaTkLength_gen_reco_genvtx_le1_t", -0.025, 0.025, "#delta PCA #\perp gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_genvtx_le1_t");
  //plotMacroIP("TkLength_gen_reco_genvtx_le1_tdp", -0.025, 0.025, "#delta PCA #\perp gen. tau decay plane [cm]", "Events", "IPLength_gen_reco_genvtx_le1_tdp");
  
  //plotMacroIP("DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_genvtx_ip50_t_tkQCut");
  //plotMacroIP("TkLength_gen_reco_genvtx_ip50_tdp_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau decay plane [cm]", "Events", "IPLength_gen_reco_genvtx_ip50_tdp_tkQCut");
  
  //plotMacroIP("DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau direction [cm]", "Events", "DeltaTkLength_gen_reco_genvtx_ip30_t_tkQCut");
  //plotMacroIP("TkLength_gen_reco_genvtx_ip30_tdp_tkQCut", -0.025, 0.025, "#delta PCA #\perp gen. tau decay plane [cm]", "Events", "IPLength_gen_reco_genvtx_ip30_tdp_tkQCut");
  

  //plotMacroIP("DeltaRfPVwrtGenX", -0.012, 0.012, "X_{vtx}^{gen} - X_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_X");
  //plotMacroIP("DeltaRfPVwrtGenY", -0.012, 0.012, "Y_{vtx}^{gen} - Y_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Y");
  //plotMacroIP("DeltaRfPVwrtGenZ", -0.045, 0.045, "Z_{vtx}^{gen} - Z_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Z");

  //plotMacroIP("DeltaRfPVwrtGenX", -0.045, 0.045, "X_{vtx}^{gen} - X_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_X_woBS");
  //plotMacroIP("DeltaRfPVwrtGenY", -0.045, 0.045, "Y_{vtx}^{gen} - Y_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Y_woBS");
  //plotMacroIP("DeltaRfPVwrtGenZ", -0.045, 0.045, "Z_{vtx}^{gen} - Z_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Z_woBS");

  //plotMacroIP("DeltaPVwrtGenX", -0.012, 0.012, "X_{vtx}^{gen} - X_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_X_opv");
  //plotMacroIP("DeltaPVwrtGenY", -0.012, 0.012, "Y_{vtx}^{gen} - Y_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Y_opv");
  //plotMacroIP("DeltaPVwrtGenZ", -0.045, 0.045, "Z_{vtx}^{gen} - Z_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Z_opv");
  
  //plotMacroIP("DeltaIanPVwrtGenX", -0.012, 0.012, "X_{vtx}^{gen} - X_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_X_ian");
  //plotMacroIP("DeltaIanPVwrtGenY", -0.012, 0.012, "Y_{vtx}^{gen} - Y_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Y_ian");
  //plotMacroIP("DeltaIanPVwrtGenZ", -0.045, 0.045, "Z_{vtx}^{gen} - Z_{vtx}^{rec}", "Events", "GGFH_PosDiffWrtGen_Z_ian");
}


void plotMacroVtx(TString histname_, float xmin_, float xmax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy();

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

  TFile *f1 = TFile::Open("GluGluToHToTauTau_M-All_MC_v18.root", "READONLY");
  TH1F* hggf = (TH1F*)f1->Get(histname_);
  hggf->Scale(1./hggf->Integral());
  hggf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  hggf->SetMaximum(hggf->GetMaximum()*100.0);
  //hggf->GetYaxis()->SetRangeUser(0.0001, 10.);
  hggf->SetLineColor(1);
  hggf->SetMarkerSize(0.5);
  hggf->SetMarkerColor(1);
  hggf->SetMarkerStyle(20);
  hggf->GetXaxis()->SetTitle(xtitle_);
  hggf->GetYaxis()->SetTitle(ytitle_);
  hggf->GetYaxis()->SetTitleOffset(1.3);
  hggf->GetXaxis()->SetTitleOffset(1.25);

  hggf->Draw();
  TF1 *func = new TF1("func", "gaus", -0.015, 0.015); //for DeltaRfPVwrtGenX,Y
  //TF1 *func = new TF1("func", "gaus(0)+gaus(3)", -0.04, 0.04);
  //hggf->Fit("gaus", "Q", "", -0.015, 0.015);
  //TF1 *fcn = hggf->GetFunction("gaus");
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.02); //for DeltaRfPVwrtGenZ
  //func->SetParameters(0.1, 0.001, 0.002, 0.001, 0.001, 0.01); //DeltaIanPVwrtGenZ, DeltaPVwrtGenZ

  //hggf->Fit("func", "R");
  //TF1 *fcn = hggf->GetFunction("func");

  TLegend* leg = new TLegend(0.5,0.75,0.8,0.85,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->AddEntry(hggf, Form("MC: MEAN = %1.2e", hggf->GetMean()), "p");
  leg->AddEntry(hggf, Form("MC: RMS = %1.2e", hggf->GetRMS()), "p");
  //leg->AddEntry(fcn, Form("Gaussian fit: mean = %1.2e", fcn->GetParameter(1)), "l");
  //leg->AddEntry(fcn, Form("Gaussian fit: #sigma = %1.2e", fcn->GetParameter(2)), "l");
  //leg->AddEntry(fcn, Form("Gaussian fit: #sigma_{2} = %1.2e", fcn->GetParameter(5)), "l");
  leg->Draw();

  c1->SaveAs("plots/"+outFileName_+".png");
  //c1->SaveAs("plots/"+outFileName_+".pdf");
}

plotMacroVtxAll()
{

  
  plotMacroVtx("DeltaPVX", -0.012, 0.012, "X_{vtx}^{orig} - X_{vtx}^{refit}", "Events", "GGFH_PosDiffWrtOpv_X");
  plotMacroVtx("DeltaPVY", -0.012, 0.012, "Y_{vtx}^{orig} - Y_{vtx}^{refit}", "Events", "GGFH_PosDiffWrtOpv_Y");
  plotMacroVtx("DeltaPVZ", -0.045, 0.045, "Z_{vtx}^{orig} - Z_{vtx}^{refit}", "Events", "GGFH_PosDiffWrtOpv_Z");
}
