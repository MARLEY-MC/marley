#include <cmath>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

void cos_plot(const std::string& file_name) {

  const int NUM_BINS = 10;
  const double COS_MIN = -1.;
  const double COS_MAX = 1.;
  const double BIN_WIDTH = ( COS_MAX - COS_MIN ) / NUM_BINS;

  TH1D* cos_theta_hist = new TH1D("cos_theta_hist",
    "scattering cosine distribution; cos#theta; #left[ d#sigma"
    "/dcos#theta #right]_{flux} (10^{-42} cm^{2})", NUM_BINS,
    COS_MIN, COS_MAX);

  TFile* f = TFile::Open(file_name.c_str());
  TTree* t = static_cast<TTree*>(f->Get("mst"));

  double pxv, pyv, pzv, pxl, pyl, pzl, xsec;
  t->SetBranchAddress("pxv", &pxv);
  t->SetBranchAddress("pyv", &pyv);
  t->SetBranchAddress("pzv", &pzv);
  t->SetBranchAddress("pxl", &pxl);
  t->SetBranchAddress("pyl", &pyl);
  t->SetBranchAddress("pzl", &pzl);
  t->SetBranchAddress("xsec", &xsec);

  Long64_t num_events = t->GetEntries();
  for (Long64_t i = 0; i < num_events; ++i) {
    t->GetEntry(i);
    if (i % 1000 == 0) std::cout << "Event " << i << '\n';

    double pnu_dot_pe = pxv * pxl + pyv * pyl + pzv * pzl;
    double norm_pnu = std::sqrt(pxv * pxv + pyv * pyv + pzv * pzv);
    double norm_pe = std::sqrt(pxl * pxl + pyl * pyl + pzl * pzl);

    double cos_theta = pnu_dot_pe / (norm_pnu * norm_pe);
    cos_theta_hist->Fill(cos_theta);
  }

  double scale_factor = xsec / BIN_WIDTH / num_events;
  cos_theta_hist->Scale(scale_factor);

  TCanvas* c = new TCanvas;
  c->cd();

  cos_theta_hist->GetXaxis()->SetTitleOffset(1.2);
  cos_theta_hist->GetYaxis()->SetTitleOffset(1.2);
  cos_theta_hist->SetStats(false);
  cos_theta_hist->SetLineColor(kBlue);
  cos_theta_hist->SetLineWidth(2);
  cos_theta_hist->Draw("hist");
}
