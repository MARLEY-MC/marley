#ifndef USE_ROOT
#error "This MARLEY validation code currently must be compiled with ROOT support."
#endif

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TTree.h"
#include "TH1D.h"

#include "marley_utils.hh"
#include "TMarleyEvent.hh"

double fermi_dirac_distribution(double C, bool e_flavor, bool anti, double nu_energy){
  double eta = 0;
  double T = 0;
  double N_nu;

  if(e_flavor && !anti){
    T = 3.5; // temperature in MeV
    N_nu = 2.8;  // total number of electron neutrinos expected (x10^57)
  }
  else if(e_flavor && anti) {
    T = 5.0; // temperature in MeV
    N_nu = 1.9;  // total number of electron anti-neutrinos expected (x10^57)
  }
  else { // !e_flavor
    T = 8.0; // temperature in MeV
    N_nu = 5.0; // total number of mu+tau neutrinos + anti-neutrinos expected (x10^57)
  }
  return (C/std::pow(T,3))*(std::pow(nu_energy,2)/(1+std::exp(nu_energy/(T-eta))))*N_nu;
}

// TODO: document this based on your rough notes
double chi_squared_p_value(double chi_squared, unsigned degrees_of_freedom,
  unsigned num_terms = 50)
{
  double chi_squared_over_two = chi_squared/2;
  double dof_over_two = static_cast<double>(degrees_of_freedom)/2.0;
  double prefactor = std::pow(chi_squared_over_two, dof_over_two)
    * std::exp(-chi_squared_over_two) / std::tgamma(dof_over_two);

  double numerator = 1.0;
  double denominator = dof_over_two;
  double sum = numerator/denominator;
  for (unsigned n = 1; n < num_terms; ++n) {
    numerator *= chi_squared_over_two;
    denominator *= (dof_over_two + n);
    sum += numerator/denominator;
  }

  return 1.0 - prefactor*sum;
}

// Numerically integrate a given function f (that takes a
// double argument to integrate over and returns a double)
// over the interval [a,b] using Clenshaw-Curtis quadrature
// at 2N sampling points.
// (see http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature)
double cc_num_integrate(const std::function<double(double)> &f,
  double a, double b, int N)
{
  double twoN = 2*N;
  double A = (b - a)/2;
  double B = (b + a)/2;
  double C = (f(b) + f(a))/twoN;
  double D = f(B)/N;

  double integral = 0;

  for(int k = 0; k <= N; ++k) {

    double npi = 0;
    int sign = 1;
    if (k % 2 == 1) sign = -1;
    double term = C + sign*D;
    for (int n = 1; n < N; ++n) {
      npi += marley_utils::pi;

      term += (f(B + A*std::cos(npi/twoN))
       + f(B - A*std::cos(npi/twoN)))*std::cos(npi*k/N)/N;
    }

    if (k != 0 && k != N) term *= 2;

    term /= 1 - 4*k*k;

    integral += term;
  }

  return A*integral;
}

// Numerically integrate a given function f (that takes a
// double argument to integrate over and returns a double)
// over the interval [a,b] using Clenshaw-Curtis quadrature
// at 2N sampling points.
// (see http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature)
double cc_num_integrate2(const std::function<double(double)> &f,
  double a, double b, int N)
{
  double twoN = 2*N;
  double A = (b - a)/2;
  double B = (b + a)/2;
  double C = (f(b) + f(a))/twoN;
  double D = f(B)/N;

  std::vector<double> eval_points;
  double npi = 0;
  for (int n = 1; n < N; ++n) {
    npi += marley_utils::pi;
    eval_points.push_back(A*std::cos(npi/twoN));
  }

  double integral = 0;

  for(int k = 0; k <= N; ++k) {

    npi = 0;
    int sign = 1;
    if (k % 2 == 1) sign = -1;
    double term = C + sign*D;
    for (int n = 1; n < N; ++n) {
      npi += marley_utils::pi;
      double epoint = eval_points[n-1];

      term += (f(B + epoint) + f(B - epoint))
        * std::cos(npi*k/N)/N;
    }

    if (k != 0 && k != N) term *= 2;

    term /= 1 - 4*k*k;

    integral += term;
  }

  return A*integral;
}

inline double PoissonPMF(double mean, int x) {
  if (mean < 50) {
    return std::exp(-mean) * std::pow(mean, x)
      / std::tgamma(x + 1);
  }
  // Approximate the PMF using a normal distribution PDF
  else {
    return 0.5 * (std::erf((x - mean + 0.5) / std::sqrt(2*mean))
      - std::erf((x - mean - 0.5) / std::sqrt(2*mean)));
  }
}

void find_prob_limits(double mean_counts, double prob_enclosed,
  double& prob_previous, int& count_low_limit, int& count_up_limit)
{
  //std::cout << "    pe = " << prob_enclosed << " pp = " << prob_previous << std::endl;
  if (prob_enclosed < 0 || prob_enclosed >= 1) throw std::runtime_error(
    std::string("Invalid enclosed probability of ") + std::to_string(prob_enclosed)
    + " passed to find_prob_limits()");
  while (prob_previous < prob_enclosed) {
    double prob_low = 0;
    if (count_low_limit > 0) prob_low = PoissonPMF(mean_counts,
      count_low_limit - 1);
    double prob_up = PoissonPMF(mean_counts, count_up_limit + 1);
    if (prob_low < prob_up) {
      prob_previous += prob_up;
      ++count_up_limit;
    }
    else if (prob_low > prob_up) {
      prob_previous += prob_low;
      --count_low_limit;
    }
    else { // prob_low == prob_up
      prob_previous += prob_up + prob_low;
      ++count_up_limit;
      --count_low_limit;
    }
  }
}

int main() {

  // Check if the ROOT file used to store the neutrino
  // energy histogram already exists
  std::string tree_file_name("n_energies.root");
  std::ifstream test_stream(tree_file_name.c_str());
  bool root_file_existed = false;
  if (test_stream.good()) {
    root_file_existed = true;
  }
  test_stream.close();

  TH1D* hist = nullptr;

  if (!root_file_existed) {

    TFile f2("n_energies.root", "CREATE");
    hist = new TH1D("ne_hist", "neutrino energies", 100, 4.36, 50);

    TFile f1("sn_1e6_events.root", "read");
    TTree* t = nullptr;
    f1.GetObject("event_tree", t);
    TMarleyEvent* e = new TMarleyEvent;
    t->GetBranch("events")->SetAddress(&e);
    int num_events = t->GetEntries();
    for (int i = 0; i < num_events; ++i) {
      t->GetEntry(i);
      hist->Fill(e->get_projectile()->get_total_energy());
      if (i % 10000 == 0) std::cout << "Read " << i << " events" << std::endl;
    }

    f2.cd();
    hist->Write();
    f1.Close();
    f2.Close();
  }

  TH1D* hist2 = nullptr;
  TFile f3("n_energies.root", "READ");
  f3.GetObject("ne_hist", hist2);

  int num_events = hist2->GetEntries();

  std::vector<double> d_bin_low_edges;
  std::vector<double> d_bin_up_edges;
  std::vector<double> d_bin_counts;

  for (int i = 0, n = hist2->GetNbinsX() + 1; i <= n; ++i) {
    d_bin_low_edges.push_back(hist2->GetXaxis()->GetBinLowEdge(i));
    d_bin_up_edges.push_back(hist2->GetXaxis()->GetBinUpEdge(i));
    d_bin_counts.push_back(hist2->GetBinContent(i));
  }

  if (d_bin_counts.front() > 0) throw std::runtime_error(std::string("Events found ")
    + "in underflow bin. Aborting analysis.");
  if (d_bin_counts.back() > 0) throw std::runtime_error(std::string("Events found ")
    + "in overflow bin. Aborting analysis.");

  // Erase the first element of each of these vectors (the underflow bin)
  d_bin_low_edges.erase(d_bin_low_edges.begin(), d_bin_low_edges.begin() + 1);
  d_bin_up_edges.erase(d_bin_up_edges.begin(), d_bin_up_edges.begin() + 1);
  d_bin_counts.erase(d_bin_counts.begin(), d_bin_counts.begin() + 1);
  // Erase the last element of each of these vectors (the overflow bin)
  d_bin_low_edges.pop_back();
  d_bin_up_edges.pop_back();
  d_bin_counts.pop_back();

  // Number of subdivisions to use when numerically integrating
  // using the trapezoid rule
  const int num_divisions = 50;
  // Minimum and maximum neutrino energy
  double minE = 4.36, maxE = 50;

  double chi_squared = 0.0;
  double expected_counts = 0.0;

  // Neutrino energies were sampled using electron flavor supernova neutrinos for C = 0.55
  std::function<double(double)> func = std::bind(fermi_dirac_distribution, 0.55,
    true, false, std::placeholders::_1);

  // Adjust the supernova neutrino spectrum by this factor so that we have a properly
  // normalized probability density function to use for our Pearson chi-squared test
  double normalization_factor = cc_num_integrate(func, minE, maxE, num_divisions);

  std::function<double(double)> pdf = [&func, &normalization_factor](double x)
    -> double { return func(x)/normalization_factor; };

  unsigned n_bins = hist2->GetNbinsX();

  std::vector<double> d_green_els, d_green_eus, d_yellow_els,
    d_yellow_eus, d_red_els, d_red_eus;

  std::vector<double> d_nbin_zeros(n_bins, 0);

  // Fill a vector with the widths of each bin
  std::vector<double> d_bin_widths;
  std::transform(d_bin_up_edges.begin(), d_bin_up_edges.end(),
    d_bin_low_edges.begin(), std::back_inserter(d_bin_widths),
    [](double up, double low) -> double { return up - low; });

  std::vector<double> d_bin_mean_counts;

  for (unsigned i = 0; i < n_bins; ++i) {

    expected_counts = num_events*(cc_num_integrate(
      pdf, d_bin_low_edges.at(i), d_bin_up_edges.at(i), num_divisions));

    d_bin_mean_counts.push_back(expected_counts);

    double observed_counts = d_bin_counts.at(i);

    chi_squared += std::pow((observed_counts - expected_counts), 2)
      / (expected_counts);

    // Calculate band limits for Brazil plot using the procedure suggested
    // in http://arxiv.org/abs/1112.2593
    const double prob_green = 0.68, prob_yellow = 0.95, prob_red = 0.99;
    int most_probable_counts = static_cast<int>(std::floor(expected_counts));
    int count_low_limit = most_probable_counts;
    int count_up_limit = most_probable_counts;
    double prob = PoissonPMF(expected_counts, most_probable_counts);

    find_prob_limits(expected_counts, prob_green, prob, count_low_limit, count_up_limit);
    d_green_els.push_back(expected_counts - (count_low_limit - 0.5));
    d_green_eus.push_back((count_up_limit + 0.5) - expected_counts);
    //std::cout << "bin = " << i << " ec = " << expected_counts << " mpc = "
    //  << most_probable_counts << " gll = "
    //  << count_low_limit << " gul = " << count_up_limit << std::endl;

    find_prob_limits(expected_counts, prob_yellow, prob, count_low_limit, count_up_limit);
    d_yellow_els.push_back(expected_counts - (count_low_limit - 0.5));
    d_yellow_eus.push_back((count_up_limit + 0.5) - expected_counts);
    //std::cout << "            yll = " << count_low_limit << " yul = "
    //  << count_up_limit;

    find_prob_limits(expected_counts, prob_red, prob, count_low_limit, count_up_limit);
    d_red_els.push_back(expected_counts - (count_low_limit - 0.5));
    d_red_eus.push_back((count_up_limit + 0.5) - expected_counts);
    //std::cout << " rll = " << count_low_limit << " rul = "
    //  << count_up_limit << std::endl;
  }

  double p_value = chi_squared_p_value(chi_squared, n_bins - 1);
  std::cout << std::scientific << std::setprecision(16);
  std::cout << "Chi-squared = " << chi_squared << std::endl;
  std::cout << "p-value = " << p_value << std::endl;

  // Prepare Brazil band plot
  TGraphAsymmErrors* green_graph = new TGraphAsymmErrors(d_bin_low_edges.size(),
    &(d_bin_low_edges.front()), &(d_bin_mean_counts.front()), &(d_nbin_zeros.front()),
    &(d_bin_widths.front()), &(d_green_els.front()), &(d_green_eus.front()));
  green_graph->SetFillColor(kGreen);

  TGraphAsymmErrors* yellow_graph = new TGraphAsymmErrors(d_bin_low_edges.size(),
    &(d_bin_low_edges.front()), &(d_bin_mean_counts.front()), &(d_nbin_zeros.front()),
    &(d_bin_widths.front()), &(d_yellow_els.front()), &(d_yellow_eus.front()));
  yellow_graph->SetFillColor(kYellow);

  TGraphAsymmErrors* red_graph = new TGraphAsymmErrors(d_bin_low_edges.size(),
    &(d_bin_low_edges.front()), &(d_bin_mean_counts.front()), &(d_nbin_zeros.front()),
    &(d_bin_widths.front()), &(d_red_els.front()), &(d_red_eus.front()));
  red_graph->SetFillColor(kRed);

  TH1D* expected_hist = static_cast<TH1D*>(hist2->Clone("expected histogram"));
  for (unsigned i = 1; i <= n_bins; ++i) {
    expected_hist->SetBinContent(i, d_bin_mean_counts.at(i - 1));
  }
  expected_hist->SetAxisColor(kBlack);

  //std::cout << "ROOT p-value = " << hist2->Chi2Test(expected_hist, "p") << std::endl;

  TCanvas canvas;
  red_graph->GetXaxis()->SetRangeUser(4.36, 14);
  red_graph->GetYaxis()->SetRangeUser(2.5e4, 4e4);
  red_graph->Draw("a2");
  yellow_graph->Draw("2 same");
  green_graph->Draw("2 same");
  expected_hist->SetLineColor(kBlack);
  expected_hist->Draw("same");
  hist2->SetMarkerStyle(7);
  hist2->Draw("p same");
  canvas.SaveAs("brazil1.pdf");

  red_graph->GetXaxis()->SetRangeUser(13.9, 20);
  red_graph->GetYaxis()->SetRangeUser(6e3, 2.7e4);
  canvas.SaveAs("brazil2.pdf");

  red_graph->GetXaxis()->SetRangeUser(20, 30);
  red_graph->GetYaxis()->SetRangeUser(0, 6e3);
  canvas.SaveAs("brazil3.pdf");

  canvas.SetLogy();
  red_graph->GetXaxis()->SetRangeUser(30, 50);
  red_graph->GetYaxis()->SetRangeUser(1, 3e3);
  canvas.Update();
  canvas.SaveAs("brazil4.pdf");

  TFile f4("brazil.root", "recreate");
  f4.cd();
  hist2->Write();
  canvas.Write();
  expected_hist->Write();
  red_graph->Write();
  yellow_graph->Write();
  green_graph->Write();

  f3.Close();
  f4.Close();
  return 0;
}
