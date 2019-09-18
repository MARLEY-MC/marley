// Standard library includes
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

// GSL includes
#include "gsl/gsl_cdf.h" // p-value for chi-squared test

// ROOT includes
#ifdef USE_ROOT
  #include "TCanvas.h"
  #include "TGraph.h"
  #include "TH1D.h"
  #include "TLatex.h"
  #include "TLegend.h"
#endif

// Catch2 includes
#include "catch.hpp"

// MARLEY includes
#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#else
  #include "marley/JSONConfig.hh"
#endif
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/Logger.hh"
#include "marley/Reaction.hh"
#include "marley/marley_kinematics.hh"
#include "marley/marley_utils.hh"

namespace {

// The number of events to generate for each sampling test
constexpr int NUM_EVENTS = 100000;

// The number of bins to use when histogramming events
constexpr int NUM_ENERGY_BINS = 50;
constexpr int NUM_COS_BINS = 5;

// Cosine range for histograms
constexpr double COS_MIN = -1.;
constexpr double COS_MAX = 1.;

#ifdef USE_ROOT
// The number of points to use when making TGraphs
// for validation plots
constexpr int NUM_TGRAPH_POINTS = 3000;
#endif

// Threshold p-value below which we'll reject the null hypothesis
// in chi-squared tests of sampling
constexpr double SIGNIFICANCE_LEVEL = 0.01;

// ROOT has better histogram classes, but we'd like the tests to be
// able to run without it. This class does just enough to get us
// what we want outside of ROOT.
class Histogram {

  public:

    Histogram(size_t num_bins, double x_min, double x_max)
      : N_bins_( num_bins ), x_min_( x_min ), x_max_( x_max ),
      bin_counts_(num_bins, 0), num_entries_( 0 ),
      x_step_( (x_max - x_min) / static_cast<double>(num_bins) )
    {}

    inline size_t num_bins() const { return N_bins_; }
    inline size_t entries() const { return num_entries_; }
    inline double bin_width() const { return x_step_; }
    inline double x_min() const { return x_min_; }
    inline double x_max() const { return x_max_; }

    // No over/underflow bins, lowest bin is bin zero
    inline int get_bin_content(size_t b) const
      { return bin_counts_.at( b ); }

    // Get the fraction of the total entries that fall in the
    // requested bin
    inline double get_bin_fraction(size_t b) const
      { return bin_counts_.at( b ) / static_cast<double>( num_entries_ ); }

    // Get the x value corresponding to the lower edge of the bin
    inline double get_bin_left_edge(size_t b) const
      { return x_min_ + b*x_step_; }

    // Get the x value corresponding to the higher edge of the bin
    inline double get_bin_right_edge(size_t b) const
      { return x_min_ + (b + 1)*x_step_; }

    // Set the contents of the bth bin
    void set_bin_content(size_t b, int counts) {
      int old_counts = get_bin_content( b );
      num_entries_ -= old_counts;
      bin_counts_.at( b ) = counts;
      num_entries_ += counts;
    }

    inline void increment_bin(size_t b) {
      ++bin_counts_.at( b );
      ++num_entries_;
    }

    inline void fill( double x ) {
      for ( size_t b = 0; b < N_bins_; ++b ) {

        double left_edge = get_bin_left_edge( b );
        double right_edge = get_bin_right_edge( b );

        if ( x >= left_edge && x < right_edge ) {
          ++bin_counts_.at( b );
          break;
        }
      }
      // Increment whether or not a suitable bin was found (correct
      // normalization is maintained even though we don't actually
      // keep the counts for underflow/overflow)
      ++num_entries_;
    }

    void print(std::ostream& os) const {
      os << N_bins_ << ' ' << x_min_ << ' ' << x_max_ << ' ' << num_entries_
        << ' ' << x_step_ << '\n';
      for (const auto& bc : bin_counts_) {
        os << ' ' << bc;
      }
    }

    friend std::istream& operator>>(std::istream& in, Histogram& h);

    #ifdef USE_ROOT
    inline TH1D th1d(const std::string& name) const {
      TH1D hist(name.c_str(), "", N_bins_, x_min_, x_max_);
      // Note that ROOT numbers the bins (except for the underflow and overflow
      // bins) from 1 to N, while we number them from 0 to N - 1 (and don't use
      // underflow and overflow)
      for ( int b = 1; b <= N_bins_; ++b ) {
        hist.SetBinContent(b, bin_counts_.at( b - 1 ));
      }
      return hist;
    }
    #endif

    /// @brief Compares this histogram's contents to a model prediction
    /// @param pdf Probability density function to use to predict bin contents
    /// @param[out] passed Boolean value set to true if the test was passed
    /// (p-value > our chosen significance level) or false otherwise
    /// @return A vector containing the bin counts predicted by the model
    std::vector<double> chi2_test(const std::function<double(double)>& pdf,
      bool& passed, double& chi2, int& degrees_of_freedom, double& p_value)
      const
    {
      std::vector<double> expected_counts_vec;
      chi2 = 0.;
      degrees_of_freedom = N_bins_ - 1;
      for ( size_t b = 0; b < N_bins_; ++b ) {
        double observed_counts = bin_counts_.at( b );
        double expected_counts = num_entries_ * marley_utils::num_integrate(
          pdf, get_bin_left_edge( b ), get_bin_right_edge( b ));

        expected_counts_vec.push_back( expected_counts );

        // To keep the chi^2 approximation valid, skip bins where the
        // expected number of counts is less than 5
        if ( expected_counts >= 5. ) {
          chi2 += std::pow(observed_counts - expected_counts, 2)
            / expected_counts;
        }
        else --degrees_of_freedom;
      }

      p_value = gsl_cdf_chisq_Q( chi2, degrees_of_freedom );

      MARLEY_LOG_INFO() << "chi2 / DOF = " << chi2
        << " / " << degrees_of_freedom;
      MARLEY_LOG_INFO() << "p-value = " << p_value;

      passed = ( p_value >= SIGNIFICANCE_LEVEL );

      return expected_counts_vec;
    }

  protected:
    size_t N_bins_;
    double x_min_;
    double x_max_;
    std::vector<int> bin_counts_;
    size_t num_entries_;
    double x_step_;
};

// Operator for writing a Histogram to a std::ostream
std::ostream& operator<<(std::ostream& out, const Histogram& h) {
  h.print( out );
  return out;
}

// Operator for reading in a Histogram from a std::istream
std::istream& operator>>(std::istream& in, Histogram& h) {
  in >> h.N_bins_;
  in >> h.x_min_;
  in >> h.x_max_;
  in >> h.num_entries_;
  in >> h.x_step_;
  h.bin_counts_.clear();

  int count;
  for (size_t b = 0; b < h.N_bins_; ++b) {
    in >> count;
    h.bin_counts_.push_back( count );
  }

  return in;
}

#ifdef USE_ROOT
/// @param hist_title Title to use when drawing the histograms in ROOT
/// @param events_hist Histogram filled with the variable of interest
/// from generated events
/// @param theory_counts Vector of expected bin counts based on the theory
/// model. This vector can be generated automatically by Histogram::chi2_test()
/// @param flux_avg_xsec Flux-averaged total cross section (10^(-40) cm^2)
/// @param chi2 The value of chi^2 from a previously-run test of agreement
/// between the model and the event distribution
/// @param ndof Number of degrees of freedom from the chi^2 test
/// @param p_value P-value obtained from the chi^2 test
/// @param output_filename Name of the output file to which the plot
/// will be saved via a call to TCanvas::SaveAs()
/// @param mod_xsec Function that takes a value within the domain of the
/// histogram and returns the model prediction for the flux-averaged
/// differential cross section
void make_plots(const std::string& hist_title, const Histogram& events_hist,
  const std::vector<double> theory_counts, double flux_avg_xsec, double chi2,
  int ndof, double p_value, const std::string& output_filename,
  const std::function<double(double)>& mod_xsec)
{
  // Prepare histograms of the energy distribution from the events
  // and bin integrals of the model PDF
  TH1D model_th1d("model", hist_title.c_str(), events_hist.num_bins(),
    events_hist.x_min(), events_hist.x_max());
  for ( size_t b = 0; b < theory_counts.size(); ++b ) {
    double counts = theory_counts.at( b );
    model_th1d.SetBinContent(b + 1, counts);
  }

  // Convert the event Histogram to a TH1D
  TH1D events_th1d = events_hist.th1d( "events" );
  events_th1d.SetTitle( hist_title.c_str() );

  double total_counts = static_cast<double>( events_hist.entries() );

  // Normalize the histograms so that they represent a flux-weighted
  // differential cross section
  double norm_factor = flux_avg_xsec / (total_counts * events_hist.bin_width());
  events_th1d.Scale( norm_factor );
  model_th1d.Scale( norm_factor );

  // Set the plotting style for the histograms
  events_th1d.SetLineWidth(2);
  events_th1d.SetLineColor(kBlack);
  events_th1d.SetStats(false);

  model_th1d.SetLineWidth(2);
  model_th1d.SetLineColor(kRed);
  model_th1d.SetLineStyle(2);
  model_th1d.SetStats(false);

  // Also plot the continuous model distribution
  std::vector<double> Xs;
  std::vector<double> model_xsecs;
  double x_min = events_hist.x_min();
  double x_max = events_hist.x_max();
  double graph_x_step = (x_max - x_min)
    / static_cast<double>(NUM_TGRAPH_POINTS);

  double max_model_diff_xs = 0.;
  for ( size_t k = 0; k <= NUM_TGRAPH_POINTS; ++k ) {
    double x = x_min + (k * graph_x_step);
    double model_diff_xs = mod_xsec( x );

    Xs.push_back( x );
    model_xsecs.push_back( model_diff_xs );

    if ( max_model_diff_xs < model_diff_xs ) max_model_diff_xs = model_diff_xs;
  }

  // Create the model TGraph and set its plotting style
  TGraph model_graph(Xs.size(), Xs.data(), model_xsecs.data());

  model_graph.SetLineColor(kBlue);
  model_graph.SetLineWidth(3);
  //model_graph.SetLineStyle(2);

  TCanvas canvas;
  events_th1d.Draw("hist e");
  //model_th1d.Draw("hist same");
  model_graph.Draw("l");

  TLegend lg(0.15, 0.65, 0.3, 0.85);
  lg.AddEntry(&events_th1d, "events", "l");
  //lg.AddEntry(&model_th1d, "model bin integrals", "l");
  lg.AddEntry(&model_graph, "model", "l");

  lg.Draw("same");

  // Add the chi2 test results to the plot
  std::ostringstream oss;
  oss << "reduced #chi^{2} = " << chi2 / ndof;

  // Pick a spot to draw the chi^2 and p-value text
  double x_ltx = (x_max - x_min) / 12.;
  double y_ltx = 0.65 * max_model_diff_xs;

  TLatex ltx(x_ltx, y_ltx, oss.str().c_str());
  ltx.SetTextColor(kBlack);
  ltx.Draw("same");
  ltx.SetTextSize(0.035);

  oss = std::ostringstream();
  oss << "p-value = " << p_value;

  TLatex ltx2(x_ltx, 0.8*y_ltx, oss.str().c_str());
  ltx2.SetTextColor(kBlack);
  ltx2.Draw("same");
  ltx2.SetTextSize(0.035);

  canvas.SaveAs( output_filename.c_str() );
}
#endif

} // anonymous namespace

// *** Tests for MC sampling of physics quantities ***
TEST_CASE( "Events match their underlying distributions", "[physics]" )
{
  // Configure the generator
  std::string config_file_name("test.js");
  #ifdef USE_ROOT
    marley::RootJSONConfig config( config_file_name );
  #else
    marley::JSONConfig config( config_file_name );
  #endif

  marley::Generator gen = config.create_generator();

  MARLEY_LOG_INFO() << "Generating " << NUM_EVENTS << " events";

  double E_min = gen.get_source().get_Emin();
  double E_max = gen.get_source().get_Emax();

  // Histograms to store the quantities of interest from the events
  // for testing
  Histogram energy_hist(NUM_ENERGY_BINS, E_min, E_max);
  Histogram cos_hist(NUM_COS_BINS, COS_MIN, COS_MAX);

  //std::ifstream ev_file("events.ascii");
  //marley::Event ev;

  //std::ofstream ev_file("events.ascii");

  for ( int e = 0; e < NUM_EVENTS; ++e ) {

    if ( e % 1000 == 0 ) MARLEY_LOG_INFO() << "Event " << e;

    // Generate a new event
    marley::Event ev = gen.create_event();
    //ev_file << ev << '\n';

    // DEBUG
    //ev_file >> ev;

    // Store the neutrino energy
    energy_hist.fill( ev.projectile().kinetic_energy() );

    // Calculate the boost parameters needed to go from the lab frame
    // to the CM frame.
    const marley::Particle& p1 = ev.projectile();
    const marley::Particle& p2 = ev.target();

    double E_tot = p1.total_energy() + p2.total_energy();

    double beta_x = ( p1.px() + p2.px() ) / E_tot;
    double beta_y = ( p1.py() + p2.py() ) / E_tot;
    double beta_z = ( p1.pz() + p2.pz() ) / E_tot;

    // Make copies of the projectile and ejectile and boost them to the CM
    // frame
    marley::Particle pr = ev.projectile();
    marley::Particle ej = ev.ejectile();

    marley_kinematics::lorentz_boost(beta_x, beta_y, beta_z, pr);
    marley_kinematics::lorentz_boost(beta_x, beta_y, beta_z, ej);

    // Store the CM frame scattering angle for the current event.
    double cos_theta_ej_cm = ( pr.px()*ej.px() + pr.py()*ej.py()
      + pr.pz()*ej.pz() ) / pr.momentum_magnitude()
      / ej.momentum_magnitude();

    cos_hist.fill( cos_theta_ej_cm );
  }

  // Reuse these variables for the different test cases
  bool passed = false;
  double chi2, p_value;
  int ndof;

  INFO("Checking sampling of neutrino energies");
  {
    std::function<double(double)> pdf = [&gen](double Ev)
      -> double { return gen.E_pdf(Ev); };

    auto expected_bin_counts = energy_hist.chi2_test(pdf, passed, chi2,
      ndof, p_value);

    // If we've built the tests against ROOT, then make a plot
    #ifdef USE_ROOT
    std::string energy_title("MARLEY reacting #nu energies;"
      "neutrino energy (MeV); #left[d#sigma/dE_{#nu}#right]_{flux}"
      " (10^{-40} cm^{2} / MeV)");

    // Flux-averaged total cross section (10^(-40) cm^2)
    double flux_avg_xsec = gen.flux_averaged_total_xs()
      * marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2;

    make_plots(energy_title, energy_hist, expected_bin_counts,
      flux_avg_xsec, chi2, ndof, p_value, "test2.pdf",
      [&pdf, &flux_avg_xsec](double Ev) -> double
      { return pdf(Ev) * flux_avg_xsec; });
    #endif

    CHECK( passed );
  }

  INFO("Checking sampling of 2->2 scattering cosines")
  {
    // First, define a couple of helper functions
    std::function<double(double, double)> total_diff_xsec
      = [&gen](double Ev, double cos_theta_c_cm) -> double
    {
      const auto& reactions = gen.get_reactions();
      auto& source = const_cast<marley::NeutrinoSource&>(gen.get_source());
      int pdg_a = source.get_pid();
      double diff_xsec = 0.;
      for ( const auto& react : reactions ) {
        diff_xsec += react->diff_xs(pdg_a, Ev, cos_theta_c_cm);
      }
      return diff_xsec;
    };

    std::function<double(double)> flux_av_total_diff_xsec
      = [&gen, &total_diff_xsec](double cos_theta_c_cm)
    {
      auto& source = const_cast<marley::NeutrinoSource&>(gen.get_source());
      double E_min = source.get_Emin();
      double E_max = source.get_Emax();

      // Compute the flux-averaged differential cross section in MeV^(-2)
      double result = marley_utils::num_integrate(
        [cos_theta_c_cm, &source, &total_diff_xsec](double Ev) -> double
        { return source.pdf(Ev) * total_diff_xsec(Ev, cos_theta_c_cm); },
        E_min, E_max);

      // Convert to (10^(-40) cm^2)
      result *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2;
      return result;
    };

    double pdf_norm = marley_utils::num_integrate(flux_av_total_diff_xsec,
      COS_MIN, COS_MAX);

    // Flux-averaged total cross section (10^(-40) cm^2)
    double flux_avg_xsec = gen.flux_averaged_total_xs()
      * marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2;

    // Integrating over energy and then angle should give the same
    // flux-averaged total cross section as integrating over angle
    // and then energy
    CHECK( pdf_norm == Approx(flux_avg_xsec) );

    std::function<double(double)> pdf = [pdf_norm, &flux_av_total_diff_xsec]
      (double cos_theta_c_cm) -> double
      { return flux_av_total_diff_xsec(cos_theta_c_cm) / pdf_norm; };

    auto expected_bin_counts = cos_hist.chi2_test(pdf, passed, chi2,
      ndof, p_value);

    // If we've built the tests against ROOT, then make a plot
    #ifdef USE_ROOT
    std::string cos_title("MARLEY CM frame scattering cosines;"
      " cos(#theta_{c}^{CM});"
      " #left[d#sigma/dcos(#theta_{c}^{CM})#right]_{flux}"
      " (10^{-40} cm^{2})");

    make_plots(cos_title, cos_hist, expected_bin_counts,
      flux_avg_xsec, chi2, ndof, p_value, "test3.pdf",
      flux_av_total_diff_xsec);
    #endif

    CHECK( passed );
  }
}
