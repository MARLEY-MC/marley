// standard library includes
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

// Catch includes
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
#include "marley/marley_utils.hh"

namespace {

constexpr int NUM_EVENTS = 100000;
constexpr int NUM_BINS = 50;

#ifdef USE_ROOT
constexpr int NUM_TGRAPH_POINTS = 10000;
#endif

// Threshold p-value below which we'll reject the null hypothesis
// in chi-squared tests of sampling
constexpr double SIGNIFICANCE_LEVEL = 0.01;

// ROOT has better histogram classes, but we'd like the tests to be
// able to run without it
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
      { return b*x_step_; }

    // Get the x value corresponding to the higher edge of the bin
    inline double get_bin_right_edge(size_t b) const
      { return (b + 1)*x_step_; }

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
        if ( x >= b*x_step_ + x_min_ && x < (b+1)*x_step_ + x_min_ ) {
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
      bool& passed, double& chi2, int& degrees_of_freedom, double& p_value) const
    {
      std::vector<double> expected_counts_vec;
      chi2 = 0.;
      degrees_of_freedom = N_bins_ - 1;
      for ( size_t b = 0; b < N_bins_; ++b ) {
        double observed_counts = bin_counts_.at( b );
        double expected_counts = num_entries_ * marley_utils::num_integrate(pdf,
          get_bin_left_edge( b ), get_bin_right_edge( b ));

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

std::ostream& operator<<(std::ostream& out, const Histogram& h) {
  h.print( out );
  return out;
}

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

} // anonymous namespace

SCENARIO( "Reacting neutrino energies in events match underlying distribution",
  "[physics]" )
{
  GIVEN( "A configured Generator object" ) {

    auto json = marley::JSON::load_file("test.js");
    #ifdef USE_ROOT
      marley::RootJSONConfig config(json);
    #else
      marley::JSONConfig config(json);
    #endif

    marley::Generator gen = config.create_generator();

    WHEN( "Events are generated" ) {

      MARLEY_LOG_INFO() << "Generating " << NUM_EVENTS << " events";

      double E_min = gen.get_source().get_Emin();
      double E_max = gen.get_source().get_Emax();

      Histogram energy_hist(NUM_BINS, E_min, E_max);

      for ( int e = 0; e < NUM_EVENTS; ++e ) {
        if ( e % 1000 == 0 ) MARLEY_LOG_INFO() << "Event " << e;
        marley::Event ev = gen.create_event();
        energy_hist.fill( ev.projectile().kinetic_energy() );
      }

      std::ofstream out_file("energy_hist.txt");
      out_file << energy_hist;
      //std::ifstream in_file("energy_hist.txt");
      //in_file >> energy_hist;

      THEN( "The energy bin probabilities match the underyling distribution" ) {

        bool passed = false;
        std::function<double(double)> pdf = [&gen](double Ev)
          -> double { return gen.E_pdf(Ev); };

        double chi2, p_value;
        int ndof;
        auto expected_bin_counts = energy_hist.chi2_test(pdf, passed, chi2,
          ndof, p_value);

        // If we've built the tests against ROOT, then make a plot
        #ifdef USE_ROOT

          // Prepare histograms of the energy distribution from the events
          // and bin integrals of the model PDF
          std::string energy_title("MARLEY reacting #nu energies;"
            "neutrino energy (MeV); #left[d#sigma/dE_{#nu}#right]_{flux}"
            " (10^{-40} cm^{2} / MeV)");

          TH1D model_th1d("energy_model", energy_title.c_str(), NUM_BINS,
            E_min, E_max);
          for ( size_t b = 0; b < expected_bin_counts.size(); ++b ) {
            double counts = expected_bin_counts.at( b );
            model_th1d.SetBinContent(b + 1, counts);
          }

          TH1D events_th1d = energy_hist.th1d( "energy_events" );
          events_th1d.SetTitle( energy_title.c_str() );

          // Retrieve the flux-averaged total cross section and
          // convert it from MeV^(-2) to 10^(-40) cm^2
          double xsec = gen.flux_averaged_total_xs();
          xsec *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2;

          double total_counts = static_cast<double>( energy_hist.entries() );

          // Factor to convert from bin probabilities to flux-weighted total
          // cross section
          double norm_factor = xsec / (total_counts * energy_hist.bin_width());

          // Normalize the histograms so that they represent flux-weighted
          // total cross section
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
          std::vector<double> energies;
          std::vector<double> model_xsecs;
          double energy_step = (energy_hist.x_max() - energy_hist.x_min())
            / static_cast<double>(NUM_TGRAPH_POINTS);
          double energy_min = energy_hist.x_min();

          double max_model_xs = 0.;
          for ( size_t k = 0; k <= NUM_TGRAPH_POINTS; ++k ) {
            double Ev = energy_min + (k * energy_step);
            double mod_xs = gen.E_pdf(Ev) * xsec;

            energies.push_back( Ev );
            model_xsecs.push_back( mod_xs );

            if ( max_model_xs < mod_xs ) max_model_xs = mod_xs;
          }

          // Create the model TGraph and set its plotting style
          TGraph model_graph(energies.size(), energies.data(),
            model_xsecs.data());

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

          double x_ltx = (E_max - E_min) / 12.;
          double y_ltx = 0.65 * max_model_xs;

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

          canvas.SaveAs("test1.pdf");
          canvas.SaveAs("test1.root");

        #endif

        if ( !passed ) {
          FAIL_CHECK("Neutrino energy distribution failed chi^2 test");
        }
      }
    }
  }
}
