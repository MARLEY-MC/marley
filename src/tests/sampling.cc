// Standard library includes
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

// Catch2 includes
#include "catch2/catch.hpp"

// ROOT includes
#ifdef USE_ROOT
  #include "TCanvas.h"
  #include "TGraph.h"
  #include "TH1D.h"
  #include "TLatex.h"
  #include "TLegend.h"
#endif

// MARLEY includes
#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#else
  #include "marley/JSONConfig.hh"
#endif
#include "marley/Error.hh"
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/JSON.hh"
#include "marley/FileManager.hh"
#include "marley/Logger.hh"
#include "marley/Reaction.hh"
#include "marley/marley_kinematics.hh"
#include "marley/marley_utils.hh"

#include "marley/tests/Histogram.hh"

namespace {

// The number of events to generate for each sampling test
constexpr int NUM_EVENTS = 100000;

// The number of bins to use when histogramming events
constexpr int NUM_ENERGY_BINS = 50;
constexpr int NUM_COS_BINS = 5;
constexpr int NUM_TF_BINS = 100;

// Cosine range for histograms
constexpr double COS_MIN = -1.;
constexpr double COS_MAX = 1.;

// Nuclear recoil kinetic energy range for histograms (MeV)
constexpr double TF_MIN = 0.;
constexpr double TF_MAX = 0.15;

// PDG code for a 40K nucleus
constexpr int PDG_40K = 1000190400;
// Excitation energy, 2*spin, and parity values to use for testing
// Hauser-Feshbach continuum decay sampling
constexpr double HF_Exi = 27.; // MeV
constexpr int HF_twoJi = 2; // Ji = 1
const auto HF_Pi = marley::Parity( true ); // positive parity

#ifdef USE_ROOT
// The number of points to use when making TGraphs
// for validation plots
constexpr int NUM_TGRAPH_POINTS = 3000;
#endif

#ifdef USE_ROOT
/// @param hist_title Title to use when drawing the histograms in ROOT
/// @param events_hist Histogram filled with the variable of interest
/// from generated events
/// @param theory_counts Pointer to a vector of expected bin counts based on
/// the theory model. This vector can be generated automatically by
/// Histogram::chi2_test(). If you don't want to plot the expected counts
/// histogram, just pass a nullptr for this argument.
/// @param flux_avg_xsec Flux-averaged total cross section (10^(-40) cm^2)
/// @param chi2 The value of chi^2 from a previously-run test of agreement
/// between the model and the event distribution
/// @param ndof Number of degrees of freedom from the chi^2 test
/// @param p_value P-value obtained from the chi^2 test
/// @param output_filename Name of the output file to which the plot
/// will be saved via a call to TCanvas::SaveAs()
/// @param mod_xsec Function that takes a value within the domain of the
/// histogram and returns the model prediction for the flux-averaged
/// differential cross section. A nullptr passed here indicates that
/// we don't want to draw a continuous model prediction.
void make_plots(const std::string& hist_title,
  const marley::tests::Histogram& events_hist,
  const std::vector<double>* theory_counts, double flux_avg_xsec, double chi2,
  int ndof, double p_value, const std::string& output_filename,
  const std::function<double(double)>* mod_xsec,
  const std::vector<std::string>* bin_labels = nullptr)
{
  // Prepare histograms of the energy distribution from the events
  // and bin integrals of the model PDF
  TH1D model_th1d("model", hist_title.c_str(), events_hist.num_bins(),
    events_hist.x_min(), events_hist.x_max());
  if ( theory_counts ) {
    for ( size_t b = 0; b < theory_counts->size(); ++b ) {
      double counts = theory_counts->at( b );
      model_th1d.SetBinContent(b + 1, counts);
    }
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

  // Get the maximum bin y-value as a reference coordinate
  // for printing the text giving the chi^2 test results
  double events_max = events_th1d.GetMaximum();

  model_th1d.SetLineWidth(2);
  model_th1d.SetLineColor(kRed);
  model_th1d.SetLineStyle(2);
  model_th1d.SetStats(false);

  // If bin labels have been supplied, then have both histograms
  // use them for corresponding bins
  if ( bin_labels ) {
    for ( size_t b = 0; b < events_hist.num_bins(); ++b ) {
      events_th1d.GetXaxis()->SetBinLabel( b + 1, bin_labels->at(b).c_str() );
      model_th1d.GetXaxis()->SetBinLabel( b + 1, bin_labels->at(b).c_str() );
    }

    events_th1d.GetXaxis()->SetLabelSize(0.075);
    events_th1d.GetXaxis()->SetTitleOffset(1.2);
  }

  // Also plot the continuous model distribution if desired
  std::unique_ptr<TGraph> model_graph;
  if ( mod_xsec ) {
    std::vector<double> Xs;
    std::vector<double> model_xsecs;
    double x_min = events_hist.x_min();
    double x_max = events_hist.x_max();
    double graph_x_step = (x_max - x_min)
      / static_cast<double>(NUM_TGRAPH_POINTS);

    for ( size_t k = 0; k <= NUM_TGRAPH_POINTS; ++k ) {
      double x = x_min + (k * graph_x_step);
      double model_diff_xs = mod_xsec->operator()( x );

      Xs.push_back( x );
      model_xsecs.push_back( model_diff_xs );
    }

    // Create the model TGraph and set its plotting style
    model_graph = std::make_unique<TGraph>(Xs.size(), Xs.data(),
      model_xsecs.data());

    model_graph->SetLineColor(kBlue);
    model_graph->SetLineWidth(3);
  }

  TCanvas canvas;
  canvas.SetLeftMargin(0.10);
  canvas.SetRightMargin(0.03);
  canvas.SetBottomMargin(0.10);

  events_th1d.Draw("hist e");
  if ( theory_counts ) model_th1d.Draw("hist same");
  if ( model_graph ) model_graph->Draw("l");

  TLegend lg(0.15, 0.65, 0.3, 0.85);
  lg.AddEntry(&events_th1d, "events", "l");
  if ( theory_counts ) lg.AddEntry(&model_th1d, "model bin integrals", "l");
  if ( model_graph ) lg.AddEntry(model_graph.get(), "model", "l");

  lg.Draw("same");

  // Add the chi2 test results to the plot
  std::ostringstream oss;
  oss << "reduced #chi^{2} = " << chi2 / ndof;

  // Pick a spot to draw the chi^2 and p-value text
  double x_min = events_hist.x_min();
  double x_max = events_hist.x_max();
  double x_ltx = (x_max - x_min) / 12.;
  double y_ltx = 0.65 * events_max;

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
  // Find the configuration file to use for these tests
  const auto& fm = marley::FileManager::Instance();
  std::string test_data_dir = fm.marley_dir() + "/data/tests/";

  std::string config_file_name = fm.find_file( "test.js", test_data_dir );

  // Configure the generator
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
  marley::tests::Histogram energy_hist( NUM_ENERGY_BINS, E_min, E_max );
  marley::tests::Histogram cos_hist( NUM_COS_BINS, COS_MIN, COS_MAX );
  marley::tests::Histogram cevns_hist( NUM_TF_BINS, TF_MIN, TF_MAX );

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

    cevns_hist.fill( ev.residue().kinetic_energy() );

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

    std::function<double(double)> diff_xs_func
      = [&pdf, &flux_avg_xsec](double Ev)
      -> double { return pdf(Ev) * flux_avg_xsec; };

    make_plots(energy_title, energy_hist, nullptr, /* &expected_bin_counts, */
      flux_avg_xsec, chi2, ndof, p_value, "test2.pdf",
      &diff_xs_func);
    #endif

    CHECK( passed );
  }

  INFO("Checking sampling of 2->2 scattering cosines");
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

    make_plots(cos_title, cos_hist, nullptr, /* expected_bin_counts, */
      flux_avg_xsec, chi2, ndof, p_value, "test3.pdf",
      &flux_av_total_diff_xsec);
    #endif

    CHECK( passed );
  }

  INFO("Checking incident neutrino energy sampling");
  {
    const auto& source = gen.get_source();
    marley::tests::Histogram source_E_hist(NUM_ENERGY_BINS, source.get_Emin(),
      source.get_Emax());

    int dummy;
    for ( int e = 0; e < NUM_EVENTS; ++e ) {
      source_E_hist.fill( source.sample_incident_neutrino(dummy, gen) );
    }

    std::function<double(double)> pdf = [&source](double Ev)
      -> double { return source.pdf(Ev); };

    auto expected_bin_counts = source_E_hist.chi2_test(pdf, passed, chi2,
      ndof, p_value);

    // If we've built the tests against ROOT, then make a plot
    #ifdef USE_ROOT
    std::string source_E_title("MARLEY incident neutrinos;"
      " E_{#nu} (MeV); PDF(E_{#nu}) = #phi(E_{#nu}) / #Phi (MeV^{-1})");

    make_plots(source_E_title, source_E_hist, nullptr,
      1., chi2, ndof, p_value, "test4.pdf", &pdf);
    #endif

    CHECK( passed );
  }

  INFO("Checking sampling of Hauser-Feshbach decays");
  {
    // Create a Particle object representing a compound 40K* ion
    // with net charge +1
    const auto& mt = marley::MassTable::Instance();
    double mass = mt.get_atomic_mass(PDG_40K)
      - mt.get_particle_mass(marley_utils::ELECTRON);
    mass += HF_Exi;

    marley::Particle compound_nuc(PDG_40K, mass, 0., 0., 0., mass, 1);

    // Now create a HauserFeshbachDecay object that can be used to decay it
    marley::HauserFeshbachDecay hfd(compound_nuc, HF_Exi, HF_twoJi, HF_Pi, gen);
    MARLEY_LOG_DEBUG() << hfd;

    // **** First, check the MC branching ratios to different final particles
    // against the underlying discrete distribution

    // Gammas aren't included in the list of "fragments," so add them
    // manually to start
    std::map<int, int> pdg_to_num_events = { { marley_utils::PHOTON, 0 } };

    // Loop over the possible fragments that can be emitted, and set their
    // counts to zero initially
    for ( const auto& f : marley::HauserFeshbachDecay::get_fragments() ) {
      pdg_to_num_events[ f.get_pid() ] = 0;
    }

    // Do a bunch of decays, and record which particle (gamma or fragment)
    // gets emitted each time
    for (int e = 0; e < NUM_EVENTS; ++e) {
      marley::Parity dummy_Pf;
      marley::Particle emitted_particle;
      marley::Particle final_nucleus;
      const auto& exit_channel = hfd.sample_exit_channel();

      int emitted_pdg = exit_channel->emitted_particle_pdg();

      pdg_to_num_events.at( emitted_pdg ) += 1;
    }

    // Sum exit channels to get partial widths. Use them to calculated
    // a predicted number of counts for each possible emitted particle.
    std::map<int, double> pdg_to_predicted_counts;
    double total_width = 0.;
    for ( const auto& ec : hfd.exit_channels() ) {
      int pdg = ec->emitted_particle_pdg();
      double width = ec->width();
      total_width += width;

      auto end = pdg_to_predicted_counts.end();
      auto iter = pdg_to_predicted_counts.find( pdg );
      if ( iter == end ) pdg_to_predicted_counts[ pdg ] = width;
      else pdg_to_predicted_counts.at( pdg ) += width;
    }

    // Normalize each of the partial widths to get branching ratios, then
    // multiply by the number of generated decays to get the expected
    // numbers of counts. These will be compared with the MC results.
    for ( auto& pair : pdg_to_predicted_counts ) {
      pair.second *= NUM_EVENTS / total_width;
    }

    // Make a histogram with the MC results and a vector of expected
    // counts. These will be used to perform a chi2 test of agreement
    // between the two. Use one histogram bin per emitted particle PDG
    // code.
    int temp_num_bins = pdg_to_predicted_counts.size();
    marley::tests::Histogram observed_hist(temp_num_bins, 0., temp_num_bins);
    std::vector<double> expected_counts;

    // Build the histogram and expected counts vector for the chi^2 test. Also
    // store the particle symbols for use as bin labels.
    std::vector<std::string> bin_labels;
    size_t bin = 0u;

    // Use these labels for each particle. The ones used in the marley_utils
    // namespace are done with unicode, which doesn't play well with ROOT's
    // plotting utilities.
    std::map<int, std::string> pdg_to_bin_label = {
      { marley_utils::PHOTON, "#gamma" },
      { marley_utils::NEUTRON, "n" },
      { marley_utils::PROTON, "p" },
      { marley_utils::DEUTERON, "d" },
      { marley_utils::HELION, "h" },
      { marley_utils::TRITON, "t" },
      { marley_utils::ALPHA, "#alpha" },
    };

    for ( const auto& pair : pdg_to_num_events ) {
      int pdg = pair.first;
      int num_observed = pair.second;
      observed_hist.set_bin_content(bin, num_observed);
      expected_counts.push_back( pdg_to_predicted_counts.at(pdg) );

      bin_labels.push_back( pdg_to_bin_label.at(pdg) );
      ++bin;
    }

    // Test for agreement
    observed_hist.chi2_test(expected_counts, passed, chi2, ndof, p_value);

    // If we've built the tests against ROOT, then make a plot
    #ifdef USE_ROOT
    std::string ec_title("MARLEY compound nucleus partial decay widths;"
      " channel; #Gamma_{channel} (MeV)");

    make_plots(ec_title, observed_hist, &expected_counts,
      total_width, chi2, ndof, p_value, "test5.pdf", nullptr, &bin_labels);
    #endif

    CHECK( passed );

    //// **** Now test differential decay widths

    // Test each continuum exit channel separately
    for ( const auto& ec : hfd.exit_channels() ) {

      // Skip discrete exit channels, which are sampled using a discrete
      // distribution (and have already been included in the previous test)
      if ( !ec->is_continuum() ) continue;

      const auto& cec = dynamic_cast<const marley::ContinuumExitChannel&>(*ec);

      // Skip sampling a final nuclear spin-parity value during calls
      // to marley::ContinuumExitChannel::do_decay(), since all we care about
      // for this test is the final nuclear excitation energy
      cec.set_skip_jpi_sampling( true );

      // Get the excitation energy bounds of the continuum
      double Ex_min = cec.Emin_;
      double Ex_max = cec.Emax_;

      // Get the total width for this exit channel
      double width = cec.width();

      std::function<double(double)> decay_pdf;

      // Mass of the initial (pre-decay) nucleus
      double mi = compound_nuc.mass();

      // Ground-state mass of the final (post-decay) nucleus
      double mf_gs = cec.gs_residue_.mass();

      // Mass and PDG code of the emitted fragment (or gamma-ray)
      double m_frag;
      int pdg_frag = cec.emitted_particle_pdg();

      if ( cec.emits_fragment() ) {
        const auto& fcec = dynamic_cast<const
          marley::FragmentContinuumExitChannel&>( cec );

        decay_pdf = [&fcec](double Ex) -> double {
          double dummy;
          return fcec.Epdf_(dummy, Ex);
        };

        m_frag = fcec.get_fragment().get_mass();
      }
      else {
        const auto& gcec = dynamic_cast<const
          marley::GammaContinuumExitChannel&>( cec );

        decay_pdf = [&gcec](double Ex) -> double { return gcec.Epdf_(Ex); };

        m_frag = 0.; // photons are massless
      }

      // Use fragment energy instead of excitation energy as the variable on
      // the x-axis. This is more intuitive, even though we're actually
      // sampling excitation energy.
      double KE_frag_min = (mi*mi - std::pow(mf_gs + Ex_max, 2) + m_frag*m_frag)
        / (2. * mi) - m_frag;
      double KE_frag_max = (mi*mi - std::pow(mf_gs + Ex_min, 2) + m_frag*m_frag)
        / (2. * mi) - m_frag;

      // Differential decay width that uses the fragment (or gamma) CM frame
      // kinetic energy as the independent variable
      std::function<double(double)> ddw = [=, &decay_pdf](double KE_fr)
        -> double
      {
        // Final nucleus excitation energy
        double Exf = marley_utils::real_sqrt(mi*mi + m_frag*m_frag
          - 2.*mi*(KE_fr + m_frag)) - mf_gs;

        // Final nucleus mass
        double mf = mf_gs + Exf;

        return width * ( mi / mf ) * decay_pdf(Exf);
      };

      double ddw_norm = marley_utils::num_integrate(ddw, KE_frag_min,
        KE_frag_max);

      // Also make a corresponding PDF normalized to unity
      std::function<double(double)> ddw_pdf = [&ddw, ddw_norm](double KE_fr)
        -> double { return ddw(KE_fr) / ddw_norm; };

      // We're ready. Do some decays and record the fragment CM frame kinetic
      // energy each time in a histogram.
      marley::tests::Histogram hf_decay_hist(50, KE_frag_min, KE_frag_max);

      for ( int e = 0; e < 50000; ++e ) {

        // Reset Ex, Jf, and Pf to the initial values for the decay
        double dummy_Ex = HF_Exi;
        int dummy_twoJf = HF_twoJi;
        marley::Parity dummy_P = HF_Pi;
        marley::Particle fragment, final_nucleus;

        // TODO: finish decay. The exit channel just gets
        // the two final masses ready. HauserFeshbachDecay
        // then samples angles and calls marley_kinematics::two_body_decay()
        // in order to set the particle momenta, etc. You could alternatively
        // just calculate the one quantity that you need (fragment CM KE)
        if ( e % 1000 == 0 ) MARLEY_LOG_INFO() << "Decay " << e;
        cec.do_decay(dummy_Ex, dummy_twoJf, dummy_P,
          fragment, final_nucleus, gen);
        double E_frag_CM = ( mi*mi - std::pow(final_nucleus.mass(), 2)
          + std::pow(fragment.mass(), 2) ) / ( 2. * mi );
        double KE_frag_CM = E_frag_CM - fragment.mass();
        hf_decay_hist.fill( KE_frag_CM );
      }

      // Test the event histogram for agreement with the differential decay
      // width
      auto expected_bin_counts = hf_decay_hist.chi2_test(ddw_pdf, passed, chi2,
        ndof, p_value);

      // If we've built the tests against ROOT, then make a plot
      #ifdef USE_ROOT
      std::string fragment_symbol = pdg_to_bin_label.at( pdg_frag );
      std::string ddw_title("MARLEY differential decay width;"
        "CM frame kinetic energy T_{" + fragment_symbol + "} (MeV);"
        "d#Gamma_{" + fragment_symbol + "}/dT_{" + fragment_symbol + "}");

      make_plots(ddw_title, hf_decay_hist, nullptr, width, chi2, ndof,
        p_value, "test" + fragment_symbol + ".pdf", &ddw);
      #endif

      CHECK( passed );

    } // loop over continuum exit channels
  }

  INFO("Checking CEvNS differential cross section");
  {
    // NOTE: right now, this test assumes that only a single
    // reaction is in use, and that it is CEvNS
    // TODO: make this treatment more general

    const auto* cevns_react = dynamic_cast<const marley::NuclearReaction*>(
      gen.get_reactions().front().get() );
    assert( cevns_react );

    constexpr double gV2 = 1.;
    double QW2 = std::pow( cevns_react->weak_nuclear_charge(), 2 );

    const auto& mt = marley::MassTable::Instance();
    // Projectile mass
    double ma = mt.get_particle_mass( cevns_react->pdg_a() );
    // Nuclear target mass
    double mb = mt.get_atomic_mass( cevns_react->pdg_b() );

    // Differential cross section with respect to lab-frame recoil kinetic
    // energy of the struck nucleus (Tf_lab) at fixed projectile kinetic
    // energy (KEa)
    std::function<double(double, double)> cevns_diff_xsec = [ma, mb, QW2, gV2](double KEa,
      double Tf_lab) -> double
    {
      // Mandelstam s
      double s = std::pow(ma + mb, 2) + 2.*mb*KEa;
      double sqrt_s = marley_utils::real_sqrt( s );

      // Target CM frame total energy
      double Eb_CM = ( s - ma*ma + mb*mb ) / (2. * marley_utils::real_sqrt(s) );
      double Eb_CM2 = std::pow(Eb_CM, 2);

      // Projectile CM frame total energy
      double Ea_CM = std::max(0., sqrt_s - Eb_CM);
      double Ea_CM2 = std::pow(Ea_CM, 2);

      // Maximum lab-frame nuclear recoil kinetic energy
      double Tf_lab_max = 2.*Ea_CM2 / mb;

      if ( Tf_lab_max <= 0. ) return 0.;
      if ( Tf_lab > Tf_lab_max ) return 0.;

      double xsec = marley_utils::GF2 * QW2 * gV2 * mb / (4. * marley_utils::pi);

      xsec *= ( Eb_CM2 / s ) * ( 1. - Tf_lab / Tf_lab_max );

      return xsec;
    };

    // Flux-averaged differential cross section
    std::function<double(double)> flux_avg_cevns_diff_xsec
      = [&gen, &cevns_diff_xsec](double Tf_lab)
    {
      auto& source = const_cast<marley::NeutrinoSource&>( gen.get_source() );
      double E_min = source.get_Emin();
      double E_max = source.get_Emax();

      // Compute the flux-averaged differential cross section in MeV^(-3)
      double result = marley_utils::num_integrate(
        [Tf_lab, &source, &cevns_diff_xsec](double KEa) -> double
        { return source.pdf(KEa) * cevns_diff_xsec(KEa, Tf_lab); },
        E_min, E_max);

      // Convert to (10^(-40) cm^2)
      result *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2;
      return result;
    };

    double cevns_pdf_norm = marley_utils::num_integrate(flux_avg_cevns_diff_xsec,
      TF_MIN, TF_MAX);

    // Flux-averaged total cross section (10^(-40) cm^2)
    double flux_avg_xsec = gen.flux_averaged_total_xs()
      * marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2;

    // Integrating over energy and then angle should give the same
    // flux-averaged total cross section as integrating over angle
    // and then energy
    CHECK( cevns_pdf_norm == Approx(flux_avg_xsec) );

    std::function<double(double)> cevns_pdf = [&flux_avg_cevns_diff_xsec,
      cevns_pdf_norm](double Tf_lab) -> double
    {
      return flux_avg_cevns_diff_xsec(Tf_lab) / cevns_pdf_norm;
    };

    auto expected_bin_counts = cevns_hist.chi2_test(cevns_pdf, passed, chi2,
      ndof, p_value);

    // If we've built the tests against ROOT, then make a plot
    #ifdef USE_ROOT
    std::string cevns_title("MARLEY CEvNS differential cross section;"
      "lab-frame nuclear recoil kinetic energy T_{f} (MeV);"
      " #left[d#sigma/dT_{f}#right]_{flux} (10^{-40} cm^{2} / MeV)");

    make_plots(cevns_title, cevns_hist, nullptr, /* &expected_bin_counts, */
      flux_avg_xsec, chi2, ndof, p_value, "test_cevns.pdf",
      &flux_avg_cevns_diff_xsec);
    #endif

    CHECK( passed );
  }
}
