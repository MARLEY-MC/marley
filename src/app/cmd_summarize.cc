// Standard library includes
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#ifdef USE_ROOT
// ROOT includes
#include "TFile.h"
#include "TTree.h"
#endif

// MARLEY includes
#include "cmd_helpers.hh"
#include "marley/CommandHandler.hh"
#include "marley/Error.hh"
#include "marley/EventFileReader.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"

#ifndef USE_ROOT

bool marley::CommandHandler::cmd_summarize(
  std::deque< std::string >& /*args*/ )
{
  std::cerr << "marley: the 'summarize' command requires linking to ROOT.";
  std::cerr << "Please rebuild MARLEY against ROOT and try again.\n";
  return false;
}

#else

bool marley::CommandHandler::cmd_summarize( std::deque< std::string >& args ) {

  std::string output_path;
  bool force = false;
  std::vector< std::string > input_files;

  while ( !args.empty() ) {
    std::string arg = args.front();
    args.pop_front();

    if ( arg == "-o" || arg == "--output" ) {
      if ( args.empty() ) {
        std::cerr << "marley summarize: missing argument after '"
          << arg << "'\n";
        return false;
      }
      output_path = args.front();
      args.pop_front();
    }
    else if ( arg == "-f" || arg == "--force" ) {
      force = true;
    }
    else if ( arg == "-h" || arg == "--help" ) {
      args.clear();
      args.push_front( "summarize" );
      return marley::CommandHandler::cmd_help( args );
    }
    else if ( arg.front() == '-' ) {
      std::cerr << "marley summarize: unrecognized option '" << arg << "'\n";
      return false;
    }
    else if ( output_path.empty() ) {
      output_path = arg;
    }
    else {
      input_files.push_back( arg );
    }
  }

  if ( output_path.empty() ) {
    std::cerr << "marley summarize: missing required output file\n";
    args.push_front( "summarize" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  if ( input_files.empty() ) {
    std::cerr << "marley summarize: no input files specified\n";
    args.push_front( "summarize" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  if ( !force ) {
    std::ifstream test( output_path );
    if ( test ) {
      bool overwrite = marley_utils::prompt_yes_no(
        "Really overwrite " + output_path + "?" );
      if ( !overwrite ) {
        std::cout << "Action aborted.\n";
        return true;
      }
    }
  }

  double flux_avg_tot_xsec;
  int proc_type;

  double Ev, KEv, pxv, pyv, pzv;
  double Mt;
  double El, KEl, pxl, pyl, pzl;
  double Er, KEr, pxr, pyr, pzr;
  int pdgv, pdgt, pdgl, pdgr;
  int np;

  std::vector< int > PDGs;
  std::vector< double > Es, KEs, pXs, pYs, pZs, Ts;

  double Ex;
  int twoJ;
  int par;

  double cv_weight;
  std::vector< double > other_weights;

  TFile out_tfile( output_path.c_str(), "recreate" );
  TTree* out_tree = new TTree( "mst", "MARLEY summary tree" );

  out_tree->Branch( "pdgv", &pdgv, "pdgv/I" );
  out_tree->Branch( "Ev", &Ev, "Ev/D" );
  out_tree->Branch( "KEv", &KEv, "KEv/D" );
  out_tree->Branch( "pxv", &pxv, "pxv/D" );
  out_tree->Branch( "pyv", &pyv, "pyv/D" );
  out_tree->Branch( "pzv", &pzv, "pzv/D" );

  out_tree->Branch( "pdgt", &pdgt, "pdgt/I" );
  out_tree->Branch( "Mt", &Mt, "Mt/D" );

  out_tree->Branch( "pdgl", &pdgl, "pdgl/I" );
  out_tree->Branch( "El", &El, "El/D" );
  out_tree->Branch( "KEl", &KEl, "KEl/D" );
  out_tree->Branch( "pxl", &pxl, "pxl/D" );
  out_tree->Branch( "pyl", &pyl, "pyl/D" );
  out_tree->Branch( "pzl", &pzl, "pzl/D" );

  out_tree->Branch( "pdgr", &pdgr, "pdgr/I" );
  out_tree->Branch( "Er", &Er, "Er/D" );
  out_tree->Branch( "KEr", &KEr, "KEr/D" );
  out_tree->Branch( "pxr", &pxr, "pxr/D" );
  out_tree->Branch( "pyr", &pyr, "pyr/D" );
  out_tree->Branch( "pzr", &pzr, "pzr/D" );

  out_tree->Branch( "Ex", &Ex, "Ex/D" );
  out_tree->Branch( "twoJ", &twoJ, "twoJ/I" );
  out_tree->Branch( "parity", &par, "parity/I" );

  out_tree->Branch( "np", &np, "np/I" );
  out_tree->Branch( "pdgp", &PDGs );
  out_tree->Branch( "Ep",  &Es );
  out_tree->Branch( "KEp", &KEs );
  out_tree->Branch( "pxp", &pXs );
  out_tree->Branch( "pyp", &pYs );
  out_tree->Branch( "pzp", &pZs );
  out_tree->Branch( "tp", &Ts );

  out_tree->Branch( "xsec", &flux_avg_tot_xsec, "xsec/D" );
  out_tree->Branch( "proc", &proc_type, "proc/I" );

  out_tree->Branch( "cv_weight", &cv_weight, "cv_weight/D" );
  out_tree->Branch( "other_weights", &other_weights );

  constexpr double XSEC_CONV = marley_utils::hbar_c2
    * marley_utils::fm2_to_minus40_cm2 * 1e2;

  bool weight_names_written = false;
  long event_count = 0;

  for_each_event( input_files,
    [ & ]( HepMC3::GenEvent& ev, bool first_event, double xsec_natural,
      const auto& first_info )
    {
      if ( first_event && !weight_names_written ) {
        auto wgt_names = first_info->weight_names();
        wgt_names.erase( wgt_names.begin() );
        out_tfile.WriteObject( &wgt_names,
          "MARLEY_other_weight_names", "WriteDelete" );
        weight_names_written = true;
      }

      if ( event_count % 1000 == 0 ) {
        std::cout << "Event " << event_count << '\n';
      }

      PDGs.clear();
      Es.clear();
      KEs.clear();
      pXs.clear();
      pYs.clear();
      pZs.clear();
      Ts.clear();

      auto projectile = marley_hepmc3::get_projectile( ev );
      pdgv = projectile->pid();

      double mv = projectile->generated_mass();
      const HepMC3::FourVector& p4v = projectile->momentum();

      Ev = p4v.e();
      KEv = std::max( 0., Ev - mv );
      pxv = p4v.px();
      pyv = p4v.py();
      pzv = p4v.pz();

      auto target = marley_hepmc3::get_target( ev );
      pdgt = target->pid();
      Mt = target->generated_mass();

      auto ejectile = marley_hepmc3::get_ejectile( ev );
      pdgl = ejectile->pid();

      int ej_id = ejectile->id();

      double ml = ejectile->generated_mass();
      const HepMC3::FourVector& p4l = ejectile->momentum();

      El = p4l.e();
      KEl = std::max( 0., El - ml );
      pxl = p4l.px();
      pyl = p4l.py();
      pzl = p4l.pz();

      auto residue = marley_hepmc3::get_residue( ev );
      pdgr = residue->pid();

      double mr = residue->generated_mass();
      const HepMC3::FourVector& p4r = residue->momentum();

      Er = p4r.e();
      KEr = std::max( 0., Er - mr );
      pxr = p4r.px();
      pyr = p4r.py();
      pzr = p4r.pz();

      auto Ex_attr = residue->attribute< HepMC3::DoubleAttribute >( "Ex" );
      Ex = Ex_attr->value();

      auto twoJ_attr = residue->attribute< HepMC3::IntAttribute >( "twoJ" );
      twoJ = twoJ_attr->value();

      auto parity_attr = residue->attribute< HepMC3::IntAttribute >( "parity" );
      par = parity_attr->value();

      flux_avg_tot_xsec = xsec_natural * XSEC_CONV;

      auto proc_type_attr = ev.attribute< HepMC3::IntAttribute >(
        "signal_process_id" );
      proc_type = proc_type_attr->value();

      np = 0;
      const auto& particles = ev.particles();
      for ( const auto& p : particles ) {
        if ( p->status() != marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS ) {
          continue;
        }

        if ( p->id() == ej_id ) continue;

        ++np;

        PDGs.push_back( p->pid() );

        double mp = p->generated_mass();
        const HepMC3::FourVector& p4p = p->momentum();

        double Ep = p4p.e();
        Es.push_back( Ep );

        double KEp = std::max( 0., Ep - mp );
        KEs.push_back( KEp );

        pXs.push_back( p4p.px() );
        pYs.push_back( p4p.py() );
        pZs.push_back( p4p.pz() );

        const HepMC3::FourVector& pos4_p = p->production_vertex()->position();
        double tp = pos4_p.t();
        tp *= marley_utils::hbar / marley_utils::hbar_c
          / marley_utils::fm_to_cm;
        Ts.push_back( tp );
      }

      other_weights = ev.weights();
      cv_weight = other_weights.front();
      other_weights.erase( other_weights.begin() );

      out_tree->Fill();
      ++event_count;
    } );

  out_tfile.cd();
  out_tree->Write();
  out_tfile.Close();
  return true;
}

#endif
