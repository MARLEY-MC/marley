#include "marley_utils.hh"
#include "RootConfigFile.hh"
#include "InterpolationGrid.hh"

// ROOT includes
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"

using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;

marley::RootConfigFile::RootConfigFile() : marley::ConfigFile(),
  root_filename("events.root"), writeroot(false),
  check_before_root_file_overwrite(true) {}

marley::RootConfigFile::RootConfigFile(const std::string& file_name)
  : marley::RootConfigFile()
{
  filename = file_name;
  parse();
}

void marley::RootConfigFile::check_pdf_pairs(const std::vector<double>& Es,
  const std::vector<double>& PDFs)
{
  size_t num_points = Es.size();
  if (PDFs.size() != num_points) throw marley::Error(
    std::string("E and PDF vectors with unequal sizes")
    + " were passed to marley::RootConfigFile::check_pdf_pairs()");

  double previous_E = marley_utils::minus_infinity;
  for(size_t i = 0; i < num_points; ++i) {
    double E = Es.at(i);
    double PDF = PDFs.at(i);

    if (E < 0.) throw marley::Error(
      std::string("Negative") + " energy "
      + std::to_string(E) + " encountered while importing a"
      + " ROOT-based neutrino source specification on line "
      + std::to_string(line_num) + " of the configuration file "
      + filename);

    if (E <= previous_E) throw marley::Error(
      std::string("Energy") + " values "
      + " from a ROOT-based neutrino source specification on line "
      + std::to_string(line_num) + " of the configuration file "
      + filename + " are not strictly increasing");

    if (PDF < 0.) throw marley::Error(
      std::string("Negative") + " probability density "
      + std::to_string(PDF) + " encountered while importing a"
      + " ROOT-based neutrino source specification on line "
      + std::to_string(line_num) + " of the configuration file "
      + filename);

    previous_E = E;
  }
}

template<typename T> T*
  marley::RootConfigFile::get_root_object(const std::string& tfile_name,
  const std::string& namecycle)
{
  // Attempt to open the TFile and complain if something goes wrong.
  // Use a std::unique_ptr so that the TFile object will be auto-deleted
  // (and therefore closed) when this function terminates.
  std::unique_ptr<TFile> file(TFile::Open(tfile_name.c_str(),
    "read"));

  if (file) {
    // Attempt to read in the ROOT object from the file using the given
    // namecycle. Return pointer to it (or nullptr if something went wrong)
    T* obj = dynamic_cast<T*>(file->Get(namecycle.c_str()));
    // Force the TFile to disown the object if it inherits from TH1 (otherwise,
    // ROOT will auto-delete it when the TFile object is deleted).
    // TODO: decide what to do if the type is TTree (disown too? throw exception?)
    // since TTrees also behave this way
    // TODO: add other checks if you discover more ROOT classes that force ownership
    TH1* th1_test = dynamic_cast<TH1*>(obj);
    if (th1_test) th1_test->SetDirectory(nullptr);
    // Return a pointer to the object
    return dynamic_cast<T*>(obj);
  }

  // couldn't open ROOT file
  else throw marley::Error(std::string("Failed to open")
    + " ROOT file '" + tfile_name + "' given on line "
    + std::to_string(line_num) + " of the configuration file " + filename);

  return nullptr;
}

bool marley::RootConfigFile::process_extra_keywords() {

  std::string arg;

  if (keyword == "writeroot") {
    next_word_from_line(arg, true, true);
    if (arg == "yes") writeroot = true;
    else if (arg == "no") writeroot = false;
    else if (arg == "overwrite") {
      writeroot = true;
      check_before_root_file_overwrite = false;
    }
    else {
      throw marley::Error(std::string("Invalid")
        + " ROOT file write flag '" + arg
        + "' encountered on line" + std::to_string(line_num)
        + " of the configuration file " + filename);
    }
    return true;
  }

  else if (keyword == "rootfile") {
    next_word_from_line(arg, true, false);
    root_filename = arg;
    return true;
  }

  return false;
}

bool marley::RootConfigFile::process_extra_source_types(
  const std::string& type, int neutrino_pid)
{
  if (type == "th1" || type == "tgraph") {

    // Get the file name and namecycle for the ROOT object
    // from the configuration file
    std::string source_tfile_name, namecycle;
    next_word_from_line(source_tfile_name, true, false);
    next_word_from_line(namecycle, true, false);

    if (type == "th1") {
      // 1D histogram from a ROOT TFile
      std::unique_ptr<TH1> source_th1(
        get_root_object<TH1>(source_tfile_name, namecycle));

      if (source_th1) {

        // The histogram was successfully recovered from the
        // TFile, so read in the (energy low edge, bin weight)
        // ordered pairs. Keep the overflow bin so that we can
        // use its left edge as the maximum energy value. Skip
        // the unneeded underflow bin.
        const TAxis* x_axis = source_th1->GetXaxis();
        if (!x_axis) throw marley::Error(std::string("Error finding x-axis")
          + " of ROOT histogram (TH1 object) with namecycle '" + namecycle
          + "' from the ROOT file '" + source_tfile_name + "' given on"
          + " line " + std::to_string(line_num)
          + " of the configuration file " + filename);
        // include the overflow bin but not the underflow bin
        size_t nbins = x_axis->GetNbins() + 1;
        std::vector<double> Es(nbins);
        std::vector<double> PDFs(nbins);

        for (size_t b = 1; b <= nbins; ++b) { // underflow bin is bin 0
          Es.at(b - 1) = x_axis->GetBinLowEdge(b);
          // the content of the overflow bin is reset to zero below
          PDFs.at(b - 1) = source_th1->GetBinContent(b);
        }

        // Convert the bin weights to probability densities (used
        // by our GridNeutrinoSource object) by dividing each bin
        // weight by the corresponding bin width.
        for (size_t c = 0; c < nbins - 1; ++c) {
          PDFs.at(c) /= Es.at(c + 1) - Es.at(c);
        }

        // Assign zero probability density to the overflow bin's
        // left edge. This ensures that neutrino energies will be
        // sampled on the half-open interval [Elow, Ehigh), where
        // Elow is the left edge of the first bin and Ehigh is the
        // left edge of the overflow bin.
        PDFs.back() = 0.;

        // Perform some standard checks to make sure that the probability
        // density function represented by these grid points is sane.
        check_pdf_pairs(Es, PDFs);

        // Now that we've processed grid points, create the grid neutrino
        // source
        set_source(std::make_unique<marley::GridNeutrinoSource>(
          Es, PDFs, neutrino_pid, InterpMethod::Constant));

        // We created the source successfully, so indicate that via the return
        // value
        return true;
      }
      else { // couldn't read the histogram
        throw marley::Error(std::string("Could not recover")
          + " ROOT histogram (TH1 object) with namecycle '" + namecycle
          + "' from the ROOT file '" + source_tfile_name + "' given on"
          + " line " + std::to_string(line_num)
          + " of the configuration file " + filename);
      }
    }
    else if (type == "tgraph") {
      // Create an interpolation grid based on the (E, PDF) pairs stored
      // in the TGraph
      std::unique_ptr<TGraph> source_tgraph(
        get_root_object<TGraph>(source_tfile_name, namecycle));

      if (source_tgraph) {

        size_t num_points = source_tgraph->GetN();
        std::vector<double> Es(num_points);
        std::vector<double> PDFs(num_points);

        // Load the energies and PDF values into our vectors
        for(size_t p = 0; p < num_points; ++p) {
          source_tgraph->GetPoint(p, Es.at(p), PDFs.at(p));
        }

        // Perform some standard checks to make sure that the probability
        // density function represented by these grid points is sane.
        check_pdf_pairs(Es, PDFs);

        // Create a neutrino source based on the grid
        set_source(std::make_unique<marley::GridNeutrinoSource>(Es,
          PDFs, neutrino_pid, InterpMethod::LinearLinear));

        return true;
      }
      else throw marley::Error(std::string("Could not recover")
          + " ROOT graph (TGraph object) with namecycle '" + namecycle
          + "' from the ROOT file '" + source_tfile_name + "' given on"
          + " line " + std::to_string(line_num)
          + " of the configuration file " + filename);
    }

  }

  // We couldn't create the source
  return false;
}
