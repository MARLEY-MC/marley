#include "marley/marley_utils.hh"
#include "marley/Logger.hh"
#include "marley/RootConfigurationFile.hh"
#include "marley/InterpolationGrid.hh"

// ROOT includes
#include "TAxis.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"

using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;

marley::RootConfigurationFile::RootConfigurationFile()
  : marley::ConfigurationFile(), check_root_overwrite_(true),
  root_filename_("events.root"), write_root_(false) {}

marley::RootConfigurationFile::RootConfigurationFile(
  const std::string& file_name) : marley::RootConfigurationFile()
{
  filename_ = file_name;
  parse();
}

void marley::RootConfigurationFile::check_pdf_pairs(
  const std::vector<double>& Es, const std::vector<double>& PDFs)
{
  size_t num_points = Es.size();
  if (PDFs.size() != num_points) throw marley::Error(
    std::string("E and PDF vectors with unequal sizes")
    + " were passed to marley::RootConfigurationFile::check_pdf_pairs()");

  double previous_E = marley_utils::minus_infinity;
  for(size_t i = 0; i < num_points; ++i) {
    double E = Es.at(i);
    double PDF = PDFs.at(i);

    if (E < 0.) throw marley::Error(
      std::string("Negative") + " energy "
      + std::to_string(E) + " encountered while importing a"
      + " ROOT-based neutrino source specification on line "
      + std::to_string(line_num_) + " of the configuration file "
      + filename_);

    if (E <= previous_E) throw marley::Error(
      std::string("Energy") + " values "
      + " from a ROOT-based neutrino source specification on line "
      + std::to_string(line_num_) + " of the configuration file "
      + filename_ + " are not strictly increasing");

    if (PDF < 0.) throw marley::Error(
      std::string("Negative") + " probability density "
      + std::to_string(PDF) + " encountered while importing a"
      + " ROOT-based neutrino source specification on line "
      + std::to_string(line_num_) + " of the configuration file "
      + filename_);

    previous_E = E;
  }
}

/// @note Since this function template is a private member of the
/// RootConfigurationFile class, we can get away with placing its definition in
/// the source file. If you derive any classes from RootConfigurationFile and
/// want to switch this function to "protected" access, or if you switch this
/// function to "public" access, you will need to put the template definition
/// in the header file.
template<typename T> T*
  marley::RootConfigurationFile::get_root_object(const std::string& tfile_name,
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
    // TODO: decide what to do if the type is TTree (disown too? throw
    // exception?) since TTrees also behave this way
    // TODO: add other checks if you discover more ROOT classes that force
    // ownership
    TH1* th1_test = dynamic_cast<TH1*>(obj);
    if (th1_test) th1_test->SetDirectory(nullptr);
    // Return a pointer to the object
    return dynamic_cast<T*>(obj);
  }

  // couldn't open ROOT file
  else throw marley::Error(std::string("Failed to open")
    + " ROOT file '" + tfile_name + "' given on line "
    + std::to_string(line_num_) + " of the configuration file " + filename_);

  return nullptr;
}

bool marley::RootConfigurationFile::process_extra_keywords() {

  std::string arg;

  if (keyword_ == "writeroot") {
    next_word_from_line(arg, true, true);
    if (arg == "yes") write_root_ = true;
    else if (arg == "no") write_root_ = false;
    else if (arg == "overwrite") {
      write_root_ = true;
      check_root_overwrite_ = false;
    }
    else {
      throw marley::Error(std::string("Invalid")
        + " ROOT file write flag '" + arg
        + "' encountered on line" + std::to_string(line_num_)
        + " of the configuration file " + filename_);
    }
    return true;
  }

  else if (keyword_ == "rootfile") {
    next_word_from_line(arg, true, false);
    root_filename_ = arg;
    return true;
  }

  return false;
}

bool marley::RootConfigurationFile::process_extra_source_types(
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
          + " line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);
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
        MARLEY_LOG_INFO() << "Created a " << type << ' '
          << neutrino_pid_to_string(neutrino_pid) << " source";

        // We created the source successfully, so indicate that via the return
        // value
        return true;
      }
      else { // couldn't read the histogram
        throw marley::Error(std::string("Could not recover")
          + " ROOT histogram (TH1 object) with namecycle '" + namecycle
          + "' from the ROOT file '" + source_tfile_name + "' given on"
          + " line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);
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
        MARLEY_LOG_INFO() << "Created a tgraph "
          << neutrino_pid_to_string(neutrino_pid) << " source";

        return true;
      }
      else throw marley::Error(std::string("Could not recover")
          + " ROOT graph (TGraph object) with namecycle '" + namecycle
          + "' from the ROOT file '" + source_tfile_name + "' given on"
          + " line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);
    }

  }

  // We couldn't create the source
  return false;
}
