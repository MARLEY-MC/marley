#include "RootConfigFile.hh"
#include "InterpolationGrid.hh"

// ROOT includes
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"

using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;

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

bool marley::RootConfigFile::process_extra_source_types(
  const std::string& type, int neutrino_pid, double weight)
{
  // 1D histogram from a ROOT TFile
  if (type == "th1" || type == "roothist") {

    // Get the file name and namecycle for the histogram
    // from the configuration file
    std::string source_tfile_name, namecycle;
    next_word_from_line(source_tfile_name, true, false);
    next_word_from_line(namecycle, true, false);

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

      // Now that we've processed grid points, create the grid neutrino
      // source
      sources.push_back(std::make_unique<marley::GridNeutrinoSource>(
        Es, PDFs, neutrino_pid, weight, InterpMethod::Constant));

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

  // We couldn't create the source
  return false;
}
