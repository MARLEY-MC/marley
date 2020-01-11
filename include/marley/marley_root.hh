/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#pragma once

// standard library includes
#include <memory>

// ROOT includes
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/NeutrinoSource.hh"

namespace marley_root {

  /// @brief Helper function for loading objects from a <a
  /// href="http://root.cern.ch">ROOT</a> <a
  /// href="https://root.cern.ch/doc/master/classTFile.html">TFile</a>
  /// @param tfile_name Name of the ROOT file to read
  /// @param namecycle Namecycle of the object to load from the file
  /// @return A pointer to the object after it has been loaded from
  /// the file, or nullptr if the load failed
  template<typename T> T* get_root_object(const std::string& tfile_name,
    const std::string& namecycle);

  std::unique_ptr<marley::GridNeutrinoSource>
    make_root_neutrino_source(int pdg_code, const TH1* th1);

  std::unique_ptr<marley::GridNeutrinoSource>
    make_root_neutrino_source(int pdg_code, const TGraph* tgraph);
}

template<typename T> T*
  marley_root::get_root_object(const std::string& tfile_name,
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
    // Force the TFile to disown the object if it inherits from TH1
    // (otherwise, ROOT will auto-delete it when the TFile object is
    // deleted).
    // TODO: decide what to do if the type is TTree (disown too? throw
    // exception?) since TTrees also behave this way
    // TODO: add other checks if you discover more ROOT classes that force
    // ownership
    TH1* th1_test = dynamic_cast<TH1*>(obj);
    if (th1_test) th1_test->SetDirectory(nullptr);
    // Return a std::unique_ptr to the object
    return dynamic_cast<T*>(obj);
  }

  // couldn't open ROOT file
  else throw marley::Error(std::string("Failed to open") + " ROOT file '"
    + tfile_name + "'");

  return nullptr;
}
