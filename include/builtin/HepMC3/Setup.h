// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
/// @file Setup.h
/// @brief Definition of \b class Setup

#ifndef HEPMC3_SETUP_H
#define HEPMC3_SETUP_H



namespace HepMC3 {

/// @brief Configuration for HepMC
///
/// Contains macro definitions for printing debug output, feature deprecation, etc.
/// Static class - configuration is shared among all HepMC events
/// and program threads
///
class Setup {

    /// Private constructor
    Setup() {}
    /// Private destructor
    ~Setup() {}


public:

    /// @name Accessors
    /// @{

    /// Get error messages printing flag
    static bool print_errors();
    /// set error messages printing flag
    static void set_print_errors(const bool flag);

    /// Get error messages printing level
    static int errors_level();
    /// set error messages printing level
    static void set_errors_level(const int flag);


    /// Get warning messages printing flag
    static bool print_warnings();
    /// Set warning messages printing flag
    static void set_print_warnings(const bool flag);

    /// Get warning messages printing level
    static int warnings_level();
    /// Set warning messages printing level
    static void set_warnings_level(const int flag);


    /// Get debug level
    static int  debug_level();
    /// Set debug level
    static void set_debug_level(const int level);
    /// @}

    /// @name Static constants
    /// @{
  
    /// Default maxUlps for AlmostEqual2sComplement function (double precision)
    static const unsigned int DEFAULT_DOUBLE_ALMOST_EQUAL_MAXULPS;

    /// Default threshold for comparing double variables
    static const double DOUBLE_EPSILON;

    /// @}


private:

    static int m_errors_level;   //!< Level for printing error messages
    static int m_warnings_level;   //!< Level for printing warnings messages
    static bool m_is_printing_errors;   //!< Flag for printing error messages
    static bool m_is_printing_warnings; //!< Flag for printing warning messages
    static int  m_debug_level;          //!< Level of debug messages printed out
};


} // namespace HepMC3

#endif
