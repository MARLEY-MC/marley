// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_CROSS_SECTION_H
#define HEPMC3_CROSS_SECTION_H
/**
 *  @file GenCrossSection.h
 *  @brief Definition of attribute \b class GenCrossSection
 *
 *  @class HepMC3::GenCrossSection
 *  @brief Stores additional information about cross-section
 *
 *  This is an example of event attribute used to store cross-section information
 *
 *  This class is meant to be used to pass, on an event by event basis,
 *  the current best guess of the total cross section.
 *  It is expected that the final cross section will be stored elsewhere.
 *
 *    - double cross_section;       // cross section in pb.
 *    - double cross_section_error; // error associated with this cross section.
 *    - long accepted_events;       ///< The number of events generated so far.
 *    - long attempted_events;      ///< The number of events attempted so far.
 *
 *  In addition, several cross sections and related info can be
 *  included in case of runs with mulltiple weights.
 *
 *  The units of cross_section and cross_section_error are expected to be pb.
 *
 *  @ingroup attributes
 *
 */
#include <iostream>
#include <algorithm>
#include "HepMC3/Attribute.h"

namespace HepMC3 {
/** Deprecated */
using namespace std;

class GenCrossSection : public Attribute {

//
// Fields
//
private:

    long accepted_events;       ///< The number of events generated so far.
    long attempted_events;      ///< The number of events attempted so far.

    std::vector<double> cross_sections;       ///< Per-weight cross-section.
    std::vector<double> cross_section_errors; ///< Per-weight errors.
//
// Functions
//
public:
    /** @brief Implementation of Attribute::from_string */
    bool from_string(const std::string &att)  override;

    /** @brief Implementation of Attribute::to_string */
    bool to_string(std::string &att) const override;
    /// @name Deprecated functionality
    /// @{
    /** @brief Set all fields
        @deprecated Use set_cross_section(const std::vector<double>& xs, const std::vector<double>& xs_err instead

        @param xs Cross section
        @param xs_err Uncertainty on cross-section
        @param n_acc Number of accepted events
        @param n_att Number of attempted events
    */
    void set_cross_section(const double& xs, const double& xs_err,const long& n_acc = -1, const long& n_att = -1);

    /// @}
    /** @brief Set all fields */
    void set_cross_section(const std::vector<double>& xs, const std::vector<double>& xs_err,const long& n_acc = -1, const long& n_att = -1);

    /** @brief Get the cross-sections
     */
    const std::vector<double>& xsecs() const { return cross_sections; }

    /** @brief Get the cross-section errors
     */
    const std::vector<double>& xsec_errs() const { return cross_section_errors; }


    /** @brief Set the number of accepted events
     */
    void set_accepted_events(const long& n_acc ) {
        accepted_events=n_acc;
    }

    /** @brief Set the number of attempted events
     */
    void set_attempted_events(const long& n_att ) {
        attempted_events=n_att;
    }

    /** @brief Get the number of accepted events
     */
    long get_accepted_events() const {
        return accepted_events;
    }

    /** @brief Get the number of attempted events
     */
    long get_attempted_events() const {
        return  attempted_events;
    }

    /** @brief Set the cross section  corresponding to the weight
        named \a wName.
     */
    void set_xsec(const std::string& wName,const double& xs) {
        int pos = windx(wName);
        if ( pos < 0 ) throw std::runtime_error("GenCrossSection::set_xsec(const std::string&,const double&): no weight with given name in this run");
        set_xsec(pos, xs);
    }

    /** @brief Set the cross section corresponding to the weight with
        index \a indx.
     */
    void set_xsec(const unsigned long& index, const double& xs) {
        if ( index >= cross_sections.size() ) {throw std::runtime_error("GenCrossSection::set_xsec(const unsigned long&): index outside of range");}
        cross_sections[index] = xs;
    }

    /** @brief Set the cross section error corresponding to the weight
        named \a wName.
     */
    void set_xsec_err(const std::string& wName, const double& xs_err) {
        int pos = windx(wName);
        if ( pos < 0 ) throw std::runtime_error("GenCrossSection::set_xsec_err(const std::string&,const double&): no weight with given name in this run");
        set_xsec_err(pos, xs_err);
    }

    /** @brief Set the cross section error corresponding to the weight
        with index \a indx.
     */
    void set_xsec_err(const unsigned long& index, const double& xs_err) {
        if ( index >= cross_section_errors.size() ) {throw std::runtime_error("GenCrossSection::set_xsec_err(const unsigned long&): index outside of range");}
        cross_section_errors[index] = xs_err;
    }

    /** @brief Get the cross section corresponding to the weight named
        \a wName.
     */
    double xsec(const std::string& wName) const {
        int pos = windx(wName);
        if ( pos < 0 ) throw std::runtime_error("GenCrossSection::xsec(const std::string&): no weight with given name in this run");
        return xsec(pos);
    }

    /** @brief Get the cross section corresponding to the weight with index
        \a indx.
     */
    double xsec(const unsigned long& index = 0) const {
        if ( index < cross_sections.size() ) { return cross_sections.at(index); }
        else  { throw std::runtime_error("GenCrossSection::xsec(const unsigned long&): index outside of range");}
        return 0.0;
    }

    /** @brief Get the cross section error corresponding to the weight
        named \a wName.
     */
    double xsec_err(const std::string& wName) const {
        int pos = windx(wName);
        if ( pos < 0 ) throw std::runtime_error("GenCrossSection::xsec_err(const std::string&): no weight with given name in this run");
        return xsec_err(pos);
    }

    /** @brief Get the cross section error corresponding to the weight
        with index \a indx.
     */
    double xsec_err(const unsigned long& index = 0) const {
        if ( index < cross_section_errors.size() ) {return cross_section_errors.at(index);}
        else  { throw std::runtime_error("GenCrossSection::xsec_err(const unsigned long&): index outside of range");}
        return 0.0;
    }

    bool operator==( const GenCrossSection& ) const; ///< Operator ==
    bool operator!=( const GenCrossSection& ) const; ///< Operator !=
    bool is_valid()                           const; ///< Verify that the instance contains non-zero information

private:

    /** @brief get the weight index given a weight name. */
    int windx(const std::string& wName) const;

};

} // namespace HepMC3

#endif
