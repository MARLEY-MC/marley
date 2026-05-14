// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2023 The HepMC collaboration (see AUTHORS for details)
//
#ifndef HEPMC3_READERASCII_H
#define HEPMC3_READERASCII_H
///
/// @file  ReaderAscii.h
/// @brief Definition of class \b ReaderAscii
///
/// @class HepMC3::ReaderAscii
/// @brief GenEvent I/O parsing for structured text files
///
/// @ingroup IO
///
#include <set>
#include <string>
#include <fstream>
#include <istream>
#include <iterator>
#include <unordered_map>
#include "HepMC3/Reader.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Data/GenEventData.h"


namespace HepMC3 {


class ReaderAscii : public Reader {
public:

    /// @brief Constructor
    ReaderAscii(const std::string& filename);
    /// The ctor to read from stream
    ReaderAscii(std::istream &);
    /// The ctor to read from stream. Useful for temp. streams
    ReaderAscii(std::shared_ptr<std::istream> s_stream);
    /// @brief Destructor
    ~ReaderAscii();

    /// @brief skip events
    bool skip(const int)  override;

    /// @brief Load event from file
    ///
    /// @param[out] evt Event to be filled
    bool read_event(GenEvent& evt)  override;

    /// @brief Return status of the stream
    bool failed()  override;

    /// @brief Close file stream
    void close()  override;

private:

    /// @brief Unsecape '\' and '\n' characters in string
    static std::string unescape(const std::string& s);

    /// @name Read helpers
    /// @{

    /// @brief Parse event
    ///
    /// Helper routine for parsing event information
    /// @param[in]  buf Line of text that needs to be parsed
    /// @return vertices count and particles count for verification
    std::pair<int,int> parse_event_information(const char *buf);

    /// @brief Parse weight value lines
    ///
    /// Helper routine for parsing weight value information
    /// @param[in]  buf Line of text that needs to be parsed
    ///
    bool parse_weight_values(const char *buf);

    /// @brief Parse units
    ///
    /// Helper routine for parsing units information
    /// @param[in]  buf Line of text that needs to be parsed
    ///
    bool parse_units(const char *buf);

    /// @brief Parse struct GenPdfInfo information
    ///
    /// Helper routine for parsing PDF information
    /// @param[in]  buf Line of text that needs to be parsed
    bool parse_pdf_info(const char *buf);

    /// @brief Parse struct GenHeavyIon information
    ///
    /// Helper routine for parsing heavy ion information
    /// @param[in]  buf Line of text that needs to be parsed
    bool parse_heavy_ion(const char *buf);

    /// @brief Parse struct GenCrossSection information
    ///
    /// Helper routine for parsing cross-section information
    /// @param[in]  buf Line of text that needs to be parsed
    bool parse_cross_section( const char *buf);

    /// @brief Parse vertex
    ///
    /// Helper routine for parsing single event information
    /// @param[in] buf Line of text that needs to be parsed
    ///
    bool parse_vertex_information(const char *buf);

    /// @brief Parse particle
    ///
    /// Helper routine for parsing single particle information
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_particle_information(const char *buf);

    /// @brief Parse attribute
    ///
    /// Helper routine for parsing single attribute information
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_attribute( const char *buf);

    /// @brief Parse run-level attribute.
    ///
    /// Helper routine for parsing single attribute information
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_run_attribute(const char *buf);

    /// @brief Parse run-level weight names.
    ///
    /// Helper routine for parsing a line with information about
    /// weight names.
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_weight_names(const char *buf);

    /// @brief Parse run-level tool information.
    ///
    /// Helper routine for parsing a line with information about
    /// tools being used.
    /// @param[in] buf Line of text that needs to be parsed
    bool parse_tool(const char *buf);
    /// @}


private:

    std::ifstream m_file; //!< Input file
    std::shared_ptr<std::istream> m_shared_stream = nullptr;///< For ctor when reading from temp. stream
    std::istream* m_stream = nullptr; ///< For ctor when reading from stream
    bool m_isstream; ///< toggles usage of m_file or m_stream

    /** @brief Temp storage for sets of incoming/outgoing ids for explicit vertices.*/
    std::map<int, std::pair< std::set<int>, std::set<int> > >  m_io_explicit;
    /** @brief Temp storage for sets of incoming/outgoing ids for implicit vertices.*/
    std::unordered_map<int, std::pair< std::set<int>, std::set<int> > >  m_io_implicit;
    /** @brief Temp storage to keep the order of implicit vertices.*/
    std::vector<int> m_io_implicit_ids;
    /** @brief Temp storage to keep the order of explicit vertices.*/
    std::set<int> m_io_explicit_ids;

    GenEventData m_data; //!< To hold event information.
};


} // namespace HepMC3

#endif
