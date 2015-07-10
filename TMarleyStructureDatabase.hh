#pragma once
#include <string>
#include <unordered_map>
//#include <unordered_set>

#include "TMarleyDecayScheme.hh"

class TMarleyStructureDatabase {

  public:
    void add_decay_scheme(const std::string nucid, TMarleyDecayScheme ds);
    // TODO: implement these
    //void add_all_decay_schemes_from_file(std::string filename,
    //  TMarleyDecayScheme::FileFormat format
    //  = TMarleyDecayScheme::FileFormat::ensdf);
    //void add_some_decay_schemes_from_file(std::string filename,
    //  std::unordered_set<std::string>& nucids,
    //  TMarleyDecayScheme::FileFormat format
    //  = TMarleyDecayScheme::FileFormat::ensdf);
    void clear();
    void remove_decay_scheme(const std::string nucid);
    TMarleyDecayScheme* get_decay_scheme(const std::string nucid);

  private:
    // Lookup table for TMarleyDecayScheme objects.
    // Keys are ENSDF-style nucIDs, values are decay schemes.
    std::unordered_map<std::string, TMarleyDecayScheme> decay_scheme_table;
};
