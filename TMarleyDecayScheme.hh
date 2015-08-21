#pragma once
#include <iostream>
#include <list>
#include <regex>
#include "TMarleyLevel.hh"
#include "TMarleyEvent.hh"

class TMarleyDecayScheme {
  public:
    // The FileFormat type is used to tell the decay scheme class
    // which format to assume when parsing the file supplied to
    // the constructor
    enum class FileFormat { ensdf, talys };

    TMarleyDecayScheme(std::string nucid, std::string filename,
      FileFormat ff = FileFormat::ensdf);
    TMarleyDecayScheme(int z, int a, std::string filename,
      FileFormat ff = FileFormat::ensdf);
    std::string get_nuc_id() const;
    void set_nuc_id(std::string id);
    TMarleyLevel* add_level(const TMarleyLevel level);
    std::list<TMarleyLevel>* get_levels();
    std::vector<TMarleyLevel*>* get_sorted_level_pointers();
    void print_report(std::ostream& ostr = std::cout) const;
    TMarleyLevel* get_pointer_to_closest_level(double E_level);
    void do_cascade(TMarleyLevel* initial_level, TMarleyEvent* p_event) const;
    void print_latex_table(std::ostream& ostr = std::cout);
    void assign_theoretical_RIs(TMarleyLevel* level_i);

  private:
    std::string nuc_id;
    std::string file_name;
    std::list<TMarleyLevel> levels;
    std::vector<TMarleyLevel*> pv_sorted_levels;
    std::vector<double> sorted_level_energies;
    static bool compare_level_energies(TMarleyLevel* first,
      TMarleyLevel* second);
    std::string process_continuation_records(std::ifstream &file_in,
      std::string &record, const std::regex &rx_cont_record) const;

    int Z; // atomic number
    int A; // mass number

    // Functions called by the constructor to parse the
    // different nuclear data formats accepted by this class
    void parse_ensdf();
    void parse_talys();

    // Default level parity
    static const TMarleyParity DEFAULT_PARITY;
    // Fermionic nuclei have default level spin 1/2
    static constexpr int DEFAULT_FERMION_TWOJ = 1;
    // Bosonic nuclei have default level spin 0
    static constexpr int DEFAULT_BOSON_TWOJ = 0;

    // ***** Constants used for ENSDF file parsing *****
    static const std::string ensdf_primary_record; 
    static const std::string ensdf_continuation_record;
    static const std::string ensdf_record_data;
    static const std::regex rx_ensdf_end_record;
    static const std::regex rx_paren_group;
    static const std::regex rx_spin;
    static const std::regex rx_jpi;
};
