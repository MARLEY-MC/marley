#pragma once
#include <iostream>
#include <map>
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
    std::string get_nuc_id() const;
    void set_nuc_id(std::string id);
    void add_level(const TMarleyLevel level);
    std::map<std::string, TMarleyLevel>* get_levels();
    TMarleyLevel* get_level(std::string energy);
    std::vector<TMarleyLevel*>* get_sorted_level_pointers();
    void print_report(std::ostream& ostr = std::cout);
    TMarleyLevel* get_pointer_to_closest_level(double E_level);
    //void do_cascade(std::string initial_energy);
    //void do_cascade(double initial_energy);
    void do_cascade(TMarleyLevel* initial_level, TMarleyEvent* p_event);
    void do_weisskopf(int i);
    void print_latex_table(std::ostream& ostr = std::cout);

  private:
    std::string nuc_id;
    std::string file_name;
    std::map<std::string, TMarleyLevel> levels;
    std::vector<TMarleyLevel*> pv_sorted_levels;
    std::vector<double> sorted_level_energies;
    static bool compare_level_energies(TMarleyLevel* first,
      TMarleyLevel* second);
    std::string process_continuation_records(std::ifstream &file_in,
      std::string &record, const std::regex &rx_cont_record) const;

    int atomic_mass_number;

    // Functions for calculating Weisskopf estimates
    double doubleFact(int i); // Calculating double factorials x!!
    double calcBE(int i);
    double calcBM(int i);
    double calcf(int i);
    double calcTE(int i, double dE, double BE_i, double f_i);
    double calcTM(int i, double dE, double BM_i, double f_i);

    // Functions called by the constructor to parse the
    // different nuclear data formats accepted by this class
    void parse_ensdf();
    void parse_talys();

    // ***** Constants used for ENSDF file parsing *****
    //static const std::regex ensdf_generic_nuc_id("^[[:alnum:] ]{5}");
    static const std::string ensdf_primary_record;
    static const std::string ensdf_continuation_record;
    static const std::string ensdf_record_data;

    static const std::regex rx_ensdf_end_record; // Matches blank lines
    //static const std::regex rx_generic_primary_identification_record;

    // ***** Constants used for TALYS file parsing *****
    //static const std::regex rx_talys_level_line;
    //static const std::regex rx_talys_gamma_line;
};
