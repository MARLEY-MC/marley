#pragma once
#include <iostream>
#include <list>
#include <regex>

#include "TMarleyEvent.hh"
#include "TMarleyLevel.hh"
#include "TMarleyParity.hh"

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
    void do_cascade(TMarleyLevel* initial_level, TMarleyEvent* p_event,
      TMarleyGenerator& gen) const;
    void print_latex_table(std::ostream& ostr = std::cout);
    void assign_theoretical_RIs(TMarleyLevel* level_i);

    friend std::ostream& operator<< (std::ostream& out,
      const TMarleyDecayScheme& ds);

    friend std::istream& operator>> (std::istream& in,
      TMarleyDecayScheme& ds);

    inline int get_Z() {
      return Z;
    }

    inline int get_A() {
      return A;
    }

  private:
    int Z; // atomic number
    int A; // mass number
    int pid; // particle ID number for this nucleus
    std::string nuc_id;
    std::list<TMarleyLevel> levels;
    std::vector<TMarleyLevel*> pv_sorted_levels;
    std::vector<double> sorted_level_energies;

    // Helper function for processing ENSDF file continuation records
    std::string process_continuation_records(std::ifstream &file_in,
      std::string &record, const std::regex &rx_cont_record) const;

    // Functions called by the constructor to parse the
    // different nuclear data formats accepted by this class
    void parse_ensdf(std::string filename);
    void parse_talys(std::string filename);

    inline void parse(std::string filename, FileFormat ff = FileFormat::ensdf) {
      // Parse the data file using the appropriate format
      switch (ff) {

        case FileFormat::ensdf:
          parse_ensdf(filename); 
          break;

        case FileFormat::talys:
          parse_talys(filename);
          break;

        default:
          throw std::runtime_error(std::string("Invalid file format ")
            + " supplied to TMarleyDecayScheme constructor.");
      }
    }

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
