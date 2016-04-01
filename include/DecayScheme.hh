#pragma once
#include <iostream>
#include <list>
#include <regex>

#include "Error.hh"
#include "Event.hh"
#include "Level.hh"
#include "Parity.hh"

// Forward-declare the class and friend operators so that we can use the
// operators in the global scope
namespace marley { class DecayScheme; }

std::istream& operator>> (std::istream& in, marley::DecayScheme& ds);
std::ostream& operator<< (std::ostream& out, const marley::DecayScheme& ds);

namespace marley {

  class DecayScheme {
    public:
      // The FileFormat type is used to tell the decay scheme class
      // which format to assume when parsing the file supplied to
      // the constructor
      enum class FileFormat { ensdf, talys };
  
      DecayScheme(std::string nucid, std::string filename,
        FileFormat ff = FileFormat::ensdf);
      DecayScheme(int z, int a, std::string filename,
        FileFormat ff = FileFormat::ensdf);
      std::string get_nuc_id() const;
      void set_nuc_id(std::string id);
      marley::Level* add_level(const marley::Level level);
      std::list<marley::Level>* get_levels();
      const std::vector<marley::Level*>* get_sorted_level_pointers() const;
      void print_report(std::ostream& ostr = std::cout) const;
      marley::Level* get_pointer_to_closest_level(double E_level);
      void do_cascade(marley::Level* initial_level, marley::Event* p_event,
        marley::Generator& gen, int qIon);
      void print_latex_table(std::ostream& ostr = std::cout);
      void assign_theoretical_RIs(marley::Level* level_i);
  
      friend std::ostream& ::operator<< (std::ostream& out,
        const marley::DecayScheme& ds);
  
      friend std::istream& ::operator>> (std::istream& in,
        marley::DecayScheme& ds);
  
      inline int get_Z() const {
        return Z;
      }
  
      inline int get_A() const {
        return A;
      }
  
    private:
      int Z; // atomic number
      int A; // mass number
      int pid; // particle ID number for this nucleus
      std::string nuc_id;
      std::list<marley::Level> levels;
      std::vector<marley::Level*> pv_sorted_levels;
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
            throw marley::Error(std::string("Invalid file format ")
              + " supplied to marley::DecayScheme constructor.");
        }
      }
  
      // Default level parity
      static const marley::Parity DEFAULT_PARITY;
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

}
