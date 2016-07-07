#pragma once
#include <iostream>
#include <memory>
#include <vector>

#include "Error.hh"
#include "Event.hh"
#include "Level.hh"
#include "Parity.hh"

namespace marley {

  /// @brief Discrete level and &gamma;-ray data for a specific nuclide
  class DecayScheme {

    public:

      // The FileFormat type is used to tell the decay scheme class
      // which format to assume when parsing the file supplied to
      // the constructor
      enum class FileFormat { talys };

      DecayScheme(int Z, int A);

      DecayScheme(int Z, int A, std::string filename,
        FileFormat ff = FileFormat::talys);

      inline const std::vector<std::unique_ptr<marley::Level> >& get_levels()
        const { return levels_; }

      marley::Level& add_level(const marley::Level& level);

      void print_report(std::ostream& ostr = std::cout) const;

      marley::Level* get_pointer_to_closest_level(double E_level);

      void do_cascade(marley::Level* initial_level, marley::Event* p_event,
        marley::Generator& gen, int qIon);

      void print_latex_table(std::ostream& ostr = std::cout);

      inline int Z() const { return Z_; }
      inline int A() const { return A_; }

      void print(std::ostream& out = std::cout) const;
      void read_from_stream(std::istream& in);

    protected:

      int Z_; ///< Atomic number
      int A_; ///< Mass number

      std::vector<std::unique_ptr<marley::Level> > levels_;

      void parse(std::string filename, FileFormat ff = FileFormat::talys);

      // Functions called by the constructor to parse the
      // different nuclear data formats accepted by this class
      void parse_talys(std::string filename);
      // Add more formats as needed


      /// @brief Get the index of the first level whose energy
      /// is not less than Ex
      /// @param Ex Excitation energy (MeV)
      /// @note Returns levels_.size() if all levels have energies below Ex
      size_t level_lower_bound_index(double Ex);
  };

}

inline std::istream& operator>>(std::istream& in, marley::DecayScheme& ds)
{
  ds.read_from_stream(in);
  return in;
}

inline std::ostream& operator<<(std::ostream& out,
  const marley::DecayScheme& ds)
{
  ds.print(out);
  return out;
}
