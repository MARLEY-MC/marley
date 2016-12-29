#pragma once
#include <iostream>
#include <memory>
#include <vector>

#include "marley/Level.hh"

namespace marley {

  class Event;
  class Generator;

  /// @brief Discrete level and &gamma;-ray data for a specific nuclide
  class DecayScheme {

    public:

      /// @brief The FileFormat type is used to tell the DecayScheme class
      /// which format to assume when parsing a discrete level data file
      /// @details Currently, only discrete level data in the format used
      /// by the <a href="http://talys.eu">TALYS</a> nuclear code is allowed.
      enum class FileFormat { talys };

      /// @brief Create a DecayScheme without any levels
      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      DecayScheme(int Z, int A);

      /// @brief Create a DecayScheme using discrete level data from a file
      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      /// @param filename Name of the file containing the discrete level data
      /// @param format FileFormat specifier that indicates which nuclear data
      /// format is used in the file
      DecayScheme(int Z, int A, std::string filename,
        FileFormat format = FileFormat::talys);

      /// @brief Get a const reference to the vector that holds the Level
      /// objects
      inline const std::vector<std::unique_ptr<marley::Level> >&
        get_levels() const;

      /// @brief Add a level to the DecayScheme
      /// @return A reference to the newly-added Level object
      marley::Level& add_level(const marley::Level& level);

      /// @brief Gets a pointer to the Level in the DecayScheme
      /// whose excitation energy is closest to E_level
      /// @param E_level excitation energy (MeV)
      /// @note Returns nullptr if the DecayScheme doesn't own any Level
      /// objects
      marley::Level* get_pointer_to_closest_level(double E_level);

      /// @brief Simulates nuclear de-excitation via &gamma;-ray emission(s)
      /// @details Gamma-rays will be randomly emitted until the nucleus
      /// reaches its ground state.
      /// @param[in] initial_level Reference to the first level that will
      /// de-excite via &gamma;-ray emission
      /// @param[in,out] event Reference to an Event object that will store
      /// the emitted &gamma;s
      /// @param[in] gen Reference to the Generator to use for random sampling
      /// @param qIon Net charge of the atom or ion whose nucleus is
      /// de-exciting
      void do_cascade(marley::Level& initial_level, marley::Event& event,
        marley::Generator& gen, int qIon);

      /// @brief Get the atomic number
      inline int Z() const;

      /// @brief Get the mass number
      inline int A() const;

      /// @brief Print this DecayScheme object to a std::ostream
      void print(std::ostream& out = std::cout) const;

      /// @brief Use a std::istream to initialize this DecayScheme object,
      /// replacing any previous data.
      /// @details The expected format for data in the std::istream
      /// is the same as the output format in DecayScheme::print()
      void read_from_stream(std::istream& in);

      /// @brief Print LaTeX source code that gives a tabular
      /// representation of the DecayScheme object
      void print_latex_table(std::ostream& ostr = std::cout);

      /// @brief Print a human-readable text representation of the
      /// DecayScheme object
      void print_report(std::ostream& ostr = std::cout) const;

    protected:

      int Z_; ///< Atomic number
      int A_; ///< Mass number

      /// @brief Level objects owned by this DecayScheme
      std::vector<std::unique_ptr<marley::Level> > levels_;

      /// @brief Get the index of the first level whose energy
      /// is not less than Ex
      /// @param Ex Excitation energy (MeV)
      /// @note Returns levels_.size() if all levels have energies below Ex
      size_t level_lower_bound_index(double Ex);

    private:

      /// @brief Helper function that selects the correct parser
      /// when constructing the DecayScheme using a data file
      void parse(std::string filename, FileFormat ff = FileFormat::talys);

      // Functions called by the constructor to parse the
      // different nuclear data formats accepted by this class

      /// @brief Initialize the Level objects in this DecayScheme using a <a
      /// href="http://talys.eu">TALYS</a>-format discrete level data file
      void parse_talys(std::string filename);
      // Add more formats as needed
  };

  // Inline function definitions
  inline int DecayScheme::Z() const { return Z_; }
  inline int DecayScheme::A() const { return A_; }

  inline const std::vector<std::unique_ptr<marley::Level> >&
    DecayScheme::get_levels() const { return levels_; }
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
