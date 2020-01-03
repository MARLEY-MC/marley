#pragma once
#include <string>
#include <vector>

#include "marley/Particle.hh"

namespace marley {

  // Status codes to use when dumping the initial and
  // final state particles to HEPEVT format
  constexpr int HEPEVT_INITIAL_STATE_STATUS_CODE = 3;
  constexpr int HEPEVT_FINAL_STATE_STATUS_CODE = 1;

  #ifndef __MAKECINT__
  // Forward declare the JSON class so that we can define a function that
  // creates a JSON representation of a Event. Hide the JSON class from
  // rootcint so that we won't have issues using marley::Event objects with
  // ROOT 5.
  class JSON;
  #endif

  /// @brief Container for ingoing and outgoing momentum 4-vectors from a
  /// reaction
  /// @details For a two-two scattering reaction a + b &rarr; c + d, this class
  /// calls the ingoing particles the projectile (a) and target (b). The
  /// outgoing particles are called the ejectile (c) and residue (d). If the
  /// residue is a composite particle (e.g., a nucleus), then the class member
  /// Ex_ stores its excitation energy just after the two-two reaction
  /// occurred.  After MARLEY simulates its de-excitation, the residue Particle
  /// object will be in its ground state, and the final_particles_ member of
  /// this class will include Particle objects representing the de-excitation
  /// products.
  /// @note This class manually manages memory (using the
  /// &ldquo;<a href="http://tinyurl.com/mohpthc">rule of 5</a>&rdquo;)
  /// for the Particle objects that it owns. A better implementation from a
  /// modern C++ perspective would be to use vectors of
  /// std::unique_ptr<marley::Particle> objects and let those manage the memory
  /// (i.e., use the
  /// &ldquo;<a href="http://tinyurl.com/mohpthc">rule of zero</a>&rdquo;).
  /// This design choice was made to ensure compatibility with ROOT 5, which
  /// cannot generate dictionaries for C++11 classes like std::unique_ptr.
  class Event {

    public:

      /// @brief Create an event with dummy particles
      Event(double Ex = 0.);

      /// @brief Create a two-two scattering event
      /// @param a the projectile
      /// @param b the target
      /// @param c the ejectile
      /// @param d the residue
      /// @param Ex residue excitation energy (MeV) just after the two-two
      /// reaction
      Event(const marley::Particle& a,
        const marley::Particle& b, const marley::Particle& c,
        const marley::Particle& d, double Ex = 0.);

      // Destructor
      ~Event();

      /// @brief Copy constructor
      Event(const Event& other_event);

      /// @brief Move constructor
      Event(Event&& other_event);

      /// @brief Copy assignment operator
      Event& operator=(const Event& other_event);

      /// @brief Move assignment operator
      Event& operator=(Event&& other_event);

      /// @brief Get a const reference to the projectile
      const marley::Particle& projectile() const;
      /// @brief Get a const reference to the target
      const marley::Particle& target() const;
      /// @brief Get a const reference to the ejectile
      const marley::Particle& ejectile() const;
      /// @brief Get a const reference to the residue
      const marley::Particle& residue() const;

      /// @brief Get a non-const reference to the projectile
      marley::Particle& projectile();
      /// @brief Get a non-const reference to the target
      marley::Particle& target();
      /// @brief Get a non-const reference to the ejectile
      marley::Particle& ejectile();
      /// @brief Get a non-const reference to the residue
      marley::Particle& residue();

      /// @brief Get a const reference to the vector of initial particles
      inline const std::vector<marley::Particle*>& get_initial_particles()
        const;

      /// @brief Get a non-const reference to the vector of initial particles
      inline std::vector<marley::Particle*>& get_initial_particles();

      /// @brief Get a const reference to the vector of final particles
      inline const std::vector<marley::Particle*>& get_final_particles() const;

      /// @brief Get a non-const reference to the vector of final particles
      inline std::vector<marley::Particle*>& get_final_particles();

      /// @brief Get the excitation energy of the residue just after the
      /// initial two-body reaction
      inline double Ex() const;

      /// @brief Add a Particle to the vector of initial particles
      void add_initial_particle(const marley::Particle& p);

      /// @brief Add a Particle to the vector of final particles
      void add_final_particle(const marley::Particle& p);

      /// @brief Write a
      /// <a href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEVT</a> record for this event to a std::ostream. Use the spacetime
      /// origin (t = 0 mm/c, x = 0 mm, y = 0 mm, z = 0 mm) as the initial
      /// position 4-vector for all particles.
      /// @param event_num The event number to use in the output HEPEVT record
      /// @param out The std::ostream to which the HEPEVT record will be written
      /// @todo Alter marley::Event::write_hepevt() so that the user can specify
      /// a vertex position 4-vector to use.
      void write_hepevt(size_t event_num, std::ostream& out) const;

      /// @brief Print this event to a std::ostream
      /// @param out The std::ostream to which this event will be written
      void print(std::ostream& out) const;

      /// @brief Read in this event from a std::istream. Any previous contents
      /// of this event will be deleted
      /// @param in The std::istream from which this event will be read
      void read(std::istream& in);

      /// @brief Read in this event from a std::istream, assuming it
      /// will appear there in HEPEVT format. Any previous contents
      /// of this event will be deleted.
      /// @details A marley::Error will be thrown if
      ///   - There are not exactly two particles in the event with the status
      ///     code HEPEVT_INITIAL_STATE_STATUS_CODE
      ///   - There is not exactly one particle with status code
      ///     HEPEVT_INITIAL_STATE_STATUS_CODE that is an ion
      ///     (a particle for which marley_utils::is_ion() returns true)
      ///   - More than one lepton (a particle for which marley_utils::is_lepton()
      ///     returns true) appears in the event with status code
      ///     HEPEVT_FINAL_STATE_STATUS_CODE
      ///   - Zero of the particles with status code
      ///     HEPEVT_FINAL_STATE_STATUS_CODE are ions
      /// When reconstructing the nuclear excitation energy Ex_ from the
      /// HEPEVT event record, it is assumed that the outgoing hadronic system
      /// immediately following the initial 2->2 scatter is on the mass shell,
      /// that (anti)neutrino-nucleus charged-current
      /// interactions produce an ion with (-1) +1 net charge, and that
      /// energy conservation between the two-particle initial state and
      /// the multi-particle final state is respected in the event record.
      /// Note also that, due to numerical roundoff, the value
      /// of Ex_ will not be preserved to full precision when making the event format
      /// conversions ASCII -> HEPEVT -> ASCII.
      /// @return True if the input stream was in a good state after attempting
      /// to read in the full HEPEVT record, or false otherwise. This behavior is
      /// designed to enable the return value of this function to be used as a while
      /// loop condition when reading in multiple HEPEVT records from an input stream.
      /// If the return value is false, the event object contents will have been cleared
      /// by this function.
      bool read_hepevt(std::istream& in);

      /// @brief Deletes all particles from the event and resets
      /// the nuclear excitation energy to zero
      void clear();

      #ifndef __MAKECINT__
      /// @brief Create a JSON representation of this event
      marley::JSON to_json() const;

      /// @brief Replace the existing event contents with those read
      /// from a JSON representation
      void from_json(const marley::JSON& json);
      #endif

    protected:

      /// @brief Vector of pointers to each of the initial state particles
      std::vector<marley::Particle*> initial_particles_;

      /// @brief Vector of pointers to each of the final state particles
      std::vector<marley::Particle*> final_particles_;

      /// @brief Excitation energy (MeV) of the residue
      /// @note The Ex_ class member is always zero for residues that have no
      /// excited states.
      double Ex_;

    private:

      /// @brief Helper function for write_hepevt()
      /// @param p Particle to write to the HEPEVT record
      /// @param os std::ostream being written to
      /// @param track whether the particle should be marked for future
      /// tracking in a simulation (true) or not (false).
      void dump_hepevt_particle(const marley::Particle& p, std::ostream& os,
        bool track = true) const;

      void delete_particles();
  };

  // Inline function definitions
  inline double Event::Ex() const { return Ex_; }

  inline const std::vector<marley::Particle*>& Event::get_initial_particles()
    const { return initial_particles_; }

  inline std::vector<marley::Particle*>& Event::get_initial_particles()
    { return initial_particles_; }

  inline const std::vector<marley::Particle*>& Event::get_final_particles()
    const { return final_particles_; }

  inline std::vector<marley::Particle*>& Event::get_final_particles()
    { return final_particles_; }
}

inline std::ostream& operator<<(std::ostream& out, const marley::Event& e) {
  e.print(out);
  return out;
}

inline std::istream& operator>>(std::istream& in, marley::Event& e) {
  e.read(in);
  return in;
}
