#pragma once
#include <string>
#include <vector>

#include "Particle.hh"

namespace marley {

  class Event {

    public:

      // Constructors
      Event(double Ex = 0.);
      Event(const marley::Particle& a,
        const marley::Particle& b, const marley::Particle& c,
        const marley::Particle& d, double Ex = 0.);

      // Destructor
      ~Event();

      // Copy constructor
      Event(const Event& other_event);

      // Move constructor
      Event(Event&& other_event);

      // Copy assignment
      Event& operator=(const Event& other_event);

      // Move assignment
      Event& operator=(Event&& other_event);

      // For a two-body reaction a(b,c)d, we may call a the target, b the
      // projectile, c the ejectile, and d the residue. When the event
      // object has been fully created, the residue will have de-excited
      // to its ground state, emitting some additional particles (e.g.,
      // gammas) which do not have special names.
      const marley::Particle& residue() const;
      const marley::Particle& ejectile() const;
      const marley::Particle& projectile() const;
      const marley::Particle& target() const;

      marley::Particle& residue();
      marley::Particle& ejectile();
      marley::Particle& projectile();
      marley::Particle& target();

      // Get const references to the vectors of initial and final particles.
      inline const std::vector<marley::Particle*>& get_initial_particles() const
        { return initial_particles_; }
      inline const std::vector<marley::Particle*>& get_final_particles() const
        { return final_particles_; }

      // Get the excitation energy of the residue just after the initial
      // two-body reaction
      inline double Ex() const;

      void add_initial_particle(const marley::Particle& p);
      void add_final_particle(const marley::Particle& p);

      // Writes a HEPEvt record for this event, using the spacetime origin (t =
      // 0 mm/c, x = 0 mm, y = 0 mm, z = 0 mm) as the initial position 4-vector
      // for all particles.
      // TODO: alter this so that the user can specify a vertex position
      // 4-vector to use.
      void write_hepevt(size_t event_num, std::ostream& out);

      void print(std::ostream& out) const;

    private:

      // Storage for the initial and final state particles
      std::vector<marley::Particle*> initial_particles_;
      std::vector<marley::Particle*> final_particles_;

      // Excitation energy (in MeV) of the residue (always zero for residues
      // that have no excited states)
      double Ex_;

      // Helper function for write_hepevt()
      void dump_hepevt_particle(const marley::Particle& p, std::ostream& os,
        bool track = true);
  };

  // Inline function definitions
  inline double Event::Ex() const { return Ex_; }
}

inline std::ostream& operator<<(std::ostream& out, const marley::Event& e) {
  e.print(out);
  return out;
}
