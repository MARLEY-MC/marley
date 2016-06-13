#pragma once
#include <list>
#include <string>

#include "Particle.hh"

// Forward-declare the class and friend operators so that we can use the
// operators in the global scope
namespace marley { class Event; }

std::ostream& operator<< (std::ostream& out, const marley::Event& e);

namespace marley {

  class Reaction;

  class Event {
    public:
      // Constructor
      Event(double E_level = 0.);

      // The ParticleRole type is used to tell the event class
      // about what role a newly added particle plays in
      // the event
      // TODO: revert this to an "enum class" statement
      // when the ROOT dictionary generator (rootcint)
      // adds support for them (more type safe)
      enum /*class*/ ParticleRole { pr_projectile, pr_target,
        pr_ejectile, pr_residue, pr_none };

      // For a two-body reaction a(b, c)d, we may call a the target, b the
      // projectile, c the ejectile, and d the residue. When the event
      // object has been fully created, the residue will have de-excited
      // to its ground state, emitting some additional particles (e.g.,
      // gammas) which do not have special names.
      marley::Particle* get_residue();
      marley::Particle* get_ejectile();
      marley::Particle* get_projectile();
      marley::Particle* get_target();

      // Get std::lists of the initial and final particles.
      inline std::list<marley::Particle>& get_initial_particles() {
        return initial_particles;
      }
      inline std::list<marley::Particle>& get_final_particles() {
        return final_particles;
      }

      // Get the excitation energy of the residue just after the initial
      // two-body reaction
      double get_E_level();

      // Print some information about this event to std::cout
      void print_event();

      void add_initial_particle(const marley::Particle& p,
        ParticleRole r = ParticleRole::pr_none);
      void add_final_particle(const marley::Particle& p,
        ParticleRole r = ParticleRole::pr_none);
      void set_reaction(marley::Reaction* r);

      friend std::ostream& ::operator<< (std::ostream& out,
        const marley::Event& e);

      // Writes a HEPEvt record for this event, using the spacetime origin (t = 0
      // mm/c, x = 0 mm, y = 0 mm, z = 0 mm) as the initial position 4-vector for
      // all particles.
      // TODO: alter this so that the user can specify a vertex position 4-vector
      // to use.
      void write_hepevt(size_t event_num, std::ostream& out);

    protected:

       // Helper function used during event creation
      void assign_particle_pointer(marley::Particle* p,
        marley::Event::ParticleRole r);

      // Helper function for write_hepevt()
      void dump_hepevt_particle(const marley::Particle& p, std::ostream& os,
        bool track = true);

      // Lists of the initial and final state particles
      std::list<marley::Particle> initial_particles;
      std::list<marley::Particle> final_particles;

      // Pointers to special elements of the list of initial particles
      marley::Particle* projectile;
      marley::Particle* target;

      // Pointers to special elements of the list of final particles
      marley::Particle* ejectile;
      marley::Particle* residue;

      // Pointer to the reaction object that created this event
      marley::Reaction* reaction; //! Don't save to ROOT file

      // Excitation energy (in MeV) of the residue (always
      // zero for residues that have no excited states)
      double E_residue_level;
  };

}
