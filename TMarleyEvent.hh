#pragma once
#include <list>
#include <string>

#include "TMarleyParticle.hh"

class TMarleyReaction;

class TMarleyEvent {
  public:
    // The ParticleRole type is used to tell the event class
    // about what role a newly added particle plays in
    // the event
    // TODO: revert this to an "enum class" statement
    // when the ROOT dictionary generator (rootcint)
    // adds support for them
    enum /*class*/ ParticleRole { pr_projectile, pr_target,
      pr_ejectile, pr_residue, pr_none };

    TMarleyEvent(double E_level = 0.0);
    //std::string create_G4_macro(); // TODO: implement this
    void print_event();
    void add_initial_particle(const TMarleyParticle& p,
      ParticleRole r = ParticleRole::pr_none);
    void add_final_particle(const TMarleyParticle& p,
      ParticleRole r = ParticleRole::pr_none);
    void set_reaction(TMarleyReaction* r);
    inline std::list<TMarleyParticle>& get_initial_particles() {
      return initial_particles;
    }
    inline std::list<TMarleyParticle>& get_final_particles() {
      return final_particles;
    }
    TMarleyParticle* get_residue();
    TMarleyParticle* get_ejectile();
    TMarleyParticle* get_projectile();
    TMarleyParticle* get_target();
    double get_E_level();

    friend std::ostream& operator<< (std::ostream& out,
      const TMarleyEvent& e);

    // Writes a HEPEvt record for this event, using the spacetime origin (t = 0
    // mm/c, x = 0 mm, y = 0 mm, z = 0 mm) as the initial position 4-vector for
    // all particles.
    // TODO: alter this so that the user can specify a vertex position 4-vector
    // to use.
    void write_hepevt(size_t event_num, std::ostream& out);

  private:
    void assign_particle_pointer(TMarleyParticle* p,
      TMarleyEvent::ParticleRole r);

    // Helper function for write_hepevt()
    void dump_hepevt_particle(const TMarleyParticle& p, std::ostream& os,
      bool track = true);

    std::list<TMarleyParticle> initial_particles;
    std::list<TMarleyParticle> final_particles;


    // Pointers to special elements of the list of initial particles
    TMarleyParticle* projectile;
    TMarleyParticle* target;

    // Pointers to special elements of the list of final particles
    TMarleyParticle* ejectile;
    TMarleyParticle* residue;

    // Pointer to the reaction object that created this event
    TMarleyReaction* reaction; //! Don't save to ROOT file


    // Excitation energy (in MeV) of the residue (always
    // zero for residues that have no excited states)
    double E_residue_level;
};
