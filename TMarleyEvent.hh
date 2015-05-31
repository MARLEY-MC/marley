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
    enum class ParticleRole { projectile, target,
      ejectile, residue, none };

    TMarleyEvent(double E_level);
    std::string create_G4_macro(); // TODO: implement this
    void print_event();
    void add_initial_particle(TMarleyParticle p,
      ParticleRole r = ParticleRole::none);
    void add_final_particle(TMarleyParticle p,
      ParticleRole r = ParticleRole::none);
    void set_reaction(TMarleyReaction* r);
    std::list<TMarleyParticle>* get_initial_particles();
    std::list<TMarleyParticle>* get_final_particles();
    TMarleyParticle* get_residue();
    TMarleyParticle* get_ejectile();
    TMarleyParticle* get_projectile();
    TMarleyParticle* get_target();
    double get_E_level();

  private:
    void assign_particle_pointer(TMarleyParticle* p,
      TMarleyEvent::ParticleRole r);

    std::list<TMarleyParticle> initial_particles;
    std::list<TMarleyParticle> final_particles;

    // Pointers to special elements of the list of initial particles
    TMarleyParticle* projectile;
    TMarleyParticle* target;

    // Pointers to special elements of the list of final particles
    TMarleyParticle* ejectile;
    TMarleyParticle* residue;

    // Pointer to the reaction object that created this event
    TMarleyReaction* reaction;

    // Excitation energy (in MeV) of the residue (always
    // zero for residues that have no excited states)
    double E_residue_level;
};
