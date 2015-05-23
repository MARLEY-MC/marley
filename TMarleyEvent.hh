#pragma once
#include <vector>
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
    std::vector<TMarleyParticle>* get_initial_particles();
    std::vector<TMarleyParticle>* get_final_particles();
    TMarleyParticle* get_residue();

  private:
    void assign_particle_pointer(TMarleyParticle* p,
      TMarleyEvent::ParticleRole r);

    std::vector<TMarleyParticle> initial_particles;
    std::vector<TMarleyParticle> final_particles;

    // Pointers to special elements of the vector of initial particles
    TMarleyParticle* projectile;
    TMarleyParticle* target;

    // Pointers to special elements of the vector of final particles
    TMarleyParticle* ejectile;
    TMarleyParticle* residue;

    // Pointer to the reaction object that created this event
    TMarleyReaction* reaction;

    // Excitation energy (in MeV) of the residue (always
    // zero for residues that have no excited states)
    double E_residue_level;
};
