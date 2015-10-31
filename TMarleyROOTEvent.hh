#pragma once
#include <list>
#include <memory>
#include <string>

#include "TMarleyParticle.hh"

class TMarleyROOTEvent {
  public:
    // The ParticleRole type is used to tell the ROOT event class
    // about what role a newly added particle plays in
    // the event
    // TODO: revert this to an "enum class" statement
    // when the ROOT dictionary generator (rootcint)
    // adds support for them
    enum /*class*/ ParticleRole { pr_projectile, pr_target,
      pr_ejectile, pr_residue, pr_none };

    TMarleyROOTEvent(double E_level = 0.0);

    void add_initial_particle(const TMarleyParticle& p,
      ParticleRole r = ParticleRole::pr_none);
    void add_final_particle(const TMarleyParticle& p,
      ParticleRole r = ParticleRole::pr_none);


    inline std::list<TMarleyParticle>& get_initial_particles() {
      return initial_particles;
    }

    inline std::list<TMarleyParticle>& get_final_particles() {
      return final_particles;
    }

    inline TMarleyParticle* get_residue() { return residue; }
    inline TMarleyParticle* get_ejectile() { return ejectile; }
    inline TMarleyParticle* get_projectile() { return projectile; }
    inline TMarleyParticle* get_target() { return target; }
    inline double get_E_level() { return E_residue_level; }

  private:
    void assign_particle_pointer(TMarleyParticle* p,
      TMarleyROOTEvent::ParticleRole r);

    std::list<TMarleyParticle> initial_particles;
    std::list<TMarleyParticle> final_particles;

    // Pointers to special elements of the list of initial particles
    TMarleyParticle* projectile;
    TMarleyParticle* target;

    // Pointers to special elements of the list of final particles
    TMarleyParticle* ejectile;
    TMarleyParticle* residue;

    // Excitation energy (in MeV) of the residue (always
    // zero for residues that have no excited states)
    double E_residue_level;
};
