#pragma once
#include <vector>
#include <string>

#include "TMarleyParticle.hh"

class TMarleyReaction;

class TMarleyEvent {
  public:
    TMarleyEvent(double E_level);
    std::string create_G4_macro(); // TODO: implement this
    void print_event();
    void add_initial_particle(TMarleyParticle p);
    void add_final_particle(TMarleyParticle p);
    void set_reaction(TMarleyReaction* r);
    std::vector<TMarleyParticle>* get_initial_particles();
    std::vector<TMarleyParticle>* get_final_particles();

  private:
    std::vector<TMarleyParticle> initial_particles;
    std::vector<TMarleyParticle> final_particles;

    // Pointers to special elements of the vector of initial particles
    //TMarleyParticle* projectile, target;

    // Pointers to special elements of the vector of final particles
    //TMarleyParticle* ejectile, residue;

    // Pointer to the reaction object that created this event
    TMarleyReaction* reaction;

    // Excitation energy (in MeV) of the residue (always
    // zero for residues that have no excited states)
    double E_residue_level;
};
