#pragma once
#include <list>
#include <string>

#include "TMarleyParticle.hh"

class TMarleyReaction;

//#ifdef USE_ROOT
//#include "TObject.h"
//class TMarleyEvent: public TObject {
//#else
class TMarleyEvent {
//#endif
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
    void add_initial_particle(TMarleyParticle p,
      ParticleRole r = ParticleRole::pr_none);
    void add_final_particle(TMarleyParticle p,
      ParticleRole r = ParticleRole::pr_none);
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
    TMarleyReaction* reaction; //! Don't save to ROOT file


    // Excitation energy (in MeV) of the residue (always
    // zero for residues that have no excited states)
    double E_residue_level;

  //#ifdef USE_ROOT
  //public:
  //  ClassDef(TMarleyEvent, 1);
  //#endif
};
