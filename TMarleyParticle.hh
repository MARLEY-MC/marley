#pragma once
#include <vector>

class TMarleyParticle {
  public:
    TMarleyParticle(int id, double E, double p_x, double p_y, double p_z, double m,
      TMarleyParticle* Parent = nullptr);
    double get_total_energy();
    double get_mass();
    void set_parent(TMarleyParticle* p);
    void add_child(TMarleyParticle* child);
    double get_px();
    double get_py();
    double get_pz();
    int get_id();
    TMarleyParticle* get_parent();
    std::vector<TMarleyParticle*>* get_children();

  private:
    double total_energy; // MeV
    double px, py, pz; // 3-momentum components
    int particle_id; // Uses the particle numbering convention given
                     // by the Particle Data Group
                     // (see http://pdg.lbl.gov/2014/reviews/rpp2014-rev-monte-carlo-numbering.pdf)
    double mass; // MeV 
    TMarleyParticle* parent; // nullptr for primary particles
    std::vector<TMarleyParticle*> children; // secondary particles created by this one
};
