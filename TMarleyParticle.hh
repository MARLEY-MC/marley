#pragma once
#include <vector>

class TMarleyParticle {
  public:
    TMarleyParticle();
    TMarleyParticle(int id, double E, double p_x, double p_y, double p_z, double m);
    double get_total_energy() const;
    double get_mass() const;
    void add_child(TMarleyParticle* child);
    double get_px() const;
    double get_py() const;
    double get_pz() const;
    int get_id() const;
    std::vector<TMarleyParticle*>* get_children();

  private:
    double total_energy; // MeV
    double px, py, pz; // 3-momentum components
    int particle_id; // Uses the particle numbering convention given
                     // by the Particle Data Group
                     // (see http://pdg.lbl.gov/2014/reviews/rpp2014-rev-monte-carlo-numbering.pdf)
    double mass; // MeV 

    // Pointers to secondary particles created by this one
    std::vector<TMarleyParticle*> children;
};
