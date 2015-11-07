#pragma once
#include <cmath>
#include <fstream>
#include <vector>

class TMarleyParticle {
  public:
    TMarleyParticle();
    TMarleyParticle(int id, double m);
    TMarleyParticle(int id, double m, int q);
    TMarleyParticle(int id, double p_x, double p_y, double p_z, double m);
    TMarleyParticle(int id, double p_x, double p_y, double p_z, double m,
      int q);
    TMarleyParticle(int id, double E, double p_x, double p_y, double p_z,
      double m);
    TMarleyParticle(int id, double E, double p_x, double p_y, double p_z,
      double m, int q);

    double get_total_energy() const;

    inline void set_total_energy(double Etot) {
      total_energy = Etot;
    }

    double get_mass() const;
    void add_child(TMarleyParticle* child);

    double get_px() const;

    inline void set_px(double p_x) {
      px = p_x;
    }

    double get_py() const;

    inline void set_py(double p_y) {
      py = p_y;
    }

    double get_pz() const;

    inline void set_pz(double p_z) {
      pz = p_z;
    }

    int get_id() const;

    inline double get_momentum_magnitude() const {
      return std::sqrt(std::pow(px, 2) + std::pow(py, 2) + std::pow(pz, 2));
    }

    inline double get_kinetic_energy() const {
      return std::max(total_energy - mass, 0.);
    }

    inline double get_charge() const {
      return charge;
    }

    inline void set_charge(int q) {
      charge = q;
    }

    inline std::vector<TMarleyParticle*>& get_children() {
      return children;
    }

    friend std::ostream& operator<< (std::ostream& out,
      const TMarleyParticle& p);

  private:
    // Helper function for the constructors
    void init(int id, double E, double p_x, double p_y, double p_z, double m,
      int q);

    double total_energy; // MeV
    double px, py, pz; // 3-momentum components
    int particle_id; // Uses the particle numbering convention given
                     // by the Particle Data Group
                     // (see http://pdg.lbl.gov/2014/reviews/rpp2014-rev-monte-carlo-numbering.pdf)
    double mass; // MeV 
    // This class member allows a TMarleyParticle to represent an atom or ion. 
    // The charge data member was added because the 2014 PDG particle ID codes do not
    // include a prescription for representing ionization states of atoms
    //
    // TODO: add exceptions for cases where the TMarleyParticle is constructed with
    // a charge that is inappropriate.
    int charge; // Electric charge of this particle (net charge in the case of atoms)
                // expressed as an integer multiple of the proton charge

    // Pointers to secondary particles created by this one
    std::vector<TMarleyParticle*> children;
};
