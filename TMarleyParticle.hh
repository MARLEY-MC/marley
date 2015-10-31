#pragma once
#include <cmath>
#include <fstream>
#include <vector>

class TMarleyParticle {
  public:
    TMarleyParticle();
    TMarleyParticle(int id, double m);
    TMarleyParticle(int id, double p_x, double p_y, double p_z, double m);
    TMarleyParticle(int id, double E, double p_x, double p_y, double p_z, double m);
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

    std::vector<TMarleyParticle*>* get_children();

    friend std::ostream& operator<< (std::ostream& out,
      const TMarleyParticle& p);

    //inline virtual TMarleyParticle* clone() const {
    //  return new TMarleyParticle(*this);
    //}

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

// Represents an atom or ion. Since the 2014 PDG particle ID codes do not
// include a prescription for representing ionization states of atoms, this
// class adds that information as a new data member. The net electric charge of
// the atom is expressed as an integer multiple of the proton charge.
//
// TODO: add exceptions for cases where the TMarleyAtom is constructed with
// a particle ID that is inappropriate. You can also add a check for a sane
// value of the charge.
class TMarleyAtom : public TMarleyParticle {
  public:

    inline TMarleyAtom(int charge = 0): TMarleyParticle() {
      net_charge = charge;
    }

    inline TMarleyAtom(int id, double m, int charge = 0): TMarleyParticle(id, m)
    {
      net_charge = charge;
    }

    TMarleyAtom(int id, double p_x, double p_y, double p_z, double m,
      int charge = 0): TMarleyParticle(id, p_x, p_y, p_z, m)
    {
      net_charge = charge;
    }

    TMarleyAtom(int id, double E, double p_x, double p_y, double p_z,
      double m, int charge = 0): TMarleyParticle(id, E, p_x, p_y, p_z, m)
    {
      net_charge = charge;
    }

    virtual inline int get_charge() const {
      return net_charge;
    }

    virtual inline void set_charge(int charge) {
      net_charge = charge;
    }

    //inline virtual TMarleyAtom* clone() const {
    //  return new TMarleyAtom(*this);
    //}

    // TODO: remove this function when ROOT's dictionaries can handle deletion
    // without explicit new and delete calls (i.e., can use smart pointers)
    inline virtual ~TMarleyAtom() {
    }

  private:
    int net_charge;
};
