#pragma once
#include <fstream>

namespace marley {

  class Particle {

    public:

      // Constructors
      Particle();
      Particle(int pdg_code, double m);
      Particle(int pdg_code, double m, int q);
      Particle(int pdg_code, double px, double py, double pz, double m);
      Particle(int pdg_code, double px, double py, double pz, double m, int q);
      Particle(int pdg_code, double E, double px, double py, double pz,
        double m);
      Particle(int pdg_code, double E, double px, double py, double pz,
        double m, int q);

      inline double mass() const;
      inline void set_mass(double m);

      inline double total_energy() const;
      inline void set_total_energy(double Etot);

      inline double px() const;
      inline void set_px(double px);

      inline double py() const;
      inline void set_py(double py);

      inline double pz() const;
      inline void set_pz(double pz);

      inline int pdg_code() const;

      inline double charge() const;
      inline void set_charge(int q);

      double momentum_magnitude() const;
      double kinetic_energy() const;

      void print(std::ostream& out) const;

    protected:

      // (E, px, py, pz) all in MeV
      // This should be replaced by the more modern std::array<double, 4>.
      // However, we use a C-style array here to allow ROOT 5, which is widely
      // used but not C++11 compliant, to generate dictionaries for this class
      // for TFile I/O. It'd be nice to include a default initializer here
      // (say, = {0., 0., 0., 0.}), but ROOT 5 can't handle that either.
      double four_momentum_[4];

      // Uses the particle numbering convention given by the Particle Data
      // Group (see
      // http://pdg.lbl.gov/2014/reviews/rpp2014-rev-monte-carlo-numbering.pdf)
      int pdg_code_ = 0;
      double mass_ = 0.; // MeV

      // This class member allows a marley::Particle to represent an atom or
      // ion.  The charge data member was added because the 2014 PDG particle
      // ID codes do not include a prescription for representing ionization
      // states of atoms.

      /// @todo add exceptions for cases where a marley::Particle is
      /// constructed with a charge that is inappropriate.
      // Electric charge of this particle (net charge in the case of atoms)
      // expressed as an integer multiple of the proton charge
      int charge_ = 0;
  };

  // Inline function definitions
  inline double Particle::mass() const { return mass_; }
  inline void Particle::set_mass(double m) { mass_ = m; }

  inline double Particle::total_energy() const { return four_momentum_[0]; }
  inline void Particle::set_total_energy(double Etot) { four_momentum_[0] = Etot; }

  inline double Particle::px() const { return four_momentum_[1]; }
  inline void Particle::set_px(double px) { four_momentum_[1] = px; }

  inline double Particle::py() const { return four_momentum_[2]; }
  inline void Particle::set_py(double py) { four_momentum_[2] = py; }

  inline double Particle::pz() const { return four_momentum_[3]; }
  inline void Particle::set_pz(double pz) { four_momentum_[3] = pz; }

  inline int Particle::pdg_code() const { return pdg_code_; }

  inline double Particle::charge() const { return charge_; }
  inline void Particle::set_charge(int q) { charge_ = q; }
}

// Operator for printing Particle objects to a std::ostream
inline std::ostream& operator<<(std::ostream& out, const marley::Particle& p)
{
  p.print(out);
  return out;
}
