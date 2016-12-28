#pragma once
#include <fstream>


namespace marley {

  // Forward declare the JSON class so that we can define a function that
  // creates a JSON representation of a Particle. Avoid including the JSON.hh
  // header so that we won't have issues using marley::Particle objects with
  // ROOT 5.
  class JSON;

  /// @brief Momentum four-vector for a simulated particle
  class Particle {

    public:

      Particle();

      /// @param pdg_code %Particle Data Group code
      /// @param m mass (MeV)
      Particle(int pdg_code, double m);

      /// @param pdg_code %Particle Data Group code
      /// @param m mass (MeV)
      /// @param q charge in units of the proton charge
      Particle(int pdg_code, double m, int q);

      /// @param pdg_code %Particle Data Group code
      /// @param px 3-momentum x component (MeV)
      /// @param py 3-momentum y component (MeV)
      /// @param pz 3-momentum z component (MeV)
      /// @param m mass (MeV)
      Particle(int pdg_code, double px, double py, double pz, double m);

      /// @param pdg_code %Particle Data Group code
      /// @param px 3-momentum x component (MeV)
      /// @param py 3-momentum y component (MeV)
      /// @param pz 3-momentum z component (MeV)
      /// @param m mass (MeV)
      /// @param q charge in units of the proton charge
      Particle(int pdg_code, double px, double py, double pz, double m, int q);

      /// @param pdg_code %Particle Data Group code
      /// @param E total energy (MeV)
      /// @param px 3-momentum x component (MeV)
      /// @param py 3-momentum y component (MeV)
      /// @param pz 3-momentum z component (MeV)
      /// @param m mass (MeV)
      Particle(int pdg_code, double E, double px, double py, double pz,
        double m);

      /// @param pdg_code %Particle Data Group code
      /// @param E total energy (MeV)
      /// @param px 3-momentum x component (MeV)
      /// @param py 3-momentum y component (MeV)
      /// @param pz 3-momentum z component (MeV)
      /// @param m mass (MeV)
      /// @param q charge in units of the proton charge
      Particle(int pdg_code, double E, double px, double py, double pz,
        double m, int q);

      /// @brief Get the particle's mass (MeV)
      inline double mass() const;
      /// @brief Set the particle's mass
      /// @param m mass (MeV)
      inline void set_mass(double m);

      /// @brief Get the particle's total energy (MeV)
      inline double total_energy() const;
      /// @brief Set the particle's total energy
      /// @param Etot total energy (MeV)
      inline void set_total_energy(double Etot);

      /// @brief Get the x component of the particle's 3-momentum (MeV)
      inline double px() const;
      /// @brief Set the x component of the particle's 3-momentum
      /// @param px 3-momentum x component (MeV)
      inline void set_px(double px);

      /// @brief Get the y component of the particle's 3-momentum (MeV)
      inline double py() const;
      /// @brief Set the y component of the particle's 3-momentum
      /// @param py 3-momentum y component (MeV)
      inline void set_py(double py);

      /// @brief Get the z component of the particle's 3-momentum (MeV)
      inline double pz() const;
      /// @brief Set the z component of the particle's 3-momentum
      /// @param pz 3-momentum z component (MeV)
      inline void set_pz(double pz);

      /// @brief Get the %Particle Data Group code for this particle
      inline int pdg_code() const;

      /// @brief Get the particle's charge in units of the proton charge
      inline double charge() const;
      /// @brief Set the particle's charge
      /// @param q charge in units of the proton charge
      inline void set_charge(int q);

      /// @brief Get the magnitude of the particle's 3-momentum (MeV)
      double momentum_magnitude() const;
      /// @brief Get the particle's kinetic energy (MeV)
      double kinetic_energy() const;

      /// @brief Print information about this particle to a std::ostream
      void print(std::ostream& out) const;

      /// @brief Read in this particle from a std::istream. Any previous
      /// contents of this particle will be deleted.
      void read(std::istream& in);

      /// @brief Create a JSON representation of this Particle
      marley::JSON to_json() const;

    protected:

      /// @brief momentum 4-vector for this particle (MeV)
      /// @details The order of the components is (E, px, py, pz)
      /// @note This should be replaced by the more modern
      /// std::array<double, 4>.  However, we use a C-style array here to allow
      /// ROOT 5, which is widely used but not C++11 compliant, to generate
      /// dictionaries for this class for TFile I/O. It'd be nice to include a
      /// default initializer here (say, = {0., 0., 0., 0.}), but ROOT 5 can't
      /// handle that either.
      double four_momentum_[4];

      /// @brief %Particle Data Group code identifying this particle
      /// @details See the <a href="http://tinyurl.com/hrvjjyl">full list</a>
      /// of PDG codes for details.
      int pdg_code_ = 0;

      /// @brief mass (MeV)
      /// @note We precompute the particle mass (as opposed to obtaining
      /// it from the momentum 4-vector) for speed and to avoid roundoff errors.
      /// @todo Refactor the Particle class to ensure that the mass,
      /// 3-momentum, and total energy are always consistent with each other.
      /// MARLEY enforces this by using the Particle class carefully, but
      /// consistency should be guaranteed in the Particle class internals.
      double mass_ = 0.;

      /// @brief Electric charge (net charge in the case of atoms) expressed as
      /// an integer multiple of the proton charge
      /// @details This class member allows a marley::Particle to represent an
      /// atom or ion.  The charge_ data member was added because the 2014 PDG
      /// codes do not include a prescription for representing ionization states
      /// of atoms.
      /// @todo add error handling for cases where a marley::Particle is
      /// constructed with a charge that is inappropriate.
      int charge_ = 0;
  };

  // Inline function definitions
  inline double Particle::mass() const { return mass_; }
  inline void Particle::set_mass(double m) { mass_ = m; }

  inline double Particle::total_energy() const { return four_momentum_[0]; }
  inline void Particle::set_total_energy(double Etot)
    { four_momentum_[0] = Etot; }

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

// Operator for reading in Particle objects from a std::istream
inline std::istream& operator>>(std::istream& in, marley::Particle& p)
{
  p.read(in);
  return in;
}
