#pragma once
#include <functional>
#include "TMarleyDecayScheme.hh"

class TMarleyReaction {
  public:
    TMarleyReaction();
    double fermi_function(double z, double t, bool electron);
    double fermi_approx(int z, double t, bool electron);
    double ejectile_energy(double E_level, double Ea, double cos_theta_c);
    double max_level_energy(double Ea);
    double get_threshold_energy();
    double differential_xs(double E_level, double Ea, double cos_theta_c);
    double total_xs(double E_level, double Ea);
    double sample_ejectile_scattering_cosine(double E_level, double Ea);
    void set_decay_scheme(TMarleyDecayScheme* scheme);
    void create_event(double Ea); // TODO: change this to return an event object

  private:
    // TODO: change this class so that these masses (and any other
    // reaction-specific data) are read in from a data file
    // or supplied during class creation. This will allow easy re-use
    // of this class for other reaction types

    // These masses are in MeV/c^2
    double ma = 0; // nu_e
    double mb = 37224.72; // 40Ar
    double mc = 0.5109989; // e-
    double md_gs = 37226.23; // 40K
    double GF = 1; //1.16637e-11; // Fermi coupling constant (MeV^(-2)) 
    double Vud = 0.97427; // abs(V_ud) (from CKM matrix)
    double Zf = 19;
    // Lab-frame total energy of the projectile at threshold
    // for this reaction (all final-state particles at rest
    // in the CM frame)
    double Ea_threshold = ((mc + md_gs)*(mc + md_gs)
      - ma*ma - mb*mb)/(2*mb);
    TMarleyDecayScheme* ds;
};
