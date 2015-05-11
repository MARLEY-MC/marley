#pragma once

class TMarleyReaction {
  public:
    TMarleyReaction();
    double fermi_function(double z, double t, bool electron);
    double fermi_approx(int z, double t, bool electron);
    double ejectile_energy(double E_level, double Ea, double cos_theta_c);

  private:
    double ma = 0; // nu_e
    double mb = 37224.72; // 40Ar
    double mc = 0.5109989; // e-
    double md_gs = 37226.23; // 40K
};
