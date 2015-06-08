#pragma once

class TMarleyEvaporationThreshold {
  public:
    TMarleyEvaporationThreshold(double E_separation, int pid);
    double get_separation_energy() const;
    int get_fragment_pid() const;
    void set_separation_energy(double E_separation);
    void set_fragment_pid(int pid);

  private:
    double separation_energy; // MeV

    // Uses the particle numbering convention given
    // by the Particle Data Group
    // (see http://pdg.lbl.gov/2014/reviews/rpp2014-rev-monte-carlo-numbering.pdf)
    int fragment_pid;
};
