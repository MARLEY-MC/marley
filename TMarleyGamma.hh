#pragma once
class TMarleyLevel;

class TMarleyGamma {
  public:
    TMarleyGamma(double energy = 0, double ri = 0, TMarleyLevel* start_level = nullptr);
    void set_start_level(TMarleyLevel* start_level);
    void set_end_level(TMarleyLevel* end_level);
    TMarleyLevel* get_end_level() const;
    TMarleyLevel* get_start_level() const;
    double get_energy() const;
    double get_ri() const;

  private:
    double fEnergy;
    double fRI;
    double fCC;
    double fTI; 
    TMarleyLevel* pStartLevel;
    TMarleyLevel* pEndLevel;
};
