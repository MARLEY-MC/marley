#pragma once
#include <iostream>

namespace marley {

  class Level;

  class Gamma {
    public:
      Gamma(double energy = 0, double ri = 0, marley::Level* start_level = nullptr);
      void set_start_level(marley::Level* start_level);
      void set_end_level(marley::Level* end_level);
      marley::Level* get_end_level() const;
      marley::Level* get_start_level() const;
      double get_energy() const;
      double get_ri() const;
  
    private:
      double fEnergy;
      double fRI;
      marley::Level* pStartLevel;
      marley::Level* pEndLevel;
  };

}
