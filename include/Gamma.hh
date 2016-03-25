#pragma once
#include <iostream>

#include "iterator_to_member.hh"

namespace marley {

  class Level;

  class Gamma {
    public:
      Gamma(double energy = 0, double ri = 0,
        marley::Level* start_level = nullptr);
      void set_start_level(marley::Level* start_level);
      void set_end_level(marley::Level* end_level);
      marley::Level* get_end_level() const;
      marley::Level* get_start_level() const;
      double get_energy() const;
      double get_ri() const;

      // Convert an iterator that points to this marley::Gamma object into an
      // iterator to the gamma's relative intensity member variable. This is
      // used to load the discrete distribution in the marley::Level object
      // with the intensities of the gammas that it owns without redundant
      // storage.
      template<typename It> static inline
        itm::iterator_to_member<It, marley::Gamma, double>
        make_intensity_iterator(It it)
      {
        return itm::iterator_to_member<It, marley::Gamma, double>(it,
          &marley::Gamma::fRI);
      }
  
    private:
      double fEnergy;
      double fRI;
      marley::Level* pStartLevel;
      marley::Level* pEndLevel;
  };

}
