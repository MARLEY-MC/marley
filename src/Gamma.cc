#include "Gamma.hh"

marley::Gamma::Gamma(double energy, double rel_intensity,
  marley::Level* start_lev, marley::Level* end_lev) : energy_(energy),
  relative_intensity_(rel_intensity), start_level_(start_lev),
  end_level_(end_lev)
{
}
