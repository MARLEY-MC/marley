#pragma once
#include <memory>
#include <random>
#include <vector>

#include "DecayChannel.hh"
#include "Fragment.hh"
#include "Generator.hh"
#include "Level.hh"
#include "Parity.hh"

namespace marley {

  class HFTable {
    public:
      inline HFTable() { total_width = 0.; }
  
      inline void add_discrete_fragment_channel(const marley::Fragment& frag,
        marley::Level& f_lev, double width)
      {
        decay_channels.push_back(std::make_unique
          <marley::DiscreteFragmentDecayChannel>(frag, f_lev));
        widths.push_back(width);
        update_dist();
        total_width += width;
      }
  
      inline void add_discrete_gamma_channel(marley::Level& f_lev, double width) {
        decay_channels.push_back(std::make_unique
          <marley::DiscreteGammaDecayChannel>(f_lev));
        widths.push_back(width);
        update_dist();
        total_width += width;
      }
  
      inline void add_continuum_fragment_channel(const marley::Fragment& frag,
        marley::Generator& gener, /*const*/ marley::SphericalOpticalModel& optmod,
        double E_min, double E_max, double mconst, double mfgs, double migs,
        double width, bool use_upper_edge = false)
      {
        decay_channels.push_back(std::make_unique
          <marley::ContinuumFragmentDecayChannel>(frag, gener, optmod, E_min,
          E_max, mconst, mfgs, migs, use_upper_edge));
        widths.push_back(width);
        update_dist();
        total_width += width;
      }
  
      inline void add_continuum_gamma_channel(marley::Generator& gener,
        int Zinitial, int Ainitial, double E_min, double E_max, double width,
        bool use_upper_edge = false)
      {
        decay_channels.push_back(std::make_unique
          <marley::ContinuumGammaDecayChannel>(gener, Zinitial, Ainitial, E_min,
          E_max, use_upper_edge));
        widths.push_back(width);
        update_dist();
        total_width += width;
      }
  
      // Update the distribution used to sample decay channels. Use the decay
      // widths as weights.
      inline void update_dist() {
        std::discrete_distribution<size_t>::param_type params(widths.begin(),
          widths.end());
        decay_dist.param(params);
      }
  
      inline marley::DecayChannel& sample_channel(marley::Generator& gen) {
        size_t c_index = gen.discrete_sample(decay_dist);
        return *decay_channels.at(c_index);
      }
  
      inline const std::vector<std::unique_ptr<marley::DecayChannel> >&
        get_channels() const
      {
        return decay_channels;
      }
  
      inline const std::vector<double>& get_widths() const {
        return widths;
      }
  
      inline double get_total_width() const {
        return total_width;
      }
  
      inline size_t get_num_channels() const {
        return decay_channels.size();
      }
  
    private:
      std::vector<std::unique_ptr<marley::DecayChannel> > decay_channels;
      std::vector<double> widths;
      std::discrete_distribution<size_t> decay_dist;
      double total_width;
  };

}
