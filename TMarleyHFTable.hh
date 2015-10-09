#pragma once
#include <random>
#include <vector>

#include "TMarleyGenerator.hh"
#include "TMarleyHFDecayChannel.hh"

class TMarleyHFTable {
  public:
    inline void add_channel(const TMarleyHFDecayChannel& chan, double width) {

      decay_channels.push_back(chan);
      widths.push_back(width);

      // Update the distribution used to sample decay channels. Use the
      // decay widths as weights.
      std::discrete_distribution<size_t>::param_type params(widths.begin(),
        widths.end());
      decay_dist.param(params);
    }

    inline const TMarleyHFDecayChannel& sample_channel(TMarleyGenerator& gen) {
      size_t c_index = gen.discrete_sample(decay_dist);
      return decay_channels.at(c_index);
    }

    inline const std::vector<TMarleyHFDecayChannel>& get_channels() const {
      return decay_channels;
    }

    inline const std::vector<double>& get_widths() const {
      return widths;
    }

  private:
    std::vector<TMarleyHFDecayChannel> decay_channels;
    std::vector<double> widths;
    std::discrete_distribution<size_t> decay_dist;
};
