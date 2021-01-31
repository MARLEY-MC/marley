/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#pragma once

// standard library includes
#include <functional>
#include <vector>

// GSL includes
#include "gsl/gsl_cdf.h" // p-value for chi-squared test

// ROOT includes
#ifdef USE_ROOT
  #include "TH1D.h"
#endif

// MARLEY includes
#include "marley/Error.hh"
#include "marley/marley_utils.hh"

namespace marley {
  namespace tests {
    class Histogram;
  }
}

// Operator for writing a Histogram to a std::ostream
std::ostream& operator<<(std::ostream& out, const marley::tests::Histogram& h);

// Operator for reading in a Histogram from a std::istream
std::istream& operator>>(std::istream& in, marley::tests::Histogram& h);

namespace marley {

  namespace tests {

    // Threshold p-value below which we'll reject the null hypothesis
    // in chi-squared tests of sampling
    constexpr double SIGNIFICANCE_LEVEL = 0.01;

    // ROOT has better histogram classes, but we'd like the tests to be
    // able to run without it. This class does just enough to get us
    // what we want outside of ROOT.
    class Histogram {

      public:

        Histogram(size_t num_bins, double x_min, double x_max)
          : N_bins_( num_bins ), x_min_( x_min ), x_max_( x_max ),
          bin_counts_(num_bins, 0), num_entries_( 0 ),
          x_step_( (x_max - x_min) / static_cast<double>(num_bins) )
        {}

        inline size_t num_bins() const { return N_bins_; }
        inline size_t entries() const { return num_entries_; }
        inline double bin_width() const { return x_step_; }
        inline double x_min() const { return x_min_; }
        inline double x_max() const { return x_max_; }

        // No over/underflow bins, lowest bin is bin zero
        inline int get_bin_content(size_t b) const
          { return bin_counts_.at( b ); }

        // Get the fraction of the total entries that fall in the
        // requested bin
        inline double get_bin_fraction(size_t b) const
          { return bin_counts_.at( b ) / static_cast<double>( num_entries_ ); }

        // Get the x value corresponding to the lower edge of the bin
        inline double get_bin_left_edge(size_t b) const
          { return x_min_ + b*x_step_; }

        // Get the x value corresponding to the higher edge of the bin
        inline double get_bin_right_edge(size_t b) const
          { return x_min_ + (b + 1)*x_step_; }

        // Set the contents of the bth bin
        void set_bin_content(size_t b, int counts) {
          int old_counts = get_bin_content( b );
          num_entries_ -= old_counts;
          bin_counts_.at( b ) = counts;
          num_entries_ += counts;
        }

        inline void increment_bin(size_t b) {
          ++bin_counts_.at( b );
          ++num_entries_;
        }

        inline void fill( double x ) {
          for ( size_t b = 0; b < N_bins_; ++b ) {

            double left_edge = get_bin_left_edge( b );
            double right_edge = get_bin_right_edge( b );

            if ( x >= left_edge && x < right_edge ) {
              ++bin_counts_.at( b );
              break;
            }
          }
          // Increment whether or not a suitable bin was found (correct
          // normalization is maintained even though we don't actually
          // keep the counts for underflow/overflow)
          ++num_entries_;
        }

        void print(std::ostream& os) const {
          os << N_bins_ << ' ' << x_min_ << ' ' << x_max_ << ' ' << num_entries_
            << ' ' << x_step_ << '\n';
          for (const auto& bc : bin_counts_) {
            os << ' ' << bc;
          }
        }

        friend std::istream& operator>>(std::istream& in, Histogram& h);

        #ifdef USE_ROOT
        inline TH1D th1d(const std::string& name) const {
          TH1D hist(name.c_str(), "", N_bins_, x_min_, x_max_);
          // Note that ROOT numbers the bins (except for the underflow and overflow
          // bins) from 1 to N, while we number them from 0 to N - 1 (and don't use
          // underflow and overflow)
          for ( int b = 1; b <= N_bins_; ++b ) {
            hist.SetBinContent(b, bin_counts_.at( b - 1 ));
          }
          return hist;
        }
        #endif

        /// @brief Compares this histogram's contents to a model prediction
        /// @param expected_counts_vec A std::vector<double> loaded with
        /// pre-computed numbers of expected counts based on a model. The
        /// number of entries must be equal to the number of histogram bins,
        /// or a marley::Error will be thrown.
        /// @param[out] passed Boolean value set to true if the test was passed
        /// (p-value > our chosen significance level) or false otherwise
        /// @param[out] chi2 The @f$ \chi^2 @f$ value obtained during the test
        /// @param[out] degrees_of_freedom The number of degrees of freedom
        /// used during the test
        /// @param[out] p_value The p-value obtained during the test
        void chi2_test(const std::vector<double>& expected_counts_vec,
          bool& passed, double& chi2, int& degrees_of_freedom, double& p_value)
          const
        {
          if ( expected_counts_vec.size() != N_bins_ ) throw marley::Error("Bin"
            " number mismatch encountered in Histogram::chi2_test");
          chi2 = 0.;
          degrees_of_freedom = N_bins_ - 1;
          for ( size_t b = 0; b < N_bins_; ++b ) {
            double observed_counts = bin_counts_.at( b );
            double expected_counts = expected_counts_vec.at( b );

            // To keep the chi^2 approximation valid, skip bins where the
            // expected number of counts is less than 5
            if ( expected_counts >= 5. ) {
              chi2 += std::pow(observed_counts - expected_counts, 2)
                / expected_counts;
            }
            else --degrees_of_freedom;
          }

          p_value = gsl_cdf_chisq_Q( chi2, degrees_of_freedom );

          MARLEY_LOG_INFO() << "chi2 / DOF = " << chi2
            << " / " << degrees_of_freedom;
          MARLEY_LOG_INFO() << "p-value = " << p_value;

          passed = ( p_value >= SIGNIFICANCE_LEVEL );
        }

        /// @brief Compares this histogram's contents to a model prediction
        /// @param pdf Probability density function to use to predict bin contents
        /// @param[out] passed Boolean value set to true if the test was passed
        /// (p-value > our chosen significance level) or false otherwise
        /// @param[out] chi2 The @f$ \chi^2 @f$ value obtained during the test
        /// @param[out] degrees_of_freedom The number of degrees of freedom
        /// used during the test
        /// @param[out] p_value The p-value obtained during the test
        /// @return A vector containing the bin counts predicted by the model
        std::vector<double> chi2_test(const std::function<double(double)>& pdf,
          bool& passed, double& chi2, int& degrees_of_freedom, double& p_value)
          const
        {
          // First, build the vector of expected counts
          std::vector<double> expected_counts_vec;
          for ( size_t b = 0; b < N_bins_; ++b ) {
            double expected_counts = num_entries_ * marley_utils::num_integrate(
              pdf, get_bin_left_edge( b ), get_bin_right_edge( b ));

            expected_counts_vec.push_back( expected_counts );
          }

          // Now delegate the rest of the work to the overloaded version of
          // chi2_test()
          this->chi2_test(expected_counts_vec, passed, chi2, degrees_of_freedom,
            p_value);

          // We're done. Return the vector of expected counts.
          return expected_counts_vec;
        }

        friend std::istream& ::operator>>(std::istream& in, marley::tests::Histogram& h);

      protected:
        size_t N_bins_;
        double x_min_;
        double x_max_;
        std::vector<int> bin_counts_;
        size_t num_entries_;
        double x_step_;
    };

  } // tests namespace

} // marley namespace

