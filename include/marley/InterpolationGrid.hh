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
#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include "marley/Error.hh"

namespace marley {

  /// @brief One-dimensional function y(x) defined using a grid of ordered
  /// pairs (x,y) and an interpolation rule
  /// @details This class implements the one-dimensional interpolation rules
  /// described in the <a href="https://www.bnl.gov/isd/documents/70393.pdf">
  /// ENDF-6 formats manual</a> for a grid of (x,y) pairs.
  /// @tparam FirstNumericType type for the x values
  /// @tparam SecondNumericType type for the y values
  template <typename FirstNumericType, typename SecondNumericType
    = FirstNumericType> class InterpolationGrid
  {
    public:

      using OrderedPair = std::pair<FirstNumericType, SecondNumericType>;
      using Grid = std::vector<OrderedPair>;
      using GridConstIterator
        = typename std::vector<OrderedPair>::const_iterator;

      /// @brief Method to use for interpolating between (x,y) grid points
      /// @details For @f$x_1 \leq x < x_2@f$, the table below describes
      /// how the corresponding y value is calculated using each of the allowed
      /// InterpolationMethod settings.
      /// <table>
      /// <caption id="interpolation_methods">Interpolation Methods</caption>
      /// <tr><th>Method<th>Description
      /// <tr><td>Constant<td> @f$y = y_1@f$
      /// <tr><td>LinearLinear<td> @f$y = y_1
      ///   + \frac{y_2 - y_1}{x_2 - x_1}(x - x_1) @f$
      /// <tr><td>LinearLog<td> @f$y = \exp\left[\ln(y_1)
      ///   + \frac{\ln(y_2) - \ln(y_1)}{x_2 - x_1}(x - x_1)\right] @f$
      /// <tr><td>LogLinear<td> @f$y = y_1
      ///   + \frac{y_2 - y_1}{\ln[x_2] - \ln[x_1]}[\ln(x) - \ln(x_1)] @f$
      /// <tr><td>LogLog<td> @f$y = \exp\left[\ln(y_1)
      ///   + \frac{\ln(y_2) - \ln(y_1)}{\ln(x_2) - \ln(x_1)}(\ln[x]
      ///   - \ln[x_1])\right] @f$
      /// </table>
      enum class InterpolationMethod { Constant = 1,
        LinearLinear = 2, LinearLog = 3, LogLinear = 4,
        LogLog = 5 };

      /// @brief Method to use for computing y(x) when the x value lies beyond
      /// the grid boundaries
      /// @details The table below describes how y(x) is calculated using each
      /// of the allowed ExtrapolationMethod settings. In the notation used
      /// below, the lowest x value is called @f$x_\text{left}@f$, the greatest
      /// x value is called @f$x_\text{right}@f$, and their corresponding y
      /// values are called @f$y_\text{left}@f$ and @f$y_\text{right}@f$,
      /// respectively.  The ExtrapolationMethod is only relevant when @f$x <
      /// x_\text{left}@f$ or @f$x > x_\text{right}@f$.
      /// <table>
      /// <caption id="extrapolation_methods">Extrapolation Methods</caption>
      /// <tr><th>Method<th>Description
      /// <tr><td>Zero<td> @f$y = 0@f$
      /// <tr><td>Endpoint<td> @f$y = \begin{cases} y_\text{left} &
      /// x < x_\text{left} \\ y_\text{right} & x > x_\text{right}
      /// \end{cases}@f$
      /// <tr><td>Continue<td>
      /// Use the usual InterpolationMethod formula with<br>
      /// @f$y = \begin{cases}
      /// y_1 = y_\text{left} \text{ and }
      /// y_2 = y_{\text{left}+1} &
      /// x < x_\text{left} \\ y_1 = y_{\text{right}-1} \text{ and }
      /// y_2 = y_\text{right} &
      /// x > x_\text{right}
      /// \end{cases}@f$<br>where
      /// @f$y_{\text{left}+1}@f$ is the y value corresponding
      /// to the second-lowest x value, and @f$y_{\text{right}-1}@f$
      /// is the y value corresponding to the second-highest x value.
      /// </table>
      enum class ExtrapolationMethod { Zero, Endpoint,
        Continue, Throw };

      /// @brief Create an InterpolationGrid without any grid points
      inline InterpolationGrid(InterpolationMethod interp_method
        = InterpolationMethod::LinearLinear, ExtrapolationMethod extrap_method
        = ExtrapolationMethod::Zero)
        : interpolation_method_(interp_method),
        extrapolation_method_(extrap_method)
      {
      }

      /// @brief Create an InterpolationGrid from a vector of ordered pairs
      inline InterpolationGrid(const Grid& grid,
        InterpolationMethod interp_method = InterpolationMethod::LinearLinear,
        ExtrapolationMethod extrap_method = ExtrapolationMethod::Zero)
        : interpolation_method_(interp_method),
        extrapolation_method_(extrap_method), ordered_pairs_(grid)
      {
        /// @todo Add error checks for the supplied grid
      }

      /// @brief Create an InterpolationGrid from vectors of x and y values
      inline InterpolationGrid(const std::vector<FirstNumericType>& xs,
        const std::vector<SecondNumericType>& ys,
        InterpolationMethod interp_method = InterpolationMethod::LinearLinear,
        ExtrapolationMethod extrap_method = ExtrapolationMethod::Zero)
        : interpolation_method_(interp_method),
        extrapolation_method_(extrap_method)
      {
        if (xs.size() != ys.size()) throw marley::Error(
          std::string("Vectors of x and y values passed to the constructor")
          + " of marley::InterpolationGrid have unequal sizes.");

        double old_x = marley_utils::minus_infinity;
        for (size_t j = 0; j < xs.size(); ++j) {
          double new_x = xs.at(j);
          if (new_x <= old_x) throw marley::Error(std::string("The grid")
            + " point x-values defined for a marley::InterpolationGrid object"
            + " are not strictly increasing");
          ordered_pairs_.push_back(OrderedPair(xs.at(j), ys.at(j)));
        }
      };

      /// @brief Compute y(x) using the current InterpolationMethod
      SecondNumericType interpolate(FirstNumericType x) const;

      /// @brief Add a new ordered pair (x, y) to the grid
      void insert(FirstNumericType x, SecondNumericType y);

      /// @brief Get the number of ordered pairs on the grid
      inline size_t size() const { return ordered_pairs_.size(); }

      /// @brief Delete all ordered pairs from the grid
      inline void clear() { ordered_pairs_.clear(); }

      /// @brief Get a reference to the jth ordered pair from the grid
      inline OrderedPair& at(size_t j) { return ordered_pairs_.at(j); }

      /// @brief Get a std::function object that represents y(x) for this
      /// InterpolationGrid
      inline std::function<SecondNumericType(FirstNumericType)>
        get_interpolating_function()
      {
        return [this](FirstNumericType x)
          -> SecondNumericType { return this->interpolate(x); };
        //return std::bind( &InterpolationGrid<FirstNumericType,
        //SecondNumericType>::interpolate, this, std::placeholders::_1);
      }

      /// @brief Returns a const_iterator to the first element of the grid for
      /// which the x value is not less than (i.e. greater than or equal to) x
      inline GridConstIterator lower_bound(const
        GridConstIterator& begin,
        const GridConstIterator& end, FirstNumericType x) const
      {
        return std::lower_bound(begin, end, x,
          [](const OrderedPair& pair, const FirstNumericType& f)
          -> bool { return pair.first < f; });
      }

      /// @brief Returns a const_iterator to the first element of the grid for
      /// which the x value is greater than x
      inline GridConstIterator upper_bound(const
        GridConstIterator& begin,
        const GridConstIterator& end, FirstNumericType x) const
      {
        return std::upper_bound(begin, end, x,
          [](const FirstNumericType& f, const OrderedPair& pair)
          -> bool { return f < pair.first; });
      }

      /// @brief Returns a reference to the first ordered pair
      inline const OrderedPair& front() const { return ordered_pairs_.front(); }

      /// @brief Returns a reference to the last ordered pair
      inline const OrderedPair& back() const { return ordered_pairs_.back(); }

      /// @brief Get the InterpolationMethod used by this InterpolationGrid
      inline InterpolationMethod interpolation_method() const
        { return interpolation_method_; }

      /// @brief Set the InterpolationMethod to use
      inline void set_interpolationMethod(InterpolationMethod method)
        { interpolation_method_ = method; }

      /// @brief Get the ExtrapolationMethod used by this InterpolationGrid
      inline ExtrapolationMethod extrapolation_method() const
        { return extrapolation_method_; }

      /// @brief Set the ExtrapolationMethod to use
      inline void set_interpolationMethod(ExtrapolationMethod method)
        { extrapolation_method_ = method; }

    private:

      /// @brief The method to use for interpolating between grid points
      InterpolationMethod interpolation_method_;

      /// @brief The method to use when extrapolating function values
      /// outside of the grid
      ExtrapolationMethod extrapolation_method_;

      /// @brief The ordered pairs to use as reference points for interpolation
      Grid ordered_pairs_;

      /// @brief Checks to make sure the grid contains at least two points. If
      /// it doesn't, throw an error.
      /// @details This function should be called by all class methods that
      /// rely on the grid making sense.
      inline void check_grid() const {
        /// @todo Improve this error message
        if (ordered_pairs_.size() < 2) throw marley::Error(std::string("A")
          + " class method was called for an InterpolationGrid object"
          + " that contains less than two grid points.");
      }

      /// @brief Returns true if the requested value of x lies on the grid
      /// (linear interpolation) or false if it is outside (linear
      /// extrapolation)
      inline bool find_bin_limits(FirstNumericType x,
        GridConstIterator& lower_point, GridConstIterator& upper_point) const
      {
        // Check to make sure that the grid contains at least two ordered pairs
        check_grid();

        // Find the first point on the grid that is not less than x
        bool extrapolate = false;
        GridConstIterator begin = ordered_pairs_.begin();
        GridConstIterator end = ordered_pairs_.end();
        GridConstIterator not_less_point = lower_bound(begin, end, x);

        // Check whether the requested grid point is within the grid limits
        if (not_less_point == begin) {
          lower_point = begin;
          upper_point = begin + 1;
          // First element of xs > x
          if (begin->first != x) extrapolate = true;
        }
        else if (not_less_point == end) {
          // last element of xs < x (extrapolate on the right)
          extrapolate = true;
          lower_point = end - 2;
          upper_point = end - 1;
        }
        else {
          // x is within the grid limits
          lower_point = not_less_point - 1;
          upper_point = not_less_point;
        }

        return extrapolate;
      }
  };

  template <typename FirstNumericType, typename SecondNumericType>
    SecondNumericType InterpolationGrid<FirstNumericType,
    SecondNumericType>::interpolate(FirstNumericType x) const
  {
    // Find grid points just below [(x1, y1)] and just above [(x2, y2)]
    InterpolationGrid::GridConstIterator lower_point, upper_point;
    // Includes a call to check_grid()
    bool extrapolate = find_bin_limits(x, lower_point, upper_point);

    // If the requested x value is outside of the grid, use the correct method
    // for dealing with this situation based on the value of
    if (extrapolate) {
      if (extrapolation_method_ == ExtrapolationMethod::Zero)
        return static_cast<SecondNumericType>(0.);
      else if (extrapolation_method_ == ExtrapolationMethod::Endpoint) {
        if (lower_point->first > x) return lower_point->second;
        else return upper_point->second;
      }
      else if (extrapolation_method_ == ExtrapolationMethod::Throw) {
        throw marley::Error(std::string("x = ") + std::to_string(x)
          + " lies outside of the current interpolation grid object"
          + " (which extends from x_min = "
          + std::to_string(ordered_pairs_.front().first)
          + " and x_max = " + std::to_string(ordered_pairs_.back().first)
          + ") and extrapolation is disabled.");
      }
    }

    // Either the requested point falls within the grid or
    // extrapolation_method == ExtrapolationMethod::Continue. If the "continue"
    // method for extrapolation is selected, then we continue to use the same
    // interpolation technique formula outside of the grid (for
    // InterpolationMethod::LinearLinear, this is linear-linear extrapolation).

    // If the constant interpolation method is selected, then use the *lower
    // bound* of each bin as the interpolated value. If we're extrapolating,
    // on the right, use the upper bound.
    if (interpolation_method_ == InterpolationMethod::Constant) {
      if (!extrapolate) return lower_point->second;
      else if (lower_point->first > x) return lower_point->second;
      else return upper_point->second;
    }

    // We'll use an interpolation formula for all other methods, so
    // get the grid point x and y values for later use.
    FirstNumericType x1 = lower_point->first;
    FirstNumericType x2 = upper_point->first;
    SecondNumericType y1 = lower_point->second;
    SecondNumericType y2 = upper_point->second;

    bool log_x = false, log_y = false;
    FirstNumericType x_to_use = x;
    if (interpolation_method_ == InterpolationMethod::LinearLog)
      log_x = true;
    else if (interpolation_method_ == InterpolationMethod::LogLinear)
      log_y = true;
    else if (interpolation_method_ == InterpolationMethod::LogLog) {
      log_x = true;
      log_y = true;
    }
    if (log_x) {
      x1 = std::log(x1);
      x2 = std::log(x2);
      x_to_use = std::log(x);
    }
    if (log_y) {
      y1 = std::log(y1);
      y2 = std::log(y2);
    }

    FirstNumericType y_interp = y1 + ((y2 - y1)/(x2 - x1))*(x_to_use - x1);
    if (log_y) y_interp = std::exp(y_interp);
    return y_interp;
  }

  template <typename FirstNumericType, typename SecondNumericType>
    void InterpolationGrid<FirstNumericType, SecondNumericType>::insert(
    FirstNumericType x, SecondNumericType y)
  {
    // Figure out where this grid point should go in the grid. Use
    // std::upper_bound so that entries with the same x value (discontinuities)
    // are inserted in the order that they are passed to
    // InterpolationGrid::insert.
    GridConstIterator insert_point = upper_bound(ordered_pairs_.begin(),
      ordered_pairs_.end(), x);

    // Insert the new grid point
    ordered_pairs_.insert(insert_point, OrderedPair(x, y));
  }

}
