#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace marley {

  // Class that implements the one-dimensional interpolation rules
  // described in the ENDF-6 formats manual for a grid of (x,y) pairs.
  template <typename FirstNumericType, typename SecondNumericType
    = FirstNumericType> class InterpolationGrid
  {
    public:
  
      using OrderedPair = std::pair<FirstNumericType, SecondNumericType>;
      using Grid = std::vector<OrderedPair>;
      using GridConstIterator = typename std::vector<OrderedPair>::const_iterator;
  
      enum class InterpolationMethod { Constant = 1,
        LinearLinear = 2, LinearLog = 3, LogLinear = 4,
        LogLog = 5 };
      enum class ExtrapolationMethod { Zero, Endpoint,
        Continue, Throw };
  
      inline InterpolationGrid(InterpolationMethod interp_method
        = InterpolationMethod::LinearLinear, ExtrapolationMethod extrap_method
        = ExtrapolationMethod::Zero)
      {
        interpolation_method = interp_method;
        extrapolation_method = extrap_method;
      };
  
      inline InterpolationGrid(const Grid& grid, InterpolationMethod interp_method
        = InterpolationMethod::LinearLinear, ExtrapolationMethod extrap_method
        = ExtrapolationMethod::Zero)
      {
        interpolation_method = interp_method;
        extrapolation_method = extrap_method;
        // TODO: add error checks for the supplied grid
        ordered_pairs = grid;
      };
  
      inline InterpolationGrid(const std::vector<FirstNumericType>& xs,
        const std::vector<SecondNumericType>& ys,
        InterpolationMethod interp_method = InterpolationMethod::LinearLinear,
        ExtrapolationMethod extrap_method = ExtrapolationMethod::Zero)
      {
        interpolation_method = interp_method;
        extrapolation_method = extrap_method;
        // TODO: add error checks for the supplied vectors
        for (size_t j = 0; j < xs.size(); ++j)
          ordered_pairs.push_back(OrderedPair(xs.at(j), ys.at(j)));
      };
  
      SecondNumericType interpolate(FirstNumericType x) const;
      void insert(FirstNumericType x, SecondNumericType y);
      inline size_t size() const { return ordered_pairs.size(); }
      inline void clear() { ordered_pairs.clear(); }
      inline OrderedPair& at(size_t j) { return ordered_pairs.at(j); }
  
      // Checks to make sure the grid contains at least two points. If it
      // doesn't, throw an error. This function is called by all class methods
      // that rely on the grid making sense.
      inline void check_grid() const {
        // TODO: improve this error message
        if (ordered_pairs.size() < 2) throw std::runtime_error(std::string("A")
          + " class method was called for an InterpolationGrid object"
          + " that contains less than two grid points.");
      }
  
      inline std::function<SecondNumericType(FirstNumericType)>
        get_interpolating_function()
      {
        return std::bind(
          &InterpolationGrid<FirstNumericType, SecondNumericType>::interpolate,
          this, std::placeholders::_1);
      }
  
      inline GridConstIterator lower_bound(const
        GridConstIterator& begin,
        const GridConstIterator& end, FirstNumericType x) const
      {
        return std::lower_bound(begin, end, x,
          [](const OrderedPair& pair, const FirstNumericType& f)
          -> bool { return pair.first < f; });
      }
  
      inline GridConstIterator upper_bound(const
        GridConstIterator& begin,
        const GridConstIterator& end, FirstNumericType x) const
      {
        return std::upper_bound(begin, end, x,
          [](const FirstNumericType& f, const OrderedPair& pair)
          -> bool { return f < pair.first; });
      }
  
      inline const OrderedPair& front() const { return ordered_pairs.front(); }
      inline const OrderedPair& back() const { return ordered_pairs.back(); }
  
      inline InterpolationMethod get_interpolation_method() const {
        return interpolation_method;
      }
  
      inline void set_interpolationMethod(InterpolationMethod method) {
        interpolation_method = method;
      }
  
      inline ExtrapolationMethod get_extrapolation_method() const {
        return extrapolation_method;
      }
  
      inline void set_interpolationMethod(ExtrapolationMethod method) {
        extrapolation_method = method;
      }
  
    private:
      InterpolationMethod interpolation_method;
      ExtrapolationMethod extrapolation_method;
      Grid ordered_pairs;
  
      // Returns true if the requested value of x lies on the grid (linear interpolation)
      // or false if it is outside (linear extrapolation)
      inline bool find_bin_limits(FirstNumericType x, GridConstIterator& lower_point,
        GridConstIterator& upper_point) const
      {
        // Check to make sure that the grid contains at least two ordered pairs
        check_grid();
  
        // Find the first point on the grid that is not less than x
        bool extrapolate = false;
        GridConstIterator begin = ordered_pairs.begin();
        GridConstIterator end = ordered_pairs.end();
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
      if (extrapolation_method == ExtrapolationMethod::Zero)
        return static_cast<SecondNumericType>(0.);
      else if (extrapolation_method == ExtrapolationMethod::Endpoint) {
        if (lower_point->first > x) return lower_point->second;
        else return upper_point->second;
      }
      else if (extrapolation_method == ExtrapolationMethod::Throw) {
        throw std::runtime_error(std::string("x = ") + std::to_string(x)
          + " lies outside of the current interpolation grid object" 
          + " (which extends from x_min = "
          + std::to_string(ordered_pairs.front().first)
          + " and x_max = " + std::to_string(ordered_pairs.back().first)
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
    if (interpolation_method == InterpolationMethod::Constant) {
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
    if (interpolation_method == InterpolationMethod::LinearLog)
      log_x = true;
    else if (interpolation_method == InterpolationMethod::LogLinear)
      log_y = true;
    else if (interpolation_method == InterpolationMethod::LogLog) {
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
    // are inserted in the order that they are passed to InterpolationGrid::insert. 
    GridConstIterator insert_point = upper_bound(ordered_pairs.begin(),
      ordered_pairs.end(), x);
  
    // Insert the new grid point
    ordered_pairs.insert(insert_point, OrderedPair(x, y));
  }

}
