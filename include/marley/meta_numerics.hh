// This code was adapted with permission from the Meta.Numerics
// (http://www.meta-numerics.net) project, by David Wright. The adapted
// code is licensed under the terms of the BSD 2-Clause License
// (see the LICENSE file included in this distribution of MARLEY or
// https://opensource.org/licenses/BSD-2-Clause for the full license text).

#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <string>

#include "marley_utils.hh"
#include "Error.hh"

namespace meta_numerics {

  // Sign function taken from http://stackoverflow.com/a/4609795/4081973
  template <typename T> inline constexpr int sign(T val)
    { return (T(0) < val) - (val < T(0)); }

  // The even Bernoulli numbers B_2n = Bernoulli[n]
  // The only nonvanishing odd Bernoulli number is B_1 = -1/2, which must be
  // handled seperately if you use these numbers in any series expansion
  constexpr std::array<double, 16> Bernoulli = { 1.0, 1.0 / 6.0,
    -1.0 / 30.0, 1.0 / 42.0, -1.0 / 30.0, 5.0 / 66.0, -691.0 / 2730.0,
    7.0 / 6.0, -3617.0 / 510.0, 43867.0 / 798.0, -174611.0 / 330.0,
    854513.0 / 138.0, -236364091.0 / 2730.0, 8553103.0 / 6.0,
    -23749461029.0 / 870.0, 8615841276005.0 / 14322.0
  };

  // maximum number of iterations of a series
  constexpr int SeriesMax = 250;

  // double dedicates 52 bits to the magnitude of the mantissa, so 2^-52 is the
  // smallest fraction difference it can detect; in order to avoid any funny
  // effects at the margin, we try for one byte less, 2^-49
  const double MaxAccuracy = std::pow(2.0, -49);

  /// <summary>
  /// The Euler constant.
  /// </summary>
  /// <remarks><para>The Euler constant &#x3B3; = 0.5772...</para></remarks>
  /// <seealso href="http://en.wikipedia.org/wiki/Euler_gamma"/>
  /// <seealso href="http://mathworld.wolfram.com/Euler-MascheroniConstant.html" />
  constexpr double EulerGamma = 0.577215664901532860606512;

  // the following infrastructure is for numerical integration of ODEs
  // eventually we should expose it, but for now it is just for computing
  // Coulomb wave functions
  class OdeStepper {

    public:
      /// <summary>
      /// The current value of the independent variable.
      /// </summary>
      double X;

      /// <summary>
      /// The current value of the dependent variable.
      /// </summary>
      double Y;

      /// <summary>
      /// The current step size.
      /// </summary>
      double DeltaX;

      inline int EvaluationCount() {
        return count;
      }

      /// <summary>
      /// The target accuracy.
      /// </summary>
      double accuracy;

      inline void set_accuracy(double value) {
        if ((value < MaxAccuracy) || (value >= 1.0))
        throw marley::Error(std::string("Invalid accuracy value ")
          + std::to_string(value) + " encountered in OdeStepper::Accuracy");
        accuracy = value;
      }

      /// <summary>
      /// The right-hand side of the differential equation.
      /// </summary>
      std::function<double(double, double)> RightHandSide;

      virtual void Step() = 0;

      void Integrate(double X1);

    protected:

      int count;

      inline double Evaluate (double x, double y) {
        count++;
        return RightHandSide(x, y);
      }
  };

  class BulrischStoerStoermerStepper : public OdeStepper {

    public:

      double YPrime;

      void Step();

    private:

      static constexpr std::array<int, 12> N = { 1, 2, 3, 4, 5, 6, 7, 8, 9,
        10, 11, 12 };

      int target_k = 0;

      // do a step consisting of n mini-steps
      void TrialStep (int n, double& Y1, double& Y1P);
  };

  double Reduce (double x, double y);

  // This class handles the Lanczos approximation to the \Gamma function and
  // the correspoding approximations to associated functions.  For basic
  // background to the Lanczos approximation, see
  // http://en.wikipedia.org/wiki/Lanczos_approximation and
  // http://mathworld.wolfram.com/LanczosApproximation.html and
  // http://www.boost.org/doc/libs/1_53_0/libs/math/doc/sf_and_dist/html/math_toolkit/backgrounders/lanczos.html.
  // The basic Lanczos formula is: \Gamma(z+1) = \sqrt{2 \pi} (z + g +
  // 1/2)^(z+1/2) e^{-(z + g + 1/2)} \left[ c_0 + \frac{c_1}{z+1} +
  // \frac{c_2}{z+2} + \cdots + \frac{c_N}{z+N} \right] Given a value of g, the
  // c-values can be computed using a complicated set of matrix equations that
  // require high precision.  We write this as: \Gamma(z) = \sqrt{2 \pi} (z + g
  // - 1/2)^(z-1/2) e^{-(z + g - 1/2)} \left[ c_0 + \frac{c_1}{z} +
  // \frac{c_2}{z+1} + \cdots + \frac{c_N}{z+N-1} \right] = \sqrt{2 \pi}
  // (\frac{z + g - 1/2}{e})^{z-1/2} e^{-g} \left[ c_0 + \frac{c_1}{z} +
  // \frac{c_2}{z+1} + \cdots + \frac{c_N}{z+N-1} \right]

  class Lanczos {

    // These are listed at http://www.mrob.com/pub/ries/lanczos-gamma.html as
    // the values used by GSL, although I don't know if they still are.
    // Measured deviations at integers 2, 2, 4, 11, 1, 17, 22, 21 X 10^(-16) so
    // this is clearly worse than Godfrey's coefficients, although it does
    // manage with slightly fewer terms.
    /*
    private const double LanczosG = 7.0; private static readonly double[]
      LanczosC = double[] { 0.99999999999980993, 676.5203681218851,
      -1259.1392167224028, 771.32342877765313, -176.61502916214059,
      12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
      1.5056327351493116e-7 };
    */

    // Godfrey's coefficients, claimed relative error < 10^(-15), documented at
    // http://my.fit.edu/~gabdo/gamma.txt and in NR 3rd edition section 6.1.
    // Measured relative deviation at integers 1, 1, 4, 1, 4, 5, 6, 3 X
    // 10^(-16) so this appears about right.  These improves to 1, 1, 2, 1, 3,
    // 3, 3 X 10^(-16) when we pull the 1/e into Math.Pow(t/e, z+1/2) instead
    // of calling Math.Exp(-t) seperately.

    private:

      static constexpr double LanczosG = 607.0 / 128.0;

      static constexpr std::array<double, 15> LanczosC =
      {
        0.99999999999999709182,
        57.156235665862923517,
        -59.597960355475491248,
        14.136097974741747174,
        -0.49191381609762019978,
        3.3994649984811888699e-5,
        4.6523628927048575665e-5,
        -9.8374475304879564677e-5,
        1.5808870322491248884e-4,
        -2.1026444172410488319e-4,
        2.1743961811521264320e-4,
        -1.6431810653676389022e-4,
        8.4418223983852743293e-5,
        -2.6190838401581408670e-5,
        3.6899182659531622704e-6
      };

      // These coefficients are given by Pugh in his thesis
      // (http://web.viu.ca/pughg/phdThesis/phdThesis.pdf) table 8.5, p. 116
      // (We record them in his form; to get the usual form, multiply by
      // Exp(-\gamma+1/2) \sqrt{2} / \pi.) He claims that they "guarantee 16
      // digit floating point accuracy in the right-half plane", and since he
      // uses fewer coefficients than Godfrey that would be fantastic. But we
      // measure relative deviations at the integers of 2, 0, 35, 23, 45, 42, 6
      // X 10^(-16), making this relatively bad.
      //
      // Unfortunately, we didn't do these measurements ourselves at first, so
      // we actually used these coefficients until version 3.  Perhaps this
      // really does give 53 bits of accuracy if you do the calculation with
      // more bits. The fact that G is relatively large makes the initial
      // coefficients relatively large, which probably leads to cancellation
      // errors.
      /*
      private const double LanczosG = 10.900511;

      private static readonly double[] LanczosC = double[] {
        +2.48574089138753565546e-5, 1.05142378581721974210,
        -3.45687097222016235469, +4.51227709466894823700,
        -2.98285225323576655721, +1.05639711577126713077,
        -1.95428773191645869583e-1, +1.70970543404441224307e-2,
        -5.71926117404305781283e-4, +4.63399473359905636708e-6,
        -2.71994908488607703910e-9, };
      */

      // From LanczosG, we derive several values that we need only compute once.
      static constexpr double LanczosGP = LanczosG - 0.5;
      static const double LanczosExpG;
      static const double LanczosExpGP;

    public:

      static double Sum(double x);
      static std::complex<double> Sum(std::complex<double> z);
      static double LogSumPrime(double x);
      static std::complex<double> LogSumPrime(std::complex<double> z);
      static double Gamma(double x);
      static double LogGamma(double x);
      static std::complex<double> LogGamma(std::complex<double> z);
      static double Psi(double x);
      static std::complex<double> Psi(std::complex<double> z);
      static double Beta(double x, double y);
      static double LogBeta(double x, double y);
  };

  std::complex<double> LogGamma_Stirling(std::complex<double> z);
  std::complex<double> LogGamma(std::complex<double> z);
  std::complex<double> Psi(std::complex<double> z);

  /// <summary>
  /// Contains a pair of solutions to a differential equation.
  /// </summary>
  /// <remarks>
  /// <para>Any linear second order differential equation has two independent
  /// solutions. For example, the Bessel differential equation (<see
  /// cref="AdvancedMath.Bessel"/>) has solutions J and Y, the Coulomb wave
  /// equation has solutions F and G, and the Airy differential equation has
  /// solutions Ai and Bi.</para>
  /// <para>A solution pair structure contains values for both solutions and
  /// for their derivatives. It is often useful to have all this information
  /// together when fitting boundary conditions.</para>
  /// <para>Which solution is considered the first and which is considered the
  /// second is a matter of convention. When one solution is regular (finite)
  /// at the origin and the other is not, we take the regular solution to be
  /// the first.</para>
  /// </remarks>
  class SolutionPair {

    private:

      double j_, jPrime_, y_, yPrime_;

      public:

      /// <summary>
      /// Gets the value of the first solution.
      /// </summary>
      inline double FirstSolutionValue() {
        return j_;
      }

      inline void FirstSolutionValue(double value) {
        j_ = value;
      }

      /// <summary>
      /// Gets the derivative of the first solution.
      /// </summary>
      inline double FirstSolutionDerivative() {
        return jPrime_;
      }

      inline void FirstSolutionDerivative(double value) {
        jPrime_ = value;
      }

      /// <summary>
      /// Gets the value of the second solution.
      /// </summary>
      inline double SecondSolutionValue() {
        return y_;
      }

      inline void SecondSolutionValue(double value) {
        y_ = value;
      }

      /// <summary>
      /// Gets the derivative of the second solution.
      /// </summary>
      inline double SecondSolutionDerivative() {
        return yPrime_;
      }

      inline void SecondSolutionDerivative(double value) {
        yPrime_ = value;
      }

      // Leaving out the Wronskian for now because it can be subject to extreme
      // cancelation error.
      /*
      /// <summary>
      /// Gets the Wronsikan of the solution pair.
      /// </summary>
      /// <remarks>
      /// <para>The Wronskian of a solution pair is the product of the first
      /// solution value and the second solution derivative minus the product
      /// of the second solution value and the first solution
      /// derivative.</para>
      /// </remarks>
      public double Wronskian {
        get {
          return (j * yPrime - y * jPrime);
        }
      }
      */

      SolutionPair() {
      }

      SolutionPair (double j, double jPrime, double y, double yPrime) {
        j_ = j;
        jPrime_ = jPrime;
        y_ = y;
        yPrime_ = yPrime;
      }
  };

  // Computes the length of a right triangle's hypotenuse.
  double Hypot (double x, double y);

  // for rho < turning point, CWF are exponential; for rho > turning point, CWF
  // are oscilatory
  // we use this in several branching calculations
  double CoulombTurningPoint (double L, double eta);

  // The Gammow factor is the coefficient of the leading power of rho in the
  // expansion of the CWF near the origin It sets the order of magnitude of the
  // function near the origin. Basically F ~ C, G ~ 1/C
  double CoulombFactorZero (double eta);
  double CoulombFactor (int L, double eta);

  // each term introduces factors of rho^2 / (L+1) and 2 eta rho / (L+1), so
  // for this to converge we need rho < sqrt(X) (1 + sqrt(L)) and 2 eta rho < X
  // (1 + L); X ~ 16 gets convergence within 30 terms
  void CoulombF_Series (int L, double eta, double rho, double& F, double& FP);

  // series for L=0 for both F and G
  // this has the same convergence properties as the L != 0 series for F above
  void Coulomb_Zero_Series (double eta, double rho, double& F, double& FP,
    double& G, double& GP);

  // gives F'/F and sgn(F)
  // converges rapidly for rho < turning point; slowly for rho > turning point,
  // but still converges
  double Coulomb_CF1 (double L, double eta, double rho, int& sign);

  // computes (G' + iF')/(G + i F)
  // converges quickly for rho > turning point; does not converge at all below
  // it
  std::complex<double> Coulomb_CF2 (double L, double eta, double rho);

  // use Steed's method to compute F and G for a given L
  // the method uses a real continued fraction (1 constraint), an imaginary
  // continued fraction (2 constraints) and the Wronskian (4 constraints) to
  // compute the 4 quantities F, F', G, G'. It is reliable past the turning
  // point, but becomes slow if used far past the turning point
  SolutionPair Coulomb_Steed (double L, double eta, double rho);

  // asymptotic region
  void Coulomb_Asymptotic (double L, double eta, double rho, double& F,
    double& G);

  void Coulomb_Recurse_Upward (int L1, int L2, double eta, double rho,
    double& U, double& UP);

  double CoulombF_Integrate (int L, double eta, double rho);

  /// <summary>
  /// Computes the regular Coulomb wave function.
  /// </summary>
  /// <param name="L">The angular momentum number, which must be
  /// non-negative.</param>
  /// <param name="eta">The charge parameter, which can be postive or
  /// negative.</param>
  /// <param name="rho">The radial distance parameter, which must be
  /// non-negative.</param>
  /// <returns>The value of F<sub>L</sub>(&#x3B7;,&#x3C1;).</returns>
  /// <remarks>
  /// <para>The Coulomb wave functions are the radial wave functions of a
  /// non-relativistic particle in a Coulomb potential.</para>
  /// <para>They satisfy the differential equation:</para>
  /// <img src="../images/CoulombODE.png" />
  /// <para>A repulsive potential is represented by &#x3B7; &gt; 0, an
  /// attractive potential by &#x3B7; &lt; 0.</para>
  /// <para>F is oscilatory in the region beyond the classical turning point.
  /// In the quantum tunneling region inside the classical turning point, F is
  /// exponentially supressed.</para>
  /// <para>Many numerical libraries compute Coulomb wave functions in the
  /// quantum tunneling region using a WKB approximation, which accurately
  /// determine only the first handfull of digits; our library computes Coulomb
  /// wave functions even in this computationaly difficult region to nearly
  /// full precision -- all but the last 3-4 digits can be trusted.</para>
  /// <para>The irregular Coulomb wave functions G<sub>L</sub>(&#x3B7;,&#x3C1;)
  /// are the complementary independent solutions of the same differential
  /// equation.</para>
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="L"/> or
  /// <paramref name="rho"/> is negative.</exception>
  /// <seealso cref="CoulombG"/>
  /// <seealso href="http://en.wikipedia.org/wiki/Coulomb_wave_function" />
  /// <seealso href="http://mathworld.wolfram.com/CoulombWaveFunction.html" />
  double CoulombF (int L, double eta, double rho);

  /// <summary>
  /// Computes the irregular Coulomb wave function.
  /// </summary>
  /// <param name="L">The angular momentum number, which must be
  /// non-negative.</param>
  /// <param name="eta">The charge parameter, which can be postive or
  /// negative.</param>
  /// <param name="rho">The radial distance parameter, which must be
  /// non-negative.</param>
  /// <returns>The value of G<sub>L</sub>(&#x3B7;,&#x3C1;).</returns>
  /// <remarks>
  /// <para>For information on the Coulomb wave functions, see the remarks on
  /// <see cref="CoulombF" />.</para>
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="L"/> or
  /// <paramref name="rho"/> is negative.</exception>
  /// <seealso cref="CoulombF"/>
  /// <seealso href="http://en.wikipedia.org/wiki/Coulomb_wave_function" />
  /// <seealso href="http://mathworld.wolfram.com/CoulombWaveFunction.html" />
  double CoulombG (int L, double eta, double rho);

  // renormalized associated legendre polynomials Pe{l,m} = sqrt( (l-m)! /
  // (l+m)! ) P{l,m} unlike the unrenormalized P{l,m}, the renormalized Pe{l,m}
  // do not get too big

  // this is not quite the same renormalization used by NR; it omits a factor
  // sqrt((2l+1)/4Pi) by omitting this factor, we avoid some unnecessary
  // factors and divisions by 4Pi

  // the l-recurrsion (l-m) P{l,m} = x (2l-1) P{l-1,m} - (l+m-1) P{l-2,m}
  // becomes sqrt((l-m)(l+m)) P{l,m} = x (2l-1) P{l-1,m} - sqrt((l-1-m)(l-1+m))
  // P{l-2,m}
  // this is stable for increasing l

  // the initial value P{m,m} = (-1)^m (2m-1)!! (1-x^2)^(m/2) becomes
  // Pe{m,m} = (-1)^m (2m-1)!! sqrt( (1-x^2)^m / (2m)! ) = (-1)^m sqrt(
  // prod_{k=1}^{m} (2k-1) (1-x^2) / (2k) )

  double LegendrePe(int l, int m, double x);

  /// <summary>
  /// Computes the value of a spherical harmonic function.
  /// </summary>
  /// <param name="l">The order, which must be non-negative.</param>
  /// <param name="m">The sub-order, which must lie between -l and l
  /// inclusive.</param>
  /// <param name="theta">The azimuthal angle &#x3B8;. This angle is usually
  /// expressed as between -&#x3C0;/2 and +&#x3C0;/2, with positive values
  /// representing the upper hemisphere and negative values representing the
  /// lower hemisphere.</param>
  /// <param name="phi">The cylindrical angle &#x3C6;. This angle is usually
  /// expressed as between 0 and 2&#x3C0;, measured counter-clockwise (as seen
  /// from above) from the positive x-axis. It is also possible to use negative
  /// values to represent clockwise movement. </param>
  /// <returns>The value of Y<sub>l,m</sub>(&#x3B8;,&#x3C6;).</returns>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="l"/> is
  /// negative, or <paramref name="m"/> lies outside the range [-l,
  /// l].</exception>
  std::complex<double> SphericalHarmonic (int l, int m, double theta,
    double phi);

  // spherical harmonics are used in spherically symmetric wave functions in
  // QM, multipole expansions in EM, expansion of any function on a sphere

  // Y{l,m} = sqrt( (2l+1)/(4Pi) (l-m)!/(l+m)! ) P{l,m}; in terms of the
  // renormalized Pe{l,m}, this is Y{l,m] = sqrt( (2l+1)/(4Pi) ) Pe{l,m}

  /// <summary>
  /// Computes the regular spherical Bessel function of integer order.
  /// </summary>
  /// <param name="n">The order parameter.</param>
  /// <param name="x">The argument.</param>
  /// <returns> The value of j<sub>n</sub>(x).</returns>
  /// <remarks>
  /// <para>The spherical Bessel functions occur in solutions to the wave
  /// equations with spherical symmetry. The regular spherical Bessel functions
  /// are finite at the origin, and thus occur in situations where the wave
  /// equation is satisfied at the origin.</para>
  /// <para>The regular spherical Bessel functions are related to the regular
  /// Bessel functions of half-integer order by j<sub>n</sub>(x) =
  /// Sqrt(&#x3C0;/2x) J<sub>n+1/2</sub>(x).</para></remarks>
  /// <seealso cref="SphericalBesselY" />
  /// <seealso cref="BesselJ(double,double)"/>
  /// <seealso href="http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html" />
  double SphericalBesselJ(int n, double x);

  /// <summary>
  /// Computes the irregular spherical Bessel function of integer order.
  /// </summary>
  /// <param name="n">The order parameter.</param>
  /// <param name="x">The argument.</param>
  /// <returns>The value of y<sub>n</sub>(x).</returns>
  /// <seealso cref="SphericalBesselJ"/>
  /// <seealso href="http://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html" />
  double SphericalBesselY(int n, double x);

  double SphericalBesselJ_Zero(double x);
  double SphericalBesselJ_SeriesOne(double x);
  double SphericalBesselJ_One(double x);


  /// <summary>
  /// Computes the factorial of an integer.
  /// </summary>
  /// <param name="n">The argument, which must be non-negative.</param>
  /// <returns>The value of n!.</returns>
  /// <remarks>
  /// <para>The factorial of an integer n is the product of all integers from 1
  /// to n. For example, 4! = 4 * 3 * 2 * 1 = 24.</para>
  /// <para>n! also has a combinatorial intrepretation as the number of
  /// permutations of n objects. For example, a set of 3 objects (abc) has 3! =
  /// 6 permutations: (abc), (bac), (cba), (acb), (cab), (bca).</para>
  /// <para>Because n! grows extremely quickly with increasing n, we return the
  /// result as a double, even though the value is always an integer. (13!
  /// would overlow an int, 21! would overflow a long, 171! overflows even a
  /// double.)</para>
  /// <para>In order to deal with factorials of larger numbers, you can use the
  /// <see cref="LogFactorial"/> method, which returns accurate values of
  /// ln(n!) even for values of n for which n! would overflow a double.</para>
  /// <para>The factorial is generalized to non-integer arguments by the
  /// &#x393; function (<see cref="AdvancedMath.Gamma(double)"/>).</para>
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is
  /// negative.</exception>
  /// <seealso cref="LogFactorial"/>
  /// <seealso cref="AdvancedMath.Gamma(double)"/>
  /// <seealso href="http://en.wikipedia.org/wiki/Factorial"/>
  inline double Factorial (int n) {
    if (n < 0) throw marley::Error(std::string("Cannot compute")
      + " n! for n = " + std::to_string(n));
    return std::round(std::tgamma(n + 1));
  }

  /// <summary>
  /// Computes the double factorial of the given integer.
  /// </summary>
  /// <param name="n">The argument, which must be positive.</param>
  /// <returns>The value of n!!.</returns>
  /// <remarks>
  /// <para>The double factorial of an integer is the product all integers of
  /// the same parity, up to and including the integer.
  /// Thus 5! = 5 * 3 * 1 = 15 and 6! = 6 * 4 * 2 = 48.</para>
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is
  /// negative.</exception>
  /// <seealso href="http://mathworld.wolfram.com/DoubleFactorial.html"/>
  double DoubleFactorial(int n);

  /// <summary>
  /// Computes the natural logarithm of the double factorial of the given
  /// number.
  /// </summary>
  /// <param name="n">The argument.</param>
  /// <returns>The value of ln(n!!).</returns>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is
  /// negative.</exception>
  /// <seealso cref="DoubleFactorial"/>
  double LogDoubleFactorial(int n);

  long DoubleFactorial_Multiply(int n);

  double LogDoubleFactorial_Gamma(int n);

  double SphericalBesselJ_Series(int n, double x);
  double SphericalBesselY_Series(int n, double x);

  double SphericalBesselY_Zero(double x);
  double SphericalBesselY_SeriesOne(double x);
  double SphericalBesselY_One(double x);

  // Miller's method assumes a value at some high N and recurs downward
  // the result is then normalized using a sum relation or a known value
  double SphericalBesselJ_Miller(int n, double x);

  // Hankel's asymptotic expansions for Bessel functions (A&S 9.2)
  //   J = \sqrt{\frac{2}{\pi x}} \left[ P \cos\phi - Q \sin\phi \right]
  //   Y = \sqrt{\frac{2}{\pi x}} \left[ P \sin\phi + Q \cos\phi \right]
  // where \phi = x - \left( \nu / 2 + 1 / 4 \right) \pi and \mu = 4 \nu^2 and
  // P = 1 - \frac{(\mu-1)(\mu-9)}{2! (8x)^2} +
  // \frac{(\mu-1)(\mu-9)(\mu-25)(\mu-49}{4! (8x)^4} + \cdots
  //   Q = \frac{(\mu-1)}{8x} - \frac{(\mu-1)(\mu-9)(\mu-25)}{3! (8x)^3} +
  //   \cdots
  // Derivatives have similiar expressions
  //   J' = - \sqrt{\frac{2}{\pi x}} \left[ R \sin\phi + S \cos\phi \right]
  //   Y' = \sqrt{\frac{2}{\pi x}} \left[ R \cos\phi - S \sin\phi \right]
  // where
  //   R = 1 - \frac{(\mu-1)(\mu+15)}{2! (8x)^2} + \cdots
  //   S = \frac{(\mu+3)}{8x} - \frac{(\mu-1)(\mu - 9)(\mu+35)}{3! (8x)^3} +
  //   \cdots

  // For nu=0, this series converges to full precision in about 10 terms at
  // x~100, and in about 25 terms even as low as x~25
  // It fails to converge at all for lower x <~ 25
  // Since the first correction term is ~ (4 \nu^2)^2 / (8 x)^2 ~ (\nu^2 / 2
  // x)^2, the minimum x should grow like \nu^2 / 2
  SolutionPair Bessel_Asymptotic(double nu, double x);

  /// <summary>
  /// Computes the natural logrithm of the Gamma function.
  /// </summary>
  /// <param name="x">The argument, which must be positive.</param>
  /// <returns>The log Gamma function ln(&#x393;(x)).</returns>
  /// <remarks>
  /// <para>Because &#x393;(x) grows rapidly for increasing positive x, it is
  /// often necessary to work with its logarithm in order to avoid overflow.
  /// This function returns accurate values of ln(&#x393;(x)) even for values
  /// of x which would cause &#x393;(x) to overflow.</para>
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is
  /// negative or zero.</exception>
  /// <seealso cref="Gamma(double)" />
  double LogGamma(double x);

  double LogGamma_Stirling(double x);
  double Sum_Stirling(double x);

  /// <summary>
  /// Computes the logarithm of the factorial of an integer.
  /// </summary>
  /// <param name="n">The argument, which must be non-negative.</param>
  /// <returns>The value of ln(n!).</returns>
  /// <remarks>
  /// <para>This function provides accurate values of ln(n!) even for values of
  /// n which would cause n! to overflow.</para>
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> is
  /// negative.</exception>
  /// <seealso cref="Factorial"/>
  inline double LogFactorial (int n) {
    if (n < 0) throw marley::Error(std::string("Cannot compute")
      + " ln(n!) for n = " + std::to_string(n));
    return LogGamma(n + 1);
  }

  // orthogonal on [0,Infinity] with weight e^{-x}
  /// <summary>
  /// Computes the value of a Laguerre polynomial.
  /// </summary>
  /// <param name="n">The order, which must be non-negative.</param>
  /// <param name="x">The argument, which must be non-negative.</param>
  /// <returns>The value L<sub>n</sub>(x).</returns>
  /// <remarks>
  /// <para>Laguerre functions are orthogonal on the interval [0,+&#8734;) with
  /// the weight e<sup>-x</sup>.</para>
  /// <img src="../images/LaguerreLOrthonormality.png" />
  /// </remarks>
  /// <exception cref="ArgumentOutOfRangeException"><paramref name="n"/> or
  /// <paramref name="x"/> is negative.</exception>
  /// <seealso href="http://en.wikipedia.org/wiki/Laguerre_polynomial" />
  /// <seealso href="http://mathworld.wolfram.com/LaguerrePolynomial.html" />
  /// <seealso cref="LaguerreL(int,double,double)"/>
  double LaguerreL(int n, double x);

  /// <summary>
  /// Computes the value of an associated Laguerre polynomial.
  /// </summary>
  /// <param name="n">The order, which must be non-negative.</param>
  /// <param name="a">The associated order, which must be greater than
  /// -1.</param>
  /// <param name="x">The argument.</param>
  /// <returns>The value L<sub>n</sub><sup>a</sup>(x).</returns>
  /// <remarks>
  /// <para>The associated Laguerre polynomials are orthonogal on the interval
  /// [0,+&#8734;) with the weight x<sup>a</sup> e<sup>-x</sup>.</para>
  /// </remarks>
  /// <seealso href="http://mathworld.wolfram.com/LaguerrePolynomial.html" />
  double LaguerreL(int n, double a, double x);
}
