// This file contains a C++ translation of code for computing the Coulomb wavefunctions
// (see, e.g, http://dlmf.nist.gov/33.2 and
// http://mathworld.wolfram.com/CoulombWaveFunction.html) excerpted from the Meta Numerics
// library (http://metanumerics.codeplex.com/), which was originally written in C#. The
// contents of this file and its accompanying header file (meta_numerics.hh), but *not*
// the rest of MARLEY, are distributed under the terms of the Microsoft Public License
// (Ms-PL). The full text of the license is reproduced at the end of this file, and it may
// also be viewed at https://metanumerics.codeplex.com/license.
//
// Original Meta Numerics Library Copyright (c) 2008-2015 David Wright
// C++ Translation and Adaptation for MARLEY Copyright (c) 2015 Steven Gardiner

// -- Begin Ms-PL licensed code
#include "meta_numerics.hh"

const std::vector<int> meta_numerics::BulrischStoerStoermerStepper::N = { 1,
  2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
};

// do a step consisting of n mini-steps
void meta_numerics::BulrischStoerStoermerStepper::TrialStep (int n, double& Y1,
  double& Y1P)
{

  // this is Stoermer's rule for 2nd order conservative equations
  double h = DeltaX / n;
  
  Y1 = Y;
  double D1 = h * (YPrime + h * Evaluate(X, Y) / 2.0);
  
  for (int k = 1; k < n; k++) {
      Y1 += D1;
      D1 += h * h * Evaluate(X + k * h, Y1);
  }
  
  Y1 += D1;
  
  Y1P = D1 / h + h * Evaluate(X + DeltaX, Y1) / 2.0;

}

void meta_numerics::BulrischStoerStoermerStepper::Step() {

    // a step consists of trial steps with different numbers of intermediate points (substep sizes)
    // the values obtained using different points are recorded and extrapolated to an infinite number
    // of points (zero substep size)

    // we store the values in a tableau whoose first row contains the measured values and
    // whoose lower rows contain values extrapolated using different degree polynomials

    // y_1    y_2    y_3    y_4
    // y_12   y_23   y_34
    // y_123  y_234
    // y_1234

    // Neville's algorithm is used to fill out this tableau

    // initialize the tableau (for both variable and derivative)
    double** T = new double* [N.size()];
    double** U = new double* [N.size()];
    //double[][] T = double[N.size()][];
    //double[][] U = double[N.size()][];

    T[0] = new double[N.size()]; U[0] = new double[N.size()];
    TrialStep(N[0], T[0][0], U[0][0]);

    // keep track of total number of evaluations
    int A = N[0];


    // set the window
    size_t kMin, kMax;
    if (target_k < 1) {
        kMin = 1;
        kMax = N.size() - 1;
    } else {
        kMin = target_k - 1;
        if (kMin < 1) kMin = 1;
        kMax = target_k + 1;
        if (kMax > N.size() - 1) kMax = N.size() - 1;
    }
    double target_work_per_step = std::numeric_limits<double>::max();
    double target_expansion_factor = 1.0;

    // try different substep sizes
    for (size_t k = 1; k <= kMax; k++) {

        // add a row to the tableau
        T[k] = new double[N.size() - k]; U[k] = new double[N.size() - k];

        // perform the substep
        TrialStep(N[k], T[0][k], U[0][k]);
        A += N[k];

        // fill out entries in the tableau
        for (size_t j = 1; j <= k; j++) {
            double x = 1.0 * N[k] / N[k - j];
            T[j][k - j] = T[j - 1][k - j + 1] + (T[j - 1][k - j + 1] - T[j - 1][k - j]) / ((x + 1.0) * (x - 1.0));
            U[j][k - j] = U[j - 1][k - j + 1] + (U[j - 1][k - j + 1] - U[j - 1][k - j]) / ((x + 1.0) * (x - 1.0));
        }

        // check for convergence and predict work in target window
        if ((k >= kMin) && (k <= kMax)) {

            double absolute_error = std::abs(T[k][0] - T[k - 1][0]);
            double relative_error = std::abs(absolute_error / T[k][0]);
            double expansion_factor = std::pow(accuracy / relative_error, 1.0 / (2 * N[k]));
            double work_per_step = A / expansion_factor;

            if (work_per_step < target_work_per_step) {
                target_k = k;
                target_work_per_step = work_per_step;
                target_expansion_factor = expansion_factor;
            }

        }

    }

    if (std::abs(T[kMax][0] - T[kMax - 1][0]) <= accuracy * std::abs(T[kMax][0])) {
        // converged
        X = X + DeltaX;
        Y = T[kMax][0];
        YPrime = U[kMax][0];
        if (target_expansion_factor > 2.0) target_expansion_factor = 2.0;
        if (target_expansion_factor < 0.25) target_expansion_factor = 0.25;
    } else {
        // didn't converge
        if (target_expansion_factor > 0.5) target_expansion_factor = 0.5;
        if (target_expansion_factor < 0.0625) target_expansion_factor = 0.0625;
    }

    DeltaX = DeltaX * target_expansion_factor;

    delete[] T; delete[] U;
}

void meta_numerics::OdeStepper::Integrate (double X1) {

  double X0 = X;

  // reverse direction, if necessary
  if (sign(DeltaX) != sign(X1 - X0)) DeltaX = -DeltaX;

  // we can't just check (X < X1) because sometimes we integrate the other way
  // so instead check that "we are on the same side of X1 as X0"
  while (sign(X-X1) == sign(X0-X1)) {

    // if we would overshoot in the next step, reduce it
    if (sign(X + DeltaX-X1) != sign(X0 - X1)) DeltaX = X1 - X;

    Step();

  }

}

double meta_numerics::Reduce(double x, double y) {

  double t = x + marley_utils::two_pi * y;
  return t;
  // This technique was used in the original MetaNumerics library. However, it relies
  // on a decimal type which might be hard to get in standard C++. For now, skip the
  // reduction. If needed, add this in later.
  //
  //
  // reduces an argument to its corresponding argument between -2 Pi < x < 2 Pi
  //if ((Math.Abs(t) < 64.0) || (Math.Abs(t) > dmax)) {
  //  // if the argument is small we don't need the high accurary reduction
  //  // if the argument is too big, we can't do the high accuracy reduction because it would overflow a decimal vairable
  //  return (t);
  //} else {
  //  // otherwise, convert to decimal, subtract a multiple of 2 Pi, and return

  //  // reduce x by factors of 2 Pi
  //  decimal dx = Convert.ToDecimal(x);
  //  decimal dn = Decimal.Truncate(dx / dPI2);
  //  dx = dx - dn * dPI2;

  //  // reduce y by factors of 1
  //  decimal dy = Convert.ToDecimal(y);
  //  decimal dm = Decimal.Truncate(dy / 1.0m);
  //  dy = dy - dm * 1.0m;

  //  // form the argument
  //  decimal dt = dx + dy * dPI2;
  //  return (Convert.ToDouble(dt));

  //}
}

double meta_numerics::Lanczos::Sum (double x) {
  double s = LanczosC[0] + LanczosC[1] / x;
  for (size_t i = 2; i < LanczosC.size(); ++i) {
    x += 1.0;
    s += LanczosC[i] / x;
  }
  return s;
}

std::complex<double> meta_numerics::Lanczos::Sum (std::complex<double> z) {
  std::complex<double> s = LanczosC[0] + LanczosC[1] / z;
  for (size_t i = 2; i < LanczosC.size(); i++) {
    z += 1.0;
    s += LanczosC[i] / z;
  }
  return s;
}

double meta_numerics::Lanczos::LogSumPrime (double x) {
  double q = LanczosC[0] + LanczosC[1] / x;
  double p = LanczosC[1] / (x * x);
  for (size_t i = 2; i < LanczosC.size(); i++) {
    x += 1.0;
    q += LanczosC[i] / x;
    p += LanczosC[i] / (x * x);
  }
  return (-p / q);
}

std::complex<double> meta_numerics::Lanczos::LogSumPrime (std::complex<double> z) {
  std::complex<double> q = LanczosC[0] + LanczosC[1] / z;
  std::complex<double> p = LanczosC[1] / (z * z);
  for (size_t i = 2; i < LanczosC.size(); i++) {
    z += 1.0;
    q += LanczosC[i] / z;
    p += LanczosC[i] / (z * z);
  }
  return (-p / q);
}

double meta_numerics::Lanczos::Gamma (double x) {
  double t = x + LanczosGP;
  return marley_utils::sqrt_two_pi *
    std::pow(t / std::exp(1), x - 0.5) * LanczosExpG *
    Sum(x);
}

double meta_numerics::Lanczos::LogGamma (double x) {
  double t = x + LanczosGP;
  return std::log(marley_utils::sqrt_two_pi * Sum(x)) +
    (x - 0.5) * std::log(t) - t;
}

std::complex<double> meta_numerics::Lanczos::LogGamma (std::complex<double> z) {
  std::complex<double> t = z + std::complex<double>(LanczosGP, 0);
  return std::log(marley_utils::sqrt_two_pi) +
    (z - 0.5) * std::log(t) - t +
    std::log(Sum(z));
}

double meta_numerics::Lanczos::Psi (double x) {
  double t = x + LanczosGP;
  return (std::log(t) - LanczosG / t + LogSumPrime(x));
}

std::complex<double> meta_numerics::Lanczos::Psi (std::complex<double> z) {
  std::complex<double> t = z + std::complex<double>(LanczosGP, 0);
  return (std::log(t) - std::complex<double>(LanczosG, 0) / t + LogSumPrime(z));
}

// If we just compute Exp( LogGamma(x) + LogGamma(y) - LogGamma(x+y) ) then several leading terms in the sum cancel,
// potentially introducing cancelation error. So we write out the ratios explicitly and take the opportunity
// to write the result in terms of some naturally occuring ratios.
double meta_numerics::Lanczos::Beta (double x, double y) {
  double tx = x + LanczosGP;
  double ty = y + LanczosGP;
  double txy = x + y + LanczosGP;
  return marley_utils::sqrt_two_pi * LanczosExpGP *
    std::pow(tx / txy, x) * std::pow(ty / txy, y) * std::sqrt(txy / tx / ty) *
    Sum(x) * Sum(y) / Sum(x + y);
}

double meta_numerics::Lanczos::LogBeta (double x, double y) {
  double tx = x + LanczosGP;
  double ty = y + LanczosGP;
  double txy = x + y + LanczosGP;
  return std::log(marley_utils::two_pi / txy) / 2.0 + (x - 0.5) * std::log( tx / txy) + (y - 0.5) * std::log(ty / txy) +
    std::log(LanczosExpGP * Sum(x) * Sum(y) / Sum(x + y));
}

const std::vector<double> meta_numerics::Lanczos::LanczosC = {
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

std::complex<double> meta_numerics::LogGamma_Stirling (std::complex<double> z) {

  // work in the upper complex plane; I think this isn't actually necessary
  if (z.imag() < 0.) return std::conj(
    meta_numerics::LogGamma_Stirling(std::conj(z)));

  std::complex<double> f = (z - 0.5) * std::log(z) - z
    + std::log(marley_utils::two_pi) / 2.0;

  // reduce f.Im modulo 2*PI
  // result is cyclic in f.Im modulo 2*PI, but if f.Im starts off too big, the corrections
  // applied below will be lost because they are being added to a big number
  f = std::complex<double>(f.real(), meta_numerics::Reduce(f.imag(), 0.0));

  std::complex<double> zz = z * z;
  std::complex<double> zp = z;
  for (size_t i = 1; i < meta_numerics::Bernoulli.size(); i++) {
    std::complex<double> f_old = f;
    f += meta_numerics::Bernoulli[i] / (2 * i) / (2 * i - 1) / zp;
    if (f == f_old) return f;
    zp *= zz;
  }
  throw std::runtime_error(std::string("LogGamma_Stirling failed to converge."));
}


/// <summary>
/// Compute the complex log Gamma function.
/// </summary>
/// <param name="z">The complex argument, which must have a non-negative real part.</param>
/// <returns>The complex value ln(&#x393;(z)).</returns>
/// <exception cref="ArgumentOutOfRangeException">The real part of <paramref name="z"/> is negative.</exception>
/// <seealso cref="AdvancedMath.LogGamma" />
std::complex<double> meta_numerics::LogGamma (std::complex<double> z) {
  if (z.real() < 0.0) {
    throw std::runtime_error(std::string("Argument z = (") + std::to_string(z.real())
      + ", " + std::to_string(z.imag()) + ") out of range in LogGamma");
  } else if (std::abs(z) < 16.0) {
    return (meta_numerics::Lanczos::LogGamma(z));
  } else {
    return (meta_numerics::LogGamma_Stirling(z));
  }
}

/// <summary>
/// Computes the complex digamma (&#x3C8;) function.
/// </summary>
/// <param name="z">The complex argument.</param>
/// <returns>The value of &#x3C8;(z).</returns>
/// <remarks>
/// <para>The image below shows the complex &#x3C8; function near the origin using domain coloring.</para>
/// <img src="../images/ComplexPsiPlot.png" />
/// </remarks>
/// <seealso cref="AdvancedMath.Psi(double)" />
std::complex<double> meta_numerics::Psi (std::complex<double> z) {
  if (z.real() < 0.5) {
    // reduce Re(z) in order to handle large real values!
    return (Psi(1.0 - z) - marley_utils::pi / std::tan(marley_utils::pi * z));
  } else {
    // add Stirling for large z
    return (meta_numerics::Lanczos::Psi(z));
  }
}

/// <summary>
/// Computes the length of a right triangle's hypotenuse.
/// </summary>
/// <param name="x">The length of one side.</param>
/// <param name="y">The length of another side.</param>
/// <returns>The length of the hypotenuse, sqrt(x<sup>2</sup> + y<sup>2</sup>).</returns>
/// <remarks>
/// <para>The length is computed accurately, even in cases where
/// x<sup>2</sup> or y<sup>2</sup> would overflow.</para>
/// </remarks>
double meta_numerics::Hypot (double x, double y) {
  if ((x == 0.0) && (y == 0.0)) {
    return (0.0);
  } else {
    double ax = std::abs(x);
    double ay = std::abs(y);
    if (ax > ay) {
      double r = y / x;
      return (ax * std::sqrt(1.0 + r * r));
    } else {
      double r = x / y;
      return (ay * std::sqrt(1.0 + r * r));
    }
  }
}

// for rho < turning point, CWF are exponential; for rho > turning point, CWF are oscilatory
// we use this in several branching calculations
double meta_numerics::CoulombTurningPoint (double L, double eta) {

  double p = L * (L + 1);
  double q = std::sqrt(p + eta * eta);

  if (eta >= 0.0) {
    return (q + eta);
  } else {
    return (p / (q - eta));
  }

}

// The Gammow factor is the coefficient of the leading power of rho in the expansion of
// the CWF near the origin It sets the order of magnitude of the function near the origin.
// Basically F ~ C, G ~ 1/C
double meta_numerics::CoulombFactorZero (double eta) {

  double x = marley_utils::two_pi * eta;

  if (std::abs(x) < 0.1) {

    // series for x / (e^x - 1)

    double x2 = x * x;
    double x4 = x2 * x2;
    double x6 = x4 * x2;
    double x8 = x6 * x2;

    return (std::sqrt(1.0 - x / 2.0 + x2 / 12.0 - x4 / 720.0 + x6 / 30240.0 - x8 / 1209600.0));

    // can we do this in some extensible way? not easily; the coefficients involve Bernoulli numbers, which
    // can be generated by a fromula, but it involves binomial coefficients and all previous Bernoulli numbers

  } else {

    return (std::sqrt(x / (std::exp(x) - 1.0)));

  }

}

double meta_numerics::CoulombFactor (int L, double eta) {

  // factor for L=0
  double C = CoulombFactorZero(eta);

  // use the recurrsion
  //   C_{L+1} = \frac{\sqrt{ L^2 + eta^2 }}{ L (2L + 1)} C_L
  // it would be better not to use this for really high L; is there another approximation we can use?
  //double abseta = std::abs(eta);
  for (int k = 1; k <= L; k++) {
    C *= Hypot(k, eta) / k / (2 * k + 1);
    /*
    if (k > abseta) {
      double x = abseta / k;
      C = std::sqrt(1.0 + x * x) / (2 * k + 1) * C;
    } else {
      double x = k / abseta;
      C = std::sqrt(1.0 + x * x) / x / (2 * k + 1) * C;
    }
    */
    //C = std::sqrt(k * k + eta_squared) / (k * (2 * k + 1)) * C;
  }

  return C;
}

// each term introduces factors of rho^2 / (L+1) and 2 eta rho / (L+1), so for this to
// converge we need rho < sqrt(X) (1 + sqrt(L)) and 2 eta rho < X (1 + L); X ~ 16 gets
// convergence within 30 terms

void meta_numerics::CoulombF_Series (int L, double eta, double rho, double& F,
  double& FP)
{

  double eta_rho = eta * rho;
  double rho_2 = rho * rho;

  double u0 = 1.0;
  double u1 = eta_rho / (L + 1);
  double u = u0 + u1;
  double v = (L + 1) * u0 + (L + 2) * u1;

  for (int k = 2; k < meta_numerics::SeriesMax; k++) {

    double u2 = (2.0 * eta_rho * u1 - rho_2 * u0) / k / (2 * L + k + 1);
    double v2 = (L + 1 + k) * u2;

    double u_old = u;
    u += u2;
    v += v2;

    if ((k % 2 == 0) && (u == u_old)) {
      double C = meta_numerics::CoulombFactor(L, eta);
      F = C * std::pow(rho, L + 1) * u;
      FP = C * std::pow(rho, L) * v;
      return;
    }

    u0 = u1; u1 = u2;

  }

  throw std::runtime_error(std::string("CoulombF_series failed to converge."));

}

// series for L=0 for both F and G
// this has the same convergence properties as the L != 0 series for F above

void meta_numerics::Coulomb_Zero_Series (double eta, double rho, double& F,
  double& FP, double& G, double& GP)
{

  double eta_rho = eta * rho;
  double rho_2 = rho * rho;

  double u0 = 0.0;
  double u1 = rho;
  double u = u0 + u1;
  double up = u1;

  double v0 = 1.0;
  double v1 = 0.0;
  double v = v0 + v1;

  for (int n = 2; n <= meta_numerics::SeriesMax; n++) {

    double u2 = (2.0 * eta_rho * u1 - rho_2 * u0) / n / (n - 1);
    double v2 = (2.0 * eta_rho * v1 - rho_2 * v0 - 2.0 * eta * (2 * n - 1) * u2) / n / (n - 1);

    double u_old = u; u += u2; up += n * u2;
    double v_old = v; v += v2;

    if ((u == u_old) && (v == v_old)) {

      double C = meta_numerics::CoulombFactorZero(eta);
      F = C * u;

      FP = C * up / rho;

      double r = Psi(std::complex<double>(1.0, eta)).real() + 2.0 * EulerGamma - 1;
      G = (v + 2.0 * eta * u * (std::log(2.0 * rho) + r)) / C;

      GP = (FP * G - 1.0) / F;

      return;
    }

    u0 = u1; u1 = u2; v0 = v1; v1 = v2;

  }

  throw std::runtime_error(std::string("Coulomb_Zero_Series failed to converge."));

}

// gives F'/F and sgn(F)
// converges rapidly for rho < turning point; slowly for rho > turning point, but still converges
double meta_numerics::Coulomb_CF1 (double L, double eta, double rho, int& sign)
{

  // maximum iterations
  int nmax = meta_numerics::SeriesMax;
  double rho0 = meta_numerics::CoulombTurningPoint(L, eta);
  if (rho > rho0) nmax += static_cast<int>(std::floor(2.0 * (rho - rho0)));

  // use Wallis method of continued fraction evalution

  double f = (L + 1.0) / rho + eta / (L + 1.0);
  sign = 1;

  double A0 = 1.0;
  double A1 = f;
  double B0 = 0.0;
  double B1 = 1.0;

  for (int n = 1; n < nmax; n++) {

    double f_old = f;

    // compute next term
    double k = L + n;
    double t = eta / k;
    double a = -(1.0 + t * t);
    double b = (2.0 * k + 1.0) * (1.0 / rho + t / (k + 1.0));

    // apply it
    double A2 = b * A1 + a * A0;
    double B2 = b * B1 + a * B0;
    // note B1 = 1 always, A1 = f always
    f = A2 / B2;
    if (B2 < 0) sign = -sign;

    // check for convergence
    if (f == f_old) {
      return f;
    }

    // renormalize by dividing by B2 and prepare for the next cycle
    A0 = A1 / B2;
    A1 = f;
    B0 = B1 / B2;
    B1 = 1.0;
  }

  throw std::runtime_error(std::string("Coulomb_CF1 failed to converge."));

}

// computes (G' + iF')/(G + i F)
// converges quickly for rho > turning point; does not converge at all below it
std::complex<double> meta_numerics::Coulomb_CF2 (double L, double eta,
  double rho)
{

  std::complex<double> a = std::complex<double>(1.0 + L, eta);
  std::complex<double> c = std::complex<double>(-L, eta);

  std::complex<double> D = 1.0 / (std::complex<double>(2.0 * (rho - eta), 2.0));
  std::complex<double> Df = a * c * D;
  std::complex<double> f = Df;

  int nmax = meta_numerics::SeriesMax;
  if (eta < 0) nmax += static_cast<int>(std::floor(-2.0 * eta));

  for (int n = 1; n < nmax; n++) {

    std::complex<double> f_old = f;

    std::complex<double> N(n, 0);
    std::complex<double> p = (a + N) * (c + N);
    std::complex<double> q = std::complex<double>(2.0 * (rho - eta), 2.0 * (n + 1));

    D = 1.0 / (q + p * D);
    Df = (q * D - 1.0) * Df;
    f += Df;

    if (f == f_old) {
      return (marley_utils::i * f / rho + std::complex<double>(0.0, 1.0 - eta / rho));
    }

  }

  throw std::runtime_error(std::string("Coulomb_CF2 failed to converge."));

  // don't use Wallis algorithm for this continued fraction! it appears that sometimes noise in what
  // should be the last few iterations prevents convergence; Steed's method appears to do better
  // an example is L=0, eta=0.1, rho=14.0
}

// use Steed's method to compute F and G for a given L
// the method uses a real continued fraction (1 constraint), an imaginary continued fraction (2 constraints)
// and the Wronskian (4 constraints) to compute the 4 quantities F, F', G, G'
// it is reliable past the turning point, but becomes slow if used far past the turning point
meta_numerics::SolutionPair meta_numerics::Coulomb_Steed (double L, double eta,
  double rho)
{

  // compute CF1 (F'/F)
  int sign;
  double f = meta_numerics::Coulomb_CF1(L, eta, rho, sign);

  // compute CF2 ((G' + iF')/(G + i F))
  std::complex<double> z = meta_numerics::Coulomb_CF2(L, eta, rho);
  double p = z.real();
  double q = z.imag();

  // use CF1, CF2, and Wronskian (FG' - GF' = 1) to solve for F, F', G, G' 
  double g = (f - p) / q;

  meta_numerics::SolutionPair result;
  result.FirstSolutionValue(sign / std::sqrt(g * g * q + q));
  result.FirstSolutionDerivative(f * result.FirstSolutionValue());
  result.SecondSolutionValue(g * result.FirstSolutionValue());
  result.SecondSolutionDerivative((p * g - q) * result.FirstSolutionValue());
  return result;

}


// asymptotic region
void meta_numerics::Coulomb_Asymptotic (double L, double eta, double rho,
  double& F, double& G)
{

  // compute phase
  // reducing the eta = 0 and eta != 0 parts seperately preserves accuracy for large rho and small eta
  double t0 = Reduce(rho, -L / 4.0);
  double t1 = Reduce(LogGamma(std::complex<double>(L + 1.0, eta)).imag()
    - eta * std::log(2.0 * rho), 0.0);
  double t = t0 + t1;
  double s = std::sin(t);
  double c = std::cos(t);

  // determine the weights of sin and cos
  double f0 = 1.0;
  double g0 = 0.0;
  double f = f0;
  double g = g0;
  for (int k = 0; true; k++) {

    // compute the next contributions to f and g
    double q = 2 * (k + 1) * rho;
    double a = (2 * k + 1) * eta / q;
    double b = ((L * (L + 1) - k * (k + 1)) + eta * eta) / q;
    double f1 = a * f0 - b * g0;
    double g1 = a * g0 + b * f0;

    // add them
    double f_old = f;
    f += f1;
    double g_old = g;
    g += g1;

    if ((f == f_old) && (g == g_old)) break;

    // check for non-convergence
    if (k > meta_numerics::SeriesMax) throw
      std::runtime_error(std::string("Coulomb_Asymptotic failed to converge."));

    // prepare for the next iteration
    f0 = f1;
    g0 = g1;

  }

  F = g * c + f * s;
  G = f * c - g * s;

}

void meta_numerics::Coulomb_Recurse_Upward (int L1, int L2, double eta,
  double rho, double& U, double& UP)
{

  if (L2 < L1) throw std::runtime_error(std::string("L2 < L1 in Coulomb_Recurse_Upward"));

  for (int K = L1 + 1; K <= L2; K++) {
    
    // compute some factors
    double S = std::sqrt(K * K + eta * eta);
    double T = K * K / rho + eta;

    // compute next higher function and derivative
    double U2 = (T * U - K * UP) / S;
    double UP2 = (S * U - T * U2) / K;

    // prepare for next iteration
    U = U2;
    UP = UP2;

  }

}

double meta_numerics::CoulombF_Integrate(int L, double eta, double rho) {

  // start at the series limit
  double rho0 = 4.0 + 2.0 * std::sqrt(L);
  if (std::abs(rho0 * eta) > 8.0 + 4.0 * L) rho0 = (8.0 + 4.0 * L) / std::abs(eta);
  double F, FP;
  CoulombF_Series(L, eta, rho0, F, FP);

  // TODO: switch so we integrate w/o the C factor, then apply it afterward
  if ((F == 0.0) && (FP == 0.0)) return (0.0);
  
  // integrate out to rho
  BulrischStoerStoermerStepper s = BulrischStoerStoermerStepper();
  s.RightHandSide = [eta, L] (double x, double U) -> double {
    return ((L * (L + 1) / x / x + 2.0 * eta / x - 1.0) * U);
  };
  s.X = rho0;
  s.Y = F;
  s.YPrime = FP;
  s.DeltaX = 0.25;
  s.set_accuracy(2.5E-13);
  s.Integrate(rho);

  // return the result
  return (s.Y);
}

double meta_numerics::CoulombF (int L, double eta, double rho) {

  if (L < 0) throw std::runtime_error(std::string("L < 0 in CoulombF"));
  if (rho < 0) throw std::runtime_error(std::string("rho < 0 in CoulombF"));

  if ((rho < 4.0 + 2.0 * std::sqrt(L)) && (std::abs(rho * eta) < 8.0  + 4.0 * L)) {
    // if rho and rho * eta are small enough, use the series expansion at the origin
    double F, FP;
    meta_numerics::CoulombF_Series(L, eta, rho, F, FP);
    return F;
  } else if (rho > 32.0 + (L * L + eta * eta) / 2.0) {
    // if rho is large enrough, use the asymptotic expansion
    double F, G;
    meta_numerics::Coulomb_Asymptotic(L, eta, rho, F, G);
    return F;
  } else {
    // transition region
    if (rho >= meta_numerics::CoulombTurningPoint(L, eta)) {
      // beyond the turning point, use Steed's method
      SolutionPair result = meta_numerics::Coulomb_Steed(L, eta, rho);
      return result.FirstSolutionValue();
    } else {
      // inside the turning point, integrate out from the series limit
      return meta_numerics::CoulombF_Integrate(L, eta, rho);
    }
  }

}

double meta_numerics::CoulombG (int L, double eta, double rho) {

  if (L < 0) throw std::runtime_error(std::string("L < 0 in CoulombG"));
  if (rho < 0) throw std::runtime_error(std::string("rho < 0 in CoulombG"));

  if ((rho < 4.0) && std::abs(rho * eta) < 8.0) {
    // for small enough rho, use the power series for L=0, then recurse upward to desired L
    double F, FP, G, GP;
    meta_numerics::Coulomb_Zero_Series(eta, rho, F, FP, G, GP);
    meta_numerics::Coulomb_Recurse_Upward(0, L, eta, rho, G, GP);
    return G;
  } else if (rho > 32.0 + (L * L + eta * eta) / 2.0) {
    // for large enough rho, use the asymptotic series
    double F, G;
    meta_numerics::Coulomb_Asymptotic(L, eta, rho, F, G);
    return G;
  } else {
    // transition region
    if (rho >= meta_numerics::CoulombTurningPoint(L, eta)) {
      // beyond the turning point, use Steed's method
      meta_numerics::SolutionPair result
        = meta_numerics::Coulomb_Steed(L, eta, rho);
      return result.SecondSolutionValue();
    } else {
      
      // we will start at L=0 (which has a smaller turning point radius) and recurse up to the desired L
      // this is okay because G increasees with increasing L

      double G, GP;

      if (rho < 2.0 * eta) {

        // if inside the turning point even for L=0, start at the turning point and integrate in
        // this is okay becaue G increases with decraseing rho

        // use Steed's method at the turning point
        // for large enough eta, we could use the turning point expansion at L=0, but it contributes
        // a lot of code for little overall performance increase so we have chosen not to
        meta_numerics::SolutionPair result;
        //if (eta > 12.0) {
        //  result = Coulomb_Zero_Turning_Expansion(eta);
        //} else {
          result = meta_numerics::Coulomb_Steed(0, eta, 2.0 * eta);
        //}

        G = result.SecondSolutionValue();
        GP = result.SecondSolutionDerivative();

        meta_numerics::BulrischStoerStoermerStepper s
          = meta_numerics::BulrischStoerStoermerStepper();
        s.RightHandSide = [eta] (double x, double U) -> double {
          return ((2.0 * eta / x - 1.0) * U);
        };
        s.X = 2.0 * eta;
        s.Y = G;
        s.YPrime = GP;
        s.DeltaX = 0.25;
        s.set_accuracy(2.5E-13);
        s.Integrate(rho);

        G = s.Y;
        GP = s.YPrime;

      } else {

        // if beyond the turning point for L=0, just use Steeds method

        meta_numerics::SolutionPair result
          = meta_numerics::Coulomb_Steed(0, eta, rho);
        G = result.SecondSolutionValue();
        GP = result.SecondSolutionDerivative();
      }

      // recurse up to desired L
      meta_numerics::Coulomb_Recurse_Upward(0, L, eta, rho, G, GP);

      return G;
    }
  }

}

double meta_numerics::LegendrePe(int l, int m, double x) {
  if (l < 0) throw std::runtime_error(std::string("Cannot compute")
    + " the associated Legendre polynomial Pe{l = " + std::to_string(l)
    + ", m = " + std::to_string(m) + "} because l must be nonnegative.");
  if (m > l || m < 0) throw std::runtime_error(std::string("Cannot compute")
    + " the associated Legendre polynomial Pe{l = " + std::to_string(l)
    + ", m = " + std::to_string(m) + "} because l and m must"
    + " satisfy the relation 0 <= m <= l");
  if (std::abs(x) > 1.0) throw std::runtime_error(std::string("Cannot compute")
    + " the associated Legendre polynomial Pe{l = " + std::to_string(l)
    + ", m = " + std::to_string(m) + "}(x = " + std::to_string(x)
    + "because x must satisfy the relation abs(x) <= 1.");
  
  double xx = (1.0 + x) * (1.0 - x);
  // determine P{m,m}
  double P0 = 1.0;
  for (int k = 1; k <= m; k++) {
    P0 *= (1.0 - 1.0 / (2 * k)) * xx;
  }
  P0 = std::sqrt(P0);
  if (m % 2 != 0) P0 = -P0;
  if (l == m) return (P0);
  // determine P{m+1,m}
  double s0 = std::sqrt(2*m + 1);
  double P1 = x * s0 * P0;
  // iterate up to P{l,m}
  for (int k = m + 2; k <= l; k++) {
    double s2 = std::sqrt((k - m) * (k + m));
    double P2 = (x * (2 * k - 1) * P1 - s0 * P0) / s2;
    // prepare for next iteration
    s0 = s2;
    P0 = P1;
    P1 = P2;
  }
  return P1;
}

std::complex<double> meta_numerics::SphericalHarmonic (int l, int m, double theta, double phi) {
  if (l < 0) throw std::runtime_error(std::string("Cannot compute")
    + " the spherical harmonic with l = " + std::to_string(l) + ", m = "
    + std::to_string(m) + ". The order l must be nonnegative.");
  if (std::abs(m) > l) throw std::runtime_error(std::string("Cannot compute")
    + " the spherical harmonic with l = " + std::to_string(l) + ", m = "
    + std::to_string(m) + ". The sub-order m must satisfy -l <= m <= l.");

  if (m < 0) {
    std::complex<double> y = SphericalHarmonic(l, -m, theta, phi);
    if ((m % 2) != 0) y = -y;
    return std::conj(y);
  }

  double LP = std::sqrt((2*l + 1) / (4.0 * marley_utils::pi))
    * LegendrePe(l, m, std::cos(theta));
  double mp = m * phi;
  return std::complex<double>(LP * std::cos(mp), LP * std::sin(mp));
}

double meta_numerics::SphericalBesselJ(int n, double x) {

  if (x < 0.0) {
    if (n % 2 == 0) return SphericalBesselJ(n, -x);
    else return -SphericalBesselJ(n, -x);
  }

  if (n < 0) {
    if ((n % 2) == 0) return SphericalBesselY(-n - 1, x);
    else return(-SphericalBesselY(-n - 1, x));
  } 
  else if (n == 0) return SphericalBesselJ_Zero(x);
  else if (n == 1) return SphericalBesselJ_One(x);
  else {
    // if close enough to the origin, use the power series
    if (x < 2.0 + 2.0 * std::sqrt(n)) return SphericalBesselJ_Series(n, x);
    // if far enough from the origin, use the asymptotic expansion
    else if (x > (32.0 + n * n / 2.0))
      return std::sqrt(marley_utils::half_pi / x)
      * Bessel_Asymptotic(n + 0.5, x).FirstSolutionValue();
    // in the transition region, use Miller's algorithm
    else return SphericalBesselJ_Miller(n, x);
  }
}

double meta_numerics::SphericalBesselY(int n, double x) {

  if (x < 0.0) {
    if (n % 2 == 0) return -SphericalBesselY(n, -x);
    else return SphericalBesselY(n, -x);
  }

  if (n < 0) {
    if ((n % 2) == 0) return -SphericalBesselJ(-n - 1, x);
    else return SphericalBesselJ(-n - 1, x);
  }
  else if (n == 0) return SphericalBesselY_Zero(x);
  else if (n == 1) return SphericalBesselY_One(x);
  else {
    if (x < (2.0 + std::sqrt(n))) return SphericalBesselY_Series(n, x);
    // if x is large enough, use asymptotic expansion
    else if (x > (30.0 + 0.5 * n * n))
      return std::sqrt(marley_utils::half_pi / x)
      * Bessel_Asymptotic(n + 0.5, x).SecondSolutionValue();
    else {
      // move up using the recursion relation
      double ym1 = SphericalBesselY_Zero(x);
      double y = SphericalBesselY_One(x);
      for (int k = 1; k < n; k++) {
        double yp1 = (2 * k + 1) / x * y - ym1;
        ym1 = y;
        y = yp1;
      }
      return y;
    }
  }
}

double meta_numerics::SphericalBesselJ_Zero(double x) {
  // It's fine to compute sin(x) / x explicitly, even for very very tiny x,
  // as long as x is not actually zero.
  if (x == 0.0) return 1.0;
  else return std::sin(x) / x;
}

double meta_numerics::SphericalBesselJ_SeriesOne(double x) {
  double xx = x * x / 2.0;
  double dj = x / 3.0;
  double j = dj;
  for (int i = 1; i < SeriesMax; i++) {
    double j_old = j;
    dj = -dj * xx / i / (2 * i + 3);
    j = j_old + dj;
    if (j == j_old) return j;
  }
  throw std::runtime_error(std::string("meta_numerics::")
    + "SphericalBesselJ_SeriesOne failed to converge.");
}

double meta_numerics::SphericalBesselJ_One(double x) {
  if (std::abs(x) < 0.5) return SphericalBesselJ_SeriesOne(x);
  else if (std::abs(x) > 100.0) return std::sqrt(marley_utils::half_pi / x)
    * Bessel_Asymptotic(1.5, x).FirstSolutionValue();
  else return (std::sin(x) / x - std::cos(x)) / x;
}

double meta_numerics::DoubleFactorial(int n) {
  if (n < 0) throw std::runtime_error(std::string("Cannot compute")
    + " n!! for n = " + std::to_string(n));
  else if (n < 32) return static_cast<double>(DoubleFactorial_Multiply(n));
  else return std::round(std::exp(LogDoubleFactorial_Gamma(n)));
}

double meta_numerics::LogDoubleFactorial(int n) {
  if (n < 0) throw std::runtime_error(std::string("Cannot compute")
    + " log(n!!) for n = " + std::to_string(n));
  else if (n < 32)
    return std::log(static_cast<double>(DoubleFactorial_Multiply(n)));
  else return LogDoubleFactorial_Gamma(n);
}

long meta_numerics::DoubleFactorial_Multiply(int n) {
  long f = 1;
  for (int k = n; k > 1; k = k - 2) f *= k;
  return f;
}

double meta_numerics::LogDoubleFactorial_Gamma(int n) {
  if (n % 2 == 0) {
    // m = n/2, n!! = 2^m Gamma(m+1)
    int m = n / 2;
    return m * marley_utils::log_2 + LogGamma(m + 1.0);
  }
  else {
    // m = (n+1)/2, n!! = 2^m Gamma(m+1/2) / Sqrt(PI)
    int m = (n + 1) / 2;
    return m * marley_utils::log_2 + LogGamma(m + 0.5)
      - std::log(marley_utils::pi) / 2.0;
  }
}

double meta_numerics::SphericalBesselJ_Series(int n, double x) {
  double xx = x * x / 2.0;
  double df = std::exp(n * std::log(x) - LogDoubleFactorial(2 * n + 1));
  double f = df;
  for (int i = 1; i < SeriesMax; i++) {
    double f_old = f;
    df = -df * xx / i / (2 * (n + i) + 1);
    f += df;
    if (f == f_old) return f;
  }
  throw std::runtime_error(std::string("meta_numerics::")
    + "SphericalBesselJ_Series failed to converge.");
}

double meta_numerics::SphericalBesselY_Series(int n, double x) {
  double xx = x * x / 2.0;
  double df = - DoubleFactorial(2 * n - 1)
    / std::pow(x, n + 1);
  double f = df;
  for (int k = 1; k < SeriesMax; k++) {
    double f_old = f;
    df = -df * xx / k / (2 * (k - n) - 1);
    f += df;
    if (f == f_old) return f;
  }
  throw std::runtime_error(std::string("meta_numerics::")
    + "SphericalBesselY_Series failed to converge.");
}

double meta_numerics::SphericalBesselY_Zero(double x) {
  if (x == 0.0) return marley_utils::minus_infinity;
  else return -std::cos(x) / x;
}

double meta_numerics::SphericalBesselY_SeriesOne(double x) {
  if (x == 0) return marley_utils::minus_infinity;
  double xx = x * x / 2.0;
  double dy = - 1.0 / (x * x);
  double y = dy;
  for (int i = 1; i < SeriesMax; i++) {
    double y_old = y;
    dy = - dy * xx / i / (2 * i - 3);
    y = y_old + dy;
    if (y == y_old) return y;
  }
  throw std::runtime_error(std::string("meta_numerics::")
    + "SphericalBesselY_SeriesOne failed to converge.");
}

double meta_numerics::SphericalBesselY_One(double x) {
  if (std::abs(x) < 1.0) return SphericalBesselY_SeriesOne(x);
  else if (std::abs(x) > 100.0) return std::sqrt(marley_utils::half_pi / x)
    * Bessel_Asymptotic(1.5, x).SecondSolutionValue();
  else return -(std::cos(x)/x + std::sin(x)) / x;
}

// Miller's method assumes a value at some high N and recurs downward
// the result is then normalized using a sum relation or a known value
double meta_numerics::SphericalBesselJ_Miller(int n, double x) {

  // pick starting value for the downward recursion that
  // takes us well into the x < nu regime, where J_nu decreases with nu
  int kmax = n;
  if (x > n) kmax = static_cast<int>(std::ceil(x));
  kmax += 50; // since J_(nu+1)/J_(nu) ~ 1/2 for x~v, taking N steps supresses by 2^(N) = 10^(16) at N ~ 50

  double jp1 = 0.0;
  double j = 1.0 / (static_cast<double>(kmax)); // look for a better guess

  // recur downward to order zero
  // the recurrence j_{k-1} = (2k+1)/x * j_k - j_{k+1} is stable in this direction
  for (int k = kmax; k > n; k--) {
    double jm1 = (2 * k + 1) / x * j - jp1;
    jp1 = j;
    j = jm1;
  }
  double jn = j;
  for (int k = n; k > 0; k--) {
    double jm1 = (2 * k + 1) / x * j - jp1;
    jp1 = j;
    j = jm1;
  }

  // compute the value we should have got and use it to normalize our result
  double j0 = SphericalBesselJ_Zero(x);
  return (j0 / j) * jn;
}

meta_numerics::SolutionPair meta_numerics::Bessel_Asymptotic(double nu, double x) {

  // pre-compute factors of nu and x as they appear in the series
  double mu = 4.0 * nu * nu;
  double xx = 8.0 * x;

  // initialize P and Q
  double P = 1.0; double R = 1.0;
  double Q = 0.0; double S = 0.0;

  // k is the current term number, k2 is (2k - 1), and t is the value of the current term
  int k = 0;
  int k2 = -1;
  double t = 1.0;

  while (true) {

    double Q_old = Q; double P_old = P;
    double R_old = R; double S_old = S;

    k++; k2 += 2;
    t /= k * xx;
    S += (mu + k2 * (k2 + 2)) * t;
    t *= (mu - k2 * k2);
    Q += t;

    k++; k2 += 2;
    t /= -k * xx;
    R += (mu + k2 * (k2 + 2)) * t;
    t *= (mu - k2 * k2);
    P += t;

    if ((P == P_old) && (Q == Q_old) && (R == R_old) && (S == S_old))
      break;

    if (k > SeriesMax) throw std::runtime_error(std::string("meta_numerics::")
      + "Bessel_Asymptotic failed to converge.");
  }

  // We attempted to move to a single trig evaluation so as to avoid errors when the two terms nearly cancel,
  // but this seemed to cause problems, perhaps because the arctan angle cannot be determined with sufficient
  // resolution. Investigate further.
  /*
  double M = N * MoreMath.Hypot(Q, P);
  double phi = Math.Atan2(Q, P);
  SolutionPair result2 = new SolutionPair(
    M * Cos(x + phi, -(nu + 0.5) / 4.0),
    -N * (R * s + S * c),
    M * Sin(x + phi, -(nu + 0.5) / 4.0),
    N * (R * c - S * s)
  );
   */

  // This technique was used in the original Meta Numerics library. To avoid having to use
  // the complicated techniques in the RangeReduction class, skip the reductions. Add them
  // back in if you experience problems.
  //
  //// Compute sin and cosine of x - (\nu + 1/2)(\pi / 2)
  //// Then we compute sine and cosine of (x1 - u1) with shift appropriate to (x0 - u0).
  //
  //// For maximum accuracy, we first reduce x = (x_0 + x_1) (\pi / 2), where x_0 is an integer and -0.5 < x1 < 0.5
  //long x0; double x1;
  //RangeReduction.ReduceByPiHalves(x, x0, x1);
  //
  //// Then we reduce (\nu + 1/2) = u_0 + u_1 where u_0 is an integer and -0.5 < u1 < 0.5
  //double u = nu + 0.5;
  //double ur = std::round(u);
  //long u0 = static_cast<long>(ur);
  //double u1 = u - ur;
  //
  //// Finally, we compute sine and cosine, having reduced the evaluation interval to -0.5 < \theta < 0.5
  //double s1 = RangeReduction.Sin(x0 - u0, x1 - u1);
  //double c1 = RangeReduction.Cos(x0 - u0, x1 - u1);
  // Assemble the solution
  //double N = std::sqrt(2.0 / marley_utils::pi / x);
  //meta_numerics::SolutionPair result = meta_numerics::SolutionPair(
  //  N * (c1 * P - s1 * Q), -N * (R * s1 + S * c1),
  //  N * (s1 * P + c1 * Q), N * (R * c1 - S * s1)
  //);

  // Assemble the solution
  double phi = x - (nu + 0.5)*marley_utils::half_pi;
  double s1 = std::sin(phi);
  double c1 = std::cos(phi);
  double N = std::sqrt(2.0 / marley_utils::pi / x);
  meta_numerics::SolutionPair result = meta_numerics::SolutionPair(
    N * (c1 * P - s1 * Q), -N * (R * s1 + S * c1),
    N * (s1 * P + c1 * Q), N * (R * c1 - S * s1)
  );

  return result;
}

/// <summary>
/// Computes the natural logrithm of the Gamma function.
/// </summary>
/// <param name="x">The argument, which must be positive.</param>
/// <returns>The log Gamma function ln(&#x393;(x)).</returns>
/// <remarks>
/// <para>Because &#x393;(x) grows rapidly for increasing positive x, it is often necessary to
/// work with its logarithm in order to avoid overflow. This function returns accurate
/// values of ln(&#x393;(x)) even for values of x which would cause &#x393;(x) to overflow.</para>
/// </remarks>
/// <exception cref="ArgumentOutOfRangeException"><paramref name="x"/> is negative or zero.</exception>
/// <seealso cref="Gamma(double)" />
double meta_numerics::LogGamma(double x) {
  if (x <= 0.0) throw std::runtime_error(std::string("Cannot compute")
    + " LogGamma(x) for x = " + std::to_string(x));
  // For small arguments, use the Lanczos approximation.
  else if (x < 16.0) return meta_numerics::Lanczos::LogGamma(x);
  // For large arguments, the asymptotic series is even faster than the Lanczos approximation.
  else return meta_numerics::LogGamma_Stirling(x);
}

double meta_numerics::LogGamma_Stirling(double x) {
  // re-write to use (x-0.5) form to eliminate one one by storing log(2\pi)?
  return (x * std::log(x) - x
    - std::log(x / (2.0 * marley_utils::pi)) / 2.0 + Sum_Stirling(x));
}

double meta_numerics::Sum_Stirling(double x) {

  double xx = x * x; // x^2 
  double xk = x; // tracks x^{2k - 1}
  double f = meta_numerics::Bernoulli[1] / 2.0 / xk; // k = 1 term
  for (size_t k = 2; k < meta_numerics::Bernoulli.size(); k++) {
    double f_old = f;
    xk *= xx;
    f += meta_numerics::Bernoulli[k] / (2 * k) / (2 * k - 1) / xk;
    if (f == f_old) return f;
  }

  throw std::runtime_error(std::string("meta_numerics::")
    + "Sum_Stirling failed to converge.");
}
// -- End Ms-PL licensed code

//Microsoft Public License (Ms-PL)
//
//This license governs use of the accompanying software. If you use the software, you accept
//this license. If you do not accept the license, do not use the software.
//
//1. Definitions
//
//The terms "reproduce," "reproduction," "derivative works," and "distribution" have the
//same meaning here as under U.S. copyright law.
//
//A "contribution" is the original software, or any additions or changes to the software.
//
//A "contributor" is any person that distributes its contribution under this license.
//
//"Licensed patents" are a contributor's patent claims that read directly on its
//contribution.
//
//2. Grant of Rights
//
//(A) Copyright Grant- Subject to the terms of this license, including the license
//conditions and limitations in section 3, each contributor grants you a non-exclusive,
//worldwide, royalty-free copyright license to reproduce its contribution, prepare
//derivative works of its contribution, and distribute its contribution or any derivative
//works that you create.
//
//(B) Patent Grant- Subject to the terms of this license, including the license conditions
//and limitations in section 3, each contributor grants you a non-exclusive, worldwide,
//royalty-free license under its licensed patents to make, have made, use, sell, offer for
//sale, import, and/or otherwise dispose of its contribution in the software or derivative
//works of the contribution in the software.
//
//3. Conditions and Limitations
//
//(A) No Trademark License- This license does not grant you rights to use any contributors'
//name, logo, or trademarks.
//
//(B) If you bring a patent claim against any contributor over patents that you claim are
//infringed by the software, your patent license from such contributor to the software ends
//automatically.
//
//(C) If you distribute any portion of the software, you must retain all copyright, patent,
//trademark, and attribution notices that are present in the software.
//
//(D) If you distribute any portion of the software in source code form, you may do so only
//under this license by including a complete copy of this license with your distribution. If
//you distribute any portion of the software in compiled or object code form, you may only
//do so under a license that complies with this license.
//
//(E) The software is licensed "as-is." You bear the risk of using it. The contributors give
//no express warranties, guarantees or conditions. You may have additional consumer rights
//under your local laws which this license cannot change. To the extent permitted under your
//local laws, the contributors exclude the implied warranties of merchantability, fitness
//for a particular purpose and non-infringement.
