// This implementation of the complex gamma function is based on the
// Lanczos approximation and its Python implementation given
// on Wikipedia (https://en.wikipedia.org/wiki/Lanczos_approximation)
// The C++ version given here is taken almost verbatim from
// http://bytes.com/topic/c/answers/576697-c-routine-complex-gamma-function
#include <complex>
#include <cmath>
#include <iomanip>
#include <iostream>

static const int g=7;
static const double pi =
3.1415926535897932384626433832795028841972;
static const double p[g+2] = {0.99999999999980993, 676.5203681218851,
-1259.1392167224028, 771.32342877765313, -176.61502916214059,
12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
1.5056327351493116e-7};

std::complex<double> gamma(std::complex<double> z)
{

  if (std::real(z) < 0.5) {
    return pi / (std::sin(pi*z)*gamma(1.0-z));
  }

  z -= 1.0;

  std::complex<double> x = p[0];

  for (int i = 1; i < g + 2; i++) {
    x += p[i]/(z+std::complex<double>(i,0));
  }

  std::complex<double> t = z + (g + 0.5);

  return std::sqrt(2*pi) * std::pow(t,z+0.5) * std::exp(-t) * x;
}

int main() {

  double re_part, im_part;

  std::cout << std::endl;

  std::cout << std::setprecision(15);

  for (int dummy = 1; dummy < 5; dummy++) {
    std::cout << "Enter real part: ";
    std::cin >> re_part;
    std::cout << "Enter imaginary part: ";
    std::cin >> im_part;
    std::cout << "The gamma function is ";
    std::cout << gamma(std::complex<double>(re_part, im_part));
    std::cout << std::endl << std::endl;
  }

}
