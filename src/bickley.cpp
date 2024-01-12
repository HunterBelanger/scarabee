#include <utils/bickley.hpp>
#include <utils/chebyshev.hpp>
#include <utils/constants.hpp>
#include <utils/gauss_kronrod.hpp>
#include <utils/scarabee_exception.hpp>

#include <array>
#include <cmath>
#include <sstream>

double Ki3(double x) {
  constexpr double a = 0.;
  constexpr double b = 10.;
  constexpr std::array<double, 60> coeffs{
      0x1.0d1d9af3e0125p-2,   -0x1.e8c3c739419b2p-3,  0x1.70ddc21b44afdp-3,
      -0x1.d75a34d619107p-4,  0x1.051673ea12c52p-4,   -0x1.017e8bdd41456p-5,
      0x1.d18dd69d23999p-7,   -0x1.8e9ee4992e7cp-8,   0x1.4f3fc1a52ebf7p-9,
      -0x1.1f6832ae58e22p-10, 0x1.030dcb842f6abp-11,  -0x1.f3ab4503c5199p-13,
      0x1.026fded3404cdp-13,  -0x1.1ca62c9820111p-14, 0x1.4a662cb73a222p-15,
      -0x1.9069780ccdddep-16, 0x1.f704bf6381111p-17,  -0x1.45c1cd8b25ddep-17,
      0x1.b13874ab24444p-18,  -0x1.26d32af122222p-18, 0x1.999ec7ead5555p-19,
      -0x1.21d733c98p-19,     0x1.a1072ea3b7777p-20,  -0x1.308dc0521999ap-20,
      0x1.c2f47f84d5555p-21,  -0x1.521224ba08888p-21, 0x1.005fd4e9ccccdp-21,
      -0x1.88fb8910ccccdp-22, 0x1.30282304aaaabp-22,  -0x1.db1d3d2eccccdp-23,
      0x1.763a93b733333p-23,  -0x1.291828ea66666p-23, 0x1.db31edacccccdp-24,
      -0x1.7ea6554622222p-24, 0x1.361d2ac755555p-24,  -0x1.f9b05f3b77777p-25,
      0x1.9ea164e155555p-25,  -0x1.55c3952088888p-25, 0x1.1b177072aaaabp-25,
      -0x1.d722635222222p-26, 0x1.89b4cf9111111p-26,  -0x1.4a4928faaaaabp-26,
      0x1.160fd58dddddep-26,  -0x1.d5aad59444444p-27, 0x1.8dbc314eeeeefp-27,
      -0x1.5194525eeeeefp-27, 0x1.1f02905444444p-27,  -0x1.e88af20222222p-28,
      0x1.9fe5498p-28,        -0x1.61c85b7p-28,       0x1.2c5466199999ap-28,
      -0x1.fbfea27bbbbbcp-29, 0x1.ab05e17bbbbbcp-29,  -0x1.639e42a666666p-29,
      0x1.23fd078dddddep-29,  -0x1.d5366b5333333p-30, 0x1.6c4da9b555555p-30,
      -0x1.0aedcd8333333p-30, 0x1.5e03d44aaaaabp-31,  -0x1.5a86574444444p-32};

  // If we are at x > 10, Ki3(x) is very small, so we just call it zero here.
  // This should be okay, as Stamm'ler and Abbate set it to 0 for x >= 9.;
  if (x > b) return 0.;

  double val;
  try {
    val = chebyshev_eval(x, a, b, coeffs);
  } catch (ScarabeeException& err) {
    std::stringstream strm;
    strm << "Ki3 is only defined on the interval [" << a << ", " << b
         << "]. Was provided with x = " << x << ".";
    err.add_to_exception(strm.str());
  }

  return val;
}

double Ki3_quad(double x) {
  auto integrand = [x](double theta) {
    const double cos_theta = std::cos(theta);
    const double cos_theta_sqrd = cos_theta * cos_theta;

    return cos_theta_sqrd * std::exp(-x / cos_theta);
  };

  GaussKronrodQuadrature<61> gk;
  auto integral = gk.integrate(integrand, 0., PI_2, 1.E-9, 10000);

  return integral.first;
}
