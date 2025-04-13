#include <data/water.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <array>
#include <cmath>

namespace scarabee {

double water_density(double temperature, double pressure) {
  /*
  This method was directly taken from OpenMC's Python API, found in
  openmc/data/data.py. It uses a polynomial fit of the data in the
  2012 IAPWS-IF97. It expects the temperature in units of Kelvin and the
  pressure in units of MPa. It returns the density of water in g/cm^3.
  */

  // Make sure the temperature and pressure are inside the min/max region 1
  // bounds.  (Relax the 273.15 bound to 273 in case a user wants 0 deg C data
  // but they only use 3 digits for their conversion to K.)
  if (pressure > 100.0) {
    const auto mssg =
        "Water density calculation not valid for pressures greater than 100 "
        "MPa.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (pressure < 0.) {
    const auto mssg =
        "Water density calculation must have a positive pressure.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (temperature < 273.) {
    const auto mssg =
        "Water density calculation not valid for temperatures less than 273 K.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  } else if (temperature > 623.15) {
    const auto mssg =
        "Water density calculation not valid for temperatures greater than "
        "623.15 K.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // IAPWS region 4 parameters
  constexpr std::array<double, 10> n4 = {
      0.11670521452767e4,  -0.72421316703206e6, -0.17073846940092e2,
      0.12020824702470e5,  -0.32325550322333e7, 0.14915108613530e2,
      -0.48232657361591e4, 0.40511340542057e6,  -0.23855557567849,
      0.65017534844798e3};

  // Compute the saturation temperature at the given pressure.
  const double beta = std::pow(pressure, 0.25);
  const double beta2 = beta * beta;
  const double E = beta2 + n4[2] * beta + n4[5];
  const double F = n4[0] * beta2 + n4[3] * beta + n4[6];
  const double G = n4[1] * beta2 + n4[4] * beta + n4[7];
  const double D = 2.0 * G / (-F - std::sqrt(F * F - 4 * E * G));
  const double T_sat =
      0.5 * (n4[9] + D -
             std::sqrt((n4[9] + D) * (n4[9] + D) - 4.0 * (n4[8] + n4[9] * D)));

  // Make sure we aren't above saturation.  (Relax this bound by .2 degrees
  // for deg C to K conversions.)
  if (temperature > T_sat + 0.2) {
    auto mssg =
        "Water density cannot be calculated for temperatures above the boiling "
        "point.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // IAPWS region 1 parameters
  constexpr double R_GAS_CONSTANT = 0.461526;  // kJ / kg / K
  constexpr double ref_p = 16.53;              // MPa
  constexpr double ref_T = 1386.;              // K
  constexpr std::array<double, 34> n1f = {
      0.14632971213167,      -0.84548187169114,     -0.37563603672040e1,
      0.33855169168385e1,    -0.95791963387872,     0.15772038513228,
      -0.16616417199501e-1,  0.81214629983568e-3,   0.28319080123804e-3,
      -0.60706301565874e-3,  -0.18990068218419e-1,  -0.32529748770505e-1,
      -0.21841717175414e-1,  -0.52838357969930e-4,  -0.47184321073267e-3,
      -0.30001780793026e-3,  0.47661393906987e-4,   -0.44141845330846e-5,
      -0.72694996297594e-15, -0.31679644845054e-4,  -0.28270797985312e-5,
      -0.85205128120103e-9,  -0.22425281908000e-5,  -0.65171222895601e-6,
      -0.14341729937924e-12, -0.40516996860117e-6,  -0.12734301741641e-8,
      -0.17424871230634e-9,  -0.68762131295531e-18, 0.14478307828521e-19,
      0.26335781662795e-22,  -0.11947622640071e-22, 0.18228094581404e-23,
      -0.93537087292458e-25};
  constexpr std::array<double, 34> I1f = {
      0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,  1,  1,  2,  2,  2,
      2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};
  constexpr std::array<double, 34> J1f = {
      -2, -1, 0,  1, 2, 3,  4,  5,  -9, -7,  -1, 0,   1,   3,   -3,  0,   1,
      3,  17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};

  // Nondimensionalize the pressure and temperature.
  const double pi = pressure / ref_p;
  const double tau = ref_T / temperature;

  // Compute the derivative of gamma (dimensionless Gibbs free energy) with
  // respect to pi.
  double gamma1_pi = 0.0;
  for (std::size_t i = 0; i < 34; i++) {
    gamma1_pi -= n1f[i] * I1f[i] * std::pow(7.1 - pi, I1f[i] - 1.) *
                 std::pow(tau - 1.222, J1f[i]);
  }

  // Compute the leading coefficient.  This sets the units at
  //   1 [MPa] * [kg K / kJ] * [1 / K]
  // = 1e6 [N / m^2] * 1e-3 [kg K / N / m] * [1 / K]
  // = 1e3 [kg / m^3]
  // = 1 [g / cm^3]
  double coeff = pressure / R_GAS_CONSTANT / temperature;

  // Compute and return the density.
  return coeff / pi / gamma1_pi;
}

std::shared_ptr<Material> borated_water(double boron_ppm, double temperature,
                                        double pressure,
                                        std::shared_ptr<NDLibrary> ndl) {
  /*
  This method was directly taken from OpenMC's Python API, found in
  openmc/model/funcs.py. It expects the temperature in units of Kelvin
  and the pressure in units of MPa. It returns the fully initialize Material
  at the desired temperature, pressure, and density.
  */

  if (boron_ppm < 0.) {
    auto mssg = "Boron concentration must be greater than zero.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  if (ndl == nullptr) {
    auto mssg = "Provided NDLibrary is nullptr.";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  // Compute the density of the solution.
  const double solution_density =
      water_density(temperature, pressure) / (1. - boron_ppm * 1.e-6);

  // Simple enough to hard-code the molar masses of H2O and B for this,
  // instead of calculating them on-the-fly like OpenMC does.
  constexpr double M_H2O = 18.01526800811539;
  constexpr double M_B = 10.8118249681472;

  // Compute the number fractions of each element.
  const double frac_H2O = (1. - boron_ppm * 1.e-6) / M_H2O;
  const double frac_B = boron_ppm * 1.e-6 / M_B;

  // Construct the material composition
  MaterialComposition water_comp(Fraction::Atoms);
  water_comp.add_nuclide("H1_H2O", 2. * frac_H2O * NATURAL_ABUNDANCES.at("H1"));
  water_comp.add_nuclide("H2", 2. * frac_H2O * NATURAL_ABUNDANCES.at("H2"));
  water_comp.add_element("O", frac_H2O);
  if (frac_B > 0.) water_comp.add_element("B", frac_B);

  return std::make_shared<Material>(water_comp, temperature, solution_density,
                                    DensityUnits::g_cm3, ndl);
}

}  // namespace scarabee