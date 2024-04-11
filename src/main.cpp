#include <cylindrical_cell.hpp>
#include <cylindrical_flux_solver.hpp>
#include <utils/constants.hpp>
#include <utils/logging.hpp>

#include <moc/pin_cell.hpp>
#include <moc/cartesian_2d.hpp>
#include <moc/moc_driver.hpp>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using namespace scarabee;

void test() {
  xt::xtensor<double, 1> Et = {1.77949E-01, 3.29805E-01, 4.80388E-01,
                               5.54367E-01, 3.11801E-01, 3.95168E-01,
                               5.64406E-01};
  xt::xtensor<double, 1> Ea = {8.02480E-03, 3.71740E-03, 2.67690E-02,
                               9.62360E-02, 3.00200E-02, 1.11260E-01,
                               2.82780E-01};
  xt::xtensor<double, 1> Ef = {7.21206E-03, 8.19301E-04, 6.45320E-03,
                               1.85648E-02, 1.78084E-02, 8.30348E-02,
                               2.16004E-01};
  xt::xtensor<double, 1> nu = {2.78145, 2.47443, 2.43383, 2.43380,
                               2.43380, 2.43380, 2.43380};
  xt::xtensor<double, 1> chi = {
      5.87910E-01, 4.11760E-01, 3.39060E-04, 1.17610E-07, 0., 0., 0.};
  xt::xtensor<double, 2> Es = {
      {1.27537E-01, 4.23780E-02, 9.43740E-06, 5.51630E-09, 0.00000E+00,
       0.00000E+00, 0.00000E+00},
      {0.00000E+00, 3.24456E-01, 1.63140E-03, 3.14270E-09, 0.00000E+00,
       0.00000E+00, 0.00000E+00},
      {0.00000E+00, 0.00000E+00, 4.50940E-01, 2.67920E-03, 0.00000E+00,
       0.00000E+00, 0.00000E+00},
      {0.00000E+00, 0.00000E+00, 0.00000E+00, 4.52565E-01, 5.56640E-03,
       0.00000E+00, 0.00000E+00},
      {0.00000E+00, 0.00000E+00, 0.00000E+00, 1.25250E-04, 2.71401E-01,
       1.02550E-02, 1.00210E-08},
      {0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.29680E-03,
       2.65802E-01, 1.68090E-02},
      {0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
       8.54580E-03, 2.73080E-01}};
  std::shared_ptr<TransportXS> UO2 =
      std::make_shared<TransportXS>(Et, Ea, Es, nu * Ef, chi);

  Et = {1.59206E-01, 4.12970E-01, 5.90310E-01, 5.84350E-01,
        7.18000E-01, 1.25445E+00, 2.65038E+00};
  Ea = {6.01050E-04, 1.57930E-05, 3.37160E-04, 1.94060E-03,
        5.74160E-03, 1.50010E-02, 3.72390E-02};
  Es = {{4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08,
         0.00000E+00, 0.00000E+00},
        {0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05,
         7.44860E-06, 1.04550E-06},
        {0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02,
         2.64430E-03, 5.03440E-04},
        {0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01,
         6.37320E-02, 1.21390E-02},
        {0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01,
         5.11820E-01, 6.12290E-02},
        {0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 2.21570E-03,
         6.99913E-01, 5.37320E-01},
        {0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
         1.32440E-01, 2.48070E+00}};
  std::shared_ptr<TransportXS> H2O = std::make_shared<TransportXS>(Et, Ea, Es);

  const double Rfuel = 0.54;
  const double Rwtr = 1.26 / std::sqrt(PI);

  // Vector of all the radii
  std::vector<double> radii;
  std::vector<std::shared_ptr<TransportXS>> mats;

  // Break fuel into 5 equal radii
  const int Nf = 5;
  const double dRfuel = Rfuel / static_cast<double>(Nf);
  double r_outer = 0.;
  for (int i = 0; i < Nf; i++) {
    r_outer += dRfuel;
    radii.push_back(r_outer);
    mats.push_back(UO2);
  }

  // Break Water into 3 equal radii
  const int Nwtr = 3;
  const double dRwtr = (Rwtr - Rfuel) / static_cast<double>(Nwtr);
  for (int i = 0; i < Nwtr; i++) {
    r_outer += dRwtr;
    radii.push_back(r_outer);
    mats.push_back(H2O);
  }

  std::shared_ptr<CylindricalCell> cell =
      std::make_shared<CylindricalCell>(radii, mats);
  cell->solve();

  CylindricalFluxSolver cell_flux(cell);
  cell_flux.set_albedo(1.);
  cell_flux.solve();
  std::cout << "\n";

  //==============================================
  // MOC
  radii = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.54, 0.57, 0.6};
  mats = {UO2, UO2, UO2, UO2, UO2, UO2, UO2, UO2, UO2, UO2, H2O, H2O, H2O};
  std::shared_ptr<PinCell> pincell1 = std::make_shared<PinCell>(radii, mats, 1.26, 1.26);
  
  std::vector<double> dx {1.26};
  std::vector<double> dy {1.26};
  std::shared_ptr<Cartesian2D> c2d = std::make_shared<Cartesian2D>(dx, dy);
  c2d->set_tiles({pincell1});

  YamamotoTabuchi<6> pq;
  PolarQuadrature pquad(pq);
  MOCDriver moc(c2d);
  moc.generate_tracks(128, 0.01, pquad);
  moc.set_keff_tolerance(1.E-5);
  moc.set_flux_tolerance(1.E-5);
  moc.solve_keff();
  std::cout << "\n";

  //=================================================
  Et = {4.52648699E-01};
  Ea = {6.9389522E-02};
  Ef = {3.97630632E-2};
  nu = {2.5};
  chi = {1.};
  Es = {{3.83259177E-01}};
  std::shared_ptr<TransportXS> fuel = std::make_shared<TransportXS>(Et, Ea, Es, nu*Ef, chi);

  Et = {8.41545641E-01};
  Ea = {3.751099E-3};
  Es = {{8.37794542E-01}};
  std::shared_ptr<TransportXS> water = std::make_shared<TransportXS>(Et, Ea, Es);

  radii = {0.1,  0.2,  0.3,  0.4,  0.5,   0.6,   std::sqrt(1.27*1.27 / PI)};
  mats  = {fuel, fuel, fuel, fuel, water, water, water};
  std::shared_ptr<CylindricalCell> cell2 =
  std::make_shared<CylindricalCell>(radii, mats); cell2->solve();
  CylindricalFluxSolver cell_flux2(cell2);
  cell_flux2.set_albedo(1.);
  cell_flux2.solve();
}

int main() {
  // set_logging_level(LogLevel::debug);
  set_logging_level(LogLevel::info);
  test();

  return 0;
}
