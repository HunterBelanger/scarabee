#include <utils/criticality_spectrum.hpp>

#include <xtensor/xbuilder.hpp>
#include <Eigen/Dense>

namespace scarabee {

void fill_A(Eigen::MatrixXd& A, std::shared_ptr<TransportXS> xs, const Eigen::MatrixXd& D, double B2) {
  const std::size_t NG = xs->ngroups();

  A.fill(0.);
  
  for (std::size_t g = 0; g < NG; g++) {
    for (std::size_t gg = 0; gg < NG; gg++) {
      A(g,gg) = B2*D(g,gg) - (xs->Es(gg, g) + xs->Es1(g, gg));
    }
    A(g,g) += xs->Et(g) + xs->Es1(g, g);
  }
}

void fill_Dinvs(Eigen::MatrixXd& Dinvs, std::shared_ptr<TransportXS> xs, const xt::xtensor<double, 1>& a) {
  const std::size_t NG = xs->ngroups();
  
  Dinvs.fill(0.);

  for (std::size_t g = 0; g < NG; g++) {
    for (std::size_t gg = 0; gg < NG; gg++) {
      Dinvs(g,gg) = -xs->Es1(gg, g);
    }
    Dinvs(g,g) += a(g) * (xs->Et(g) + xs->Es1(g, g));
  }

  Dinvs *= 3.;
}

xt::xtensor<double, 2> P1_spectrum(std::shared_ptr<TransportXS> xs) {
  const std::size_t NG = xs->ngroups();

  // First, we create and fill the Dinvs matrix for the current
  Eigen::MatrixXd Dinvs(NG, NG);
  fill_Dinvs(Dinvs, xs, xt::ones<double>({NG}));
  
  // Get the D matrix
  Eigen::MatrixXd D = Dinvs.inverse();

  // Make the chi and vEf vectors
  Eigen::MatrixXd chi(NG), vEf(NG);
  for (std::size_t g = 0; g < NG; g++) {
    chi(g) = xs->chi(g);
    vEf(g) = xs->vEf(g);
  }

  // For P1 approximation, alphas are all 1
  xt::xtensor<double, 1> a = xt::ones<double>({NG});

  // Make vectors for flux and current
  Eigen::MatrixXd flx(NG);
  Eigen::MatrixXd cur(NG);

  // Create the A matrix
  Eigen::MatrixXd A(NG, NG);

  // Get A for k_inf
  double B2 = 0.;
  fill_A(A, xs, D, B2);
  auto A_solver = A.colPivHouseholderQr();
  flx = A_solver.solve(chi);
  const double k_inf = vEf.dot(flx);

  // Get A for small B2, from Stamm'ler and Abbate
  B2 = 0.001;
  if (k_inf < 1.) B2 = -B2;
  fill_A(A, xs, D, B2);
  A_solver = A.colPivHouseholderQr();
  flx = A_solver.solve(chi);
  const double k_1 = vEf.dot(flx);

  // Calculate the slope constant
  const double k_inf_M2 = B2 / ((1. / k_1) - (1. / k_inf));

  double k = k_1;

  while (std::abs(k - 1.) > 1.E-6) {
    B2 += k_inf_M2 * (1. - (1. / k));
    fill_A(A, xs, D, B2);
    A_solver = A.colPivHouseholderQr();
    flx = A_solver.solve(chi);
    k = vEf.dot(flx);
  }

  // We have converged on k = 1.
  const double sqrt_abs_B2 = std::sqrt(std::abs(B2));
  const double B = B2 > 0. ? sqrt_abs_B2 : -sqrt_abs_B2;

  // Get the current
  cur = B * D * flx;

} 

} // namespace scarabee
