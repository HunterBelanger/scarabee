#include <utils/criticality_spectrum.hpp>

#include <xtensor/xbuilder.hpp>
#include <Eigen/Dense>

namespace scarabee {

void fill_A(Eigen::MatrixXd& A, std::shared_ptr<CrossSection> xs,
            const Eigen::MatrixXd& D, double B2) {
  const std::size_t NG = xs->ngroups();

  A.fill(0.);

  for (std::size_t g = 0; g < NG; g++) {
    for (std::size_t gg = 0; gg < NG; gg++) {
      A(g, gg) = B2 * D(g, gg) - xs->Es(gg, g);
    }
    A(g, g) += xs->Et(g);
  }
}

void fill_Dinvs(Eigen::MatrixXd& Dinvs, std::shared_ptr<CrossSection> xs,
                const xt::xtensor<double, 1>& a) {
  const std::size_t NG = xs->ngroups();

  Dinvs.fill(0.);

  for (std::size_t g = 0; g < NG; g++) {
    for (std::size_t gg = 0; gg < NG; gg++) {
      Dinvs(g, gg) = -xs->Es1(gg, g);
    }
    Dinvs(g, g) += a(g) * xs->Et(g);
  }

  Dinvs *= 3.;
}

xt::xtensor<double, 2> P1_spectrum(std::shared_ptr<CrossSection> xs) {
  const std::size_t NG = xs->ngroups();

  // For P1 approximation, alphas are all 1
  const xt::xtensor<double, 1> a = xt::ones<double>({NG});

  // First, we create and fill the Dinvs matrix for the current
  Eigen::MatrixXd Dinvs(NG, NG);
  fill_Dinvs(Dinvs, xs, xt::ones<double>({NG}));

  // Get the D matrix
  const Eigen::MatrixXd D = Dinvs.inverse();

  // Make the chi and vEf vectors
  Eigen::VectorXd chi(NG), vEf(NG);
  for (std::size_t g = 0; g < NG; g++) {
    chi(g) = xs->chi(g);
    vEf(g) = xs->vEf(g);
  }

  // Make vectors for flux and current
  Eigen::VectorXd flx(NG);
  Eigen::VectorXd cur(NG);

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

  // The output info is in the format (2,NG) where the first line has
  // the flux spectrum, and the second has the diffusion coefficients.
  xt::xtensor<double, 2> out;
  out.resize({2,NG});
  for (std::size_t g = 0; g < NG; g++) {
    out(0,g) = flx(g);
    out(1,g) = cur(g) / (B*flx(g));
  }

  return out;
}

}  // namespace scarabee
