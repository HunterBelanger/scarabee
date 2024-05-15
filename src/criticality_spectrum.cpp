#include <utils/criticality_spectrum.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

#include <xtensor/xbuilder.hpp>
#include <Eigen/Dense>

#include <sstream>

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

void fill_alphas(xt::xtensor<double, 1>& a,
                 const std::shared_ptr<CrossSection>& xs, const double B2) {
  const std::size_t NG = xs->ngroups();

  for (std::size_t g = 0; g < NG; g++) {
    const double Et2 = xs->Et(g) * xs->Et(g);
    const double x2 = std::abs(B2 / Et2);
    const double x = std::sqrt(x2);

    // Et2 will always be > 0, so we only need to check the sign of B2
    if (x2 < 0.1) {
      const double y = B2 / Et2;
      const double x = (1. / 3.) - y * (1. / 5. - y / 7.);
      a(g) = (1. / 3.) * (1. / (x - y));
    } else if (B2 > 0.) {
      a(g) = (1. / 3.) * x2 * ((std::atan(x)) / (x - std::atan(x)));
    } else {
      const double xp1 = 1. + x;
      const double xm1 = 1. - x;
      const double arg = xp1 / xm1;
      a(g) = (1. / 3.) * x2 * ((std::log(arg)) / (std::log(arg) - 2. * x));
    }
  }
}

P1CriticalitySpectrum::P1CriticalitySpectrum(std::shared_ptr<CrossSection> xs) {
  if (xs->fissile() == false) {
    std::stringstream mssg;
    mssg << "Cannot compute P1 spectrum of homogenized material that is not "
            "fissile.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  const std::size_t NG = xs->ngroups();

  // For P1 approximation, alphas are all 1
  const xt::xtensor<double, 1> a = xt::ones<double>({NG});

  // Make vectors for flux and current
  Eigen::VectorXd flx(NG);
  Eigen::VectorXd cur(NG);
  Eigen::VectorXd chi(NG), vEf(NG);
  for (std::size_t g = 0; g < NG; g++) {
    chi(g) = xs->chi(g);
    vEf(g) = xs->vEf(g);
  }

  // Create the A matrix
  Eigen::MatrixXd A(NG, NG);
  Eigen::MatrixXd Dinvs(NG, NG);
  fill_Dinvs(Dinvs, xs, a);
  const Eigen::MatrixXd D = Dinvs.inverse();

  // Solve B2 = 0 for k_inf
  B2_ = 0.;
  fill_A(A, xs, D, B2_);
  auto A_solver = A.colPivHouseholderQr();
  flx = A_solver.solve(chi);
  k_inf_ = vEf.dot(flx);

  // Solve for small B2 for k1
  B2_ = 0.001;
  if (k_inf_ < 1.) B2_ = -B2_;
  fill_A(A, xs, D, B2_);
  A_solver = A.colPivHouseholderQr();
  flx = A_solver.solve(chi);
  const double k_1 = vEf.dot(flx);

  // Calculate the slope constant
  const double k_inf_M2 = B2_ / ((1. / k_1) - (1. / k_inf_));

  double k = k_1;
  while (std::abs(k - 1.) > 1.E-6) {
    B2_ += k_inf_M2 * (1. - (1. / k));
    fill_A(A, xs, D, B2_);
    A_solver = A.colPivHouseholderQr();
    flx = A_solver.solve(chi);
    k = vEf.dot(flx);
  }

  // We have converged on k = 1.
  const double sqrt_abs_B2 = std::sqrt(std::abs(B2_));
  const double B = B2_ > 0. ? sqrt_abs_B2 : -sqrt_abs_B2;

  // Get the current
  cur = B * D * flx;

  // The output info is in the format (2,NG) where the first line has
  // the flux spectrum, and the second has the diffusion coefficients.
  flux_.resize({NG});
  current_.resize({NG});
  diff_coeff_.resize({NG});
  for (std::size_t g = 0; g < NG; g++) {
    flux_(g) = flx(g);
    current_(g) = cur(g);
    diff_coeff_(g) = cur(g) / (B * flx(g));
  }
}

B1CriticalitySpectrum::B1CriticalitySpectrum(std::shared_ptr<CrossSection> xs) {
  if (xs->fissile() == false) {
    std::stringstream mssg;
    mssg << "Cannot compute B1 spectrum of homogenized material that is not "
            "fissile.";
    spdlog::error(mssg.str());
    throw ScarabeeException(mssg.str());
  }

  const std::size_t NG = xs->ngroups();

  // For P1 approximation, alphas are all 1
  xt::xtensor<double, 1> a = xt::zeros<double>({NG});

  // First, we create and fill the Dinvs matrix for the current
  Eigen::MatrixXd Dinvs(NG, NG);

  // Declare the D matrix
  Eigen::MatrixXd D;

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
  B2_ = 0.;
  fill_alphas(a, xs, B2_);
  fill_Dinvs(Dinvs, xs, a);
  D = Dinvs.inverse();
  fill_A(A, xs, D, B2_);
  auto A_solver = A.colPivHouseholderQr();
  flx = A_solver.solve(chi);
  k_inf_ = vEf.dot(flx);

  // Get A for small B2, from Stamm'ler and Abbate
  B2_ = 0.001;
  if (k_inf_ < 1.) B2_ = -B2_;
  fill_alphas(a, xs, B2_);
  fill_Dinvs(Dinvs, xs, a);
  D = Dinvs.inverse();
  fill_A(A, xs, D, B2_);
  A_solver = A.colPivHouseholderQr();
  flx = A_solver.solve(chi);
  const double k_1 = vEf.dot(flx);

  // Calculate the slope constant
  const double k_inf_M2 = B2_ / ((1. / k_1) - (1. / k_inf_));

  double k = k_1;
  while (std::abs(k - 1.) > 1.E-6) {
    B2_ += k_inf_M2 * (1. - (1. / k));
    fill_alphas(a, xs, B2_);
    fill_Dinvs(Dinvs, xs, a);
    D = Dinvs.inverse();
    fill_A(A, xs, D, B2_);
    A_solver = A.colPivHouseholderQr();
    flx = A_solver.solve(chi);
    k = vEf.dot(flx);
  }

  // We have converged on k = 1.
  const double sqrt_abs_B2 = std::sqrt(std::abs(B2_));
  const double B = B2_ > 0. ? sqrt_abs_B2 : -sqrt_abs_B2;

  // Get the current
  cur = B * D * flx;

  // The output info is in the format (2,NG) where the first line has
  // the flux spectrum, and the second has the diffusion coefficients.
  flux_.resize({NG});
  current_.resize({NG});
  diff_coeff_.resize({NG});
  for (std::size_t g = 0; g < NG; g++) {
    flux_(g) = flx(g);
    current_(g) = cur(g);
    diff_coeff_(g) = cur(g) / (B * flx(g));
  }
}

}  // namespace scarabee
