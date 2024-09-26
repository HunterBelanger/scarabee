#include <utils/spherical_harmonics.hpp>
#include <utils/constants.hpp>

#include <cmath>

namespace scarabee {

SphericalHarmonics::SphericalHarmonics(
    const std::size_t& L, const std::vector<double>& azimuthal_angle,
    const std::vector<double>& polar_angle)
    : L_(L), Nlj_((L_ + 1) * (L_ + 1)) {
  // calculate associated legendre for all cosine of polar angles,
  // the polar_angle contains half of the angle between [0, PI/2]
  // therefor, evaluate the remaining half between [PI/2 , PI]
  // if the indexing of angle between [0, PI/2] is i, then
  // index of next half will be (no of angle ) + i
  std::size_t n_polar_angle = polar_angle.size();
  xt::xtensor<double, 2> associated_legendre_;
  // associated_legendre_ = std::vector<std::vector<double>>(2 * n_polar_angle,
  // std::vector<double>((L_+1)*(L_+2) / 2, 0.));
  associated_legendre_.resize({2 * n_polar_angle, (L_ + 1) * (L_ + 2) / 2});
  associated_legendre_.fill(0.);
  const double sqrt_inv_4PI = 0.5 / std::sqrt(PI);
  std::size_t p = 0;
  double theta;
  for (std::size_t pp = 0; pp < 2 * n_polar_angle; pp++) {
    if (pp < n_polar_angle) {
      // first half angles between [0, PI/2]
      p = pp;
      theta = polar_angle[p];
    } else {
      // remaining half angles between [PI/2, PI]
      p = pp - n_polar_angle;
      theta = PI - polar_angle[p];
    }
    const double cos_theta = std::cos(theta);
    std::size_t it_lj = 0;
    for (std::size_t l = 0; l <= L_; l++) {
      const double factor =
          std::sqrt(static_cast<double>(2 * l + 1)) * sqrt_inv_4PI;
      for (std::size_t j = 0; j <= l; j++) {
        const double P_j_l =
            eval_associated_legendre(l, static_cast<int>(j), cos_theta);
        if (j == 0) {
          // associated_legendre_[pp][it_lj] = factor * P_j_l;
          associated_legendre_(pp, it_lj) = factor * P_j_l;
        } else {
          const double sqrt_lj =
              std::sqrt(2. * factorial(l - j) / factorial(l + j));
          // associated_legendre_[pp][it_lj] = factor * sqrt_lj * P_j_l;
          associated_legendre_(pp, it_lj) = factor * sqrt_lj * P_j_l;
        }
        it_lj++;
      }
    }
  }

  // calculate the cos and sin of azimuthal angle phi for
  // angles corresponds to forward and backward direction entry
  const std::size_t n_azimuthal_angle = azimuthal_angle.size();
  xt::xtensor<double, 2> cos_phi_;
  cos_phi_.resize({n_azimuthal_angle * 2, L_});
  cos_phi_.fill(0.);
  xt::xtensor<double, 2> sin_phi_;
  sin_phi_.resize({n_azimuthal_angle * 2, L_});
  sin_phi_.fill(0.);
  // cos_phi_ = std::vector<std::vector<double>>(n_azimuthal_angle * 2,
  // std::vector<double>(L_, 0.)); sin_phi_ =
  // std::vector<std::vector<double>>(n_azimuthal_angle * 2,
  // std::vector<double>(L_, 0.));
  std::size_t m = 0;
  double phi;
  for (std::size_t mm = 0; mm < 2 * n_azimuthal_angle; mm++) {
    if (mm < n_azimuthal_angle) {
      // azimuthal angle in the forward direction of tracks
      m = mm;
      phi = azimuthal_angle[m];
    } else {
      // azimuthal angle in the backward direction of tracks
      m = mm - n_azimuthal_angle;
      phi = PI + azimuthal_angle[m];
    }

    for (std::size_t j = 1; j <= L_; j++) {
      // cos_phi_[mm][j-1] = std::cos(j * phi);
      // sin_phi_[mm][j-1] = std::sin(j * phi);
      cos_phi_(mm, j - 1) = std::cos(static_cast<double>(j) * phi);
      sin_phi_(mm, j - 1) = std::sin(static_cast<double>(j) * phi);
    }
  }

  // pre-calculate the spherical harmonics for all given azimuthal and polar
  // angles
  all_harmonics_.resize(
      {2 * n_azimuthal_angle, 2 * n_polar_angle, (L_ + 1) * (L_ + 1)});
  all_harmonics_.fill(0.);
  for (std::size_t azm = 0; azm < 2 * n_azimuthal_angle; azm++) {
    for (std::size_t p = 0; p < 2 * n_polar_angle; p++) {
      std::size_t it_lj = 0;
      for (std::size_t l = 0; l <= L_; l++) {
        for (int j = -static_cast<int>(l); j <= static_cast<int>(l); j++) {
          if (j == 0) {
            all_harmonics_(azm, p, it_lj) =
                associated_legendre_(p, l * (l + 1) / 2);
          } else if (j > 0) {
            const std::size_t abs_j = std::abs(j);
            all_harmonics_(azm, p, it_lj) =
                associated_legendre_(p, l * (l + 1) / 2 + abs_j) *
                cos_phi_(azm, abs_j - 1);
          } else {
            // j < 0
            const std::size_t abs_j = std::abs(j);
            all_harmonics_(azm, p, it_lj) =
                associated_legendre_(p, l * (l + 1) / 2 + abs_j) *
                sin_phi_(azm, abs_j - 1);
          }
          it_lj++;
        }

      }  // all scattering moments
    }    // all polar angles
  }      // all azimuthal angles
}

double SphericalHarmonics::eval_spherical_harmonics(const std::size_t& l,
                                                    const int& j,
                                                    const double& theta,
                                                    const double& phi) const {
  // theta is the polar angle, get the cosine of it.
  const double u = std::cos(theta);
  if (std::abs(j) > l) {
    return 0.;
  } else if (j > 0) {
    const double factor =
        std::sqrt((2. * static_cast<double>(l) + 1.) * factorial(l - j) /
                  (2. * PI * factorial(l + j)));
    return factor * eval_associated_legendre(l, j, u) * std::cos(j * phi);
  } else if (j == 0) {
    const double factor =
        std::sqrt((2. * static_cast<double>(l) + 1.) / (4. * PI));
    return factor * eval_associated_legendre(l, j, u);
  } else {
    // j < 1
    const auto abs_j = std::abs(j);
    const double factor =
        std::sqrt((2. * static_cast<double>(l) + 1.) * factorial(l - abs_j) /
                  (2. * PI * factorial(l + abs_j)));
    return factor * eval_associated_legendre(l, abs_j, u) *
           std::sin(abs_j * phi);
  }

  // it should not get here.
  return 0.;
}

double SphericalHarmonics::eval_associated_legendre(const std::size_t& order,
                                                    const int& j,
                                                    const double& x) const {
  std::size_t abs_j = std::abs(j);
  const double associated_l = eval_legendre_differentiation(order, abs_j, x);
  double factor = std::pow(-1., abs_j) *
                  std::pow(1. - x * x, static_cast<double>(abs_j) / 2);

  if (j < 0) {
    factor *=
        std::pow(-1., j) * factorial(order - abs_j) / factorial(order + abs_j);
  }
  return factor * associated_l;
}

double SphericalHarmonics::eval_legendre_differentiation(
    const std::size_t& order, const std::size_t& n, const double& x) const {
  if (n > 0) {
    if (order > 1) {  // a recursion relation for order > 1 and n > 0
      const double term1 = eval_legendre_differentiation(order - 1, n, x) * x;
      const double term2 = static_cast<double>(n) *
                           eval_legendre_differentiation(order - 1, n - 1, x);
      const double term3 = eval_legendre_differentiation(order - 2, n, x);
      return (static_cast<double>(2 * order - 1) * (term1 + term2) -
              static_cast<double>(order - 1) * term3) /
             static_cast<double>(order);
    } else if (order == 1 && n == 1) {
      return 1.;
    }
  } else if (n == 0) {  // special case
    return eval_legendre(order, x);
  }

  // for any n > 0, at order 0 and 1 return 0.
  return 0.;
}

}  // namespace scarabee
