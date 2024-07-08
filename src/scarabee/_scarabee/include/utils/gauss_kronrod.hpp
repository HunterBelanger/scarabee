#ifndef SCARABEE_GUASSKRONROD_H
#define SCARABEE_GAUSSKRONROD_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace scarabee {

template <std::size_t NK>
struct GaussKronrodQuadrature {
  static std::pair<double, double> integrate(
      const std::function<double(double)>& f, double x_low, double x_hi) {
    // First calculate integral using the Gauss-Legendre quadrature
    double glIntegral = 0.;
    double gkIntegral = 0;
    double f_x = 0.;
    double xi = 0.;
    double x = 0.;
    for (std::size_t i = 0; i < glWeights.size(); i++) {
      f_x = 0.;

      // Get the positive point
      xi = abscissae[i];
      x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
      f_x = f(x);
      glIntegral += f_x * glWeights[i];
      gkIntegral += f_x * weights[i];

      // Get the negative point
      if (xi != 0.) {
        xi = -xi;
        x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
        f_x = f(x);
        glIntegral += f_x * glWeights[i];
        gkIntegral += f_x * weights[i];
      }
    }

    // Now calculate integral using the Gauss-Kronrod quadrature
    // Get contribution from Gauss-Kronrod points
    for (std::size_t i = glWeights.size(); i < weights.size(); i++) {
      // Get the positive point
      xi = abscissae[i];
      x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
      gkIntegral += f(x) * weights[i];

      // Get the negative point
      if (xi != 0.) {
        xi = -xi;
        x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
        gkIntegral += f(x) * weights[i];
      }
    }

    glIntegral *= 0.5 * (x_hi - x_low);
    gkIntegral *= 0.5 * (x_hi - x_low);

    // Calculate the error
    const double err = std::abs((glIntegral - gkIntegral) / gkIntegral);

    return {gkIntegral, err};
  }

  static std::pair<double, double> integrate(
      const std::function<double(double)>& f, double x_low, double x_hi,
      double max_rel_err, std::size_t max_splits) {
    // First calculate integral using the Gauss-Legendre quadrature
    double glIntegral = 0.;
    double gkIntegral = 0;
    double f_x = 0.;
    double xi = 0.;
    double x = 0.;
    for (std::size_t i = 0; i < glWeights.size(); i++) {
      f_x = 0.;

      // Get the positive point
      xi = abscissae[i];
      x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
      f_x = f(x);
      glIntegral += f_x * glWeights[i];
      gkIntegral += f_x * weights[i];

      // Get the negative point
      if (xi != 0.) {
        xi = -xi;
        x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
        f_x = f(x);
        glIntegral += f_x * glWeights[i];
        gkIntegral += f_x * weights[i];
      }
    }

    // Now calculate integral using the Gauss-Kronrod quadrature
    // Get contribution from Gauss-Kronrod points
    for (std::size_t i = glWeights.size(); i < weights.size(); i++) {
      // Get the positive point
      xi = abscissae[i];
      x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
      gkIntegral += f(x) * weights[i];

      // Get the negative point
      if (xi != 0.) {
        xi = -xi;
        x = 0.5 * (x_hi - x_low) * xi + 0.5 * (x_low + x_hi);
        gkIntegral += f(x) * weights[i];
      }
    }

    glIntegral *= 0.5 * (x_hi - x_low);
    gkIntegral *= 0.5 * (x_hi - x_low);

    // Calculate the error
    double err = std::abs((glIntegral - gkIntegral) / gkIntegral);

    if (err > max_rel_err && max_splits > 0) {
      // We split into two intervals
      const double x_mid = 0.5 * (x_low + x_hi);
      const double new_max_err = 0.5 * max_rel_err;
      const std::size_t new_max_splits = max_splits - 1;

      auto left = integrate(f, x_low, x_mid, new_max_err, new_max_splits);
      auto right = integrate(f, x_mid, x_hi, new_max_err, new_max_splits);

      gkIntegral = left.first + right.first;
      err = left.second + right.second;
    }

    return {gkIntegral, err};
  }

  static constexpr std::size_t order() { return 3 * NK + 1; }

  static const std::vector<double> abscissae;
  static const std::vector<double> weights;
  static const std::vector<double> glWeights;
};

}  // namespace scarabee

#endif
