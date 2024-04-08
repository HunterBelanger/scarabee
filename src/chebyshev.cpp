#include <utils/chebyshev.hpp>
#include <utils/constants.hpp>
#include <utils/scarabee_exception.hpp>
#include <utils/logging.hpp>

#include <cmath>

namespace scarabee {

std::vector<double> chebyshev_fit(const std::function<double(double)>& func,
                                  double a, double b, std::size_t n) {
  const double bma = 0.5 * (b - a);
  const double bpa = 0.5 * (b + a);
  const double dn = static_cast<double>(n);
  const double fac = 2. / dn;

  std::vector<double> f(n, 0.);
  std::vector<double> c(n, 0.);

  for (std::size_t k = 0; k < n; k++) {
    const double dk = static_cast<double>(k);
    const double y = std::cos(PI * (dk + 0.5) / dn);
    f[k] = func(y * bma + bpa);
  }

  for (std::size_t j = 0; j < n; j++) {
    const double dj = static_cast<double>(j);
    double sum = 0.;

    for (std::size_t k = 0; k < n; k++) {
      const double dk = static_cast<double>(k);
      sum += f[k] * std::cos(PI * dj * (dk + 0.5) / dn);
    }

    c[j] = fac * sum;
  }

  return c;
}

double chebyshev_eval(double x, double a, double b, std::span<const double> c) {
  if (x < a || x > b) {
    auto mssg = "Argument x must be in the interval [a,b].";
    spdlog::error(mssg);
    throw ScarabeeException(mssg);
  }

  const double y = (2. * x - a - b) / (b - a);
  const double y2 = 2 * y;

  double d = 0.;
  double sv = 0.;
  double dd = 0.;
  for (std::size_t j = c.size() - 1; j >= 1; j--) {
    sv = d;
    d = y2 * d - dd + c[j];
    dd = sv;
  }

  return y * d - dd + 0.5 * c[0];
}

}  // namespace scarabee
