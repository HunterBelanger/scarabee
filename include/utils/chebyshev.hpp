#ifndef SCARABEE_CHEBYSHEV_H
#define SCARABEE_CHEBYSHEV_H

#include <functional>
#include <span>
#include <vector>

namespace scarabee {

std::vector<double> chebyshev_fit(const std::function<double(double)>& func,
                                  double a, double b, std::size_t n);

double chebyshev_eval(double x, double a, double b, std::span<const double> c);

}  // namespace scarabee

#endif
