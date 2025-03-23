#ifndef SCARABEE_CONSTANTS_H
#define SCARABEE_CONSTANTS_H

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace scarabee {

constexpr double SQRT_2{1.41421356237309504880168872420969808};
constexpr double PI{3.14159265358979323846264338327950288};
constexpr double PI_2{0.5 * PI};
constexpr double INF{std::numeric_limits<double>::max()};
constexpr double SURFACE_COINCIDENT{1E-11};
constexpr double VEC_FP_TOL{1.E-10};
constexpr double N_MASS_AMU{1.00866491595};
constexpr double N_AVAGADRO{0.602214076};  // [10^24 / mol]
constexpr std::size_t MAX_SURFS{5};
#define NDL_ENV_VAR "SCARABEE_ND_LIBRARY"

extern const std::map<std::string, double> NATURAL_ABUNDANCES;
extern const std::map<std::string, std::vector<std::string>> ELEMENT_ISOTOPES;
extern const std::map<std::string, double> ISOTOPE_MASSES;

struct ElementInfo {
  std::string name;
  std::string symbol;
};

extern const std::array<ElementInfo, 119> ELEMENTS;



}  // namespace scarabee

#endif
