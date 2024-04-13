#include <pybind11/pybind11.h>

#include <moc/boundary_condition.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_BoundaryCondition(py::module& m) {
  py::enum_<BoundaryCondition>(m, "BoundaryCondition")
  .value("Reflective", BoundaryCondition::Reflective)
  .value("Vacuum", BoundaryCondition::Vacuum);
}

