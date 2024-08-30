#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <xtensor-python/pytensor.hpp>

#include <reflector_sn.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_ReflectorSN(py::module& m) {
  py::class_<ReflectorSN>(m, "ReflectorSN")
  .def(py::init<const std::vector<std::shared_ptr<CrossSection>>& /*xs*/,
                const xt::xtensor<double, 1>& /*dx*/>())
  .def("solve", &ReflectorSN::solve);
}