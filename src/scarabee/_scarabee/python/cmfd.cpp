#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/cmfd.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CMFD(py::module& m) {
  py::class_<CMFD, std::shared_ptr<CMFD>>(m, "CMFD")
    .def(py::init<const std::vector<double>& /*dx*/,
                  const std::vector<double>& /*dy*/>(),
        "A Cartesian mesh for accelerating MOC convergence.\n\n"
        "Parameters\n"
        "----------\n"
        "dx : list of float\n"
        "     List of all x widths.\n"
        "dy : list of float\n"
        "     List of all y heights.\n",
        py::arg("dx"), py::arg("dy"));
}