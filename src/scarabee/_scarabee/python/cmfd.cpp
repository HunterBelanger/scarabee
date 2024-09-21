#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <moc/cmfd.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_CMFD(py::module& m) {
  py::class_<CMFD, std::shared_ptr<CMFD>>(m, "CMFD")
    .def(py::init<const std::vector<double>& /*dx*/,
                  const std::vector<double>& /*dy*/,
                  const std::vector<std::pair<std::size_t, std::size_t>>& /*groups*/>(),
        "A Cartesian mesh for accelerating MOC convergence.\n\n"
        "Parameters\n"
        "----------\n"
        "dx : list of float\n"
        "     List of all x widths.\n"
        "dy : list of float\n"
        "     List of all y heights.\n"
        "groups : list of 2D tuples of ints.\n"
        "         The scheme for condensing energy groups.\n",
        py::arg("dx"), py::arg("dy"), py::arg("groups"));
}