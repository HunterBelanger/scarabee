#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <data/depletion_chain.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_DepletionChain(py::module& m) {
  py::class_<NoTarget>(m, "NoTarget",
                       "Represents the end of a decay / transmutation chain "
                       "with no resulting nuclide.")
      .def(py::init<>());

  py::class_<SingleTarget>(m, "SingleTarget",
                           "Represents a single possible resulting nuclide "
                           "from a decay or transmutation.")
      .def(py::init<const std::string& /*target*/>(),
           "Creates a single nuclide depletion target.\n\n"
           "Parameters\n"
           "----------\n"
           "target : string\n"
           "    Name of the target nuclide.\n",
           py::arg("target"))

      .def_property_readonly("target", &SingleTarget::target,
                             "Name of the target nuclide.");

  py::class_<BranchingTargets::Branch>(
      m, "Branch",
      "Represents a single branch (or target) when multiple nuclide targets "
      "are possible from a decay or nuclear reaction. This should not be used "
      "for fission yields ! Instead use :class:`FissionYields`.")
      .def(py::init<const std::vector<Branch>& /*branches*/>(),
           "Creates a BranchingTarget object.\n\n"
           "Parameters\n"
           "----------\n"
           "branches : list of :class:`Branch`\n"
           "    Possible branches for the decay or transmutation.\n",
           py::arg("branches"))
      .def_property_readonly("branches", &BranchingTarget::branches,
        "List of all branches with the target name and branch ratio.");
}