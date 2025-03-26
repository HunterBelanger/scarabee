#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <xtensor-python/pytensor.hpp>

#include <data/depletion_chain.hpp>
#include <utils/logging.hpp>
#include <utils/scarabee_exception.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_DepletionChain(py::module& m) {
  //===========================================================================
  // TARGETS
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
      .def(py::init<>())
      .def_readwrite("target", &BranchingTargets::Branch::target,
                     "Target nuclide.")
      .def_readwrite("branch_ratio", &BranchingTargets::Branch::branch_ratio,
                     "Branch ratio.");

  py::class_<BranchingTargets>(
      m, "BranchingTargets",
      "Represents a set of possible targets for a sinble decay or reaction.")
      .def(
          py::init<const std::vector<BranchingTargets::Branch>& /*branches*/>(),
          "Creates a BranchingTargets object.\n\n"
          "Parameters\n"
          "----------\n"
          "branches : list of :class:`Branch`\n"
          "    Possible branches for the decay or transmutation.\n",
          py::arg("branches"))
      .def_property_readonly(
          "branches", &BranchingTargets::branches,
          "List of all branches with the target name and branch ratio.");
  
  //===========================================================================
  // Fission Yields
  py::class_<FissionYields>(
      m, "FissionYields",
      "Represents a set of multiple targets with different yields (not "
      "probabilities) resulting from fission. The yields can be a function of "
      "the incident neutron energy.")
      .def(py::init<const std::vector<std::string>& /*targets*/,
                    const std::vector<double>& /*incident_energies*/,
                    const xt::xtensor<double, 2>& /*yields*/>(),
           "Constructs a tabulated set of fission yields.\n\n"
           "Parameters\n"
           "----------\n"
           "targets : list of string\n"
           "    List of all possible targets.\n"
           "incident_energies : list of float\n"
           "    Sorted list of tabulated incident energies for which yields "
           "are provided.\n"
           "yields : ndarray\n"
           "    2D Numpy array containing the fission yields where first axis "
           "is incident energy and second is target.\n")
      .def_property_readonly("size", &FissionYields::size, "Number of targets.")
      .def_property_readonly("targets", &FissionYields::targets,
                             "List of possible target nuclides.")
      .def_property_readonly(
          "incident_energies", &FissionYields::incident_energies,
          "List of incident energies at which yields are tabulated.")
      .def("yield", &FissionYields::yield,
           "Computes the fission yield for a specified target and incident "
           "neutron energy.\n\n"
           "Parameters\n"
           "----------\n"
           "t : int\n"
           "    Index of the target nuclide.\n"
           "E : float\n"
           "    Incident energy of the neutron inducing fission.\n\n"
           "Returns\n"
           "-------\n"
           "float\n"
           "    Yield of the specified target at the specified incident "
           "energy.\n",
           py::arg("t"), py::arg("E"));
  
  //===========================================================================
  // Chain Entry
  py::class_<ChainEntry>(
      m, "ChainEntry",
      "Represents all possible means by which a nuclide can be removed from a "
      "material. This includes radioactive decay and possible transmutations. "
      "For all possible loss paths, the possible resulting nuclide are also "
      "provided as different types of targets.")
      .def(py::init<>())

      .def_property(
          "half_life", [](const ChainEntry& c) { return c.half_life(); },
          [](ChainEntry& c, std::optional<double> hl) {
            if (hl.has_value() && hl.value() <= 0.) {
              const auto mssg = "Half life must be > 0.";
              spdlog::error(mssg);
              throw ScarabeeException(mssg);
            }
            c.half_life() = hl;
          },
          "Half life of the nuclide.")

      .def_property(
          "decay_targets",
          [](const ChainEntry& c) { return c.decay_targets(); },
          [](ChainEntry& c, std::optional<Target> t) { c.decay_targets() = t; },
          "Targets for the radioactive decay of the nuclide.")
          
      .def_property(
          "n_gamma",
          [](const ChainEntry& c) { return c.n_gamma(); },
          [](ChainEntry& c, std::optional<Target> t) { c.n_gamma() = t; },
          "Targets for the (n,gamma) transmutation of the nuclide.")
          
      .def_property(
          "n_2n",
          [](const ChainEntry& c) { return c.n_2n(); },
          [](ChainEntry& c, std::optional<Target> t) { c.n_2n() = t; },
          "Targets for the (n,2n) transmutation of the nuclide.")
      
      .def_property(
          "n_3n",
          [](const ChainEntry& c) { return c.n_3n(); },
          [](ChainEntry& c, std::optional<Target> t) { c.n_3n() = t; },
          "Targets for the (n,3n) transmutation of the nuclide.")
      
      .def_property(
          "n_p",
          [](const ChainEntry& c) { return c.n_p(); },
          [](ChainEntry& c, std::optional<Target> t) { c.n_p() = t; },
          "Targets for the (n,p) transmutation of the nuclide.")

      .def_property(
          "n_alpha",
          [](const ChainEntry& c) { return c.n_alpha(); },
          [](ChainEntry& c, std::optional<Target> t) { c.n_alpha() = t; },
          "Targets for the (n,alpha) transmutation of the nuclide.")
          
      .def_property(
          "n_fission",
          [](const ChainEntry& c) { return c.n_fission(); },
          [](ChainEntry& c, std::optional<FissionYields> f) { c.n_fission() = f; },
          "Energy dependent fission yields for the nuclide.");

  //===========================================================================
  // Depletion Chain
}