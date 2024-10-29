#include <pybind11/pybind11.h>

#include <utils/simulation_mode.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_SimulationMode(py::module& m) {
  py::enum_<SimulationMode>(m, "SimulationMode")
      .value("FixedSource", SimulationMode::FixedSource)
      .value("Keff", SimulationMode::Keff);
}
