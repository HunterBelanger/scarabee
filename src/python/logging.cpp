#include <pybind11/pybind11.h>
#include <utils/logging.hpp>

namespace py = pybind11;

using namespace scarabee;

void init_Logging(py::module& m) {
  py::enum_<LogLevel>(m, "LogLevel")
      .value("Critical", LogLevel::critical)
      .value("Debug", LogLevel::debug)
      .value("Error", LogLevel::err)
      .value("Info", LogLevel::info)
      .value("Off", LogLevel::off)
      .value("Trace", LogLevel::trace)
      .value("Warning", LogLevel::warn);

  m.def("set_logging_level", &set_logging_level);
  m.def("set_output_file", &set_output_file);
}
