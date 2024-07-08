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

  m.def("set_logging_level", &set_logging_level,
        R"(Sets the verbosity of logging output.
     
Parameters
----------
level : LogLevel
        Minimmum logging information level written to console/file.
     )",
        py::arg("level"));

  m.def("set_output_file", &set_output_file,
        R"(Sets the name of an optional output file for logging.
    
Parameters
----------
fname : str
        Name of the output file.
    )",
        py::arg("fname"));
}
