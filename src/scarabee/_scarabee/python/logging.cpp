#include <pybind11/pybind11.h>
#include <utils/logging.hpp>

#include <string>

namespace py = pybind11;

using namespace scarabee;

void scarabee_log(LogLevel lvl, const std::string& mssg) {
  switch (lvl) {
    case LogLevel::critical:
      spdlog::critical(mssg);
      break;

    case LogLevel::debug:
      spdlog::debug(mssg);
      break;

    case LogLevel::err:
      spdlog::error(mssg);
      break;

    case LogLevel::info:
      spdlog::info(mssg);
      break;

    case LogLevel::warn:
      spdlog::warn(mssg);
      break;

    default:
      spdlog::info(mssg);
      break;
  }
}

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

  m.def("scarabee_log", &scarabee_log,
        "Logs program information.\n\n"
        "Parameters\n"
        "----------\n"
        "level : LogLevel\n"
        "        Type of information to be logged.\n"
        "mssg : str\n"
        "       Message to be logged.\n",
        py::arg("level"), py::arg("mssg"));

  m.def("set_output_file", &set_output_file,
      "Sets the output file for the Scarabee log.\n\n"
      "Parameters\n"
      "----------\n"
      "fname : str\n"
      "        Name of log file.\n");
}
