#include <utils/logging.hpp>

#include <spdlog/sinks/basic_file_sink.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <memory>

namespace scarabee {

struct _private_logger_init {
  _private_logger_init() {
    // First, we need to remove the default sink that goes to stdout
    spdlog::default_logger()->sinks().clear();

    // Create a new custom sink
    auto python_sink = std::make_shared<PythonSinkMT>();

    // Set the patern for logging output
    python_sink->set_pattern("[%^%l%$] %v");

    // Save the sink to the logger
    spdlog::default_logger()->sinks().push_back(python_sink);
  }
};

static _private_logger_init _pli;

void set_logging_level(LogLevel level) { spdlog::set_level(level); }

void set_output_file(const std::string& fname) {
  auto file_sink =
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(fname, true);
  spdlog::default_logger()->sinks().push_back(file_sink);
}

}  // namespace scarabee
