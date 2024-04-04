#include <utils/logging.hpp>

#include <spdlog/sinks/basic_file_sink.h>

#include <memory>

struct _private_logger_init {
  _private_logger_init() {
    spdlog::set_pattern("[%^%l%$] %v");
  }
};

static _private_logger_init _pli;

void set_logging_level(LogLevel level) {
  spdlog::set_level(level);
}

void set_output_file(const std::string& fname) {
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(fname, true);
  spdlog::default_logger()->sinks().push_back(file_sink);
}

