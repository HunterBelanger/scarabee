#ifndef SCARABEE_LOGGING_H
#define SCARABEE_LOGGING_H

#include <spdlog/spdlog.h>

#include <string>

using LogLevel = spdlog::level::level_enum;

void set_logging_level(LogLevel level);

void set_output_file(const std::string& fname);

#endif
