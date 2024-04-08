#ifndef SCARABEE_EXCEPTION_H
#define SCARABEE_EXCEPTION_H

#include <exception>
#include <filesystem>
#include <source_location>
#include <string>

namespace scarabee {

/**
 * @brief Class used for exceptions by the library.
 */
class ScarabeeException : public std::exception {
 public:
  ScarabeeException() : message("\n") {}
  /**
   * @param mssg Error message.
   * @param file Location where the error occurred;
   */
  ScarabeeException(
      const std::string& mssg,
      std::source_location location = std::source_location::current())
      : message("\n") {
    add_to_error_message(mssg, location);
  }
  ~ScarabeeException() = default;

  /**
   * @brief Adds details to the exception message as it is passed up the stack.
   * @param mssg Error message.
   * @param location Location where the error was thrown.
   */
  void add_to_exception(
      const std::string& mssg,
      std::source_location location = std::source_location::current()) {
    add_to_error_message(mssg, location);
  }

  const char* what() const noexcept override { return message.c_str(); }

 private:
  std::string message;

  void add_to_error_message(const std::string& mssg,
                            const std::source_location& location) {
    // Go through original string and determine line breaks
    std::string mssg_tmp = mssg;
    int nbreaks = static_cast<int>(mssg_tmp.size()) / 80;
    if ((mssg_tmp.size() % 80) == 0) nbreaks--;
    if (nbreaks > 0) {
      for (size_t b = 0; b < static_cast<size_t>(nbreaks); b++) {
        // Get index of break.
        size_t ind = 80 * (b + 1);

        // Work backwards to first space
        while (mssg_tmp[ind] != ' ') {
          if (ind > 0)
            ind--;
          else {
            ind = 80 * (b + 1);
            break;
          }
        }

        mssg_tmp[ind] = '\n';
      }
    }

    std::filesystem::path src_file(location.file_name());
    std::string fname;
    fname.reserve(src_file.string().size());
    bool save = false;
    for (const auto& comp : src_file) {
      if (save) {
        fname += "/" + comp.string();
      } else if (comp.string() == "src" || comp.string() == "include") {
        save = true;
        fname += comp.string();
      }
    }

    std::string tmp = "\n";
    tmp +=
        " #--------------------------------------------------------------------"
        "-------------\n";
    // tmp += " # File: " + std::string(location.file_name()) + "\n";
    tmp += " # File: " + fname + "\n";
    tmp += " # Function: " + std::string(location.function_name()) + "\n";
    tmp += " # Line: " + std::to_string(location.line()) + "\n";
    tmp += " # \n";
    tmp += " # ";

    for (const auto& c : mssg_tmp) {
      if (c == '\n') {
        tmp += "\n # ";
      } else {
        tmp += c;
      }
    }
    tmp += "\n";
    tmp +=
        " #--------------------------------------------------------------------"
        "-------------";

    message = tmp + message;
  }
};

}  // namespace scarabee

#endif
