#include <map>
#include <string>

#include <docopt.h>

static const std::string help =
"Scarabée\n\n"

"Usage:\n"
"   scarabée <input>\n"
"   scarabée (-h | --help)\n"
"   scarabée --version\n\n"

" Options:\n"
"   <input>       Name of the input file.\n"
"   --version     Show the version information.\n"
"   -h --help     Show this help message.\n";

int main(int argc, char** argv) {
  // Initialize docopt
  std::map<std::string, docopt::value> args = docopt::docopt(help, {argv + 1, argv + argc}, true, "Scarabée 0.1.0");

  return 0;
}