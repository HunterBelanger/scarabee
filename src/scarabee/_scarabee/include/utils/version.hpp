#ifndef SCARABEE_VERSION_H
#define SCARABEE_VERSION_H

#define XSTR(x) #x
#define STR(x) XSTR(x)

namespace scarabee {

const unsigned int VERSION_MAJOR{SCARABEE_MAJOR_VERSION};
const unsigned int VERSION_MINOR{SCARABEE_MINOR_VERSION};
const unsigned int VERSION_PATCH{SCARABEE_PATCH_VERSION};
const char* VERSION_STRING{STR(SCARABEE_MAJOR_VERSION) "." STR(
    SCARABEE_MINOR_VERSION) "." STR(SCARABEE_PATCH_VERSION)};

}  // namespace scarabee

#endif
