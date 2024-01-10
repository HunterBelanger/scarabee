cmake_minimum_required(VERSION 3.11)
project(scarabee CXX)

# Get FetchContent for downloading dependencies
include(FetchContent)

#===============================================================================
# Get YAML-CPP version 0.8.0 message(STATUS "Downloading yaml-cpp 0.8.0")
FetchContent_Declare(yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG        0.8.0
)
set(YAML_CPP_BUILD_SHARED_LIBS OFF CACHE BOOL "Build yaml-cpp shared library")
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Enable testing")
set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "Enable parse tools")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Enable contrib stuff in library")
set(YAML_CPP_INSTALL OFF CACHE BOOL "Enable generation of install target")
FetchContent_MakeAvailable(yaml-cpp)

#===============================================================================
# Get PapillonNDL develop
message(STATUS "Downloading PapillonNDL feature/sab branch")
FetchContent_Declare(PapillonNDL
  GIT_REPOSITORY https://github.com/HunterBelanger/papillon-ndl.git
  GIT_TAG        develop
)
set(PNDL_SHARED OFF CACHE BOOL "Build PapillonNDL as a shared library")
set(PNDL_INSTALL OFF CACHE BOOL "Install the PapillonNDL library and header files")
set(PNDL_PYTHON OFF CACHE BOOL "Enable Python interface to PapillonNDL")
FetchContent_MakeAvailable(PapillonNDL)

#===============================================================================
# Get Eigen version 3.4.0
message(STATUS "Downloading Eigen")
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/HunterBelanger/eigen.git
  GIT_TAG        fix/uninstall-target
)
set(BUILD_TESTING OFF)
set(EIGEN_TEST_NOQT ON)
set(EIGEN_BUILD_CMAKE_PACKAGE OFF)
FetchContent_MakeAvailable(eigen)

#===============================================================================
# Get docopt
message(STATUS "Downloading docopt.cpp v0.6.3")
FetchContent_Declare(docopt.cpp
  GIT_REPOSITORY https://github.com/docopt/docopt.cpp.git
  GIT_TAG        v0.6.3
)
if(NOT docopt.cpp_POPULATED)
  FetchContent_Populate(docopt.cpp)
  add_subdirectory(${docopt.cpp_SOURCE_DIR} ${docopt.cpp_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

#===============================================================================
# Get NDArray
message(STATUS "Downloading NDArray")
FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray.git
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

#===============================================================================
add_executable(scarabee src/main.cpp src/gauss_kronrod.cpp src/chebyshev.cpp src/bickley.cpp src/mg_cross_sections.cpp src/cylindrical_cell.cpp src/cylindrical_flux_solver.cpp)
set_target_properties(scarabee PROPERTIES OUTPUT_NAME "scarabée")
target_include_directories(scarabee PRIVATE include)
target_compile_features(scarabee PRIVATE cxx_std_20)
target_link_libraries(scarabee PUBLIC docopt_s yaml-cpp NDArray::NDArray Eigen3::Eigen PapillonNDL::PapillonNDL)

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Comile options for Windows
  target_compile_options(scarabee PRIVATE /W3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
  target_compile_options(scarabee PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
  target_compile_options(scarabee PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
  target_compile_options(scarabee PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
endif()

#===============================================================================
# Allow for install of executable
install(TARGETS scarabee
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
