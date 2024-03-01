cmake_minimum_required(VERSION 3.11)
project(scarabee CXX)

# Get FetchContent for downloading dependencies
include(FetchContent)

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
# Get NDArray
message(STATUS "Downloading NDArray")
FetchContent_Declare(NDArray
  GIT_REPOSITORY https://github.com/HunterBelanger/ndarray.git
  GIT_TAG        master
)
set(NDARRAY_INSTALL OFF CACHE BOOL "Install NDArray")
FetchContent_MakeAvailable(NDArray)

#===============================================================================
add_executable(scarabee src/main.cpp src/gauss_kronrod.cpp src/chebyshev.cpp src/bickley.cpp src/transport_xs.cpp src/cylindrical_cell.cpp src/cylindrical_flux_solver.cpp)
set_target_properties(scarabee PROPERTIES OUTPUT_NAME "scarabée")
target_include_directories(scarabee PRIVATE include)
target_compile_features(scarabee PRIVATE cxx_std_20)
target_link_libraries(scarabee PUBLIC NDArray::NDArray Eigen3::Eigen)

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
