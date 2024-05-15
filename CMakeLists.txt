cmake_minimum_required(VERSION 3.11)
project(scarabee CXX)

option(SCARABEE_USE_OMP "Compile Scarabée with OpenMP for shared memory parallelism" ON)

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
# Get SPDLOG
message(STATUS "Downloading spdlog v1.13.0")
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG        v1.13.0
)
set(SPDLOG_BUILD_PIC ON CACHE BOOL "Build position independent code (-fPIC)")
FetchContent_MakeAvailable(spdlog)

#===============================================================================
# Get HTL
message(STATUS "Downloading HTL")
FetchContent_Declare(htl
  GIT_REPOSITORY https://github.com/HunterBelanger/htl.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(htl)

#===============================================================================
# Get XTL
message(STATUS "Downloading xtl v0.7.7")
FetchContent_Declare(xtl
  GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
  GIT_TAG        0.7.7
)
FetchContent_MakeAvailable(xtl)

#===============================================================================
# Get XSIMD
message(STATUS "Downloading xsimd v12.1.1")
FetchContent_Declare(xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG        12.1.1
)
FetchContent_MakeAvailable(xsimd)

#===============================================================================
# Get XTENSOR
message(STATUS "Downloading xtensor v0.25.0")
FetchContent_Declare(xtensor
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
  GIT_TAG        0.25.0
)
set(XTENSOR_USE_XSIMD ON CACHE BOOL "Enable SIMD acceleration.")
FetchContent_MakeAvailable(xtensor)

#===============================================================================
# Get Pybind11
set(PYBIND11_FINDPYTHON ON)
message(STATUS "Downloading pybind11 v2.12.0")
FetchContent_Declare(pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11
  GIT_TAG        v2.12.0
)
FetchContent_MakeAvailable(pybind11)

#===============================================================================
# Get XTENSOR-PYTHON
message(STATUS "Downloading xtensor-python v0.27.0")
FetchContent_Declare(xtensor-python
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-python.git
  GIT_TAG        0.27.0
)
FetchContent_MakeAvailable(xtensor-python)

#===============================================================================
# Get HighFive
message(STATUS "Downloading HighFive")
set(HDF5_USE_STATIC_LIBRARIES ON)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_USE_XTENSOR OFF)
set(HighFive_FIND_QUIETLY ON)
FetchContent_Declare(HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
  GIT_TAG        v2.9.0
)
FetchContent_MakeAvailable(HighFive)

#===============================================================================
# Make python library
pybind11_add_module(scarabee  src/logging.cpp
                              src/gauss_kronrod.cpp
                              src/chebyshev.cpp
                              src/math.cpp
                              src/cross_section.cpp
                              src/material.cpp
                              src/nd_library.cpp
                              src/cylindrical_cell.cpp
                              src/cylindrical_flux_solver.cpp
                              src/surface.cpp
                              src/flat_source_region.cpp
                              src/cell.cpp
                              src/empty_cell.cpp
                              src/simple_pin_cell.cpp
                              src/pin_cell.cpp
                              src/cartesian_2d.cpp
                              src/track.cpp
                              src/legendre.cpp
                              src/yamamoto_tabuchi.cpp
                              src/moc_driver.cpp
                              src/criticality_spectrum.cpp
                              #=================================================
                              src/python/scarabee.cpp
                              src/python/vector.cpp
                              src/python/direction.cpp
                              src/python/logging.cpp
                              src/python/cross_section.cpp
                              src/python/material.cpp
                              src/python/nd_library.cpp
                              src/python/cylindrical_cell.cpp
                              src/python/cylindrical_flux_solver.cpp
                              src/python/polar_quadrature.cpp
                              src/python/boundary_condition.cpp
                              src/python/simulation_mode.cpp
                              src/python/track.cpp 
                              src/python/cell.cpp
                              src/python/empty_cell.cpp
                              src/python/simple_pin_cell.cpp
                              src/python/pin_cell.cpp
                              src/python/cartesian_2d.cpp
                              src/python/moc_driver.cpp
                            )

target_include_directories(scarabee PRIVATE include)
target_compile_features(scarabee PRIVATE cxx_std_20)
target_link_libraries(scarabee PUBLIC xtl xtensor xtensor-python htl HighFive Eigen3::Eigen spdlog::spdlog)
target_include_directories(scarabee PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
)

# Find OpenMP if desired
if(SCARABEE_USE_OMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(scarabee PUBLIC OpenMP::OpenMP_CXX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Compile options for Windows
      target_compile_options(scarabee PRIVATE /openmp:llvm) 
    endif()
    target_compile_definitions(scarabee PUBLIC SCARABEE_USE_OMP)
  endif()
endif()
