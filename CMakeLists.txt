cmake_minimum_required(VERSION 3.11)
project(${SKBUILD_PROJECT_NAME}
        VERSION ${SKBUILD_PROJECT_VERSION}
        LANGUAGES CXX)

option(SCARABEE_USE_OMP "Compile Scarabée with OpenMP for shared memory parallelism" ON)

# Get FetchContent for downloading dependencies
include(FetchContent)

#===============================================================================
# Get Eigen version 3.4.0
message(STATUS "Downloading Eigen")
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/HunterBelanger/eigen.git
  GIT_TAG        fix/uninstall-target
  EXCLUDE_FROM_ALL
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
  EXCLUDE_FROM_ALL
)
set(SPDLOG_BUILD_PIC ON CACHE BOOL "Build position independent code (-fPIC)")
FetchContent_MakeAvailable(spdlog)

#===============================================================================
# Get HTL
message(STATUS "Downloading HTL")
FetchContent_Declare(htl
  GIT_REPOSITORY https://github.com/HunterBelanger/htl.git
  GIT_TAG        master
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(htl)

#===============================================================================
# Get XTL
message(STATUS "Downloading xtl v0.7.7")
FetchContent_Declare(xtl
  GIT_REPOSITORY https://github.com/xtensor-stack/xtl.git
  GIT_TAG        0.7.7
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(xtl)

#===============================================================================
# Get XSIMD
message(STATUS "Downloading xsimd v12.1.1")
FetchContent_Declare(xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG        12.1.1
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(xsimd)

#===============================================================================
# Get XTENSOR
message(STATUS "Downloading xtensor v0.25.0")
FetchContent_Declare(xtensor
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
  GIT_TAG        0.25.0
  EXCLUDE_FROM_ALL
)
set(XTENSOR_USE_XSIMD ON CACHE BOOL "Enable SIMD acceleration.")
FetchContent_MakeAvailable(xtensor)

#===============================================================================
# Find Pybind11 from Scikit-Build-Core
set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG QUIET)
if (NOT pybind11_FOUND)
  set(PYBIND11_FINDPYTHON ON)
  message(STATUS "Downloading pybind11 v2.12.0")
  FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG        v2.12.0
  )
  FetchContent_MakeAvailable(pybind11)
endif()

execute_process(COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/get_numpy_includes.py
                OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
                OUTPUT_STRIP_TRAILING_WHITESPACE
               )

#===============================================================================
# Get XTENSOR-PYTHON
message(STATUS "Downloading xtensor-python v0.27.0")
FetchContent_Declare(xtensor-python
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-python.git
  GIT_TAG        0.27.0
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(xtensor-python)
message(STATUS "Numpy includes: ${NUMPY_INCLUDE_DIRS}")

#===============================================================================
# Get HDF5
message(STATUS "Downloading HDF5 v1.14.3")
set(HDF5_ENABLE_Z_LIB_SUPPORT OFF)
set(HDF5_EXTERNALLY_CONFIGURED  1)
set(HDF5_ENABLE_SZIP_SUPPORT  OFF)
set(HDF5_BUILD_EXAMPLES       OFF)
set(HDF5_BUILD_TOOLS          OFF)
set(HDF5_BUILD_UTILS          OFF)
set(HDF5_BUILD_HL_LIB         OFF)
set(BUILD_SHARED_LIBS         OFF)
FetchContent_Declare(hdf5
  GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
  GIT_TAG        hdf5-1_14_3
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(hdf5)

#===============================================================================
# Get HighFive
message(STATUS "Downloading HighFive")
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_TEST_SPAN OFF)
set(HIGHFIVE_EXAMPLES OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_FIND_HDF5 OFF)
FetchContent_Declare(HighFive
  GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
  GIT_TAG        master
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(HighFive)

#===============================================================================
# Get ImApp
message(STATUS "Downloading ImApp")
FetchContent_Declare(ImApp
  GIT_REPOSITORY https://github.com/HunterBelanger/ImApp.git
  GIT_TAG        master
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(ImApp)

#===============================================================================
# Get Cereal
message(STATUS "Downloading Cereal")
set(BUILD_DOC OFF)
set(BUILD_SANDBOX OFF)
set(SKIP_PERFORMANCE_COMPARISON ON)
FetchContent_Declare(cereal
  GIT_REPOSITORY https://github.com/USCiLab/cereal.git
  GIT_TAG        master
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(cereal)

#===============================================================================
# Make python library
pybind11_add_module(_scarabee src/scarabee/_scarabee/constants.cpp
                              src/scarabee/_scarabee/water.cpp
                              src/scarabee/_scarabee/logging.cpp
                              src/scarabee/_scarabee/gauss_kronrod.cpp
                              src/scarabee/_scarabee/chebyshev.cpp
                              src/scarabee/_scarabee/math.cpp
                              src/scarabee/_scarabee/cross_section.cpp
                              src/scarabee/_scarabee/diffusion_cross_section.cpp
                              src/scarabee/_scarabee/material.cpp
                              src/scarabee/_scarabee/nd_library.cpp
                              src/scarabee/_scarabee/flux_calculator.cpp
                              src/scarabee/_scarabee/cylindrical_cell.cpp
                              src/scarabee/_scarabee/cylindrical_flux_solver.cpp
                              src/scarabee/_scarabee/surface.cpp
                              src/scarabee/_scarabee/flat_source_region.cpp
                              src/scarabee/_scarabee/cell.cpp
                              src/scarabee/_scarabee/empty_cell.cpp
                              src/scarabee/_scarabee/simple_pin_cell.cpp
                              src/scarabee/_scarabee/pin_cell.cpp
                              src/scarabee/_scarabee/cartesian_2d.cpp
                              src/scarabee/_scarabee/cmfd.cpp
                              src/scarabee/_scarabee/track.cpp
                              src/scarabee/_scarabee/legendre.cpp
                              src/scarabee/_scarabee/yamamoto_tabuchi.cpp
                              src/scarabee/_scarabee/moc_driver.cpp
                              src/scarabee/_scarabee/moc_plotter.cpp
                              src/scarabee/_scarabee/criticality_spectrum.cpp
                              src/scarabee/_scarabee/diffusion_data.cpp
                              src/scarabee/_scarabee/diffusion_geometry.cpp
                              src/scarabee/_scarabee/fd_diffusion_driver.cpp
                              src/scarabee/_scarabee/nem_diffusion_driver.cpp
                              src/scarabee/_scarabee/fuel_pin.cpp
                              src/scarabee/_scarabee/guide_tube.cpp
                              src/scarabee/_scarabee/burnable_poison_pin.cpp
                              src/scarabee/_scarabee/pwr_assembly.cpp
                              src/scarabee/_scarabee/reflector_sn.cpp
                              src/scarabee/_scarabee/spherical_harmonics.cpp
                              #=================================================
                              src/scarabee/_scarabee/python/scarabee.cpp
                              src/scarabee/_scarabee/python/vector.cpp
                              src/scarabee/_scarabee/python/direction.cpp
                              src/scarabee/_scarabee/python/logging.cpp
                              src/scarabee/_scarabee/python/xs1d.cpp
                              src/scarabee/_scarabee/python/xs2d.cpp
                              src/scarabee/_scarabee/python/cross_section.cpp
                              src/scarabee/_scarabee/python/diffusion_cross_section.cpp
                              src/scarabee/_scarabee/python/material.cpp
                              src/scarabee/_scarabee/python/nd_library.cpp
                              src/scarabee/_scarabee/python/flux_calculator.cpp
                              src/scarabee/_scarabee/python/cylindrical_cell.cpp
                              src/scarabee/_scarabee/python/cylindrical_flux_solver.cpp
                              src/scarabee/_scarabee/python/polar_quadrature.cpp
                              src/scarabee/_scarabee/python/boundary_condition.cpp
                              src/scarabee/_scarabee/python/simulation_mode.cpp
                              src/scarabee/_scarabee/python/track.cpp 
                              src/scarabee/_scarabee/python/cell.cpp
                              src/scarabee/_scarabee/python/empty_cell.cpp
                              src/scarabee/_scarabee/python/pin_cell_type.cpp
                              src/scarabee/_scarabee/python/simple_pin_cell.cpp
                              src/scarabee/_scarabee/python/pin_cell.cpp
                              src/scarabee/_scarabee/python/cartesian_2d.cpp
                              src/scarabee/_scarabee/python/cmfd.cpp
                              src/scarabee/_scarabee/python/moc_driver.cpp
                              src/scarabee/_scarabee/python/criticality_spectrum.cpp
                              src/scarabee/_scarabee/python/diffusion_data.cpp
                              src/scarabee/_scarabee/python/diffusion_geometry.cpp
                              src/scarabee/_scarabee/python/fd_diffusion_driver.cpp
                              src/scarabee/_scarabee/python/nem_diffusion_driver.cpp
                              src/scarabee/_scarabee/python/pins.cpp
                              src/scarabee/_scarabee/python/pwr_assembly.cpp
                              src/scarabee/_scarabee/python/reflector_sn.cpp
                              src/scarabee/_scarabee/python/water.cpp
                            )

target_include_directories(_scarabee PRIVATE include)
target_compile_features(_scarabee PRIVATE cxx_std_20)
target_link_libraries(_scarabee PUBLIC xtl xtensor xtensor-python htl HighFive hdf5-static Eigen3::Eigen spdlog::spdlog ImApp::ImApp cereal::cereal)
target_include_directories(_scarabee PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/scarabee/_scarabee/include>
  $<BUILD_INTERFACE:${NUMPY_INCLUDE_DIRS}>
)
target_compile_definitions(_scarabee PRIVATE SCARABEE_MAJOR_VERSION=${PROJECT_VERSION_MAJOR})
target_compile_definitions(_scarabee PRIVATE SCARABEE_MINOR_VERSION=${PROJECT_VERSION_MINOR})
target_compile_definitions(_scarabee PRIVATE SCARABEE_PATCH_VERSION=${PROJECT_VERSION_PATCH})

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Comile options for Windows
  target_compile_options(_scarabee PRIVATE /W3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
  target_compile_options(_scarabee PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
  target_compile_options(_scarabee PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
  target_compile_options(_scarabee PRIVATE -W -Wall -Wextra -Wconversion -Wpedantic)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Comile options for Windows
    target_compile_options(_scarabee PRIVATE /O2)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # Compile options for GCC
    target_compile_options(_scarabee PRIVATE -O3)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") # Compile options for Clang
    target_compile_options(_scarabee PRIVATE -O3)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") # Compile options for Intel
    target_compile_options(_scarabee PRIVATE -O3)
  endif()
endif()

# Find OpenMP if desired
if(SCARABEE_USE_OMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(_scarabee PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(_scarabee PUBLIC SCARABEE_USE_OMP)
  endif()
endif()

if (SKBUILD_PROJECT_NAME)
  # Generate stub file for type completion
  add_custom_command(TARGET _scarabee POST_BUILD
      COMMAND stubgen --include-docstrings -m _scarabee -o ${SKBUILD_PLATLIB_DIR}/${SKBUILD_PROJECT_NAME}
      WORKING_DIRECTORY $<TARGET_FILE_DIR:_scarabee> 
      COMMENT "Generating Python stub file for _scarabee"
  )

  # The install directory is the output (wheel) directory
  install(TARGETS _scarabee LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME} ) 
endif()
