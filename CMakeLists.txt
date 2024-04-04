cmake_minimum_required(VERSION 3.11)
project(scarabee CXX)

option(SCARABEE_USE_OMP "Compile Scarabée with OpenMP for shared memory parallelism" OFF)

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
add_executable(scarabee src/main.cpp
                        src/logging.cpp
                        src/gauss_kronrod.cpp
                        src/chebyshev.cpp
                        src/bickley.cpp
                        src/transport_xs.cpp
                        src/cylindrical_cell.cpp
                        src/cylindrical_flux_solver.cpp
                        src/surface.cpp
                        src/cell.cpp
                        src/pin_cell.cpp
                        src/cartesian_2d.cpp
                        src/track.cpp
                        src/legendre.cpp
                        src/yamamoto_tabuchi.cpp
                        src/moc_driver.cpp)
set_target_properties(scarabee PROPERTIES OUTPUT_NAME "scarabée")
target_include_directories(scarabee PRIVATE include)
target_compile_features(scarabee PRIVATE cxx_std_20)
target_link_libraries(scarabee PUBLIC xtl xtensor htl Eigen3::Eigen spdlog::spdlog)

# Find OpenMP if desired
if(SCARABEE_USE_OMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(scarabee PUBLIC OpenMP::OpenMP_CXX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") # Compile options for Windows
      target_compile_options(scarabee PRIVATE /openmp:llvm) 
    endif()
    target_compile_definitions(scarabee PUBLIC SCARABEE_USE_OMP)
  endif()
endif()

#===============================================================================
# Allow for install of executable
install(TARGETS scarabee
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
