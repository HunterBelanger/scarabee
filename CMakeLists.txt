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
# Get XTENSOR
message(STATUS "Downloading xtensor v0.25.0")
FetchContent_Declare(xtensor
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor.git
  GIT_TAG        0.25.0
)
FetchContent_MakeAvailable(xtensor)

#===============================================================================
add_executable(scarabee src/main.cpp
                        src/gauss_kronrod.cpp
                        src/chebyshev.cpp
                        src/bickley.cpp
                        src/transport_xs.cpp
                        src/cylindrical_cell.cpp
                        src/cylindrical_flux_solver.cpp
                        src/surface.cpp
                        src/flat_source_region.cpp
                        src/cell.cpp
                        src/pin_cell.cpp
                        src/cartesian_2d.cpp
                        src/track.cpp
                        src/moc_driver.cpp)
set_target_properties(scarabee PROPERTIES OUTPUT_NAME "scarabée")
target_include_directories(scarabee PRIVATE include)
target_compile_features(scarabee PRIVATE cxx_std_20)
target_link_libraries(scarabee PUBLIC xtl xtensor htl Eigen3::Eigen)

#===============================================================================
# Allow for install of executable
install(TARGETS scarabee
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
