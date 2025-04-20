cmake_minimum_required(VERSION 3.10)
project(WaterSimulation VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(include)

set(SOURCES
    src/SpatialGrid.cpp
    src/SPHKernels.cpp
    src/Solver.cpp
    src/Renderer.cpp
    src/Simulation.cpp
)

add_library(WaterSim STATIC ${SOURCES})

install(TARGETS WaterSim DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)

option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

message(STATUS "Configuration Summary:")
message(STATUS "  C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install Prefix: ${CMAKE_INSTALL_PREFIX}")
