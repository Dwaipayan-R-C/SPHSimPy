cmake_minimum_required(VERSION 3.11)

project(sphproject)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

include(tools)

add_external_package(Eigen3 VERSION 3.4.0 CONFIG)
set(SPHProject_HDRS
        ${CMAKE_CURRENT_SOURCE_DIR}/headers/atom_structure.h
        ${CMAKE_CURRENT_SOURCE_DIR}/headers/simulation.h
        ${CMAKE_CURRENT_SOURCE_DIR}/headers/xyz.h
        ${CMAKE_CURRENT_SOURCE_DIR}/headers/neighbors.h
        ${CMAKE_CURRENT_SOURCE_DIR}/headers/smoothing_kernel.h
        ${CMAKE_CURRENT_SOURCE_DIR}/headers/force_utils.h
        )

set(SPHProject_SRCS
        # Add source files here
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/xyz.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/simulation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/neighbors.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/smoothing_kernel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/srcs/force_utils.cpp

)

set(SPHProject_executable_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
       )

add_subdirectory(tests)
# Simulation code
add_executable(sphproject ${SPHProject_SRCS} ${SPHProject_HDRS} ${SPHProject_executable_SRCS})
target_link_libraries(sphproject PUBLIC Eigen3::Eigen )
