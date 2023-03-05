include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/eaf9a3fd77869cf95befb87455a2e2a2e85044ff.zip
)
FetchContent_MakeAvailable(googletest)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(sphproject_tests_SRCS
        # ${CMAKE_CURRENT_SOURCE_DIR}/test_hello_world.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/neighbor_list_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/smoothing_kernel_test.cpp
        )

# Tests
include_directories(sphproject_tests ${CMAKE_SOURCE_DIR})
add_executable(sphproject_tests ${SPHProject_SRCS} ${SPHProject_HDRS} ${sphproject_tests_SRCS})
target_link_libraries(sphproject_tests PUBLIC Eigen3::Eigen gtest gtest_main)

include(GoogleTest)
gtest_discover_tests(sphproject_tests)

add_test(NAME sphproject_tests COMMAND sphproject_tests)
