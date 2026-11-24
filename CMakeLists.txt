cmake_minimum_required(VERSION 3.16)
project(sobol_burley_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Library target (header-only; the target just declares the include dir).
add_library(sobol_burley INTERFACE)
target_include_directories(sobol_burley INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Example: reproduces examples/plots.rs -> writes PBM scatter plots.
add_executable(plots plots.cpp)
target_link_libraries(plots PRIVATE sobol_burley)

# Self-verification tests.
add_executable(tests tests.cpp)
target_link_libraries(tests PRIVATE sobol_burley)

# Generator: emits sobol_4096spp_256d.h.
add_executable(gen_table gen_sobol_table.cpp)
target_link_libraries(gen_table PRIVATE sobol_burley)

# Verifier: checks the generated table against the live library.
add_executable(verify_table verify_table.cpp)
target_link_libraries(verify_table PRIVATE sobol_burley)

enable_testing()
add_test(NAME sobol_burley_tests COMMAND tests)
