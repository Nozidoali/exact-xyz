# set required cmake version
cmake_minimum_required(VERSION 3.20)

project(
    xyz
    LANGUAGES CXX
    VERSION 0.1.0
    DESCRIPTION "quantum computing library"
)

option(BUILD_XYZ_TESTS "Build the XYZ tests" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(XYZ_SOURCE_DIR ${CMAKE_SOURCE_DIR})

set(PACKAGE_NAME cxyz)

add_library(xyz_headers INTERFACE)
target_include_directories(xyz_headers INTERFACE ${XYZ_SOURCE_DIR}/include)
target_compile_features(xyz_headers INTERFACE cxx_std_17)

add_library(xyz_third_party INTERFACE)
target_include_directories(xyz_third_party INTERFACE ${XYZ_SOURCE_DIR}/thrid-party)

add_subdirectory(${XYZ_SOURCE_DIR}/lib)
add_subdirectory(${XYZ_SOURCE_DIR}/tools)

if (BUILD_XYZ_TESTS)
    add_subdirectory(${XYZ_SOURCE_DIR}/tests)
endif()