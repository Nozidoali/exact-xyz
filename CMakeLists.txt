cmake_minimum_required(VERSION 3.20)

project(xyz
    LANGUAGES CXX
    VERSION 0.1.0
    DESCRIPTION "quantum computing library"
)

# Options
option(BUILD_XYZ_TESTS "Build the XYZ tests" ON)

# Global settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Interface libraries for headers
add_library(xyz_headers INTERFACE)
target_include_directories(xyz_headers INTERFACE include)
target_compile_features(xyz_headers INTERFACE cxx_std_17)

add_library(xyz_third_party INTERFACE)
target_include_directories(xyz_third_party INTERFACE thrid-party)

# Main library
file(GLOB LIB_SOURCES CONFIGURE_DEPENDS "lib/*.cpp")
add_library(cxyz SHARED ${LIB_SOURCES})
target_link_libraries(cxyz PUBLIC xyz_headers)

# Tools executables
file(GLOB TOOL_SOURCES CONFIGURE_DEPENDS "tools/*.cpp")
foreach(TOOL_SOURCE ${TOOL_SOURCES})
    get_filename_component(TOOL_NAME ${TOOL_SOURCE} NAME_WE)
    add_executable(${TOOL_NAME} ${TOOL_SOURCE})
    target_link_libraries(${TOOL_NAME} PRIVATE cxyz xyz_third_party)
endforeach()

# Tests
if(BUILD_XYZ_TESTS)
    add_subdirectory(tests)
endif()