file( GLOB LIBRARY_SOURCE_FILES "*.cpp" )

add_library( cxyz SHARED ${LIBRARY_SOURCE_FILES} )
add_subdirectory( qgate )

target_include_directories( cxyz PUBLIC ${XYZ_SOURCE_DIR} / include )

find_package(fmt REQUIRED)
target_link_libraries(cxyz fmt::fmt)