file( GLOB LIBRARY_SOURCE_FILES "*.cpp" )
add_library( ${PACKAGE_NAME} SHARED ${LIBRARY_SOURCE_FILES} )
target_include_directories( ${PACKAGE_NAME} PUBLIC ${XYZ_SOURCE_DIR}/include )

