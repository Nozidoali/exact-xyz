file(GLOB CPP_FILES CONFIGURE_DEPENDS "*.cpp")

foreach(CPP_FILE ${CPP_FILES})
  get_filename_component(EXECUTABLE_NAME ${CPP_FILE} NAME_WE)
  add_executable(${EXECUTABLE_NAME} ${CPP_FILE})
  target_link_libraries(${EXECUTABLE_NAME} cxyz xyz_third_party)
endforeach()
