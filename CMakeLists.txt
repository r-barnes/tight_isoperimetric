# Created by the script cgal_create_cmake_script (and then adapted manually).
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.9...3.15)
project( tight_isoperimetric )

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

add_subdirectory(submodules/graphlite)

find_package(CGAL REQUIRED)

if ( CGAL_FOUND )
  include(${CGAL_USE_FILE})

  message(STATUS "Found CGAL: ${CGAL_LIBRARIES}")

  add_executable(tightiso src/main.cpp)
  target_link_libraries(tightiso ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} graphlite)
  target_compile_options(tightiso PRIVATE -Wall -Wextra -pedantic -march=native)
  target_compile_features(tightiso PRIVATE cxx_std_20)
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS tightiso)

  # include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  # cgal_add_compilation_test(tightiso)

else()

  message(STATUS "NOTICE: This requires CGAL and will not be compiled.")

endif()
