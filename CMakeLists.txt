cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_DISABLE_SOURCE_CHANGES on)

project(marley)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories(${INCLUDE_DIR})

file(GLOB SOURCES "src/[A-Z]*.cc" "src/marley_utils.cc" "src/meta_numerics.cc")

find_program(ROOTCINT rootcint PATHS ENV PATH)
find_program(ROOTCONFIG root-config PATHS ENV PATH)
find_program(ROOT root PATHS ENV PATH)

if(ROOT AND ROOTCINT AND ROOTCONFIG)
  exec_program(${ROOTCONFIG} ARGS --version OUTPUT_VARIABLE ROOT_VERSION)
  message("Found ROOT version ${ROOT_VERSION} in ${ROOT}")
  message("MARLEY will be built with ROOT support.")
  exec_program(${ROOTCONFIG} ARGS --cflags OUTPUT_VARIABLE ROOT_CFLAGS)
  exec_program(${ROOTCONFIG} ARGS --ldflags --libs OUTPUT_VARIABLE ROOT_LDFLAGS)

  set(ROOT_DICT_INCLUDES ${INCLUDE_DIR}/Particle.hh
    ${INCLUDE_DIR}/Event.hh ${INCLUDE_DIR}/marley_linkdef.hh)
  
  # Create a ROOT dictionary target for our analysis classes based on their latest header files
  add_custom_command(OUTPUT marley_root_dict.cc marley_root_dict.h
                      COMMAND rm -f marley_root_dict.*
                      COMMAND echo "Building MARLEY ROOT dictionaries..."
                      COMMAND ${ROOTCINT} marley_root_dict.cc -c ${ROOT_DICT_INCLUDES})
  add_custom_target(root_dictionaries DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc)

  # Add the dictionary to the source file list
  set(SOURCES ${SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc)

  # Tell the MARLEY sources to use optional ROOT-dependent features
  add_definitions(-DUSE_ROOT)
else()
  message("WARNING: Could not find a valid ROOT installation. MARLEY will be built without ROOT support.")
endif()

# Build shared and static libraries with the same name (other than the suffix)
# on systems other than Windows (doing this on Windows causes problems).
add_library(MARLEY SHARED ${SOURCES})
add_library(MARLEY_static STATIC ${SOURCES})
if (NOT WIN32)
  set_target_properties(MARLEY_static PROPERTIES OUTPUT_NAME MARLEY)
endif()

if(ROOT AND ROOTCINT AND ROOTCONFIG)
  # Now that we've created the library targets, add the dictionary to their
  # dependencies
  add_dependencies(MARLEY root_dictionaries)
  add_dependencies(MARLEY_static root_dictionaries)
endif()

add_executable(marley "src/marley.cc")
  #${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc)
target_link_libraries(marley MARLEY)#_static)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 \
  -std=c++14 -Wall -Wextra -Wpedantic -Werror \
  ${ROOT_CFLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ROOT_LDFLAGS}")
