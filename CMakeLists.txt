cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_DISABLE_SOURCE_CHANGES on)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

project(marley)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories(${INCLUDE_DIR})

file(GLOB SOURCES "src/[A-Z]*.cc" "src/marley_*.cc" "src/meta_numerics.cc")
list(REMOVE_ITEM SOURCES "${SOURCE_DIR}/RootConfigFile.cc")

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

add_executable(marley "src/marley.cc")

if(ROOT AND ROOTCINT AND ROOTCONFIG)
  # Create a shared library to allow the ROOT dictionary to be easily loaded on demand
  add_library(MARLEY_ROOT SHARED ${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc
    ${SOURCE_DIR}/RootConfigFile.cc)
  target_link_libraries(MARLEY_ROOT MARLEY)
  set_target_properties(MARLEY_ROOT PROPERTIES LINK_FLAGS ${ROOT_LDFLAGS})

  # Create the example executable using ROOT
  target_link_libraries(marley MARLEY_ROOT)
  set_target_properties(marley PROPERTIES LINK_FLAGS ${ROOT_LDFLAGS})
else()
  # Create the example executable without using ROOT
  target_link_libraries(marley MARLEY)
endif()

# ROOT dictionaries generated with rootcint trigger clang++'s keyword-macro
# warning
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS} -O3 \
  -Wall -Wextra -Wshadow -pedantic \
  -Wcast-align -Wunused -Woverloaded-virtual \
  -Werror -Wno-keyword-macro")
