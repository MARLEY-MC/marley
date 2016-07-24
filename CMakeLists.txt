cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_DISABLE_SOURCE_CHANGES on)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

project(marley)

set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories(${INCLUDE_DIR})

file(GLOB SOURCES "src/*.cc")
list(REMOVE_ITEM SOURCES "${SOURCE_DIR}/RootConfigurationFile.cc")
list(REMOVE_ITEM SOURCES "${SOURCE_DIR}/marley.cc")

find_program(ROOTCINT rootcint rootcling PATHS ENV PATH)
find_program(ROOTCONFIG root-config PATHS ENV PATH)
find_program(ROOT root PATHS ENV PATH)

if(ROOT AND ROOTCINT AND ROOTCONFIG)
  exec_program(${ROOTCONFIG} ARGS --version OUTPUT_VARIABLE ROOT_VERSION)
  message("Found ROOT version ${ROOT_VERSION} in ${ROOT}")
  message("MARLEY will be built with ROOT support.")
  exec_program(${ROOTCONFIG} ARGS --cflags OUTPUT_VARIABLE ROOT_CFLAGS)
  exec_program(${ROOTCONFIG} ARGS --ldflags OUTPUT_VARIABLE ROOT_LDFLAGS)
  exec_program(${ROOTCONFIG} ARGS --libdir OUTPUT_VARIABLE ROOT_LIBDIR)
  set(ROOT_LDFLAGS "${ROOT_LDFLAGS} -L${ROOT_LIBDIR} -lCore -lRIO \
    -lHist -lTree -rdynamic")

  set(ROOT_DICT_INCLUDES -I${INCLUDE_DIR} marley/Particle.hh
    marley/Event.hh marley/marley_linkdef.hh)

  # Create a ROOT dictionary target for our analysis classes based on their latest header files
  add_custom_command(OUTPUT marley_root_dict.cc marley_root_dict.h
                      COMMAND rm -f marley_root_dict.*
                      COMMAND echo "Building MARLEY ROOT dictionaries..."
                      COMMAND ${ROOTCINT} -f marley_root_dict.cc -c ${ROOT_DICT_INCLUDES})
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


# Create the example executable. Use the name "marley.exe" on OS X
# to avoid problems with the case-insensitive file system. Also use
# the name marley.exe on Windows to match the standard convention.
if (APPLE OR WIN32)
  set(EXECUTABLE_NAME "marley.exe")
else()
  set(EXECUTABLE_NAME "marley")
endif()

add_executable(${EXECUTABLE_NAME} "src/marley.cc")

if(ROOT AND ROOTCINT AND ROOTCONFIG)
  # Create a shared library to allow the ROOT dictionary to be easily loaded on demand
  add_library(MARLEY_ROOT SHARED ${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc
    ${SOURCE_DIR}/RootConfigurationFile.cc)
  target_link_libraries(MARLEY_ROOT MARLEY)
  set_target_properties(MARLEY_ROOT PROPERTIES LINK_FLAGS ${ROOT_LDFLAGS})

  # Create the example executable using ROOT.
  target_link_libraries(${EXECUTABLE_NAME} MARLEY_ROOT)
  set_target_properties(${EXECUTABLE_NAME} PROPERTIES LINK_FLAGS ${ROOT_LDFLAGS})
else()
  # Create the example executable without using ROOT
  target_link_libraries(${EXECUTABLE_NAME} MARLEY)
endif()

# ROOT dictionaries generated with rootcint trigger clang++'s keyword-macro
# warning, so disable it. Since g++ 4.9.3 gives many false positives for
# -Wshadow, also disable that warning.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CFLAGS} -O3 \
  -Wall -Wextra -pedantic \
  -Wcast-align -Wunused -Woverloaded-virtual \
  -Wno-keyword-macro -Wno-missing-braces -Wno-shadow ")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # Linking to ROOT libraries can be problematic on distributions (e.g., Ubuntu)
  # that set the g++ flag -Wl,--as-needed by default (see
  # http://www.bnikolic.co.uk/blog/gnu-ld-as-needed.html for details),
  # so disable this behavior
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
endif()

# Set installation directories for MARLEY components
get_property(LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

# Allow the user to set -DMARLEY_LIB_DIR_SUFFIX in case CMake gets confused
# about whether the libraries should be installed in lib, lib64, or lib32.
# This can happen when you're cross-compiling, for example.
if (DEFINED MARLEY_LIB_DIR_SUFFIX)
  set(LIBSUFFIX "${MARLEY_LIB_DIR_SUFFIX}")
elseif ("${LIB64}" STREQUAL "TRUE")
  set(LIBSUFFIX 64)
else()
  set(LIBSUFFIX "")
endif()

set(INSTALL_LIB_DIR lib${LIBSUFFIX} CACHE PATH "Installation
directory for libraries")
mark_as_advanced(INSTALL_LIB_DIR)

install(TARGETS marley MARLEY MARLEY_static
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

install(DIRECTORY include/ DESTINATION include/marley
  FILES_MATCHING PATTERN "*.hh")

install(DIRECTORY react/ DESTINATION share/marley/react
  FILES_MATCHING PATTERN "*.react")

install(FILES config.txt DESTINATION share/marley)
