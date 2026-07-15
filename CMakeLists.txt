# CMakeLists.txt for MARLEY neutrino event generator
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Compute "source revision" (runtime version string) from .VERSION or Git
# ------------------------------------------------------------------------------
if(EXISTS "${CMAKE_SOURCE_DIR}/.VERSION")
  file(STRINGS "${CMAKE_SOURCE_DIR}/.VERSION" MARLEY_COMPUTED_VERSION)
  set(MARLEY_COMPUTED_GIT_REVISION "${MARLEY_COMPUTED_VERSION}")
  set(VERSION_PREFIX "v")
else()
  find_package(Git QUIET)
  if(Git_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_REVISION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    execute_process(
      COMMAND ${GIT_EXECUTABLE} diff-index --quiet HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      RESULT_VARIABLE GIT_DIRTY
      ERROR_QUIET
    )
    if(NOT GIT_DIRTY EQUAL 0)
      set(GIT_REVISION "${GIT_REVISION}-dirty")
    endif()
  else()
    set(GIT_REVISION "unknown version")
  endif()
  set(MARLEY_COMPUTED_VERSION "${GIT_REVISION}")
  set(MARLEY_COMPUTED_GIT_REVISION "${GIT_REVISION}")
  set(VERSION_PREFIX "")
endif()
if(NOT DEFINED MARLEY_COMPUTED_GIT_REVISION)
  set(MARLEY_COMPUTED_GIT_REVISION "${MARLEY_COMPUTED_VERSION}")
endif()

# ------------------------------------------------------------------------------
# Semantic project version (for packaging + find_package version checks)
# ------------------------------------------------------------------------------
project(MARLEY VERSION 1.2.0 LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Build configuration / warnings
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_definitions(MARLEY_COMPILED_LOG_LEVEL=TRACE)
else()
  add_compile_definitions(MARLEY_COMPILED_LOG_LEVEL=INFO)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Embed runtime version info
add_compile_definitions(
  MARLEY_VERSION="${MARLEY_COMPUTED_VERSION}"
  MARLEY_GIT_REVISION="${MARLEY_COMPUTED_GIT_REVISION}"
)

# Stricter warnings by default
add_compile_options(-Wall -Wextra -Wpedantic -Wno-error=unused-parameter -Wcast-align)

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(STATUS "Compiling with Clang")
  add_compile_options(-Wno-keyword-macro -Wno-missing-braces)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  message(STATUS "Compiling with GCC ${CMAKE_CXX_COMPILER_VERSION}")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    add_compile_options(-Wno-shadow)
  endif()
  if(UNIX AND NOT APPLE)
    add_link_options(-Wl,--no-as-needed)
  endif()
endif()

# ------------------------------------------------------------------------------
# Discovery modes (auto | cmake | config) and toggles
# ------------------------------------------------------------------------------
set(ROOT_FIND_MODE   "auto" CACHE STRING "How to find ROOT: auto\ncmake\nconfig")
set(GSL_FIND_MODE    "auto" CACHE STRING "How to find GSL: auto\ncmake\nconfig")
set(HEPMC3_FIND_MODE "auto" CACHE STRING "How to find HepMC3: auto\ncmake\nconfig")
set_property(CACHE ROOT_FIND_MODE   PROPERTY STRINGS auto cmake config)
set_property(CACHE GSL_FIND_MODE    PROPERTY STRINGS auto cmake config)
set_property(CACHE HEPMC3_FIND_MODE PROPERTY STRINGS auto cmake config)

option(IGNORE_ROOT   "Build without ROOT support even if available" OFF)
option(IGNORE_HEPMC3 "Build without external HepMC3 even if available" OFF)
option(IGNORE_GSL    "Use the built-in GSL subset even if a system GSL is available" OFF)

# ------------------------------------------------------------------------------
# Helper: query config-based tools and return flags
# ------------------------------------------------------------------------------
function(find_library_with_config lib_name config_exe required_var cflags_var libs_var)
  find_program(${lib_name}_CONFIG_EXECUTABLE ${config_exe})
  if(${lib_name}_CONFIG_EXECUTABLE)
    execute_process(COMMAND ${${lib_name}_CONFIG_EXECUTABLE} --cflags
                    OUTPUT_VARIABLE ${lib_name}_CFLAGS
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${${lib_name}_CONFIG_EXECUTABLE} --libs
                    OUTPUT_VARIABLE ${lib_name}_LIBS
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${${lib_name}_CONFIG_EXECUTABLE} --version
                    OUTPUT_VARIABLE ${lib_name}_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    separate_arguments(${lib_name}_CFLAGS)
    separate_arguments(${lib_name}_LIBS)
    set(${required_var} TRUE PARENT_SCOPE)
    set(${cflags_var} ${${lib_name}_CFLAGS} PARENT_SCOPE)
    set(${libs_var}   ${${lib_name}_LIBS}   PARENT_SCOPE)
    if(${lib_name}_VERSION)
      message(STATUS "Found ${lib_name} ${${lib_name}_VERSION} via ${config_exe}")
    else()
      message(STATUS "Found ${lib_name} via ${config_exe}")
    endif()
  else()
    set(${required_var} FALSE PARENT_SCOPE)
  endif()
endfunction()

# ------------------------------------------------------------------------------
# Helper: prefer imported targets, else apply flags
# ------------------------------------------------------------------------------
# prefer_imported_or_flags(<target> <scope> "<imported_targets_list>" <cflags_var> <libs_var>)
function(prefer_imported_or_flags tgt scope imported_list cflags_var libs_var)
  set(_linked FALSE)
  if(imported_list)
    foreach(_it IN LISTS imported_list)
      if(TARGET ${_it})
        target_link_libraries(${tgt} ${scope} ${_it})
        set(_linked TRUE)
      endif()
    endforeach()
  endif()
  if(NOT _linked)
    if(DEFINED ${cflags_var} AND NOT "${${cflags_var}}" STREQUAL "")
      target_compile_options(${tgt} ${scope} ${${cflags_var}})
    endif()
    if(DEFINED ${libs_var} AND NOT "${${libs_var}}" STREQUAL "")
      target_link_libraries(${tgt} ${scope} ${${libs_var}})
    endif()
  endif()
endfunction()

# ------------------------------------------------------------------------------
# GSL (required unless IGNORE_GSL=ON)
# ------------------------------------------------------------------------------
set(GSL_FOUND FALSE)
if(NOT IGNORE_GSL)
  if(GSL_FIND_MODE STREQUAL "cmake" OR GSL_FIND_MODE STREQUAL "auto")
    find_package(GSL QUIET)
    if(GSL_FOUND)
      message(STATUS "Found GSL via CMake")
      set(MARLEY_FOUND_GSL_DEFINE TRUE)
    endif()
  endif()
  if(NOT GSL_FOUND AND (GSL_FIND_MODE STREQUAL "config" OR GSL_FIND_MODE STREQUAL "auto"))
    find_library_with_config(GSL gsl-config GSL_FOUND GSL_CFLAGS GSL_LIBS)
    if(GSL_FOUND)
      set(MARLEY_FOUND_GSL_DEFINE TRUE)
    endif()
  endif()
endif()

if(NOT GSL_FOUND)
  # Built-in GSL fallback (IGNORE_GSL=ON or no system GSL detected)
  message(STATUS "Using built-in GSL subset (src/marley_gsl.cc)")
  # No external includes or link flags needed — the builtin sources
  # (src/marley_gsl.cc) are already included via the MARLEY_SOURCES glob.
  # Link against libm, and hide the built-in GSL symbols to avoid clashes
  # with any system GSL loaded by other libraries.
  set(GSL_CFLAGS "")
  set(GSL_LIBS   "-lm")
endif()

# ------------------------------------------------------------------------------
# HepMC3 (with built-in fallback)
# ------------------------------------------------------------------------------
set(USE_BUILTIN_HEPMC3 FALSE)
set(HEPMC3_FOUND FALSE)

if(NOT IGNORE_HEPMC3)
  if(HEPMC3_FIND_MODE STREQUAL "cmake" OR HEPMC3_FIND_MODE STREQUAL "auto")
    find_package(HepMC3 QUIET)
    if(HepMC3_FOUND)
      message(STATUS "Found HepMC3 via CMake config")
      set(HEPMC3_FOUND TRUE)
      get_target_property(HEPMC3_INCLUDE_DIR HepMC3::HepMC3 INTERFACE_INCLUDE_DIRECTORIES)
    endif()
  endif()

  if(NOT HEPMC3_FOUND AND (HEPMC3_FIND_MODE STREQUAL "config" OR HEPMC3_FIND_MODE STREQUAL "auto"))
    find_library_with_config(HepMC3 HepMC3-config HEPMC3_FOUND HEPMC3_CFLAGS HEPMC3_LIBS)
    if(HEPMC3_FOUND)
      execute_process(
        COMMAND ${HepMC3_CONFIG_EXECUTABLE} --includedir
        OUTPUT_VARIABLE HEPMC3_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      set(HEPMC3_CFLAGS ${HEPMC3_CFLAGS} -DMARLEY_FOUND_HEPMC3)
    endif()
  endif()
endif()

if(NOT HEPMC3_FOUND)
  message(STATUS "Using built-in HepMC3 library (src/marley_hepmc3.cc)")
  set(USE_BUILTIN_HEPMC3 TRUE)
  # Headers now live under include/builtin/ (moved in copilot/refactor-builtin-hepmc3-library)
  set(HEPMC3_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/builtin)
  add_library(HepMC3_builtin SHARED ${CMAKE_SOURCE_DIR}/src/marley_hepmc3.cc)
  target_include_directories(HepMC3_builtin PUBLIC
    $<BUILD_INTERFACE:${HEPMC3_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
  )
  # Suppress deprecation warnings from sprintf usage in the amalgamated source
  target_compile_options(HepMC3_builtin PRIVATE -Wno-deprecated-declarations -Wno-strict-aliasing)
  prefer_imported_or_flags(HepMC3_builtin PUBLIC "GSL::gsl;GSL::gslcblas" GSL_CFLAGS GSL_LIBS)
  set_target_properties(HepMC3_builtin PROPERTIES
    OUTPUT_NAME HepMC3
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  )
  set(HEPMC3_CFLAGS -I${HEPMC3_INCLUDE_DIR})
  set(HEPMC3_LIBS   -L${CMAKE_BINARY_DIR}/lib -lHepMC3)
endif()

# ------------------------------------------------------------------------------
# ROOT (optional)
# ------------------------------------------------------------------------------
set(USE_ROOT FALSE)
if(NOT IGNORE_ROOT)
  set(ROOT_FOUND FALSE)
  if(ROOT_FIND_MODE STREQUAL "cmake" OR ROOT_FIND_MODE STREQUAL "auto")
    set(_cmake_message_log_level ${CMAKE_MESSAGE_LOG_LEVEL})
    set(CMAKE_MESSAGE_LOG_LEVEL ERROR)
    find_package(ROOT QUIET COMPONENTS Core RIO Hist Tree Graf Gpad)
    set(CMAKE_MESSAGE_LOG_LEVEL ${_cmake_message_log_level})
    if(ROOT_FOUND)
      set(ROOT_HAS_ALL_COMPONENTS TRUE)
      foreach(comp Core RIO Hist Tree Graf Gpad)
        if(NOT TARGET ROOT::${comp})
          set(ROOT_HAS_ALL_COMPONENTS FALSE)
          break()
        endif()
      endforeach()
      if(ROOT_HAS_ALL_COMPONENTS)
        set(USE_ROOT TRUE)
        message(STATUS "MARLEY will be built with ROOT support (CMake)")
        add_compile_definitions(USE_ROOT)
      else()
        set(ROOT_FOUND FALSE)
      endif()
    endif()
  endif()

  if(NOT USE_ROOT AND (ROOT_FIND_MODE STREQUAL "config" OR ROOT_FIND_MODE STREQUAL "auto"))
    find_library_with_config(ROOT root-config ROOT_FOUND ROOT_CFLAGS ROOT_LIBS)
    if(ROOT_FOUND)
      string(FIND "${ROOT_CFLAGS}" "c++23" HAS_CXX23)
      string(FIND "${ROOT_CFLAGS}" "c++2b" HAS_CXX2B)
      string(FIND "${ROOT_CFLAGS}" "c++20" HAS_CXX20)
      string(FIND "${ROOT_CFLAGS}" "c++2a" HAS_CXX2A)
      if(NOT HAS_CXX23 EQUAL -1 OR NOT HAS_CXX2B EQUAL -1)
        set(CMAKE_CXX_STANDARD 23)
        message(STATUS "Upgrading to C++23 to match ROOT")
      elseif(NOT HAS_CXX20 EQUAL -1 OR NOT HAS_CXX2A EQUAL -1)
        set(CMAKE_CXX_STANDARD 20)
        message(STATUS "Upgrading to C++20 to match ROOT")
      endif()
      find_program(ROOTCLING_EXECUTABLE NAMES rootcling rootcint)
      if(ROOTCLING_EXECUTABLE)
        set(USE_ROOT TRUE)
        message(STATUS "MARLEY will be built with ROOT support (root-config)")
        add_compile_definitions(USE_ROOT)
      else()
        message(WARNING "Found ROOT via root-config but not rootcling/rootcint; disabling ROOT support")
        set(ROOT_FOUND FALSE)
      endif()
    endif()
  endif()
endif()
if(NOT USE_ROOT)
  message(STATUS "MARLEY will be built without ROOT support")
endif()

# ------------------------------------------------------------------------------
# Main MARLEY library
# ------------------------------------------------------------------------------
file(GLOB MARLEY_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/app/*.cc)
list(FILTER MARLEY_SOURCES EXCLUDE REGEX "marley\\.cc$")
list(FILTER MARLEY_SOURCES EXCLUDE REGEX "marley_hepmc3\\.cc$")
if(NOT USE_ROOT)
  list(FILTER MARLEY_SOURCES EXCLUDE REGEX "marley_root\\.cc$")
  list(FILTER MARLEY_SOURCES EXCLUDE REGEX "OutputFileRoot\\.cc$")
endif()

add_library(marley_lib SHARED ${MARLEY_SOURCES})

if(MARLEY_FOUND_GSL_DEFINE)
  target_compile_definitions(marley_lib PRIVATE MARLEY_FOUND_GSL)
else()
  # When we're using the built-in GSL code, tell the compiler to hide the
  # symbols it provides from external users of the MARLEY shared library. This
  # is intended to prevent clashes with the system GSL library later. The extra
  # compile flag is only set for this source file when compiling the marley_lib
  # shared library target
  set_property(SOURCE ${CMAKE_SOURCE_DIR}/src/marley_gsl.cc TARGET marley_lib
    APPEND PROPERTY COMPILE_OPTIONS "-fvisibility=hidden")
endif()

# Public includes and transitive deps
target_include_directories(marley_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${HEPMC3_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:include>
)

prefer_imported_or_flags(marley_lib PUBLIC "GSL::gsl;GSL::gslcblas" GSL_CFLAGS GSL_LIBS)

if(USE_BUILTIN_HEPMC3)
  target_link_libraries(marley_lib PUBLIC HepMC3_builtin)
else()
  target_compile_definitions(marley_lib PRIVATE MARLEY_FOUND_HEPMC3)
  prefer_imported_or_flags(marley_lib PUBLIC "HepMC3::HepMC3" HEPMC3_CFLAGS HEPMC3_LIBS)
endif()

if(USE_ROOT)
  prefer_imported_or_flags(marley_lib PUBLIC "ROOT::Core;ROOT::RIO;ROOT::Hist;ROOT::Tree" ROOT_CFLAGS ROOT_LIBS)
endif()

set_target_properties(marley_lib PROPERTIES
  OUTPUT_NAME MARLEY           # produces libMARLEY.so (matches marley-config --libs)
  EXPORT_NAME MARLEY           # for MARLEY::MARLEY namespace
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  PUBLIC_HEADER ""             # optional
)

# ------------------------------------------------------------------------------
# ROOT dictionary
# ------------------------------------------------------------------------------
if(USE_ROOT)
  set(DICT_HEADER ${CMAKE_SOURCE_DIR}/include/marley/marley_linkdef.hh)
  file(GLOB HEPMC3_DATA_HEADERS "${HEPMC3_INCLUDE_DIR}/HepMC3/Data/*.h")

  include(ROOTGenerateDictionary OPTIONAL)
  if(COMMAND ROOT_GENERATE_DICTIONARY)
    set(DICT_NAME G__marley_root_dict)
    ROOT_GENERATE_DICTIONARY(${DICT_NAME}
      ${HEPMC3_DATA_HEADERS}
      LINKDEF ${DICT_HEADER}
      OPTIONS -I${CMAKE_SOURCE_DIR}/include
              -I${HEPMC3_INCLUDE_DIR}
              -I${HEPMC3_INCLUDE_DIR}/HepMC3/Data
    )
    target_sources(marley_lib PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/${DICT_NAME}.cxx)
    target_include_directories(marley_lib PRIVATE ${HEPMC3_INCLUDE_DIR}/HepMC3/Data)
    add_custom_command(TARGET marley_lib POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/lib
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict_rdict.pcm
        ${CMAKE_BINARY_DIR}/lib/libmarley_root_dict_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict.rootmap
        ${CMAKE_BINARY_DIR}/lib/libmarley_root_dict.rootmap
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict.rootmap
    )
  else()
    find_program(ROOTCLING_EXECUTABLE NAMES rootcling rootcint)
    if(NOT ROOTCLING_EXECUTABLE)
      message(FATAL_ERROR "ROOT dictionary helper not found and rootcling is missing.")
    endif()
    set(DICT_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc)
    add_custom_command(
      OUTPUT ${DICT_OUTPUT}
      COMMAND ${ROOTCLING_EXECUTABLE} -f ${DICT_OUTPUT} -c
        -I${CMAKE_SOURCE_DIR}/include
        -I${HEPMC3_INCLUDE_DIR}
        -I${HEPMC3_INCLUDE_DIR}/HepMC3/Data
        ${HEPMC3_DATA_HEADERS}
        ${DICT_HEADER}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/lib
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict_rdict.pcm
        ${CMAKE_BINARY_DIR}/lib/libmarley_root_dict_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict.rootmap
        ${CMAKE_BINARY_DIR}/lib/libmarley_root_dict.rootmap
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libmarley_root_dict.rootmap
      DEPENDS ${DICT_HEADER}
      COMMENT "Generating ROOT dictionary"
    )
    target_sources(marley_lib PRIVATE ${DICT_OUTPUT})
    target_include_directories(marley_lib PRIVATE ${HEPMC3_INCLUDE_DIR}/HepMC3/Data)
  endif()
endif()

# ------------------------------------------------------------------------------
# Executables (keep light; depend on MARLEY which carries transitive deps)
# ------------------------------------------------------------------------------
add_executable(marley ${CMAKE_SOURCE_DIR}/src/marley.cc)
target_link_libraries(marley PRIVATE marley_lib)
if(USE_ROOT)
  if(TARGET ROOT::Graf AND TARGET ROOT::Gpad)
    target_link_libraries(marley PRIVATE ROOT::Graf ROOT::Gpad)
  else()
    if(DEFINED ROOT_CFLAGS AND NOT "${ROOT_CFLAGS}" STREQUAL "")
      target_compile_options(marley PRIVATE ${ROOT_CFLAGS})
    endif()
    if(DEFINED ROOT_LIBS AND NOT "${ROOT_LIBS}" STREQUAL "")
      target_link_libraries(marley PRIVATE ${ROOT_LIBS} -lGraf -lGpad)
    endif()
  endif()
endif()
set_target_properties(marley PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
  INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_BINARY_DIR}/lib"
  BUILD_RPATH   "${CMAKE_BINARY_DIR}/lib"
)

if(USE_ROOT)
  configure_file(${CMAKE_SOURCE_DIR}/src/scripts/mroot ${CMAKE_BINARY_DIR}/bin/mroot COPYONLY)
  file(CHMOD ${CMAKE_BINARY_DIR}/bin/mroot
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE
  )
endif()

# ------------------------------------------------------------------------------
# marley-config script generation (classic Make-style CLI)
# ------------------------------------------------------------------------------
set(MARLEY_CONFIG_INPUT  ${CMAKE_SOURCE_DIR}/src/scripts/marley-config.in)
set(MARLEY_CONFIG_OUTPUT ${CMAKE_BINARY_DIR}/bin/marley-config)

set(CONFIG_VERSION        "${MARLEY_COMPUTED_VERSION}")
set(CONFIG_GIT_REVISION   "${MARLEY_COMPUTED_GIT_REVISION}")
set(CONFIG_CXX_STD        "c++${CMAKE_CXX_STANDARD}")
set(CONFIG_HEPMC3_CFLAGS  "${HEPMC3_CFLAGS}")
string(REPLACE ";" " " CONFIG_HEPMC3_CFLAGS "${CONFIG_HEPMC3_CFLAGS}")
set(CONFIG_HEPMC3_LIBS    "${HEPMC3_LIBS}")
string(REPLACE ";" " " CONFIG_HEPMC3_LIBS   "${CONFIG_HEPMC3_LIBS}")
set(CONFIG_GSL_CFLAGS     "${GSL_CFLAGS}")
string(REPLACE ";" " " CONFIG_GSL_CFLAGS    "${CONFIG_GSL_CFLAGS}")
set(CONFIG_GSL_LIBS       "${GSL_LIBS}")
string(REPLACE ";" " " CONFIG_GSL_LIBS      "${CONFIG_GSL_LIBS}")
if(USE_ROOT)
  set(CONFIG_USE_ROOT "yes")
else()
  set(CONFIG_USE_ROOT "no")
endif()

file(READ ${MARLEY_CONFIG_INPUT} MARLEY_CONFIG_CONTENT)
string(REPLACE "@@VERSION@@"       "\"${CONFIG_VERSION}\""       MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@GIT_REVISION@@"  "\"${CONFIG_GIT_REVISION}\""  MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@CXX_STD@@"       "\"${CONFIG_CXX_STD}\""       MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@HEPMC3_CFLAGS@@" "\"${CONFIG_HEPMC3_CFLAGS}\"" MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@HEPMC3_LIBS@@"   "\"${CONFIG_HEPMC3_LIBS}\""   MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@GSL_CFLAGS@@"    "\"${CONFIG_GSL_CFLAGS}\""    MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@GSL_LIBS@@"      "\"${CONFIG_GSL_LIBS}\""      MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REPLACE "@@USE_ROOT@@"      "\"${CONFIG_USE_ROOT}\""      MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")
string(REGEX REPLACE "\n##[^\n]*" "" MARLEY_CONFIG_CONTENT "${MARLEY_CONFIG_CONTENT}")

file(WRITE ${MARLEY_CONFIG_OUTPUT} "${MARLEY_CONFIG_CONTENT}")
file(CHMOD ${MARLEY_CONFIG_OUTPUT}
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
              GROUP_READ GROUP_EXECUTE
              WORLD_READ WORLD_EXECUTE
)

# ------------------------------------------------------------------------------
# Optional examples (requires Geant4)
# ------------------------------------------------------------------------------
find_package(Geant4 QUIET)
if(Geant4_FOUND)
  add_subdirectory(examples/marg4 EXCLUDE_FROM_ALL)
endif()

# ------------------------------------------------------------------------------
# Package export (build-tree + install-tree)
# ------------------------------------------------------------------------------
include(CMakePackageConfigHelpers)

# Build-tree package config directory
set(MARLEY_BUILD_CONFIG_DIR "${CMAKE_BINARY_DIR}/lib/cmake/MARLEY")
file(MAKE_DIRECTORY "${MARLEY_BUILD_CONFIG_DIR}")

# Collect targets to export from the build tree
set(_marley_export_targets marley_lib)
if(USE_BUILTIN_HEPMC3)
  list(APPEND _marley_export_targets HepMC3_builtin)
endif()

# Export targets from the build tree
export(TARGETS ${_marley_export_targets}
  FILE "${MARLEY_BUILD_CONFIG_DIR}/MARLEYTargets.cmake"
  NAMESPACE MARLEY::
)

# Build-tree package config files
configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/cmake/MARLEYConfig.cmake.in"
  "${MARLEY_BUILD_CONFIG_DIR}/MARLEYConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/MARLEY"
)
write_basic_package_version_file(
  "${MARLEY_BUILD_CONFIG_DIR}/MARLEYConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)

# Install-tree export
install(TARGETS marley_lib
  EXPORT MARLEYTargets
  LIBRARY DESTINATION lib
)
install(TARGETS marley
  RUNTIME DESTINATION bin
)
if(USE_ROOT)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/mroot
    DESTINATION bin
  )
  install(FILES
    ${CMAKE_BINARY_DIR}/lib/libmarley_root_dict_rdict.pcm
    ${CMAKE_BINARY_DIR}/lib/libmarley_root_dict.rootmap
    DESTINATION lib
    OPTIONAL
  )
endif()

# If built-in HepMC3 is used, install & export it (and its headers) for a clean install package
if(USE_BUILTIN_HEPMC3)
  install(TARGETS HepMC3_builtin
    EXPORT MARLEYTargets
    LIBRARY DESTINATION lib
  )
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/builtin/HepMC3
    DESTINATION include
  )
endif()

install(PROGRAMS ${MARLEY_CONFIG_OUTPUT}
  DESTINATION bin
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.hh"
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/
  DESTINATION share/marley/data
)

# Install the export set and package config files
install(EXPORT MARLEYTargets
  NAMESPACE MARLEY::
  DESTINATION lib/cmake/MARLEY
)
configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/cmake/MARLEYConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/MARLEYConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/MARLEY"
)
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/MARLEYConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(FILES
  "${CMAKE_BINARY_DIR}/MARLEYConfig.cmake"
  "${CMAKE_BINARY_DIR}/MARLEYConfigVersion.cmake"
  DESTINATION "lib/cmake/MARLEY"
)

# ------------------------------------------------------------------------------
# Documentation targets (optional)
# ------------------------------------------------------------------------------
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  add_custom_target(doxygen
    COMMAND ${CMAKE_COMMAND} -E env
      MARLEY_VERSION=${VERSION_PREFIX}${MARLEY_COMPUTED_VERSION}
      ${CMAKE_MAKE_PROGRAM} doxygen
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
    COMMENT "Generating Doxygen documentation"
  )
endif()
find_program(SPHINX_BUILD sphinx-build)
if(SPHINX_BUILD)
  add_custom_target(docs
    COMMAND ${CMAKE_COMMAND} -E env
      MARLEY_VERSION=${VERSION_PREFIX}${MARLEY_COMPUTED_VERSION}
      ${CMAKE_MAKE_PROGRAM} html
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/docs
    COMMENT "Generating Sphinx documentation"
  )
endif()

# ------------------------------------------------------------------------------
# Uninstall target
# Reads the install_manifest.txt generated by cmake --install and removes
# every listed file.  This is the standard CMake pattern.
# ------------------------------------------------------------------------------
configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
  @ONLY
)
add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_BINARY_DIR}/cmake_uninstall.cmake"
  COMMENT "Uninstalling MARLEY (reading install_manifest.txt)"
)

# ------------------------------------------------------------------------------
# CTest integration
# Builds the Catch2-based test executable (martest) and registers it with CTest.
# Run with: ctest --test-dir build
# The MARLEY env-var must point at the source root so that FileManager can find
# data files (same requirement as the GNU Make 'make test' path).
# ------------------------------------------------------------------------------
enable_testing()

file(GLOB MARLEY_TEST_SOURCES
  "${CMAKE_SOURCE_DIR}/src/tests/*.cc"
)
if(NOT GSL_FOUND)
  # If the system GSL library was not found (or it was ignored), then we need
  # to include the built-in GSL code when building the test executable
  list(APPEND MARLEY_TEST_SOURCES "${CMAKE_SOURCE_DIR}/src/marley_gsl.cc")
endif()

add_executable(martest EXCLUDE_FROM_ALL ${MARLEY_TEST_SOURCES})
target_link_libraries(martest PRIVATE marley_lib)
prefer_imported_or_flags(martest PRIVATE "${GSL_IMPORTED_TARGETS}"
  GSL_CFLAGS GSL_LIBS)
target_include_directories(martest PRIVATE
  "${CMAKE_SOURCE_DIR}/include"
)
set_target_properties(martest PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
# Inherit ROOT compiler and linker settings when building with ROOT support
if(USE_ROOT)
  if(TARGET ROOT::Graf AND TARGET ROOT::Gpad)
    target_link_libraries(martest PRIVATE ROOT::Graf ROOT::Gpad)
  else()
    if(DEFINED ROOT_CFLAGS AND NOT "${ROOT_CFLAGS}" STREQUAL "")
      target_compile_options(martest PRIVATE ${ROOT_CFLAGS})
    endif()
    if(DEFINED ROOT_LIBS AND NOT "${ROOT_LIBS}" STREQUAL "")
      target_link_libraries(martest PRIVATE ${ROOT_LIBS} -lGraf -lGpad)
    endif()
  endif()
endif()

add_test(
  NAME martest
  COMMAND martest
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
# Point the test at the source root so FileManager finds data/ and data/tests/
set_tests_properties(martest PROPERTIES
  ENVIRONMENT "MARLEY=${CMAKE_SOURCE_DIR}"
)
