# CMakeLists.txt for MARLEY neutrino event generator
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Compute version string (from .VERSION or Git) and git revision
# ------------------------------------------------------------------------------
if(EXISTS "${CMAKE_SOURCE_DIR}/.VERSION")
  file(STRINGS "${CMAKE_SOURCE_DIR}/.VERSION" MARLEY_COMPUTED_VERSION)
  set(VERSION_PREFIX "v")
else()
  set(VERSION_PREFIX "")
endif()

# Always try to get the actual git revision (used for --git-revision)
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

set(MARLEY_COMPUTED_GIT_REVISION "${GIT_REVISION}")

# Fall back to git revision if no .VERSION file
if(NOT DEFINED MARLEY_COMPUTED_VERSION)
  set(MARLEY_COMPUTED_VERSION "${GIT_REVISION}")
endif()

# ------------------------------------------------------------------------------
# Semantic project version (for packaging + find_package version checks)
# ------------------------------------------------------------------------------
project(MARLEY VERSION 2.0.0 LANGUAGES CXX)

# ------------------------------------------------------------------------------
# Install directories
#
# libdir, bindir and includedir must resolve to exactly ONE directory below
# the prefix, because the data root is found by a single ".." ascent from the
# library's own location (design spec 2.2, 4.3). The GNUInstallDirs names are
# reused because those are the names Spack, conda-forge and RPM recipes
# already pass; a MARLEY_* namespace would silently ignore them, which is the
# very defect this check exists to prevent.
#
# CMAKE_INSTALL_LIBDIR is genuinely ambiguous in the wild -- relative on some
# systems, absolute on others (Fedora's %cmake macro historically set it
# absolute, producing /usr/usr/lib64 in projects that concatenated it). Both
# spellings are accepted. The prefix comparison is done with LENGTH/SUBSTRING
# rather than a regex because the prefix would otherwise be interpolated into
# a pattern, and Spack prefixes routinely contain '.', '+' and '['.
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_LIBDIR     "lib"     CACHE STRING "Library directory name")
set(CMAKE_INSTALL_BINDIR     "bin"     CACHE STRING "Executable directory name")
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "Header directory name")

function(marley_relative_install_dir var_name out_var)
  set(_value "${${var_name}}")
  set(_prefix "${CMAKE_INSTALL_PREFIX}/")
  string(LENGTH "${_prefix}" _plen)
  string(LENGTH "${_value}" _vlen)
  if(_vlen GREATER_EQUAL _plen)
    string(SUBSTRING "${_value}" 0 ${_plen} _head)
    if(_head STREQUAL "${_prefix}")
      string(SUBSTRING "${_value}" ${_plen} -1 _value)
    endif()
  endif()
  if(_value STREQUAL "" OR _value MATCHES "/")
    message(FATAL_ERROR
      "${var_name}=\"${${var_name}}\" must be exactly one directory below "
      "CMAKE_INSTALL_PREFIX=\"${CMAKE_INSTALL_PREFIX}\" (resolved to "
      "\"${_value}\"); use e.g. -D${var_name}=lib64")
  endif()
  set(${out_var} "${_value}" PARENT_SCOPE)
endfunction()

marley_relative_install_dir(CMAKE_INSTALL_LIBDIR     MARLEY_LIBDIR_REL)
marley_relative_install_dir(CMAKE_INSTALL_BINDIR     MARLEY_BINDIR_REL)
marley_relative_install_dir(CMAKE_INSTALL_INCLUDEDIR MARLEY_INCLUDEDIR_REL)

# datadir is fixed: unlike the others it is used for DESCENT, so its name must
# be known at runtime rather than merely counted. Setting it is therefore an
# error rather than something to ignore silently.
foreach(_fixed CMAKE_INSTALL_DATADIR CMAKE_INSTALL_DATAROOTDIR)
  if(DEFINED ${_fixed} AND NOT "${${_fixed}}" STREQUAL "share")
    message(FATAL_ERROR
      "${_fixed}=\"${${_fixed}}\" is not supported: MARLEY's data directory "
      "name is compiled in and must remain \"share\".")
  endif()
endforeach()

# The build tree mirrors the install layout (design spec 2.2), so that the
# single hop taken by self-location works before installation.
set(MARLEY_BUILD_LIBDIR "${CMAKE_BINARY_DIR}/${MARLEY_LIBDIR_REL}")

# Token the loader expands to the containing object's own directory. Every
# linked artifact locates its dependencies relative to itself, so no absolute
# path is baked into anything and the install tree can be moved (spec 4.2).
if(APPLE)
  set(MARLEY_ORIGIN "@loader_path")
else()
  set(MARLEY_ORIGIN "$ORIGIN")
  # Ask for DT_RUNPATH rather than DT_RPATH. Which one -Wl,-rpath produces is
  # a binutils configure-time default, and the reason for preferring RUNPATH
  # -- that a user's LD_LIBRARY_PATH still wins -- depends on getting it.
  add_link_options(-Wl,--enable-new-dtags)
endif()

# ------------------------------------------------------------------------------
# Build configuration / warnings
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MARLEY_LOG_LEVEL TRACE)
else()
  set(MARLEY_LOG_LEVEL INFO)
endif()
add_compile_definitions(MARLEY_COMPILED_LOG_LEVEL=${MARLEY_LOG_LEVEL})

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
      get_target_property(_gsl_incdirs GSL::gsl INTERFACE_INCLUDE_DIRECTORIES)
      if(_gsl_incdirs)
        set(GSL_CFLAGS "")
        foreach(_d IN LISTS _gsl_incdirs)
          list(APPEND GSL_CFLAGS "-I${_d}")
        endforeach()
        list(APPEND GSL_CFLAGS "-DMARLEY_FOUND_GSL")
      endif()
      set(GSL_LIBS "")
      foreach(_t gsl gslcblas)
        get_target_property(_loc GSL::${_t} IMPORTED_LOCATION)
        if(NOT _loc)
          get_target_property(_loc GSL::${_t} IMPORTED_LOCATION_RELEASE)
        endif()
        if(_loc)
          get_filename_component(_libdir "${_loc}" DIRECTORY)
          get_filename_component(_libfile "${_loc}" NAME)
          string(REGEX REPLACE "^lib(.+)\\.(so|dylib|a)(\\..*)?$" "\\1" _libname "${_libfile}")
          if(NOT GSL_LIBS)
            list(APPEND GSL_LIBS "-L${_libdir}")
          endif()
          list(APPEND GSL_LIBS "-l${_libname}")
        endif()
      endforeach()
    endif()
  endif()
  if(NOT GSL_FOUND AND (GSL_FIND_MODE STREQUAL "config" OR GSL_FIND_MODE STREQUAL "auto"))
    find_library_with_config(GSL gsl-config GSL_FOUND GSL_CFLAGS GSL_LIBS)
    if(GSL_FOUND)
      set(MARLEY_FOUND_GSL_DEFINE TRUE)
      list(APPEND GSL_CFLAGS "-DMARLEY_FOUND_GSL")
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
      set(HEPMC3_CFLAGS "-I${HEPMC3_INCLUDE_DIR}")
      get_target_property(_hepMC3_loc HepMC3::HepMC3 IMPORTED_LOCATION)
      if(NOT _hepMC3_loc)
        get_target_property(_hepMC3_loc HepMC3::HepMC3 IMPORTED_LOCATION_RELEASE)
      endif()
      if(_hepMC3_loc)
        get_filename_component(_hepMC3_dir "${_hepMC3_loc}" DIRECTORY)
        get_filename_component(_hepMC3_file "${_hepMC3_loc}" NAME)
        string(REGEX REPLACE "^lib(.+)\\.(so|dylib)(\\..*)?$" "\\1" _hepMC3_name "${_hepMC3_file}")
        set(HEPMC3_LIBS "-L${_hepMC3_dir} -l${_hepMC3_name}")
      endif()
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
    # Must track the install destination below, or the exported target hands
    # consumers an include path with no HepMC3/ beneath it
    $<INSTALL_INTERFACE:${MARLEY_INCLUDEDIR_REL}/marley/builtin>
  )
  # Suppress deprecation warnings from sprintf usage in the amalgamated source
  target_compile_options(HepMC3_builtin PRIVATE -Wno-deprecated-declarations -Wno-strict-aliasing)
  prefer_imported_or_flags(HepMC3_builtin PUBLIC "GSL::gsl;GSL::gslcblas" GSL_CFLAGS GSL_LIBS)
  set_target_properties(HepMC3_builtin PROPERTIES
    # Renamed from "HepMC3": a libHepMC3 of ours in a shared prefix collides
    # with a real one, and the amalgamated copy does not implement all of it.
    # The exported *target* keeps the name HepMC3_builtin, which says why it
    # exists; nobody types it, and target names differing from file names is
    # ordinary in CMake.
    OUTPUT_NAME MARLEYHepMC3
    LIBRARY_OUTPUT_DIRECTORY ${MARLEY_BUILD_LIBDIR}
    INSTALL_RPATH "${MARLEY_ORIGIN}"
    BUILD_RPATH   "${MARLEY_ORIGIN}"
    MACOSX_RPATH ON
    INSTALL_NAME_DIR "@rpath"
  )
  set(HEPMC3_CFLAGS -I${HEPMC3_INCLUDE_DIR})
  set(HEPMC3_LIBS   -L${MARLEY_BUILD_LIBDIR} -lMARLEYHepMC3)
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

# dladdr() is used for self-location (design spec 2.1). CMAKE_DL_LIBS is
# "dl" on glibc and empty on macOS, where dladdr lives in libSystem --
# hence the guard: target_link_libraries() with nothing after PRIVATE is an
# error, not a no-op.
if(CMAKE_DL_LIBS)
  target_link_libraries(marley_lib PRIVATE ${CMAKE_DL_LIBS})
endif()

# Make the build tree mirror the install layout, so that the single hop taken
# by self-location works before installation. Created at configure time;
# safe as a one-time step because the top-level Makefile's `clean` removes
# the whole build directory and so forces a reconfigure. The links are
# removed first because file(CREATE_LINK) fails if the destination exists,
# which it does on every reconfigure after the first.
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/share/marley")
file(REMOVE "${CMAKE_BINARY_DIR}/include"
            "${CMAKE_BINARY_DIR}/share/marley/data")
file(CREATE_LINK "${CMAKE_SOURCE_DIR}/include" "${CMAKE_BINARY_DIR}/include"
  SYMBOLIC)
file(CREATE_LINK "${CMAKE_SOURCE_DIR}/data"
  "${CMAKE_BINARY_DIR}/share/marley/data" SYMBOLIC)

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
  LIBRARY_OUTPUT_DIRECTORY ${MARLEY_BUILD_LIBDIR}
  PUBLIC_HEADER ""             # optional
  # libMARLEY needs its own RUNPATH, and this is a correctness requirement
  # rather than a convenience. ld.so(8): DT_RUNPATH directories "do not
  # apply to those objects' children", so the executable's RUNPATH cannot
  # help libMARLEY find libMARLEYHepMC3. Neither build system set this
  # before -- CMakeLists only ever set rpath on the *executable*.
  INSTALL_RPATH "${MARLEY_ORIGIN}"
  BUILD_RPATH   "${MARLEY_ORIGIN}"
  # Without this the recorded install name is the build-tree path and the
  # library is not relocatable at all
  MACOSX_RPATH ON
  INSTALL_NAME_DIR "@rpath"
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
    # The dictionary source generated here is compiled directly into
    # marley_lib (which is built as libMARLEY.{so,dylib}, see the
    # OUTPUT_NAME property set above) rather than into a standalone
    # library of its own. Without an explicit MODULE argument, the
    # generated .rootmap file would incorrectly advertise a
    # "libmarley_root_dict" library that is never actually built,
    # causing ROOT's TClass autoloading to fail for the HepMC3 Data
    # classes (e.g. HepMC3::GenEventData) in ordinary ROOT/PyROOT
    # sessions. Setting MODULE to the real output library name fixes
    # the autoload entries to point at libMARLEY instead.
    ROOT_GENERATE_DICTIONARY(${DICT_NAME}
      ${HEPMC3_DATA_HEADERS}
      MODULE MARLEY
      LINKDEF ${DICT_HEADER}
      OPTIONS -I${CMAKE_SOURCE_DIR}/include
              -I${HEPMC3_INCLUDE_DIR}
    )
    target_sources(marley_lib PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/${DICT_NAME}.cxx)
    # NOTE: because MODULE MARLEY was given above, ROOT names the generated
    # pcm/rootmap files after the module (libMARLEY), not after DICT_NAME.
    add_custom_command(TARGET marley_lib POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E make_directory ${MARLEY_BUILD_LIBDIR}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY_rdict.pcm
        ${MARLEY_BUILD_LIBDIR}/libMARLEY_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY.rootmap
        ${MARLEY_BUILD_LIBDIR}/libMARLEY.rootmap
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY.rootmap
    )
  else()
    find_program(ROOTCLING_EXECUTABLE NAMES rootcling rootcint)
    if(NOT ROOTCLING_EXECUTABLE)
      message(FATAL_ERROR "ROOT dictionary helper not found and rootcling is missing.")
    endif()
    set(DICT_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/marley_root_dict.cc)
    # See the comment above the ROOT_GENERATE_DICTIONARY() call for why the
    # "-s" flag (naming the actual output library, libMARLEY) is needed here:
    # without it, the generated .rootmap advertises a nonexistent
    # "libmarley_root_dict" library and ROOT's autoloader fails to find the
    # HepMC3 Data class dictionaries in ordinary ROOT/PyROOT sessions.
    # "-rml"/"-rmf" are additionally required on some ROOT/rootcling versions
    # to actually force generation of the .rootmap file itself ("-s" alone
    # only embeds the library name in the .pcm on those versions). These are
    # the same flags used internally by ROOT's own ROOT_GENERATE_DICTIONARY()
    # CMake macro (see upstream cmake/modules/RootMacros.cmake), so they
    # should be robust across ROOT6 versions.
    set(DICT_LIB_NAME ${CMAKE_SHARED_LIBRARY_PREFIX}MARLEY${CMAKE_SHARED_LIBRARY_SUFFIX})
    add_custom_command(
      OUTPUT ${DICT_OUTPUT}
      COMMAND ${ROOTCLING_EXECUTABLE} -f ${DICT_OUTPUT} -c
        -s ${DICT_LIB_NAME} -rml ${DICT_LIB_NAME} -rmf libMARLEY.rootmap
        -I${CMAKE_SOURCE_DIR}/include
        -I${HEPMC3_INCLUDE_DIR}
        ${HEPMC3_DATA_HEADERS}
        ${DICT_HEADER}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${MARLEY_BUILD_LIBDIR}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY_rdict.pcm
        ${MARLEY_BUILD_LIBDIR}/libMARLEY_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY.rootmap
        ${MARLEY_BUILD_LIBDIR}/libMARLEY.rootmap
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY_rdict.pcm
      COMMAND ${CMAKE_COMMAND} -E rm -f
        ${CMAKE_CURRENT_BINARY_DIR}/libMARLEY.rootmap
      DEPENDS ${DICT_HEADER}
      COMMENT "Generating ROOT dictionary"
    )
    target_sources(marley_lib PRIVATE ${DICT_OUTPUT})
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
  # Relative to the executable's own location, so the install is relocatable.
  # The previous value baked in both CMAKE_INSTALL_PREFIX and the *build
  # directory*, which meant an installed binary kept working only while the
  # build tree survived -- see spec 10.3.
  INSTALL_RPATH "${MARLEY_ORIGIN}/../${MARLEY_LIBDIR_REL}"
  BUILD_RPATH   "${MARLEY_ORIGIN}/../${MARLEY_LIBDIR_REL}"
)

# ------------------------------------------------------------------------------
# Generated shell scripts (design spec 6.4)
#
# mroot and setup_marley.sh are each generated TWICE from one template: once
# for the build tree and once for the install tree. The two trees differ by
# more than a prefix -- a vendored HepMC3 sits at include/builtin before
# installation and at include/marley/builtin after it -- so one self-locating
# script applying one set of relative hops cannot serve both.
#
# The two copies share a basename, so the install-tree copies are staged in
# their own directory and installed from there.
# ------------------------------------------------------------------------------
set(MARLEY_INSTALL_SCRIPT_DIR "${CMAKE_BINARY_DIR}/install-scripts")

set(lib_rel "${MARLEY_LIBDIR_REL}")

# The build tree's include symlink is named "include" unconditionally, while
# the install tree honours CMAKE_INSTALL_INCLUDEDIR. Likewise bin.
set(MARLEY_BUILD_BIN_REL     "bin")
set(MARLEY_BUILD_INCLUDE_REL "include")
if(USE_BUILTIN_HEPMC3)
  set(MARLEY_BUILD_BUILTIN_REL   "include/builtin")
  set(MARLEY_INSTALL_BUILTIN_REL "${MARLEY_INCLUDEDIR_REL}/marley/builtin")
else()
  set(MARLEY_BUILD_BUILTIN_REL   "")
  set(MARLEY_INSTALL_BUILTIN_REL "")
endif()

function(marley_generate_script template out_path bin_rel include_rel
                                builtin_rel)
  set(lib_rel "${MARLEY_LIBDIR_REL}")
  configure_file("${CMAKE_SOURCE_DIR}/src/scripts/${template}" "${out_path}"
    @ONLY)
  file(CHMOD "${out_path}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                GROUP_READ GROUP_EXECUTE
                WORLD_READ WORLD_EXECUTE
  )
endfunction()

foreach(_script setup_marley.sh)
  marley_generate_script("${_script}.in" "${CMAKE_BINARY_DIR}/bin/${_script}"
    "${MARLEY_BUILD_BIN_REL}" "${MARLEY_BUILD_INCLUDE_REL}"
    "${MARLEY_BUILD_BUILTIN_REL}")
  marley_generate_script("${_script}.in"
    "${MARLEY_INSTALL_SCRIPT_DIR}/${_script}"
    "${MARLEY_BINDIR_REL}" "${MARLEY_INCLUDEDIR_REL}"
    "${MARLEY_INSTALL_BUILTIN_REL}")
endforeach()

if(USE_ROOT)
  marley_generate_script("mroot.in" "${CMAKE_BINARY_DIR}/bin/mroot"
    "${MARLEY_BUILD_BIN_REL}" "${MARLEY_BUILD_INCLUDE_REL}"
    "${MARLEY_BUILD_BUILTIN_REL}")
  marley_generate_script("mroot.in" "${MARLEY_INSTALL_SCRIPT_DIR}/mroot"
    "${MARLEY_BINDIR_REL}" "${MARLEY_INCLUDEDIR_REL}"
    "${MARLEY_INSTALL_BUILTIN_REL}")
endif()

# ------------------------------------------------------------------------------
# marley-config and marley.pc (design spec 6.2, 6.3)
#
# marley-config is concatenated, not substituted: a generated header of
# KEY='value' assignments followed by the static body in src/scripts/.
#
# Note what this section does NOT emit: any MARLEY-owned path. Those are
# computed by the body from the script's own location, so an installed
# marley-config structurally cannot report a build-tree -L or a source-tree
# -I. That pair of defects is what forced a patch in larsim.
# ------------------------------------------------------------------------------
set(MARLEY_CONFIG_BODY   ${CMAKE_SOURCE_DIR}/src/scripts/marley-config-body.sh)
set(MARLEY_CONFIG_OUTPUT ${CMAKE_BINARY_DIR}/bin/marley-config)

# ROOT's flags are captured here and emitted verbatim thereafter. The script
# used to shell out to root-config at run time, which was the last remaining
# way for it to describe a ROOT other than the one MARLEY was compiled
# against -- a consumer would then compile successfully and misbehave. Frozen
# flags can go stale instead, which fails loudly by naming a path that is no
# longer there.
# Both discovery paths ask root-config, so the published flags are the same
# parse in both build systems rather than two reconstructions that drift. On
# the config-package path there are no flat flags to reuse -- only imported
# targets -- so root-config is located explicitly; it lives in ROOT's bindir
# and is present in every ROOT installation that provides the CMake package.
set(MARLEY_ROOT_CFLAGS "")
set(MARLEY_ROOT_LIBS   "")
if(USE_ROOT)
  if(DEFINED ROOT_CONFIG_EXECUTABLE AND ROOT_CONFIG_EXECUTABLE)
    set(_rootconfig "${ROOT_CONFIG_EXECUTABLE}")
  else()
    find_program(MARLEY_ROOTCONFIG root-config
      HINTS ${ROOT_BINDIR} ${ROOT_bindir} $ENV{ROOTSYS}/bin)
    set(_rootconfig "${MARLEY_ROOTCONFIG}")
  endif()

  if(_rootconfig)
    execute_process(COMMAND ${_rootconfig} --cflags
      OUTPUT_VARIABLE _root_cflags OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${_rootconfig} --libs
      OUTPUT_VARIABLE MARLEY_ROOT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
    # cxx_std is authoritative -- it is the maximum of MARLEY's floor, any
    # requested standard and ROOT's -- so ROOT's own -std can only duplicate
    # it or contradict it.
    string(REGEX REPLACE "-std=[^ \t]*" "" MARLEY_ROOT_CFLAGS "${_root_cflags}")
    string(REGEX REPLACE "  +" " " MARLEY_ROOT_CFLAGS "${MARLEY_ROOT_CFLAGS}")
    string(STRIP "${MARLEY_ROOT_CFLAGS}" MARLEY_ROOT_CFLAGS)
  else()
    # Fall back to the include directories alone. This is incomplete -- it
    # omits -pthread and anything else ROOT considers necessary -- so say so
    # rather than publishing quietly degraded flags.
    foreach(_d IN LISTS ROOT_INCLUDE_DIRS)
      list(APPEND MARLEY_ROOT_CFLAGS "-I${_d}")
    endforeach()
    if(DEFINED ROOT_LIBRARY_DIR AND NOT "${ROOT_LIBRARY_DIR}" STREQUAL "")
      list(APPEND MARLEY_ROOT_LIBS "-L${ROOT_LIBRARY_DIR}")
    endif()
    list(APPEND MARLEY_ROOT_LIBS -lCore -lRIO -lHist -lTree)
    message(WARNING "root-config was not found, so marley-config and "
      "marley.pc will publish approximate ROOT flags (include and library "
      "directories only).")
  endif()
endif()

# CMake lists are semicolon-separated; the shell wants spaces.
function(marley_flatten_flags in_var out_var)
  string(REPLACE ";" " " _flat "${${in_var}}")
  set(${out_var} "${_flat}" PARENT_SCOPE)
endfunction()

marley_flatten_flags(MARLEY_ROOT_CFLAGS MARLEY_ROOT_CFLAGS_STR)
marley_flatten_flags(MARLEY_ROOT_LIBS   MARLEY_ROOT_LIBS_STR)

if(USE_BUILTIN_HEPMC3)
  set(MARLEY_HEPMC3_MODE "builtin")
  set(MARLEY_EXT_HEPMC3_CFLAGS_STR "")
  set(MARLEY_EXT_HEPMC3_LIBS_STR   "")
else()
  set(MARLEY_HEPMC3_MODE "external")
  # MARLEY_FOUND_HEPMC3 is an internal build switch -- no installed header
  # tests it -- so it must not reach a consumer's compiler (design spec 6.2).
  set(_hepmc3_public ${HEPMC3_CFLAGS})
  list(REMOVE_ITEM _hepmc3_public "-DMARLEY_FOUND_HEPMC3")
  marley_flatten_flags(_hepmc3_public MARLEY_EXT_HEPMC3_CFLAGS_STR)
  marley_flatten_flags(HEPMC3_LIBS    MARLEY_EXT_HEPMC3_LIBS_STR)
endif()

if(GSL_FOUND)
  set(MARLEY_GSL_MODE "external")
  marley_flatten_flags(GSL_CFLAGS MARLEY_EXT_GSL_CFLAGS_STR)
  marley_flatten_flags(GSL_LIBS   MARLEY_EXT_GSL_LIBS_STR)
else()
  set(MARLEY_GSL_MODE "builtin")
  set(MARLEY_EXT_GSL_CFLAGS_STR "")
  set(MARLEY_EXT_GSL_LIBS_STR   "")
endif()

if(USE_ROOT)
  set(MARLEY_USE_ROOT_STR "yes")
else()
  set(MARLEY_USE_ROOT_STR "no")
endif()

if(APPLE)
  set(MARLEY_NEW_DTAGS "")
else()
  set(MARLEY_NEW_DTAGS "-Wl,--enable-new-dtags")
endif()

set(MARLEY_CXX_STD "c++${CMAKE_CXX_STANDARD}")

string(CONCAT MARLEY_CONFIG_HEADER
  "#!/bin/bash\n"
  "# Generated at configure time by CMakeLists.txt -- do not edit.\n"
  "# The static half of this script is src/scripts/marley-config-body.sh.\n"
  "version='${MARLEY_COMPUTED_VERSION}'\n"
  "git_revision='${MARLEY_COMPUTED_GIT_REVISION}'\n"
  "cxx_std='${MARLEY_CXX_STD}'\n"
  "compiled_log_level='${MARLEY_LOG_LEVEL}'\n"
  "use_root='${MARLEY_USE_ROOT_STR}'\n"
  "hepmc3_mode='${MARLEY_HEPMC3_MODE}'\n"
  "gsl_mode='${MARLEY_GSL_MODE}'\n"
  "root_cflags='${MARLEY_ROOT_CFLAGS_STR}'\n"
  "root_libs='${MARLEY_ROOT_LIBS_STR}'\n"
  "ext_hepmc3_cflags='${MARLEY_EXT_HEPMC3_CFLAGS_STR}'\n"
  "ext_hepmc3_libs='${MARLEY_EXT_HEPMC3_LIBS_STR}'\n"
  "ext_gsl_cflags='${MARLEY_EXT_GSL_CFLAGS_STR}'\n"
  "ext_gsl_libs='${MARLEY_EXT_GSL_LIBS_STR}'\n"
  "new_dtags='${MARLEY_NEW_DTAGS}'\n"
  "libdir_rel='${MARLEY_LIBDIR_REL}'\n"
  "bindir_rel='${MARLEY_BINDIR_REL}'\n"
  "includedir_rel='${MARLEY_INCLUDEDIR_REL}'\n"
)

file(READ ${MARLEY_CONFIG_BODY} MARLEY_CONFIG_BODY_CONTENT)
file(WRITE ${MARLEY_CONFIG_OUTPUT}
  "${MARLEY_CONFIG_HEADER}${MARLEY_CONFIG_BODY_CONTENT}")
file(CHMOD ${MARLEY_CONFIG_OUTPUT}
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
              GROUP_READ GROUP_EXECUTE
              WORLD_READ WORLD_EXECUTE
)

# marley.pc carries the same flags marley-config does, minus the rpath. The
# conditional fields are whole-string substitution keys rather than a bracket
# notation, because configure_file cannot express a conditional and both
# build systems must produce the same file.
#
# There is no Requires: line at all -- see the template for why. The short
# version: neither HepMC3 nor ROOT ships a pkg-config file, and listing GSL
# as Requires.private would push GSL's -I back into "pkg-config --cflags",
# since pkg-config includes private Cflags unconditionally.
#
# -DUSE_ROOT travels with ROOT's include path and cannot be separated from
# it: EventFileReader.hh includes TFile.h and TTree.h under that macro, so
# publishing the macro alone yields a public header that will not parse.
set(pc_cflags_extra "")
set(pc_libs_extra "")
if(USE_ROOT)
  string(APPEND pc_cflags_extra " -DUSE_ROOT ${MARLEY_ROOT_CFLAGS_STR}")
  # ROOT's own -Wl,-rpath is stripped here and nowhere else. marley-config
  # may carry it, because a hand-written Makefile has no other way to make
  # its binary runnable; a .pc file must not, because a packager's toolchain
  # manages rpath itself and this one would bake the build machine's ROOT
  # location into a redistributed binary (design spec 6.3).
  string(REGEX REPLACE "-Wl,-rpath,[^ \t]*" "" _pc_root_libs
    "${MARLEY_ROOT_LIBS_STR}")
  string(REGEX REPLACE "  +" " " _pc_root_libs "${_pc_root_libs}")
  string(STRIP "${_pc_root_libs}" _pc_root_libs)
  string(APPEND pc_libs_extra " ${_pc_root_libs}")
endif()
if(USE_BUILTIN_HEPMC3)
  string(PREPEND pc_cflags_extra " -I\${includedir}/marley/builtin")
  string(PREPEND pc_libs_extra " -lMARLEYHepMC3")
else()
  string(PREPEND pc_cflags_extra " ${MARLEY_EXT_HEPMC3_CFLAGS_STR}")
  string(PREPEND pc_libs_extra " ${MARLEY_EXT_HEPMC3_LIBS_STR}")
endif()
string(STRIP "${pc_cflags_extra}" pc_cflags_extra)
string(STRIP "${pc_libs_extra}" pc_libs_extra)

set(libdir_rel     "${MARLEY_LIBDIR_REL}")
set(includedir_rel "${MARLEY_INCLUDEDIR_REL}")
set(version        "${MARLEY_COMPUTED_VERSION}")
set(cxx_std        "${MARLEY_CXX_STD}")
configure_file(
  ${CMAKE_SOURCE_DIR}/src/scripts/marley.pc.in
  ${MARLEY_BUILD_LIBDIR}/pkgconfig/marley.pc
  @ONLY
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
set(MARLEY_BUILD_CONFIG_DIR "${MARLEY_BUILD_LIBDIR}/cmake/MARLEY")
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
  INSTALL_DESTINATION "${MARLEY_LIBDIR_REL}/cmake/MARLEY"
)
write_basic_package_version_file(
  "${MARLEY_BUILD_CONFIG_DIR}/MARLEYConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)

# Install-tree export
install(TARGETS marley_lib
  EXPORT MARLEYTargets
  LIBRARY DESTINATION ${MARLEY_LIBDIR_REL}
)
install(TARGETS marley
  RUNTIME DESTINATION ${MARLEY_BINDIR_REL}
)
# From the staging directory, not bin/: the two copies share a basename and
# differ in the hops they take (design spec 6.4).
install(PROGRAMS ${MARLEY_INSTALL_SCRIPT_DIR}/setup_marley.sh
  DESTINATION ${MARLEY_BINDIR_REL}
)
if(USE_ROOT)
  install(PROGRAMS ${MARLEY_INSTALL_SCRIPT_DIR}/mroot
    DESTINATION ${MARLEY_BINDIR_REL}
  )
  install(FILES
    ${MARLEY_BUILD_LIBDIR}/libMARLEY_rdict.pcm
    ${MARLEY_BUILD_LIBDIR}/libMARLEY.rootmap
    DESTINATION ${MARLEY_LIBDIR_REL}
    OPTIONAL
  )
endif()

# If built-in HepMC3 is used, install & export it (and its headers).
#
# Both the library name and the header location are deliberately MARLEY's,
# not HepMC3's. This is an amalgamated, partial HepMC3 built from
# src/marley_hepmc3.cc, produced only when no system HepMC3 was found -- but
# nothing stops a real HepMC3 arriving in the same prefix later, and
# /usr/local and Spack views make that ordinary. Two different libHepMC3, or
# two different include/HepMC3, fail *silently*: whichever is reached first
# wins, and this one does not implement all of HepMC3.
#
# Nesting the headers also means no #include directive changes anywhere --
# the -I root travels with them, so "HepMC3/FourVector.h" is spelled
# identically before and after.
if(USE_BUILTIN_HEPMC3)
  install(TARGETS HepMC3_builtin
    EXPORT MARLEYTargets
    LIBRARY DESTINATION ${MARLEY_LIBDIR_REL}
  )
  install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/builtin/HepMC3
    DESTINATION ${MARLEY_INCLUDEDIR_REL}/marley/builtin
  )
endif()

install(PROGRAMS ${MARLEY_CONFIG_OUTPUT}
  DESTINATION ${MARLEY_BINDIR_REL}
)
# pkgconfig follows libdir, so a lib64 build gets lib64/pkgconfig. The file's
# own ${pcfiledir}/../.. still reaches the prefix, because libdir_rel is one
# segment by construction (design spec 4.3).
install(FILES ${MARLEY_BUILD_LIBDIR}/pkgconfig/marley.pc
  DESTINATION ${MARLEY_LIBDIR_REL}/pkgconfig
)
# Only MARLEY's own headers. The previous rule swept include/ for "*.hh",
# which also installed include/builtin/marley_gsl.hh to <prefix>/include/
# builtin/ -- a top-level entry MARLEY does not own, and a divergence from
# the Make path, which copies only include/marley. marley_gsl.hh is an
# internal GSL shim; the one test that needs it includes it by relative path.
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/marley
  DESTINATION ${MARLEY_INCLUDEDIR_REL}
  FILES_MATCHING PATTERN "*.hh"
)
# data/tests/ is used only by martest, which runs from the build tree
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/
  DESTINATION share/marley/data
  PATTERN "tests" EXCLUDE
)
# examples/ ships: the documentation references it 15 times across 5 files,
# including commands users are told to type verbatim, and
# examples/config/annotated.js is the primary job-configuration reference.
# None of it worked on an installed MARLEY.
install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples
  DESTINATION share/marley
)
# COPYING and GUIDELINES ship because the licence header of *every installed
# header* points at them, and MARLEY is GPLv3, so distributing binaries
# without the licence text is a compliance gap. CITATION.bib joins them
# because the documentation asks users to cite.
install(FILES
  ${CMAKE_SOURCE_DIR}/COPYING
  ${CMAKE_SOURCE_DIR}/GUIDELINES
  ${CMAKE_SOURCE_DIR}/CITATION.bib
  DESTINATION share/marley
)

# Install the export set and package config files
install(EXPORT MARLEYTargets
  NAMESPACE MARLEY::
  DESTINATION ${MARLEY_LIBDIR_REL}/cmake/MARLEY
)
configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/cmake/MARLEYConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/MARLEYConfig.cmake"
  INSTALL_DESTINATION "${MARLEY_LIBDIR_REL}/cmake/MARLEY"
)
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/MARLEYConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(FILES
  "${CMAKE_BINARY_DIR}/MARLEYConfig.cmake"
  "${CMAKE_BINARY_DIR}/MARLEYConfigVersion.cmake"
  DESTINATION "${MARLEY_LIBDIR_REL}/cmake/MARLEY"
)

# ------------------------------------------------------------------------------
# Documentation targets (optional)
# ------------------------------------------------------------------------------
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  add_custom_target(doxygen
    COMMAND ${CMAKE_COMMAND} -E env
      MARLEY_VERSION=${VERSION_PREFIX}${MARLEY_COMPUTED_VERSION}
      DOXYGEN_OUTPUT_DIR=${CMAKE_BINARY_DIR}/docs/doxygen
      BUILDDIR=${CMAKE_BINARY_DIR}/docs
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
      DOXYGEN_OUTPUT_DIR=${CMAKE_BINARY_DIR}/docs/doxygen
      BUILDDIR=${CMAKE_BINARY_DIR}/docs
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
# Test executable
#
# 'make test' (via the custom target below) builds martest without running it,
# mirroring the GNU Make path in make/build.mk where 'test' depends only on
# the test executable.  Run the tests manually with:
#   MARLEY=<source_root> ./build/bin/martest
# ------------------------------------------------------------------------------
file(GLOB MARLEY_TEST_SOURCES
  "${CMAKE_SOURCE_DIR}/src/tests/*.cc"
)
if(NOT GSL_FOUND)
  # If the system GSL library was not found (or it was ignored), then we need
  # to include the built-in GSL code when building the test executable.
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
  # martest is EXCLUDE_FROM_ALL and never installed, but it still has to find
  # libMARLEY when run from the build tree
  BUILD_RPATH "${MARLEY_ORIGIN}/../${MARLEY_LIBDIR_REL}"
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

# 'make test' builds the test executable without running it, consistent with
# the GNU Make path.  CMake's built-in 'test' target (from enable_testing()) is
# intentionally not used so that both build systems share the same semantics.
add_custom_target(test DEPENDS martest)
