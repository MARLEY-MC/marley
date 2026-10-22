# A minimal external consumer of an installed MARLEY, and the only thing in
# the tree that calls find_package(MARLEY).
#
# It lives under .github/ because it is CI infrastructure rather than
# something MARLEY ships.
#
# Run it by hand against any prefix:
#
#   cmake -S .github/consumer -B /tmp/c -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build /tmp/c && /tmp/c/consumer

cmake_minimum_required(VERSION 3.18)
project(marley_consumer LANGUAGES CXX)

# A version is given so that a package reporting the wrong one fails here.
find_package(MARLEY 2.0 REQUIRED)

# Report what the package defined, for the CI log.
get_property(_tgts DIRECTORY PROPERTY IMPORTED_TARGETS)
list(FILTER _tgts INCLUDE REGEX "^MARLEY::")
message(STATUS "MARLEY imported targets: ${_tgts}")
message(STATUS "MARLEY_VERSION=${MARLEY_VERSION}")

add_executable(consumer main.cc)

# No target_include_directories and no target_compile_features on purpose:
# both must arrive through the imported target, or this fails to compile.
target_link_libraries(consumer PRIVATE MARLEY::MARLEY)
