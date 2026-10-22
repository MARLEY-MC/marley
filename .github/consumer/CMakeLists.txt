# A minimal external consumer of an INSTALLED MARLEY.
#
# Nothing else in the tree calls find_package(MARLEY), which is exactly why
# the config package's commented-out find_dependency block went unnoticed:
# find_package() resolved, handed back imported targets referring to HepMC3
# and ROOT without either having been looked for, and the consumer failed
# later at link time.
#
# It lives under .github/ rather than examples/ on purpose. It is CI
# infrastructure, not something MARLEY ships -- examples/ is installed, and
# a consumer project inside an installation would be odd.
#
# Run it by hand against any prefix:
#
#   cmake -S .github/consumer -B /tmp/c -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build /tmp/c && /tmp/c/consumer

cmake_minimum_required(VERSION 3.18)
project(marley_consumer LANGUAGES CXX)

# REQUIRED, and with a version: the package once declared itself version 6
# by matching the leading digit of a git hash, which no unversioned call
# would have caught.
find_package(MARLEY 2.0 REQUIRED)

# Report what the package actually defined. A wrong guess at the target name
# and a package that exported nothing are otherwise indistinguishable.
get_property(_tgts DIRECTORY PROPERTY IMPORTED_TARGETS)
list(FILTER _tgts INCLUDE REGEX "^MARLEY::")
message(STATUS "MARLEY imported targets: ${_tgts}")
message(STATUS "MARLEY_VERSION=${MARLEY_VERSION}")

add_executable(consumer main.cc)

# No target_include_directories and no target_compile_features here, on
# purpose: both must arrive through the imported target. If they do not,
# this fails to compile, which is the point.
target_link_libraries(consumer PRIVATE MARLEY::MARLEY)
