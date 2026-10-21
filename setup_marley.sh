# Wrapper. The real script is generated into the build directory, because it
# has to know the library directory name and whether HepMC3 was vendored --
# neither of which is knowable before the project is configured.
#
# Every documented "source setup_marley.sh" invocation keeps working through
# this file; it just forwards to the generated one.
#
# The failure mode this replaces is worth stating: before, sourcing this
# script in an unbuilt tree quietly exported paths into a build/ directory
# that did not exist, and the user found out later as "command not found".

_marley_wrapper_dir="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"

if [ -r "${_marley_wrapper_dir}/build/bin/setup_marley.sh" ]; then
  . "${_marley_wrapper_dir}/build/bin/setup_marley.sh"
else
  echo "setup_marley.sh: ${_marley_wrapper_dir}/build/bin/setup_marley.sh" \
    "does not exist." >&2
  echo "  Build MARLEY first (run 'make'), then source this script again." >&2
  echo "  For an out-of-source CMake build, source the copy in that build" >&2
  echo "  directory instead: source <builddir>/bin/setup_marley.sh" >&2
  echo "  For an installed MARLEY, source <prefix>/bin/setup_marley.sh" >&2
fi

unset _marley_wrapper_dir
