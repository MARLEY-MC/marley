# Prepends the 'build/bin' subfolder onto the system PATH, allowing MARLEY
# to be run easily from an in-source build.
#
# Usage: source set_marley_path.sh
#
# For most use cases, setting the PATH is sufficient. In situations where other
# environment variables (e.g., ROOT_INCLUDE_PATH) are needed, one may source the
# generated setup_marley.sh script from the build/bin/ folder or from
# <prefix>/bin/ in any MARLEY installation.

_marley_src_dir="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"

if [ -x "${_marley_src_dir}/build/bin/marley" ]; then
  export PATH="${_marley_src_dir}/build/bin${PATH:+:${PATH}}"
else
  # A missing build/bin would otherwise fail later as "command not found".
  echo "set_marley_path.sh: ${_marley_src_dir}/build/bin/marley does not" \
    "exist." >&2
  echo "  Build MARLEY first (run 'make'), then source this script again." >&2
  echo "  For an out-of-source CMake build, add <builddir>/bin to PATH" >&2
  echo "  yourself. For an installed MARLEY, add <prefix>/bin to PATH, or" >&2
  echo "  source <prefix>/bin/setup_marley.sh for the full environment." >&2
fi

unset _marley_src_dir
