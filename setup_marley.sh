#!/bin/bash

# Finds the directory in which this script is located. This method isn't
# foolproof. See https://stackoverflow.com/a/246128/4081973 if you need
# something more robust for edge cases (e.g., you're calling the script using
# symlinks). To make this technique work in zsh as well as bash, use the
# approach described here: https://stackoverflow.com/a/54755784
THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]:-${(%):-%x}}" )" && pwd )"

# The root MARLEY folder
export MARLEY=${THIS_DIRECTORY}

# For running MARLEY
export PATH=${PATH}:${THIS_DIRECTORY}/build

if [ "$(uname)" = "Darwin" ]; then
  # macOS platform
  export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THIS_DIRECTORY}/build
else
  # Assume a GNU/Linux platform
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THIS_DIRECTORY}/build
fi

# For using MARLEY classes in ROOT 6
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${THIS_DIRECTORY}/include
