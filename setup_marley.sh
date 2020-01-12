#!/bin/bash

# Finds the directory in which this script is located. This method isn't
# foolproof. See https://stackoverflow.com/a/246128/4081973 if you need
# something more robust for edge cases (e.g., you're calling the script using
# symlinks).
THIS_DIRECTORY="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# The root MARLEY folder
export MARLEY=${THIS_DIRECTORY}

# For running MARLEY
export PATH=$PATH:${THIS_DIRECTORY}/build

if [ "$(uname)" == "Darwin" ]; then
  # Mac OS X platform
  export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${THIS_DIRECTORY}/build
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
  # GNU/Linux platform
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${THIS_DIRECTORY}/build
fi

# For using MARLEY classes in ROOT 6
export ROOT_INCLUDE_PATH=${THIS_DIRECTORY}/include

# For building external software against MARLEY
export CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:${THIS_DIRECTORY}/include
export LIBRARY_PATH=${LIBRARY_PATH}:${THIS_DIRECTORY}/build
