#!/usr/bin/env bash

# This script builds the CAPE analysis package in the `packages` directory,
# the installs it into the `Rlib` directory. From there it is accessible,
# and can be imported into rpy2, as a proper R package.
# This script exists to ensure a consistent process that can also be run
# by a CI server.
# Run this from the root of the package directory

# clear out old bits
# rm cape/package/*.gz

echo "${BASH_SOURCE[0]}"
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $DIR
echo 'Compiling and installing the CAPE package...'
PKG=$(Rscript "${DIR}/compile_package.R")

