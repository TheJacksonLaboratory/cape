#!/usr/bin/env bash

# Run this command from the top level directory, e.g.:

# > cd cape/config/..
# > ./config/install_prerequisites.sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

Rscript "${DIR}/install_prerequisites.R"