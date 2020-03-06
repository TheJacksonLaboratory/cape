#!/usr/bin/env bash

# This script runs the tests of the CAPE analysis package 
# Run this from the root of the package's source directory

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo 'Running the CAPE package tests...'
#echo 'Running the CAPE package tests in ${DIR}'
#echo 'direcotry contents...'
#ls -d */
#echo 'end directory contents.'
Rscript --vanilla "${DIR}/ci-tests.R"
