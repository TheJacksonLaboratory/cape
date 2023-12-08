## New in version 3.1.1

## New in version 3.1.2

### Bug Fixes

* Edited vignette to better explain how to set up the code
    and directories so that the demo code runs.

* We removed examples from function documentation that do not 
    work outside of run_cape

### New Features

* `pracma::pinv()` replaced `pseudoinverse()`.

* Added param_file argument to `load_input_and_run_cape()`.

* Added ability for `cape2mpp()` to recognize and format multi-parent crosses

* Added ability to read in multi-parent crosses in qtl format.

### Bug Fixes

* Fixed display of marker effects in `plot_singlescan()`.

* Updated missing parameters from `load_input_and_run_cape()` documentation

* Changed binary file extensions of saved files to .RData.

* Additional minor bug fixes.

## New in version 3.1.0

### New Features

* `run_cape()` now runs a full cape analysis from a yml parameter file.

* Dots were replaced with underscores in all function names.

* All code was refactored to S6.

* `qtl2cape()` now accepts a covariate matrix as an argument.

* Added demos starting from different file formats. 

* Added formatting capabilities for PLINK files.

### Bug Fixes

* `bin_curve()` was updated to handle very short vectors.

* Additional minor bug fixes.

* Fixed covariate handling in `pairscan()`.