# Introduction
#### An R package for the Combined Analysis of Epistasis and Pleiotropy  
The CAPE2.0 R package implements a method, originally described in Carter et al. (2012), that infers directed interaction networks between genetic variants for predicting the inﬂuence of genetic perturbations on phenotypes. This method takes advantage of complementary information in partially pleiotropic genetic variants to resolve directional inﬂuences between variants that interact epistatically. **cape** can be applied to a variety of genetic variants, such as single nucleotide polymorphisms (SNPs), copy number variations (CNVs) or structural variations (SVs).

CAPE 2.0 has a vignette, which can be found [here](). 

## New Features!
- new run_cape function runs the entire cape pipeline with one command
- read in data in multiple formats (R/qtl, R/qtl2, and PLINK)
- performs kinship correction using linear mixed models as described in Kang et al. (2008)
- Handles multi-parent populations
- R6 reformatting improves speed and handling of large data

Kang, H. M. et al. Efficient control of population structure in model organism association mapping. Genetics 178, 1709–1723 (2008).

## Installation
CAPE 2.0 requires R 3.6+ to run. 

Once you download CAPE 2.0, navigate to the the project's root directory and in it create a new directory called **/cape_pkg**, so that the path **/cape/cape_pkg** exists.

To set your R working directory, start R and enter:
```sh
setwd('path/to/cape2')
```

To download and install all **cape** dependencies, start R and enter:
```sh
source('config/install_prerequisites.R')
```

To compile the package, start R and enter:
**Note**: if you are using Windows, you will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 
```sh
source('config/compile_package.R')
```


## Tests
CAPE2.0 provides demo scripts, which you can run to verify that the installation was successful.

```sh
source('demo/a.example.R')
```

```sh
source('demo/a.example.DO.R')
```

## To-Do:
- enable CAPE 2.0 run in parallel

## License
CAPE2.0 is licensed under [GPL-3](https://www.r-project.org/Licenses/GPL-3)