# CAPE

<!-- badges: start -->
[![License: 
GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

# Introduction
#### An R package for the Combined Analysis of Epistasis and Pleiotropy  
The CAPE R package implements a method, originally described in Carter et al. (2012), that 
infers directed interaction networks between genetic variants for predicting the inﬂuence of 
genetic perturbations on phenotypes. This method takes advantage of complementary information 
in partially pleiotropic genetic variants to resolve directional inﬂuences between variants 
that interact epistatically. **CAPE** can be applied to a variety of genetic variants, such 
as single nucleotide polymorphisms (SNPs), copy number variations (CNVs) or structural variations (SVs).

For detailed documentation about how to format data, load data, and analyze data, please see 
the CAPE vignette. 

## New Features!
- new *run_cape()* function runs the entire cape pipeline with one command
- read in data in multiple formats (R/qtl, R/qtl2, and PLINK)
- performs kinship correction using linear mixed models as described in Kang et al. (2008)
- Handles multi-parent populations
- R6 reformatting improves speed and handling of large data

## Installation
CAPE requires R 3.6+ to run.

```sh
# Install the released version from CRAN
install.packages("cape")

# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("TheJacksonLaboratory/cape")
```

## Demos
CAPE provides demo scripts, which you can run to verify that the installation was successful.

```r
demo(package = "cape")
```

```r
demo(demo_plink)
```

```r
demo(demo_qtl)
```

```r
demo(demo_qtl2)
```

## To-Do:
- enable CAPE run in parallel

## License
CAPE is licensed under [GPL-3](https://www.r-project.org/Licenses/GPL-3)

## References

Tyler, A. L., Lu, W., Hendrick, J. J., Philip, V. M. & Carter, G. W. CAPE: an R package for 
combined analysis of pleiotropy and epistasis. PLoS Comput. Biol. 9, e1003270 (2013).

Kang, H. M. et al. Efficient control of population structure in model organism association 
mapping. Genetics 178, 1709–1723 (2008).

## Related Publications
Carter, G. W. Inferring gene function and network organization in Drosophila signaling by 
combined analysis of pleiotropy and epistasis. G3 (Bethesda) 3, 807–814 (2013).

Carter, G. W., Hays, M., Sherman, A. & Galitski, T. Use of pleiotropy to model genetic interactions 
in a population. PLoS Genet. 8, e1003010 (2012).

Tyler, A. L., McGarr, T. C., Beyer, B. J., Frankel, W. N. & Carter, G. W. A genetic interaction 
network model of a complex neurological disease. Genes, Brain and Behavior 13, 831–840 (2014).

Tyler, A. L. et al. Epistatic Networks Jointly Influence Phenotypes Related to Metabolic Disease 
and Gene Expression in Diversity Outbred Mice. Genetics 206, 621–639 (2017).

Tyler, A. L. et al. Epistatic networks jointly influence phenotypes related to metabolic disease 
and gene expression in diversity outbred mice. Genetics 206, 621–639 (2017).

Tyler, A. L., Donahue, L. R., Churchill, G. A. & Carter, G. W. Weak Epistasis Generally 
Stabilizes Phenotypes in a Mouse Intercross. PLoS Genet. 12, e1005805–22 (2016).