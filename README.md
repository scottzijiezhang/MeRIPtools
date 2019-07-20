# MeRIPtools
Tool sets to analyze high throughput data for RNA modifications

### Install the R package from Github

Depends: GenomicFeatures, Rsamtools, ggplot2, doParallel, foreach,grid,rtracklayer,GenomicAlignments,reshape2,Rcpp,RcppArmadillo,
Guitar, stringr,vcfR,gamlss, broom, DESeq2

	install.packages("devtools")
	library(devtools)
	install_github("scottzijiezhang/MeRIPtools")
	library("MeRIPtools")

## Manual page

Please refer to [manual page](https://scottzijiezhang.github.io/MeRIPtoolsManual/) for detailed instructions.  

## Citation 
MeRIPtools is a tool sets that implemented functions for peak calling, QTL calling, differential methylation analysis, visualization. 

If you used MeRIPtools in your publication, please cite:
*Something to be published...*

**Note** MeRIPtools also have wrapper functions to call functions from other R packages to do specific analysis.  
If you used the `plotMetaGene` or `MetaGene` function, please cite the original R package [`Guitar`](https://bioconductor.org/packages/release/bioc/html/Guitar.html)  
Cui X, Wei Z, Zhang L, Liu H, Sun L, Zhang s, Huang Y, Meng J (2016). “Guitar: an R/Bioconductor package for gene annotation guided transcriptomic analysis of RNA related genomic features.” BioMed Research International. 
