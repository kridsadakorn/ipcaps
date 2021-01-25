
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IPCAPS <img src="man/figures/ipcaps_logo.png" align="right" />

<!-- badges: start -->

<!-- badges: end -->

`IPCAPS` is an unsupervised clustering algorithm based on iterative
pruning to capture population structure. This version supports ordinal
data which can be applied directly to SNP data to identify fine-level
population structure and it is built on the iterative pruning Principal
Component Analysis (ipPCA) algorithm by Intarapanich et al. (2009)
\<doi: 10.1186/1471-2105-10-382\> and Limpiti et al. (2011)\<doi:
10.1186/1471-2105-12-255\>. The IPCAPS involves an iterative process
using multiple splits based on multivariate Gaussian mixture modeling of
principal components and Clustering EM estimation as in Lebret et
al. (2015) \<doi: 10.18637/jss.v067.i06\>. In each iteration, rough
clusters and outliers are also identified using the function
`rubikclust()` from the R package `KRIS`.

## Installation

You can install the released version of IPCAPS from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("IPCAPS")
```

Alternatively, you can install the dev version of IPCAPS from
[Gitlab](https://gitlab.com/chaichoompu/ipcaps) with

``` r
install.packages("remotes")
remotes::install_gitlab("chaichoompu/ipcaps", dependencies = TRUE)
```

## Document

You can see the reference manual from:
<https://chaichoompu.gitlab.io/ipcaps_doc/index.html>

## Example

This is a basic example which shows you how to use the packages:

``` r
library(IPCAPS)

BED.file <- system.file("extdata", "ipcaps_example.bed", package = "IPCAPS")
LABEL.file <- system.file("extdata", "ipcaps_example_individuals.txt.gz",
                          package = "IPCAPS")
my.cluster1 <- ipcaps(bed = BED.file, label.file = LABEL.file, lab.col = 2,
out = tempdir())
#> Running ... IPCAPS 
#>  output: /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4 
#>  label file: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IPCAPS/extdata/ipcaps_example_individuals.txt.gz
#>  label column: 2
#>  threshold: 0.18
#>  minimum Fst: 8e-04
#>  maximum thread: 1
#>  method: mix
#>  bed: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/IPCAPS/extdata/ipcaps_example.bed
#>  minimum in group: 20
#>  data type: snp
#>  missing: NA
#> In preprocessing step
#> Note: the directory '/var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4' is existed. 
#> The result files will be saved to this directory: /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output 
#> Input data: 103 individuals, 2000 markers
#> Start calculating
#> Node 1: Start the process
#> Node 1: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 1: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/condition.RData
#> Node 1: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/rawdata.RData
#> Node 1: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node1.RData
#> Node 1: Check for splitting
#> Node 1: Reducing matrix
#> Node 1: EigenFit = 0.652486925995028, Threshold = 0.18, no. significant PCs = 3
#> Node 1: Start clustering
#> Node 1: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node1.RData
#> Node 1: done for clustering
#> Time difference of 0.11854 secs
#> Node 1: Checking Fst of group 1 and group 2
#> Node 1: Checking Fst of group 1 and group 3
#> Node 1: Checking Fst of group 1 and group 4
#> Node 1: Checking Fst of group 2 and group 3
#> Node 1: Checking Fst of group 2 and group 4
#> Node 1: Checking Fst of group 3 and group 4
#> Node 1 times took 0.346548080444336 secs
#> Node 1: 103 individuals were splitted into: 2 groups
#> Node 1: Return status 0
#> Node 1: Split to sub-nodes
#> Node 1: Saving /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node2.RData
#> Node 1: Saving /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node3.RData
#> Node 1: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 1: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 1: Done!
#> Time difference of 0.617434 secs
#> Node 2: Start the process
#> Node 2: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 2: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/condition.RData
#> Node 2: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/rawdata.RData
#> Node 2: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node2.RData
#> Node 2: Check for splitting
#> Node 2: A number of node is lower than the minimum number (20), therefore split was not performed
#> Node 2: Return status 1
#> Node 2: No split was performed, Status=1
#> Node 2: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/leafnode.RData
#> Node 2: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 2: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 2: Done!
#> Time difference of 0.03329301 secs
#> Node 3: Start the process
#> Node 3: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 3: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/condition.RData
#> Node 3: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/rawdata.RData
#> Node 3: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node3.RData
#> Node 3: Check for splitting
#> Node 3: Reducing matrix
#> Node 3: EigenFit = 0.343228072167861, Threshold = 0.18, no. significant PCs = 3
#> Node 3: Start clustering
#> Node 3: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node3.RData
#> Node 3: done for clustering
#> Time difference of 0.1172371 secs
#> Node 3: Checking Fst of group 1 and group 2
#> Node 3 times took 0.15182089805603 secs
#> Node 3: 100 individuals were splitted into: 2 groups
#> Node 3: Return status 0
#> Node 3: Split to sub-nodes
#> Node 3: Saving /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node4.RData
#> Node 3: Saving /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node5.RData
#> Node 3: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 3: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 3: Done!
#> Time difference of 0.3071558 secs
#> Node 4: Start the process
#> Node 4: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 4: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/condition.RData
#> Node 4: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/rawdata.RData
#> Node 4: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node4.RData
#> Node 4: Check for splitting
#> Node 4: Reducing matrix
#> Node 4: EigenFit = 0.0221170820186458, Threshold = 0.18, no. significant PCs = 
#> Node 4: No split was performed because EigenFit is lower than threshold
#> Node 4: Return status 1
#> Node 4: No split was performed, Status=1
#> Node 4: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/leafnode.RData
#> Node 4: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 4: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 4: Done!
#> Time difference of 0.130151 secs
#> Node 5: Start the process
#> Node 5: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 5: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/condition.RData
#> Node 5: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/rawdata.RData
#> Node 5: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/node5.RData
#> Node 5: Check for splitting
#> Node 5: Reducing matrix
#> Node 5: EigenFit = 0.0367738545252401, Threshold = 0.18, no. significant PCs = 
#> Node 5: No split was performed because EigenFit is lower than threshold
#> Node 5: Return status 1
#> Node 5: No split was performed, Status=1
#> Node 5: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/leafnode.RData
#> Node 5: Loading /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 5: Updating /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/RData/tree.RData
#> Node 5: Done!
#> Time difference of 0.1171911 secs
#> In post process step
#> Exporting node 2 as group 1
#> Exporting node 4 as group 2
#> Exporting node 5 as group 3
#> Note: save as /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output/groups.txt
#> The result files were saved at:  /var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output 
#> Total runtime is 2.92840123176575 sec
```

The function `ipcaps` does unsupervised clusering, and here is the
result:

``` r
table(my.cluster1$cluster$label, my.cluster1$cluster$group)
#>           
#>             1  2  3
#>   outlier3  3  0  0
#>   pop1      0 50  0
#>   pop2      0  0 50
```

The output directory will be indicated in the console or in
`my.cluster1$output.dir`. All result files are saved at:You can naviage
to check the `html` visualizations in the output directory.

``` r
print(my.cluster1$output.dir)
#> [1] "/var/folders/sp/hhmj9xvx53z4g4dktf5f503r0000gp/T//RtmpVw7ED4/cluster_output"

list.files(my.cluster1$output.dir)
#> [1] "groups.txt"                "images"                   
#> [3] "RData"                     "tree_scatter_cluster.html"
#> [5] "tree_scatter_label.html"   "tree_scree.html"          
#> [7] "tree_text.html"
```

## About

  - Prof. Kristel Van Steen, visit
    <a href="http://bio3.giga.ulg.ac.be/" border=0 style="border:0; text-decoration:none; outline:none"><img width="40px" src="man/figures/bio3_logo.png" align="center" /></a><br />
  - Kridsadakorn Chaichoompu, visit
    <a href="http://www.biostatgen.org/" border=0 style="border:0; text-decoration:none; outline:none"><img width="110px" src="man/figures/biostatgen_logo.png" align="center" /></a><br />
