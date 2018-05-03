#' IPCAPS - Iterative Pruning to CApture Population Structure
#'
#' The R package IPCAPS (Iterative Pruning to CApture Population Structure) is 
#' an unsupervised clustering algorithm based on iterative pruning to capture 
#' population structure. This version supports ordinal data which can be applied 
#' directly to SNP data to identify fine-level population structure and it is 
#' built on the iterative pruning Principal Component Analysis (ipPCA) algorithm 
#' (Intarapanich et al., 2009; Limpiti et al., 2011). The IPCAPS involves an 
#' iterative process using multiple splits based on multivariate Gaussian 
#' mixture modeling of principal components and Clustering EM estimation as in 
#' Lebret et al. (2015). In each iteration, rough clusters and outliers are also 
#' identified using our own method called RubikClust (R package \pkg{kris}).
#'
#' The R package \pkg{IPCAPS} requires the package \pkg{kris}.
#'
#' Here is the list of functions in the R package \pkg{IPCAPS}:
#' \itemize{
#' \item \code{\link{export.groups}}
#' \item \code{\link{get.node.info}}
#' \item \code{\link{ipcaps}}
#' \item \code{\link{save.html}}
#' \item \code{\link{save.plots}}
#' \item \code{\link{top.discriminator}}
#' }
#'
#' Moreover, here is the list of example datasets in the R package \pkg{IPCAPS}:
#' \itemize{
#' \item \code{\link{sim5pop}}
#' \item \code{\link{pop_labels}}
#' \item \code{\link{PCs}}
#' }
#' @keywords internal
#' @import kris
#' @references 
#' 
#' Intarapanich, A., Shaw, P.J., Assawamakin, A., Wangkumhang, P., Ngamphiw, C., 
#' Chaichoompu, K., Piriyapongsa, J., and Tongsima, S. (2009). Iterative pruning 
#' PCA improves resolution of highly structured populations. BMC Bioinformatics 
#' 10, 382.
#' 
#' Lebret, R., Iovleff, S., Langrognet, F., Biernacki, C., Celeux, G., and 
#' Govaert, G. (2015). Rmixmod: TheRPackage of the Model-Based Unsupervised, 
#' Supervised, and Semi-Supervised ClassificationMixmodLibrary. J. Stat. Softw. 
#' 67.
#' 
#' Limpiti, T., Intarapanich, A., Assawamakin, A., Shaw, P.J., Wangkumhang, P., 
#' Piriyapongsa, J., Ngamphiw, C., and Tongsima, S. (2011). Study of large and 
#' highly stratified population datasets by combining iterative pruning 
#' principal component analysis and structure. BMC Bioinformatics 12, 255.

"_PACKAGE"
#> [1] "_PACKAGE"

NULL