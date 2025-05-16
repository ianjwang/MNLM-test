#' Multiple neutral landscape models using Gaussian random fields
#'
#' Generates a set of neutral landscape models (NLMs) using spatially correlated random fields (Gaussian random fields) with defined correlations between NLMs
#'
#' @param nlayers The number of NLMs to generate
#' @param r The correlation coefficient between the first NLM and each successive NLM
#' @param ncol The number of columns in the rasters
#' @param nrow The number of rows in the rasters
#' @param resolution The resolution of the rasters (default = 1)
#' @param autocorr_range Maximum range (in raster units) of spatial autocorrelation (default = 10)
#' @param mag_var Magnitude of variation across the raster (default = 5)
#' @param nug Magnitude of variation in the scale of autocorr_range; smaller values generate more homogeneous rasters (default = 0.2)
#' @param mean Mean value over the field (default = 0.5)
#' @param user_seed Set seed for simulation (optional)
#' @param rescale If TRUE (default), raster values are scaled from 0 to 1
#' @details
#' Gaussian random fields are a collection of random numbers on a set of discrete coordinates (the NLM raster). The -----
#' @examples
#' NLMs <- mnlm_gaussianfield(nlayers = 3, r = c(0.3, 0.6), ncol = 20, nrow = 20)
#' @references
#' Sciaini, M., Fritsch, M., Scherer, C., & Simpkins, C. E. (2018). NLMR and landscapetools: An integrated environment for simulating and modifying neutral landscape models in R. Methods in Ecology and Evolution, 9(11), 2240–2248. doi:10.1111/2041-210X.13076
#' Kéry & Royle (2016) Applied Hierarchical Modeling in Ecology Chapter 20
#' @importFrom NLMR nlm_gaussianfield
#' @importFrom terra ext rast
#' @import RandomFields
#' @export
mnlm_gaussianfield <- function(nlayers = 2, r, ncol, nrow, resolution = 1, autocorr_range = 10, mag_var = 5, nug = 0.2, mean = 0.5, user_seed = NULL, rescale = TRUE){
  if(length(r) == 1) r <- rep(r, nlayers - 1)

  for(i in 1:nlayers){
    nlm <- terra::rast(NLMR::nlm_gaussianfield(ncol, nrow, resolution = resolution, autocorr_range = autocorr_range, mag_var = mag_var,
                                         nug = nug, mean = mean, user_seed = user_seed, rescale = rescale))
    if(i == 1) nlm.s <- nlm
    else nlm.s <- c(nlm.s, nlm)
  }

  for(j in 2:nlayers){
    s <- c(nlm.s[[1]], nlm.s[[j]])
    newNLM <- rasterQR(s, r[j-1])
    names(newNLM) <- paste0("nlm.", j)
    terra::ext(newNLM) <- terra::ext(nlm.s[[1]]) # Match extents
    nlm.s[[j]] <- newNLM
  }

  names(nlm.s[[1]]) <- "nlm.1"
  return(nlm.s)

}

