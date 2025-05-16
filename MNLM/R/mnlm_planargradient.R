#' Multiple neutral landscape models using planar gradients
#'
#' @param nlayers The number of NLMs to generate.
#' @param r The correlation coefficient between the first NLM and each successive NLM.
#' @param nrow The number of rows in the rasters.
#' @param ncol The number of columns in the rasters .
#' @param resolution The resolution of the rasters (default = 1).
#' @param direction Direction of the gradient (between 0 and 360 degrees). If unspecified the direction is determined randomly.
#' @param rescale If TRUE (default), raster values are scaled from 0 to 1.
#' @details
#' Generates multiple neutral landscape models with linear gradients.
#' The r argument can accept either a single value, in which case all NLMs produced will have the same correlation with the first layer, or a vector containing the desired correlation coefficients for each layer.
#' @examples
#' NLMs <- mnlm_planargradient(nlayers = 3, r = c(0.3, 0.6), ncol = 20, nrow = 20)
#' @references
#' Sciaini, M., Fritsch, M., Scherer, C., & Simpkins, C. E. (2018). NLMR and landscapetools: An integrated environment for simulating and modifying neutral landscape models in R. Methods in Ecology and Evolution, 9(11), 2240–2248. doi:10.1111/2041-210X.13076
#' Palmer, M.W. (1992) The coexistence of species in fractal landscapes. The American Naturalist, 139, 375 - 397.
#' @importFrom NLMR nlm_planargradient
#' @importFrom raster stack layerStats
#' @export
mnlm_planargradient <- function(nlayers = 2, r, ncol, nrow, resolution = 1, direction = NA, rescale = TRUE){
  if(length(r) == 1) r <- rep(r, nlayers - 1)

  for(i in 1:nlayers){
    nlm <- terra::rast(NLMR::nlm_planargradient(ncol = ncol, nrow = nrow, resolution = resolution, direction = direction, rescale = rescale))
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


