#' Raster QR decomposition
#'
#' @param s A RasterStack
#' @param r A correlation coefficient
#' @details
#' Uses QR decomposition to generate a new neutral landscape model (NLM) raster with a specified correlation to the original raster.
#' @examples
#' r1 <- terra::rast(x = matrix(runif(25), nrow = 5, ncol = 5))
#' r2 <- terra::rast(x = matrix(runif(25), nrow = 5, ncol = 5))
#' s <- c(r1, r2)
#' qr <- rasterQR(s, r = 0.5)
#' @import terra
#' @export
rasterQR <- function(s, r){
  # Extract values and bind
  vals.1 <- terra::values(s[[1]])
  vals.2 <- terra::values(s[[2]])
  m <- cbind(vals.1, vals.2)
  m.ctr  <- scale(m, center = TRUE, scale = FALSE)

  # Prepare matrix projection
  Id.m <- diag(length(m.ctr[, 1]))  # Identity matrix
  Q <- qr.Q(qr(m.ctr[, 1, drop = FALSE]))  # QR decomposition
  P <- tcrossprod(Q)  # Cross-product of matrix for projection
  xvals.o <- (Id.m - P) %*% m.ctr[ , 2]  # Make xvals orthogonal to m.ctr
  m.2  <- cbind(m.ctr[, 1], xvals.o)
  Y <- m.2 %*% diag(1 / sqrt(colSums(m.2^2)))  # Scale columns to length 1

  # Generate new NLM raster
  theta <- acos(r) # Calculate corresponding angle for desired correlation
  newNLM.v <- Y[, 2] + (1 / tan(theta)) * Y[, 1]  # Vector of new values
  newNLM.v <- (newNLM.v - min(newNLM.v)) / (max(newNLM.v) - min(newNLM.v))  # Scale new values from 0 to 1
  newNLM.m <- matrix(data = newNLM.v, nrow = nrow(s[[2]]), ncol = ncol(s[[2]]), byrow = TRUE)
  newNLM <- terra::rast(newNLM.m)
  return(newNLM)
}
