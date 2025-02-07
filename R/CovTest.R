#' @title Test for structure of data's covariance matrix
#'
#' @description This function conducts the test for the covariance matrix of
#' data regarding structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the
#' Anova-type-statistic(ATS) based on a specified number of runs.
#' @param X a matrix containing the observation vectors as columns
#' (one group only)
#' @param structure a character specifying the structure regarding the
#' covariance matrix should be checked. Options are "autoregressive" ("ar"),
#' "FO-autoregressive" ("FO-ar"), "diagonal" ("diag"), "sphericity" ("spher"),
#' "compoundsymmetry" ("cs") and "toeplitz" ("toep").
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the
#' predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen
#' method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed a seed, if it should be set for reproducibility. Predefined value
#' is NULL, which means no seed is set.
#' @return an object of the class \code{\link{CovTest}}
#'
#'
#' @references
#' Sattler, P. & Dobler, D. (2024) Testing for patterns and structures in covariance and correlation matrices
#'
#'
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' # Select only the males with the diagnosis AD
#' X <- as.matrix(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",
#'                           c("brainrate_temporal", "brainrate_frontal",
#'                           "brainrate_central","complexity_temporal",
#'                           "complexity_frontal", "complexity_central")])
#'
#' TestCovariance_structure(X = X, structure = "diagonal", method = "MC")
#'
#' @export
TestCovariance_structure <- function(X, structure, method = "BT",
                                     repetitions = 1000, seed = NULL){

  structure <- tolower(structure)
  method <- toupper(method)
  if(!(method == "MC" | method == "BT")){
    stop("method must be bootstrap ('BT') or Monte-Carlo-technique('MC')")
  }

  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }
  if(is.list(X)){
    if(length(X) > 1){
      warning("The input X must be a matrix but is a list. Only the first
              element of the list is used.")
      X <- X[[1]]
    }
    if(length(X) == 1){
      X <- X[[1]]
    }

  }

  n1 <- dim(X)[2]
  d <- dim(X)[1]
  if(d==1){ stop("Structures can be only investigated for more than one
                 dimension") }
  if(!(structure %in% c("autoregressive", "ar", "fo-autoregressive", "fo-ar",
                        "diagonal", "diag", "sphericity", "spher",
                        "compoundsymmetry", "cs", "toeplitz", "toep") )){
    stop("no predefined hypothesis")
  }

  if(d > 1){
    p <- d * (d + 1) / 2
    a <- cumsum(c(1, (d):2))

    vX <- dvech(stats::var(t(X)), a, d, p, inc_diag = TRUE)
    Xq <- apply(X - rowMeans(X), 2, vdtcrossprod, a, d, p)
    HatCov <- stats::var(t(Xq))

    if(structure == "autoregressive" | structure == "ar"){
      C <- diag(1,d,d)
      for(l in 2:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
      C <- matrixcalc::direct.sum(C, Pd(d - 1))
      Xi <- c(rep(1,d),rep(0, times = p - 1))
      Jacobi <- Jacobian(vX, a, d, p, 'subdiagonal_mean_ratio_fct')
      HatCovg <- QF(Jacobi, HatCov)

      if(method == "MC"){
        ResamplingResult <- ATSwS(QF(C, HatCovg), repetitions)
      }
      if(method == "BT"){
        ResamplingResult <- sapply(1:repetitions, Bootstrap_trans, n1, a, d, p,
                                   C, MSroot(HatCov), vX, 'subdiagonal_mean_ratio_fct')
      }

      Teststatistic <- ATS(n1, subdiagonal_mean_ratio_fct(vX, a, d), C, HatCovg, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)
    }


    if(structure == "fo-autoregressive" | structure == "fo-ar"){
      C <- Pd(d)
      for(l in 2:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
      C <- matrixcalc::direct.sum(C, Pd(d - 1))
      Xi <- rep(0, times = p + d - 1)
      Jacobi <- Jacobian(vX, a, d, p, 'subdiagonal_mean_ratio_fct')
      HatCovg <- QF(Jacobi, HatCov)

      if(method == "MC"){
        ResamplingResult <- ATSwS(QF(C, HatCovg), repetitions)
      }
      if(method == "BT"){
        ResamplingResult <- vapply(1:repetitions, Bootstrap_trans,
                                   FUN.VALUE = numeric(1),
                                   n1, a, d, p,
                                   C, MSroot(HatCov), vX,
                                   'subdiagonal_mean_ratio_fct')
      }
      Teststatistic <- ATS(n1, subdiagonal_mean_ratio_fct(vX, a, d), C,
                           HatCovg, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)
    }

    if(structure %in% c("diagonal", "diag", "sphericity", "spher",
                        "compoundsymmetry", "cs", "toeplitz", "toep")){

      Xi <- rep(0, p)

      if(structure == "diagonal" | structure == "diag"){
        C <- matrixcalc::direct.sum(matrix(0, d, d), diag(1, p - d, p - d))
      }
      if(structure == "sphericity" | structure == "spher"){
        C <- matrixcalc::direct.sum(Pd(d), diag(1, p - d, p - d))
      }
      if(structure == "compoundsymmetry" | structure == "cs"){
        C <- matrixcalc::direct.sum(Pd(d), Pd(p - d))
      }
      if(structure == "toeplitz" | structure == "toep"){
        C <- Pd(d)
        for(l in 2:d){
          C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
        }
      }

      if(method == "MC"){
        ResamplingResult <- ATSwS(QF(C, HatCov), repetitions)
      }
      if(method == "BT"){
        ResamplingResult <- vapply(1:repetitions, Bootstrap,
                                   FUN.VALUE = numeric(1),
                                   n1,
                                   C, MSroot(HatCov))
      }

      Teststatistic <- ATS(n1, vX, C, HatCov, Xi)
      pvalue <- mean(ResamplingResult > Teststatistic)

    }
  }
  CovTest <- list("method" = "Covariance",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = HatCov,
                  "C" = C,
                  "Xi" = Xi,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = structure,
                  "nv" = 1)

  class(CovTest) <- "CovTest"

  return(CovTest)
}


