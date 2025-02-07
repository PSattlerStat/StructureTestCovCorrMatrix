#' @title Test for structure of data's correlation matrix
#'
#' @description With this function the correlation matrix of data can be checked
#' for one of the predefined structures. Depending on the chosen method a
#' bootstrap, the Taylor-based Monte-Carlo approach or  Monte-Carlo-technique
#' is used to calculate the p-value of the Anova-type-statistic(ATS) based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns (one group)
#' @param structure a character specifying the structure regarding them the
#' correlation matrix should be checked. Options are "Hautoregressive" ("Har"),
#' "diagonal" ("diag"), "Hcompoundsymmetry" ("Hcs") and "Htoeplitz" ("Hteop").
#' @param method a character, to chose whether bootstrap("BT") or Taylor-based
#' Monte-Carlo-approach("TAY") or Monte-Carlo-technique("MC") is used, while
#' bootstrap is the predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen
#' method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed a seed, if it should be set for reproducibility. Predefined value
#' is NULL, which means no seed is set.
#' @return an object of the class \code{\link{CovTest}}
#'
#' @references
#' Sattler, P. & Dobler, D. (2024) Testing for patterns and structures in covariance and correlation matrices
#'
#' @import MANOVA.RM
#' @examples
#' # Load the data
#' data("EEGwide", package = "MANOVA.RM")
#'
#' # Select only the males with the diagnosis AD
#' X <- as.matrix(EEGwide[EEGwide$sex == "W" & EEGwide$diagnosis == "AD",
#'              c("brainrate_temporal", "brainrate_frontal","brainrate_central",
#'              "complexity_temporal","complexity_frontal",
#'              "complexity_central")])
#'
#' TestCorrelation_structure(X = X, structure = "diagonal", method = "MC")
#'
#' @export
TestCorrelation_structure <- function(X, structure, method = "BT",
                                      repetitions = 1000, seed = NULL){
  if(!is.null(seed)){
    old_seed <- .Random.seed
    on.exit({ .Random.seed <<- old_seed })
    set.seed(seed)
  }
  structure <- tolower(structure)
  method <- toupper(method)
  if(!(method == "MC" | method == "BT" | method == "TAY")){
    stop("method must be bootstrap ('BT'), Monte-Carlo-technique('MC') or
         Taylor-based Monte-Carlos-approach('Tay')")
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

  if(d == 1){
    stop("Correlation is only defined for dimension higher than one")
  }

  if(!(structure %in% c("hautoregressive", "har", "diagonal", "diag" ,
                        "hcompoundsymmetry", "hcs", "htoeplitz", "htoep"))){
    stop("no predefined hypothesis")
  }

  p <- d * (d + 1) / 2
  pu <- d * (d - 1) / 2
  a <- cumsum(c(1, (d):2))
  H <- matrix(rep(a, d), d, d)
  L <- diag(1, p, p)[-a, ]
  M <- matrix(0, p, p)
  Q <- diag(as.vector(matrixcalc::vech(diag(1, d, d))), p, p)
  for(i in 1:p){
    M[i, c(matrixcalc::vech(t(H))[i], matrixcalc::vech(H)[i])] <- 1
  }
  M1 <- L %*% M
  VarData <- stats::var(t(X))
  CorData <- stats::cov2cor(VarData)
  vCorData <- dvech(CorData, a, d, p, inc_diag = FALSE)
  Xq <- matrix(apply(X - rowMeans(X), 2, vtcrossprod), nrow = p, ncol = n1)
  HatCov <- stats::var(t(Xq))
  MvrH1 <- (L - 1 / 2 * vechp(CorData) * M1)

  MvrH2 <- sqrt(diag(as.vector(1 / vtcrossprod(matrix(
    matrixcalc::vech(VarData)[a]))), p, p))
  Atilde <- matrix(0, pu, pu)
  for(l in 1:(d - 1)){
    for(k in 1:(d - l)){
      Atilde[a[l] + k - l, a[k] - (k - l)] <- 1
    }
  }

  Upsidv <- QF(Atilde %*% MvrH1 %*% MvrH2, HatCov)
  Xi <- rep(0, pu)

  if(structure == "hautoregressive" | structure == "har"){
    Jacobi <- Jacobian(vCorData, a, d, p, fun = "subdiagonal_mean_ratio_cor")
    Upsidvhtilde <- QF(Jacobi, Upsidv)
    C <- Pd(d - 1)
    for(l in 3:d){
      C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
    }
    C <- matrixcalc::direct.sum(C, Pd(d-2))
    Xi=rep(0,p-2)

    Teststatistic <- ATS(n1, subdiagonal_mean_ratio_cor(vCorData, a, d), C,
                         Upsidvhtilde, Xi)
    if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsidvhtilde),
                                                  repetitions) }
    if(method == "BT"){ ResamplingResult <- vapply(X = 1:repetitions,
                                                   FUN = Bootstrap_trans,
                                                   FUN.VALUE = numeric(1),
                                                   n1, a, d,
                                                   p, C, MSroot(Upsidv),
                                                   vCorData,
                                                   fun = "subdiagonal_mean_ratio_cor") }
    if(method == "TAY"){
      P <- diag(1, p, p)[a, ]
      StUpsi <- QF(MvrH2, HatCov)
      Trace <- sum(diag(QF(C, Upsidvhtilde)))

      ResamplingResult <- Tayapp1G(repetitions, C, MSroot(StUpsi), CorData,
                                   MvrH1, Trace, M, L, P, Q, n1, Atilde, Jacobi)
    }
  }
  else{
    if(structure == "diagonal" | structure == "diag"){
      C <- diag(1, pu, pu)
    }
    if(structure == "hcompoundsymmetry" | structure == "hcs"){
      C <- Pd(pu)
    }
    if(structure == "htoeplitz" | structure == "htoep"){
      C <- Pd(d - 1)
      for(l in 3:d){
        C <- matrixcalc::direct.sum(C, Pd(d - l + 1))
      }
    }
    Teststatistic <- ATS(n1, vCorData, C, Upsidv, Xi)

    if(method == "MC"){ ResamplingResult <- ATSwS(QF(C, Upsidv), repetitions) }
    if(method == "BT"){ ResamplingResult <- vapply(X = 1:repetitions,
                                                   FUN = Bootstrap,
                                                   FUN.VALUE = numeric(1),
                                                   n1,
                                                   C, MSroot(Upsidv)) }
    if(method == "TAY"){
      P <- diag(1, p, p)[a, ]
      StUpsi <- QF(MvrH2, HatCov)
      Trace <- sum(diag(QF(C, Upsidv)))
      ResamplingResult <- Tayapp1G(repetitions, C, MSroot(StUpsi), CorData,
                                   MvrH1, Trace, M, L, P, Q,n1)
    }
  }

  pvalue <- mean(ResamplingResult > Teststatistic)


  CovTest <- list("method" = "Correlation",
                  "pvalue" = pvalue,
                  "Teststatistic" = Teststatistic,
                  "CovarianceMatrix" = Upsidv,
                  "C" = C,
                  "Xi" = Xi,
                  "resampling_method" = method,
                  "repetitions" = repetitions,
                  "hypothesis" = structure,
                  "nv" = n1)

  class(CovTest) <- "CovTest"

  return(CovTest)

}
