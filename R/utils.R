#' @title Centering matrix
#'
#' @description matrix Pd for testing equality of the d components of a vector
#' @param d a scalar, characterizing the matrix and set its dimension
#' @return a matrix
#'
#'
#' @noRd
Pd <- function(d){
  return( diag(1,d,d) - matrix(1/d,d,d) )
}


#' @title Quadratic form for vectors and matrices
#'
#' @param A,B matrices or vectors
#' @return a metrix or vector
#'
#' @noRd
QF <- function(A, B){
  return( A %*% B %*% t(A) )
}

#' @title Square root of a matrix
#'
#' @param X matrix
#' @return matrix
#'
#' @noRd
MSroot <- function(X){
  if(length(X) == 1){
    MSroot <- matrix(sqrt(X),1,1)
  }
  else{
    SVD <- svd(X)
    MSroot <- SVD$u %*% ( tcrossprod(sqrt(diag(SVD$d)), (SVD$v)) )
  }
  return(MSroot)
}



#' @title Diagonal vectorisation
#' @description  Diagonal vectorisation of the upper triangular matrix
#' @param X quadratic matrix which should be diagonalized
#' @param a vector containing the indices which belong to the diagonal of the
#' matrix
#' @param d dimension of the matrix which should be vectorised
#' @param p dimension of the vectorised matrix
#' @param inc_diag TRUE or FALSE: should the diagonal be included?
#'
#' @return vector
#'
#' @keywords internal
#' @export
dvech <- function(X, a, d, p, inc_diag){
  if(!matrixcalc::is.square.matrix(X)){
    stop("argument X is not a square numeric matrix")
  }
  else{
    E <- rep(X[1,d],p)
    for(i in 1:(d-1)){
      E[a[i]:(a[i+1]-1)] <- diag(X[1:(d-i+1),i:d])
    }
    # without the diagonal
    if(!inc_diag){
      E <- E[-(1:d)]
    }

    return(E)
  }
}

#' Vectorization of the upper triangular part of the matrix
#'
#' @param X
#'
#' @return vector
#'
#' @keywords internal
#' @export
vechp <- function(X){
  if(!matrixcalc::is.square.matrix(X)){
    stop("argument X is not a square numeric matrix")
  }

  return(as.vector(t(X)[!upper.tri(X,TRUE)]))
}


#' @title Weighted direct sums for lists
#' @description Hereby the matrices which are part of a list are multiplied with
#' the corresponding components of a matrix w, containing the weights. These,
#' now weighted matrices are put together to one larger block-diagonal matrix.
#'
#' @param X matrix
#' @param w weight matrix
#'
#' @return matrix
#'
#' @keywords internal
#' @export
WDirect.sumL <- function(X, w){
  groups <- length(X)
  if(groups == 1){
    Result <- X*w
  }
  else{
    Result <- matrixcalc::direct.sum(w[1]*X[[1]], w[2]*X[[2]])
    if(groups > 2){
      for(i in 3:groups){
        Result <-  matrixcalc::direct.sum(Result, w[i]*X[[i]])
      }
    }
  }
  return(Result)
}



#' @title Function to calculate vech(X t(X))
#'
#' @param X matrix
#' @return vector
#'
#' @export
#' @keywords internal
vtcrossprod <- function(X){
  return(matrixcalc::vech(tcrossprod(X,X)))
}

#' @title Function to calculate dvech(X t(X))
#'
#' @param X matrix
#' @param a indices that belong to the diagonal of the matrix
#' @param d dimension of the matrix
#' @param p dimension of the vectorised matrix
#'
#' @return vector
#'
#' @export
#' @keywords internal
vdtcrossprod <- function(X,a,d,p){
  return(dvech(tcrossprod(X,X),a,d,p, inc_diag = TRUE))
}



#' @title Auxiliary function to calculate the covariance of the vectorized
#' correlation matrix
#'
#' @param X matrix
#' @param n number of columns
#'
#' @return matrix
#'
#' @export
#' @keywords internal
Qvech <- function(X, n){
  return(matrix(apply(X,2,vtcrossprod), ncol=n))
}


#' @title Root transformation of the vectorised correlation matrix
#'
#' @description A function calculating the roots of a vectorised covariance
#' matrix. The roots increasing, so square root for the first secondary
#' diagonals, third root for the second secondary diagonal and so on. For roots
#' with even order the absolute value of the argument is used, since the
#' arguments can be negative.
#'
#' @param v vectorised correlation matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' correlation matrix
#' @param d dimension of the correlation matrix
#' @return a transformed vector
#'
#' @keywords internal
#' @export
subdiagonal_mean_ratio_cor <- function(v, a, d){
  ratio <- rep(0, d - 2)
  ae <-  c(a, a[d] + 1)
  for(l in 3:(d)){
    ratio[l - 2] <-  mean(v[(ae[l]-d):(ae[l + 1] - 1-d)]) / mean(v[(ae[l - 1]-d):(ae[l] - 1-d)])
  }
  return(c(v,ratio))
}







#' Function for the Taylor-based Monte-Carlo-approximation for one group
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation for one group. After receiving some
#' auxiliary matrices and data, the Monte-Carlo observations are generated and
#' different parts of the final sum are defined. Based on this a number of the
#' Taylor-based ATS are calculated, where the number can be chosen.
#' @param repetitions a number specifying the number of runs for the
#' approximation
#' @param C the used hypothesis matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param CorData the calculated correlation matrix
#' @param MvrH an auxiliary matrix for the transformation from vectorised
#' covariances to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M4 a auxiliary matrix for the transformation from vectorised
#' covariances to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorised
#' covariances to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorised
#' covariances to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorised
#' covariances to vectorized correlations
#' @param n1 the total sample size, a scalar
#' @param Atilde an auxiliary matrix for the transformation from row-wise
#' vectorisation
#' @param Jacobi the Jacobian matrix of the transformation function applied
#' for the diagonal vectorised correlation  to diagonalwise vectorisation. used
#' for the transformation function 'subdiagonal_mean_ratio_cor', else NULL
#' @return a matrix containing the values of the Taylor ATS for a number of
#' repetitions
#'
#' @keywords internal
#' @export
Tayapp1G <- function(repetitions, C, MSrootStUpsi, CorData, MvrH, Trace, M4, L,
                     P, Q, n1, Atilde = NULL, Jacobi = NULL){
  vechCorData <- matrixcalc::vech(CorData)
  DvechCorDataM <- as.vector(vechCorData) * M4
  XPB <- generateData(MSrootStUpsi, repetitions)

  PUX <- P %*% XPB
  MUX0 <- M4 %*% Q %*% XPB

  HX <- Qvech(PUX, repetitions)

  Part1 <- MvrH %*% XPB
  Part2 <- 1 / (4) * as.vector(vechCorData) * HX
  Part3 <- 1 / (2) * XPB * MUX0
  Part4 <- 3 / (8) * DvechCorDataM %*% HX

  # without the transformation
  if(is.null(Atilde)){
    CXTay <- C %*% (Part1 + L %*% (Part2 - Part3 + Part4) / sqrt(n1)) / sqrt(n1)
    Result <- n1 * apply(CXTay, 2, crossprod) / Trace
  }
  # with the transformation
  else{
    XTaydv <- Atilde %*% (Part1 + L%*%(Part2 - Part3 + Part4)/sqrt(n1))/sqrt(n1)
    CXTaydv <- C %*% Jacobi %*% XTaydv
    Result <- n1 * apply(CXTaydv, 2, crossprod) / Trace
  }
  return(Result)
}


#' @title  Transformation of the vectorised covariance matrix by
#' quotients of means
#'
#' @description A function which calculates the mean of the secondary diagonals
#' and divide them through the next one. Since the elements can be negative, for
#' the denominator absolute values are used.
#' @param v vectorised covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#'
#' @keywords internal
#' @export
subdiagonal_mean_ratio_fct <- function(v, a, d){
  ratio <- rep(0, d - 1)
  ae <-  c(a, a[d] + 1)
  for(l in 2:(d)){
    ratio[l - 1] <-  mean(v[ae[l]:(ae[l + 1] - 1)]) /
      mean(v[ae[l - 1]:(ae[l] - 1)])
  }
  return(c(v, ratio))
}


#' @title Jacobian matrix for transformation functions
#'
#' @description A function which calculates the Jacobian matrix for a given
#' transformation function \code{\link{subdiagonal_mean_ratio_cor}} or
#' \code{\link{subdiagonal_mean_ratio_fct}}
#' @param X vectorised covariance matrix for which the Jacobian matrix is
#' applied
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param fun transformation function, that should be used.
#' \code{\link{subdiagonal_mean_ratio_fct}} or \code{\link{subdiagonal_mean_ratio_cor}}
#' @return the Jacobian matrix applied for the given vector
#'
#' @keywords internal
#' @export
Jacobian <- function(X, a, d, p, fun){
  if(fun == "subdiagonal_mean_ratio_fct"){
    J <-  matrix(0, d - 1, p)
    for(l in 1:(d - 1)){
      S1 <- sum((X[a[l]:(a[l] + d - l)]))
      S2 <- sum(X[a[l + 1]:(a[l + 1] + d - l - 1)])

      J[l, a[l + 1] + 0:(d - l - 1)] <- (d - l + 1) / (d - l) / S1
      J[l, a[l] + 0:(d - l)] <- (d - l + 1) / (d - l)  * (-S2 / (S1^2))
    }
    return(rbind(diag(1, p, p), J))
  }
  else{
    if(fun == "subdiagonal_mean_ratio_cor"){
      J <-  matrix(0, d - 2, p-d)
      for (l in 2:(d - 1)){
        S1 <- sum((X[(a[l]-d):(a[l] - l)]))
        S2 <- sum(X[(a[l + 1]-d):(a[l + 1]- l - 1)])

        J[l-1, (a[l + 1] + 0:(d - l - 1))-d] <- (d - l + 1) / (d - l) / S1
        J[l-1, (a[l] + 0:(d - l))-d] <- (d - l + 1) / (d - l)  * (-S2 / (S1^2))
      }
      return(rbind(diag(1, p-d, p-d), J))

    }

    else{
      stop("fun must be 'subdiagonal_mean_ratio_fct', 'subdiagonal_mean_ratio_cor'")
    }
  }

}


#' @title ATS for transformed vectors
#'
#' @description A function which calculates the Anova-type-statistic based on
#' a transformation function
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param v vectorised empirical covariance matrix of the original data
#' @param a vector containing the indices which belongs to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param fun transformation function, that should be used.
#' \code{\link{subdiagonal_mean_ratio_fct}} or \code{\link{subdiagonal_mean_ratio_cor}}
#' @return a scalar, the value of the ATS
#'
#' @keywords internal
#' @export
ATS_fun <- function(N, X, C, v, a, d, p, fun){
  Xmean <-  rowMeans(X)
  CDiff <- C %*% (do.call(fun, list(Xmean, a, d)) - do.call(fun, list(v, a, d)))
  Jacobi <-  Jacobian(Xmean, a, d, p, fun)
  HatCov <-  stats::var(t(X))
  Trace <-  sum(diag(QF(C %*% Jacobi, HatCov)))
  return(c(N * crossprod(CDiff) / Trace))
}


#' @title Bootstrap using transformation for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix HatCov, which matrix root MSrootHatCov is given and
#' expectation vector vX. For the generated bootstrap sample the value of the
#' ATS based on a transformation is calculated
#' @param N.sim control variable for using sapply
#' @param n1 a scalar, declaring the sample size for the bootstrap sample
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param C a hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix root of the covariance matrix HatCov, to generate
#' the bootstrap sample
#' @param vX the expectation vector for the bootstrap sample
#' @param fun transformation function, that should be used.
#' \code{\link{subdiagonal_mean_ratio_fct}} or
#' \code{\link{subdiagonal_mean_ratio_cor}}
#' @return a scalar, the value of the ATS
#'
#' @keywords internal
#' @export
Bootstrap_trans <- function(N.sim, n1, a, d, p, C, MSrootHatCov, vX, fun){
  XPB <- generateData(MSrootHatCov, n1) + vX
  return(ATS_fun(n1, XPB, C, vX, a, d, p, fun))
}


#' @title Bootstrap for one and multiple groups
#'
#' @description This function generates normal distributed random vectors.
#' For one group, nv random vectors with covariance matrix HatCov are generated
#' and the corresponding value of the ATS is generated. For multiple groups the
#' corresponding sample sizes from nv are used.
#' The weighted sum of covariance matrices is calculated and used to calculate
#' the value of the ATS.
#' @param N.sim control variable for using sapply
#' @param nv scalar (one group) or vector (multiple groups) of sample sizes for
#' the bootstrap samples
#' @param C hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix (one group) or list of matrices (multiple groups)
#' of roots of the covariance matrices, to generate
#' the bootstrap sample
#' @return a scalar, the value of the ATS
#'
#' @keywords internal
#' @export
Bootstrap <- function(N.sim, nv, C, MSrootHatCov){
  # one group
  if(length(nv) == 1){
    XPB <- generateData(MSrootHatCov, nv)
    PBHatCov <- stats::var(t(XPB))
    return(ATS(nv, rowMeans(XPB), C, PBHatCov))
  }
  # multiple groups
  else{
    N <- sum(nv)
    kappainvv <- N / nv

    DataPB <- mapply(generateData, MSrootHatCov, nv, SIMPLIFY = FALSE)
    PBHatCov <- WDirect.sumL(lapply(DataPB,
                                    function(X) stats::var(t(X))), kappainvv)
    return(ATS(N, unlist(lapply(DataPB, rowMeans)), C, PBHatCov))
  }
}




#' Anova-Type-Statistic with weighted sum
#'
#' @description
#' Calculation of a Anova-Type-Statistic using
#'
#' @param A a matrix
#' @param repetitions a scalar, number of runs
#' @return a vector of the length of repetitions
#'
#' @export
#' @keywords internal
ATSwS <- function(A, repetitions){
  Chi <- matrix(stats::rchisq(dim(A)[1] * repetitions, df = 1),
                ncol = repetitions)
  return(colSums(crossprod(eigen(A, only.values = 1)$value, Chi))/sum(diag(A)))
}


#' Anova-Type-statistic
#'
#' @param N number of observations
#' @param vVarData a matrix of vectorized covariance/correlation data
#' @param C hypothesis matrix for calculating the ATS
#' @param HatCov covariance matrix
#' @param Xi a vector defining together with C the investigated hypothesis
#'
#' @return a vector
#'
#' @export
#' @keywords internal
ATS <- function(N, vVarData, C, HatCov, Xi = 0){
  CDiff <-  C %*% vVarData - Xi
  statisticATS <-  N * crossprod(CDiff) / (sum(diag(QF(C, HatCov))))
  return(as.numeric(statisticATS))
}

#' Function to generate bootstrap observations
#'
#' @param WSigma weight matrix
#' @param nv number of observations in the groups
#'
#' @return a matrix
#'
#' @export
#' @keywords internal
generateData <- function(WSigma, nv){
  data <- WSigma %*% matrix(stats::rnorm(dim(WSigma)[1] * nv), ncol = nv)
  return(data)
}

