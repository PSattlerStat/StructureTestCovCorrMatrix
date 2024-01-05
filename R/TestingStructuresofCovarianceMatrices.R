#' @title Centering matrix
#'
#' @description matrix Pd for testing equality of the d components of a vector
#' @param d a scalar, characterizing the matrix and set its dimension
#' @return the centering matrix
#' @export
Pd <- function(d) {
  return(diag(1, d, d) - matrix(1 / d, d, d))
}

#' @title Function to calculate vech(X t(X))
vtcrossprod <- function(X) {
  return(vech(tcrossprod(X, X)))
}

#' @title Auxiliary function to calculate the covariance of the vectorized correlation matrix
Qvech <- function(X, n) {
  return(matrix(apply(X, 2, vtcrossprod), ncol = n))
}

#' @title Quadratic form for vectors and matrices
#' @export
QF <- function(A, B) {
  return(A %*% B %*% t(A))
}


#' @title Weighted sum to calculate the asymptotic distribution of the ATS
#' @description Generating a number of realisations of weighted sums of chi-square
#' distributed random variables which is the asymptotic distribution of the ATS
#' @param A the matrix which is used to calculate the ATS
#' @param repetitions the number of realisations of the weighted sum
#' @return a vector containing the values of the weighted sum for a number of
#' repetitions
#' @export
ATSwS <- function(A, repetitions)
{
  Chi = matrix(rchisq(dim(A)[1] * repetitions, df = 1), ncol = repetitions)
  return(colSums(crossprod(eigen(A, only.values = 1)$value, Chi)) / sum(diag(A)))
}


#' @title Function to calculate the variance of transposed observations
tvar <- function(X) {
  return(var(t(X)))
}


#' @title Square root of a matrix
#' @export
MSroot <- function(X)
{
  if (length(X) == 1) {
    MSroot = matrix(sqrt(X), 1, 1)
  }
  else{
    SVD = svd(X)
    MSroot = SVD$u %*% (tcrossprod(sqrt(diag(SVD$d)), (SVD$v)))
  }
  return(MSroot)
}


#' @title Diagonal vectorisation including the diagonal
#' @description Diagonal vectorisation of the upper triangular part of a matrix
#' containing the diagonal elements.
#' @param X quadratic matrix which should be vectorized
#' @param a vector containing the indices which belong to the diagonal of the
#' matrix
#' @param d dimension of the matrix which should be vectorised
#' @param p dimension of the vectorised matrix
#' @export
dvech <- function(X, a, d, p)
{
  if (!is.square.matrix(X))
    stop("argument X is not a square numeric matrix")
  else{
    E = rep(X[1, d], p)
    for (i in 1:(d - 1))
    {
      E[a[i]:(a[i + 1] - 1)] = diag(X[1:(d - i + 1), i:d])
    }
    return(E)
  }
}



#' @title Anova-Type-statistic
ATS <- function(N, vVarData, C, HatCov, Xi = 0)
{
  CDiff = C %*% vVarData - Xi
  statisticATS = N * crossprod(CDiff) / (sum(diag(QF(C, HatCov))))
  Trace = (sum(diag(QF(C, HatCov))))
  return(as.numeric(statisticATS))
}


#' @title Function to generate normal distributed bootstrap observations
#' A function to generate nv normal distributed random vectors with expectation
#' zero and covariance matrix Sigma, while the matrix root of Sigma is given.
gData <- function(WSigma, nv)
{
  Data = WSigma %*% matrix(rnorm(dim(WSigma)[1] * nv), ncol = nv)
  return(Data)
}


#' @title Root transformation of the vectorised covariance matrix
#'
#' @description A function calculating the roots of a vectorised covariance matrix.
#' The roots increasing, so square root for the first secondary diagonals, third
#' root for the second secondary diagonal and so on. For roots with even order the
#' absolute value of the argument is used, since the arguments can be negative.
#' @param x vectorised covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#' @export
h <- function(x, a, d)
{
  for (i in 3:d)
  {
    if ((i %% 2 == 1)) {
      x[(0:(d - i)) + a[i]] = abs(x[(0:(d - i)) + a[i]]) ^ (1 / (i - 1))
    }
    if ((i %% 2 == 0)) {
      x[(0:(d - i)) + a[i]] = (x[(0:(d - i)) + a[i]] <= 0) * (-abs(x[(0:(d - i)) +
                                                                       a[i]]) ^ (1 / (i - 1))) + (x[(0:(d - i)) + a[i]] > 0) * (abs(x[(0:(d - i)) +
                                                                                                                                        a[i]]) ^ (1 / (i - 1)))
    }
  }
  return(x)
}




#' @title  Transformation of the vectorised covariance matrix by quotients of means
#'
#' @description A function which calculates the mean of the secondary diagonals
#' and divide them through the next one.Since the elements can be negative, for
#' the denominator absolute valued are used.
#' @param x vectorised covariance matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#' @export
g <- function(v, a, d)
{
  Ratio = rep(0, d - 1)
  ae = c(a, a[d] + 1)
  for (l in 2:(d))
    Ratio[l - 1] = mean(v[ae[l]:(ae[l + 1] - 1)]) / mean(abs(v[ae[l - 1]:(ae[l] -
                                                                            1)]))
  return(c(v, Ratio))
}

#' @title Jacobian matrix for the function h
#'
#' @description A function which calculates the Jacobian matrix for the root
#' transformation h  applied for a given vector
#' @param x vectorised covariance matrix  for which the Jacobian matrix is applied
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return the Jacobian matrix applied for the given vector
#' @export
Jacobianh <- function(x, a, d, p)
{
  E = rep(1, p)
  for (i in 3:d)
  {
    if ((i %% 2 == 1)) {
      E[(0:(d - i)) + a[i]] = x[(0:(d - i)) + a[i]] / ((i - 1) * abs(x[(0:(d -
                                                                             i)) + a[i]]) ^ (2 - 1 / (i - 1)))
    }
    if ((i %% 2 == 0)) {
      E[(0:(d - i)) + a[i]] = 1 / ((i - 1) * abs(x[(0:(d - i)) + a[i]]) ^ (1 -
                                                                             (1 / (i - 1))))
    }
  }
  return(diag(E, p, p))
}


#' @title Jacobian matrix for the function g
#'
#' @description A function which calculates the Jacobian matrix for
#' transformation function based on quotients of means. This Jacobian matrix is
#' applied for a given vector
#' @param x vectorised covariance matrix for which the Jacobian matrix is applied
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return the Jacobian matrix applied for the given vector
#' @export
Jacobiang <- function(X, a, d, p)
{
  J = matrix(0, d - 1, p)
  for (l in 1:(d - 1))
  {
    S1 = sum(abs(X[a[l]:(a[l] + d - l)]))
    S2 = sum(X[a[l + 1]:(a[l + 1] + d - l - 1)])
    
    J[l, a[l + 1] + 0:(d - l - 1)] = (d - l + 1) / (d - l) / S1
    J[l, a[l] + 0:(d - l)] = (d - l + 1) / (d - l) * sign(X[a[l] + 0:(d - l)]) *
      (-S2 / (S1) ^ 2)
  }
  return(rbind(diag(1, p, p), J))
}


#' @title ATS for vectors transformed with the function h
#'
#' @description A function which calculates the Anova type statistic based on
#' the transformation h
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param v vectorised empirical covariance matrix of the original data
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return a scalar, the value of the ATS
#' @export
ATSh <- function(N, X, C, v, a, d, p)
{
  Xmean = rowMeans(X)
  CDiff = C %*% (h(Xmean, a, d) - h(v, a, d))
  HatCov = tvar(X)
  Jacobi = Jacobianh(Xmean, a, d, p)
  JHatCov = QF(Jacobi, HatCov)
  Trace = sum(diag(QF(C, JHatCov)))
  return(c(N * crossprod(CDiff) / Trace))
}



#' @title ATS for vectors transformed with the function g
#'
#' @description A function which calculates the Anova type statistic based on
#' the transformation g
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param v vectorised empirical covariance matrix of the original data
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return a scalar, the value of the ATS
#' @export
ATSg <- function(N, X, C, v, a, d, p)
{
  Xmean = rowMeans(X)
  CDiff = C %*% (g(Xmean, a, d) - g(v, a, d))
  Jacobi = Jacobiang(Xmean, a, d, p)
  HatCov = tvar(X)
  Trace = sum(diag(QF(C %*% Jacobi, HatCov)))
  return(c(N * crossprod(CDiff) / Trace))
}


#' @title Bootstrap for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix HatCov, which matrix root MSrootHatCov is given. For
#' the generated bootstrap sample the value of the ATS is calculated
#' @param N.sim control variable for using sapply
#' @param n1 a scalar, declaring the sample size for the bootstrap sample
#' @param C a hypothesis matrix for calculating the ATS
#' @param MSrootHatCov matrix root of the covariance matrix HatCov, to generate
#' the bootstrap sample
#' @return a scalar, the value of the ATS
#' @export
Bootstrap <- function(N.sim, n1, C, MSrootHatCov)
{
  XPB = gData(MSrootHatCov, n1)
  PBHatCov = tvar(XPB)
  return(ATS(n1, rowMeans(XPB), C, PBHatCov))
}



#' @title Bootstrap using the transformation h for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix HatCov, which matrix root MSrootHatCov is given and
#' expectation vector vX. For the generated bootstrap sample the value of the
#' ATS based on transformation h is calculated
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
#' @return a scalar, the value of the ATS
#' @export
Bootstraph <- function(N.sim, n1, a, d, p, C, MSrootHatCov, vX)
{
  XPB = gData(MSrootHatCov, n1) + vX
  return(ATSh(n1, XPB, C, vX, a, d, p))
}


#' @title Bootstrap using the transformation g for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix HatCov, which matrix root MSrootHatCov is given and
#' expectation vector vX. For the generated bootstrap sample the value of the
#' ATS based on transformation g is calculated
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
#' @return a scalar, the value of the ATS
#' @export
Bootstrapg <- function(N.sim, n1, a, d, p, C, MSrootHatCov, vX)
{
  XPB = gData(MSrootHatCov, n1) + vX
  return(ATSg(n1, XPB, C, vX, a, d, p))
}



#' @title Test for structure of data's covariance matrix
#'
#' @description With this function the covariance matrix of data can be checked
#' for one of the usual structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns
#' @param structure a character specifying the structure regarding them the
#' covariance matrix should be checked. Options are "autoregressive", "FO-autoregressive"
#' "diagonal", "sphericity", "compoundsymmetry" and "toeplitz".
#' @param method a character, to chose whether bootstrap("BT") or
#' Monte-Carlo-technique("MC") is used, while bootstrap is the predefined method.
#' @param repetitions a scalar,  indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X=matrix(rnorm(5*20),5,20)
#' TestCovStructure(X,structure="toeplitz",method="MC")
#' @export
TestCovStructure <-
  function(X,
           structure,
           method = "BT",
           repetitions = 1000,
           seed = 0)
  {
    if (seed != 0)
      (set.seed(seed))
    n1=dim(X)[2]
    d=dim(X)[1]
    if(d==1){stop("Structures can be only investigated for more than one dimension")}
    if((1-(structure %in% c("autoregressive", "FO-autoregressive", "diagonal", "sphericity", "compoundsymmetry", "toeplitz") ))==1){stop("no predefined hypothesis")}
    else{if(d>1){
      p=d*(d+1)/2#dimension vectorized covariance matrix
      a=cumsum(c(1,(d):2))
      
      vX=dvech(tvar(X),a,d,p)
      Xq=apply(X-rowMeans(X),2,function(X,a,d,p){ return(dvech(tcrossprod(X,X),a,d,p))},a,d,p)
      HatCov=tvar(Xq)
      if(structure %in% c("autoregressive","FO-autoregressive")){
        if(structure=="autoregressive")
        {C=direct.sum(diag(1,d,d),Pd(p-d))
        Xi=c(rep(1,times=d),rep(0,times=p-d))
        Jacobi=Jacobianh(vX,a,d,p)
        HatCovh=QF(Jacobi,HatCov)
        
        if(method=="MC"){ResamplingResult=ATSwS(QF(C,HatCovh),repetitions)}
        if(method=="BT"){
          ResamplingResult=sapply(1:repetitions,Bootstraph,n1,a,d,p,C,MSroot(HatCov),vX)}
        Teststatistic=ATS(n1,h(vX,a,d),C,HatCovh,Xi)}
        if(structure=="FO-autoregressive")
        {C=Pd(d)
        for (l in 2:d){C=direct.sum(C,Pd(d-l+1))}
        C=direct.sum(C,Pd(d-1));  Xi=rep(0,times=p+d-1)
        Jacobi=Jacobiang(vX,a,d,p)
        HatCovg=QF(Jacobi,HatCov)
        
        if(method=="MC"){ResamplingResult=ATSwS(QF(C,HatCovg),repetitions)}
        if(method=="BT"){
          ResamplingResult=sapply(1:repetitions,Bootstrapg,n1,a,d,p,C,MSroot(HatCov),vX)}
        Teststatistic=ATS(n1,g(vX,a,d),C,HatCovg,Xi)}}
      else
      {Xi=rep(0,p)
      if(structure=="diagonal")
      {C=direct.sum(matrix(0,d,d),diag(1,p-d,p-d)) }
      if(structure=="sphericity"){C=direct.sum(Pd(d),diag(1,p-d,p-d))}
      if(structure=="compoundsymmetry"){C=direct.sum(Pd(d),Pd(p-d))}
      if(structure=="toeplitz"){C=Pd(d)
      for (l in 2:d){C=direct.sum(C,Pd(d-l+1))}}
      
      if (method=="MC"){  ResamplingResult=ATSwS(QF(C,HatCov),repetitions)}
      if (method=="BT"){  ResamplingResult=sapply(1:repetitions,Bootstrap,n1,C,MSroot(HatCov))}
      
      Teststatistic=ATS(n1,vX,C,HatCov,Xi)}
      pvalue=mean(ResamplingResult<Teststatistic)
      if(seed!=0){set.seed(NULL)}
      return(list("pvalue"=pvalue,"Teststatistic"=Teststatistic,"CovarianceMatrix"=HatCov))}}}
