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


#' @title Upper triangular vectorisation without the diagonal
#' @description Row-wise vectorisation of the upper triangular part of a matrix without the
#' diagonal elements.
#' @export
vechp <- function(x)
{
  if (!is.square.matrix(x))
    stop("argument x is not a square numeric matrix")
  return(as.vector(t(x)[!upper.tri(x, TRUE)]))
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


#' @title Diagonal vectorisation excluding the diagonal
#' @description  Diagonal vectorisation of the upper triangular part of a matrix
#' without the diagonal elements.
#' @param X quadratic matrix which should be diagonalized
#' @param a vector containing the indices which belong to the diagonal of the
#' matrix
#' @param d dimension of the matrix which should be vectorised
#' @param pu dimension of the vectorised matrix
#' @export
dvechp <- function(X, a, d, pu)
{
  if (!is.square.matrix(X))
    stop("argument X is not a square numeric matrix")
  else {
    E = rep(X[1, d], pu)
    for (i in 2:(d - 1))
    {
      E[(a[i] - d):(a[i + 1] - 1 - d)] = diag(X[1:(d - i + 1), i:d])
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


#' Function for the Taylor-based Monte-Carlo-approximation for one group
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation for one group.  After receiving some
#' auxiliary matrices and data, the Monte-Carlo observations are generated and
#' different parts of the final sum are defined. Based on this a number of the
#' Taylor-based ATS are calculated, where the number can be chosen.
#' @param repetitions a number specifying the number of runs for the approximation
#' @param C the used hypothesis matrix
#' @param p the dimension of the vectorized matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param Cordata the calculated correlation matrix
#' @param MvrH1 an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Atilde an auxiliary matrix for the transformation from row-wise vectorisation
#' to diagonalwise vectorisation
#' @param n1 the sample size, a scalar
#' @return a matrix containing the values of the Taylor ATS for a number of repetitions
Tayapp <-
  function(repetitions,
           C,
           MSrootStUpsi,
           CorData,
           MvrH1,
           Trace,
           M,
           L,
           P,
           Q,
           Atilde,
           n1)
  {
    vechCorData = vech(CorData)
    DvechCorDataM = as.vector(vechCorData) * M
    XPB = gData(MSrootStUpsi, repetitions)

    PUX = P %*% XPB
    MUX0 = M %*% Q %*% XPB

    HX = Qvech(PUX, repetitions)

    Part1 = MvrH1 %*% XPB
    Part2 = 1 / (4) * as.vector(vechCorData) * HX
    Part3 = 1 / (2) * XPB * MUX0
    Part4 = 3 / (8) * DvechCorDataM %*% HX
    CXTaydv = C %*% Atilde %*% (Part1 + L %*% (Part2 - Part3 + Part4) / sqrt(n1)) /
      sqrt(n1)

    Result = n1 * apply(CXTaydv, 2, crossprod) / Trace
    return(Result)
  }


#' @title ATS for vectors transformed with the function htilde
#'
#' @description A function which calculates the Anova type statistic based on
#' the transformation htilde
#' @param N sample size
#' @param X matrix containing the bootstrap observations as columns
#' @param C the hypothesis matrix
#' @param r vectorised empirical covariance matrix of the original data
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return a scalar, the value of the ATS
#' @export
ATShtilde <- function(N, X, C, r, a, d, p)
{
  Xmean = rowMeans(X)
  CDiff = C %*% (htilde(Xmean, a, d) - htilde(r, a, d))

  Jacobi = Jacobianhtilde(Xmean, a, d, p)
  HatCov = tvar(X)
  JHatCov = QF(Jacobi, HatCov)
  Trace = sum(diag(QF(C, JHatCov)))
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


#' @title Bootstrap using the transformation htilde for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix Upsidv, which matrix root is given and expectation
#' vector vCorData. For the generated bootstrap sample the value of the
#' ATS based on transformation h is calculated.
#' @param N.sim control variable for using sapply
#' @param n1 a scalar, declaring the sample size for the bootstrap sample
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param C a hypothesis matrix for calculating the ATS
#' @param MSUpsidv matrix root of the covariance matrix Upsidv, to generate
#' the bootstrap sample
#' @param vCorData the expectation vector for the bootstrap sample
#' @return a scalar, the value of the ATS
#' @export
Bootstraphtilde <- function(N.sim, n1, a, d, p, C, MSUpsidv, vCorData)
{
  XPB = gData(MSUpsidv, n1) + vCorData
  return(ATShtilde(n1, XPB, C, vCorData, a, d, p))
}





#' @title Root transformation of the vectorised correlation matrix
#'
#' @description A function which calculates the roots of a vectorised covariance matrix, while
#' the roots increasing, so square root for the first secondary diagonals, third
#' root for the second secondary diagonal and so on. For roots with even order the
#' absolute value of the argument is used, since the arguments can be negative.
#' @param x diagonal vectorised correlation matrix which should be transformed
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @return a transformed vector
#' @export
htilde <- function(x, a, d)
{
  xt = x
  for (i in 3:d)
  {
    if ((i %% 2 == 1)) {
      xt[(0:(d - i)) + a[i - 1] - (i - 2)] = abs(x[(0:(d - i)) + a[i - 1] - (i -
                                                                               2)]) ^ (1 / (i - 1))
    }
    if ((i %% 2 == 0)) {
      xt[(0:(d - i)) + a[i - 1] - (i - 2)] = (x[(0:(d - i)) + a[i - 1] - (i -
                                                                            2)] <= 0) * (-abs(x[(0:(d - i)) + a[i - 1] - (i - 2)]) ^ (1 / (i - 1))) + (x[(0:(d -
                                                                                                                                                               i)) + a[i - 1] - (i - 2)] > 0) * (abs(x[(0:(d - i)) + a[i - 1] - (i - 2)]) ^
                                                                                                                                                                                                   (1 / (i - 1)))
    }
  }
  return(xt)
}

#' @title Jacobian matrix for the function htilde
#'
#' @description A function which calculates the Jacobian matrix for the root
#' transformation htilde applied for a given vector
#' @param x vectorised covariance matrix for which the Jacobian matrix is applied
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @return the Jacobian matrix applied for the given vector
#' @export
Jacobianhtilde <- function(x, a, d, p)
{
  E = rep(1, p - d)
  for (i in 3:d)
  {
    if ((i %% 2 == 1)) {
      E[(0:(d - i)) + a[i - 1] - (i - 2)] = x[(0:(d - i)) + a[i - 1] - (i - 2)] /
        ((i - 1) * abs(x[(0:(d - i)) + a[i - 1] - (i - 2)]) ^ (2 - 1 / (i - 1)))
    }
    if ((i %% 2 == 0)) {
      E[(0:(d - i)) + a[i - 1] - (i - 2)] = 1 / ((i - 1) * abs(x[(0:(d - i)) +
                                                                   a[i - 1] - (i - 2)]) ^ (1 - (1 / (i - 1))))
    }
  }
  return(diag(E, p - d, p - d))
}

#' @title Bootstrap using the transformation htilde for one group
#'
#' @description This function generates n1 normal distributed random vectors
#' with covariance matrix Upsidv, which matrix root is given and expectation
#' vector vCorData. For the generated bootstrap sample the value of the
#' ATS based on transformation h is calculated.
#' @param N.sim control variable for using sapply
#' @param n1 a scalar, declaring the sample size for the bootstrap sample
#' @param a vector containing the indices which belong to the diagonal of the
#' covariance matrix
#' @param d dimension of the covariance matrix
#' @param p dimension of the vectorised matrix
#' @param C a hypothesis matrix for calculating the ATS
#' @param MSUpsidv matrix root of the covariance matrix Upsidv, to generate
#' the bootstrap sample
#' @param vCorData the expectation vector for the bootstrap sample
#' @return a scalar, the value of the ATS
#' @export
Bootstraphtilde <- function(N.sim, n1, a, d, p, C, MSUpsidv, vCorData)
{
  XPB = gData(MSUpsidv, n1) + vCorData
  return(ATShtilde(n1, XPB, C, vCorData, a, d, p))
}




#' @title Test for the correlation matrix of data regarding their structure
#'
#' @description With this function the covariance matrix of data can be checked
#' for one of the usual structures. Depending on the chosen method a bootstrap or
#' Monte-Carlo-technique is used to calculate p-value of the ATS based on a
#' specified number of runs.
#' @param X  a matrix containing the observation vectors as columns
#' @param structure a character specifying the structure regarding them the
#' covariance matrix should be checked. Options are "Hautoregressive",
#' "diagonal", "Hcompoundsymmetry" and "Htoeplitz", where H stands for Heterogenous.
#' @param method a character, to chose whether bootstrap("BT"),
#'  Monte-Carlo-technique("MC") or a Taylor-based Monte-Carlos-approach("Tay")
#'  is used, while bootstrap is the predefined method.
#' @param repetitions a scalar, indicate the number of runs for the chosen method.
#' The predefined value is 1,000, and the number should not be below 500.
#' @param seed A seed, if it should be set for reproducibility. Predefined values
#' is 0, which means no seed is set. A chosen seed is deleted at the end.
#' @return a list containing the p-value, the value of the test statistic and the
#' value of the estimated covariance matrix used in the test
#' @examples
#' X=matrix(rnorm(5*20),5,20)
#' TestCorStructure(X,structure="Htoeplitz",method="MC")
#' @export
TestCorStructure <-
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
    if(d==1){stop("Correlation is only defined for dimension higher than one")}
    if((1-(structure %in% c("Hautoregressive", "diagonal", "Hcompoundsymmetry", "Htoeplitz") ))==1){stop("no predefined hypothesis")}
    else{if(d>1) {
      p=d*(d+1)/2
      pu=d*(d-1)/2
      a=cumsum(c(1,(d):2))
      H=matrix(rep(a,d),d,d,byrow=TRUE)
      L=diag(1,p,p)[-a,]
      M=matrix(0,p,p)
      Q=diag(as.vector(vech(diag(1,d,d))),p,p)
      for (i in 1:p)
      {M[i,c(vech(t(H))[i],vech(H)[i])]=1}
      M1=L%*%(M+Q)
      VarData=tvar(X)
      CorData=cov2cor(VarData)
      vCorData=dvechp(CorData,a,d,pu)
      Xq=matrix(apply(X-rowMeans(X),2,vtcrossprod),nrow=p,ncol=n1)
      HatCov=tvar(Xq)
      MvrH1=(L-1/2*vechp(CorData)*M1)
      MvrH2=sqrt(diag(as.vector(1/vtcrossprod(matrix(vech(VarData)[a]))),p,p))
      Atilde=matrix(0,pu,pu)
      for (l in 1:(d-1))
      {for (k in 1:(d-l))
      {Atilde[a[l]+k-l,a[k]-(k-l)]=1}}
      Upsidv=QF(Atilde%*%MvrH1%*%MvrH2,HatCov)

      Xi=rep(0,pu)
      if(structure=="Hautoregressive")
      { Jacobi=Jacobianhtilde(vCorData,a,d,p)
      Upsidvhtilde=QF(Jacobi,Upsidv)
      C=Pd(pu)
      Teststatistic=ATS(n1,htilde(vCorData,a,d),C,Upsidvhtilde,Xi)
      if(method=="MC"){ResamplingResult=ATSwS(QF(C,Upsidvhtilde),repetitions)}
      if(method=="BT"){ResamplingResult=sapply(1:repetitions,Bootstraphtilde,n1,a,d,p,C,MSroot(Upsidv),vCorData)}
      if(method=="Tay"){P=diag(1,p,p)[a,]
      StUpsi=QF(MvrH2,HatCov)
      Trace=sum(diag(QF(C,Upsidvhtilde)))

      ResamplingResult=Tayapphtilde(repetitions,C,MSroot(StUpsi),CorData,Jacobi,MvrH1,Trace,M,L,P,Q,Atilde,n1)}}


      else{
        if(structure=="diagonal"){C=diag(1,pu,pu)}
        if(structure=="Hcompoundsymmetry"){C=Pd(pu)}
        if(structure=="Htoeplitz"){C=Pd(d-1)
        for (l in 3:d)
        {C=direct.sum(C,Pd(d-l+1))}}
        Teststatistic=ATS(n1,vCorData,C,Upsidv,Xi)
        if (method=="MC"){ResamplingResult=ATSwS(QF(C,Upsidv),repetitions)}
        if (method=="BT"){ResamplingResult=sapply(1:repetitions,Bootstrap,n1,C,MSroot(Upsidv))}
        if (method=="Tay"){P=diag(1,p,p)[a,]
        StUpsi=QF(MvrH2,HatCov)
        Trace=sum(diag(QF(C,Upsidv)))
        ResamplingResult=Tayapp(repetitions,C,MSroot(StUpsi),CorData,MvrH1,Trace,M,L,P,Q,Atilde,n1)}    }
      pvalue=mean(ResamplingResult<Teststatistic)
      if(seed!=0){set.seed(NULL)}
      return(list("pvalue"=pvalue,"Teststatistic"=Teststatistic,"CovarianceMatrix"=Upsidv))}}}


#' Taylor-based Monte-Carlo-approximation with transformation htilde
#'
#' @description An auxiliary function to calculate the values for the
#' Taylor-based Monte-Carlo-approximation with the correlation transformation
#' htilde .  After receiving some auxiliary matrices and data, the Monte-Carlo
#' observations are generated and different parts of the final sum are defined.
#' Based on this a number of the  Taylor-based ATS are calculated, where the
#' number can be chosen.
#' @param repetitions a number specifying the number of runs for the approximation
#' @param C the used hypothesis matrix
#' @param MSrootStUpsi the matrix root of the covariance matrix for the Taylor
#' observations
#' @param Cordata the calculated correlation matrix
#' @param Jacobi the Jacobian matrix of the function hthilde applied for the
#' diagonal vectorised correlation
#' @param MvrH1 an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Trace a trace used in the ATS for the test statistic
#' @param M an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param L an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param P an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Q an auxiliary matrix for the transformation from vectorised covariances
#' to vectorized correlations
#' @param Atilde an auxiliary matrix for the transformation from row-wise vectorisation
#' to diagonalwise vectorisation
#' @param n1 the sample size, a scalar
#' @return a matrix containing the values of the Taylor ATS for a number of repetitions
Tayapphtilde <-
  function(repetitions,
           C,
           MSrootStUpsi,
           CorData,
           Jacobi,
           MvrH1,
           Trace,
           M,
           L,
           P,
           Q,
           Atilde,
           n1)
  {
    vechCorData = vech(CorData)
    DvechCorDataM = as.vector(vechCorData) * M
    XPB = gData(MSrootStUpsi, repetitions)

    PUX = P %*% XPB
    MUX0 = M %*% Q %*% XPB

    HX = Qvech(PUX, repetitions)

    Part1 = MvrH1 %*% XPB
    Part2 = 1 / (4) * as.vector(vechCorData) * HX
    Part3 = 1 / (2) * XPB * MUX0
    Part4 = 3 / (8) * DvechCorDataM %*% HX
    XTaydv = Atilde %*% (Part1 + L %*% (Part2 - Part3 + Part4) / sqrt(n1))
    CXTaydv = C %*% Jacobi %*% XTaydv
    Result = apply(CXTaydv, 2, crossprod) / Trace
    return(Result)
  }
