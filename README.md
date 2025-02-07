This package is to conduct tests regarding the structure of covariance and correlation matrices. The test as well as the used methods and possible hypotheses are described in the article „Testing for patterns and structures in covariance and correlation matrices“.

The key functions here are TestCovStructure(X,structure,method = "BT",repetitions = 1000, seed = 0) and TestCorStructure(X,structure, method = "BT", repetitions = 1000, seed = 0). The first argument X are the data, which has to be structured as a d x n matrix, while d is the dimension of the observation vectors and n is the number of observations. The second argument structure specifies regarding which structure the matrix should be tested. For covariance matrices possible choices are
"autoregressive", "FO-autoregressive", "diagonal", "sphericity", "compoundsymmetry" and "toeplitz"
while for correlation matrices they are "Hautoregressive", "diagonal", "Hcompoundsymmetry" and "Htoeplitz".
Possible methods are bootstrap("BT"), a Monte-Carlo approach("MC") and for correlations also a Taylor-based Monte-Carlo approach ("Tay"), while the predefined method is the bootstrap. It is also possible to choose the number of repetitions for the resampling approaches and to set a seed to allow reproducibility.

The functions return a p-value, the value of the used test statistic as well as the estimated covariance matrix used in the test.
