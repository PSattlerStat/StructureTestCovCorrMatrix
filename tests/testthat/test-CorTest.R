# Loading the dataset
data("EEGwide", package = "MANOVA.RM")

vars <- colnames(EEGwide)[1:6]
d <- 6
p <- d * (d + 1) / 2


X_list <- list(t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "AD", vars]),
               t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "MCI", vars]),
               t(EEGwide[EEGwide$sex == "M" &
                           EEGwide$diagnosis == "SCC", vars]),
               t(EEGwide[EEGwide$sex == "W" &
                           EEGwide$diagnosis == "AD", vars]),
               t(EEGwide[EEGwide$sex == "W" &
                           EEGwide$diagnosis == "MCI", vars]),
               t(EEGwide[EEGwide$sex == "W" &
                           EEGwide$diagnosis == "SCC", vars]))
X_matrix <- matrix(unlist(X_list), nrow = 6)
X <- X_list[[1]]
nv <- c(12, 27, 20, 24, 30, 47)



## Correlation Structure
test_that("TestCorrelation_structure teststatic", {
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "BT"
    )$Teststatistic,
    3.95677381790978
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "TAY"
    )$Teststatistic,
    3.95677381790978
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "MC"
    )$Teststatistic,
    3.95677381790978
  )

  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "BT"
    )$Teststatistic,
    68.9068261789833
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "TAY"
    )$Teststatistic,
    68.9068261789833
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "MC"
    )$Teststatistic,
    68.9068261789833
  )

  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "BT"
    )$Teststatistic,
    5.26105347405293
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "TAY"
    )$Teststatistic,
    5.26105347405293
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "MC"
    )$Teststatistic,
    5.26105347405293
  )

  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "BT"
    )$Teststatistic,
    4.88304727273834
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "TAY"
    )$Teststatistic,
    4.88304727273834
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "MC"
    )$Teststatistic,
    4.88304727273834
  )

})

test_that("TestCorrelation_structure pvalue", {
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.021
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "TAY",
      seed = 31415
    )$pvalue,
    0.104
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Har",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.02
  )

  expect_equal(
    TestCorrelation_structure(
      X = X_list[[2]],
      structure = "diag",
      method = "BT",
      seed = 31415
    )$pvalue,
    0
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "TAY",
      seed = 31415
    )$pvalue,
    0
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "diag",
      method = "MC",
      seed = 31415
    )$pvalue,
    0
  )

  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.007
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "TAY",
      seed = 31415
    )$pvalue,
    0.065
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Hcs",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.006
  )

  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.009
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "TAY",
      seed = 31415
    )$pvalue,
    0.124
  )
  expect_equal(
    TestCorrelation_structure(
      X = X_list[[1]],
      structure = "Htoep",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.006
  )

})


test_that("TestCorrelation_structure wrong method/hypothesis", {
  expect_error(
    TestCorrelation_structure(
      X = X,
      structure = "hcs",
      method = "abc",
      repetitions = 1000,
      seed = 31415
    )
  )
  expect_equal(
    TestCorrelation_structure(
      X = X,
      structure = "hcs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$Teststatistic,
    5.26105347405293
  )
  expect_error(TestCorrelation_structure(X = X, structure = "a"))
})

test_that("TestCorrelation_structure input list", {
  expect_warning(expect_equal(
    TestCorrelation_structure(
      X = X_list,
      structure = "hcs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.006
  ))
  expect_equal(
    TestCorrelation_structure(
      X = list(X),
      structure = "hcs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.006
  )
})

test_that("TestCorrelation_structure d=1", {
  expect_error(TestCorrelation_structure(X = X_list[[1]][1, 1:12, drop = FALSE],
                                         structure = "Har"))
})


## Covariance Structure
test_that("TestCovariance_structure teststatistics", {
  expect_equal(
    TestCovariance_structure(X, structure = "autoregressive",
                             method = "MC")$Teststatistic,
    2.14534594320388
  )
  expect_equal(
    TestCovariance_structure(X, structure = "ar", method = "BT")$Teststatistic,
    2.14534594320388
  )
  expect_equal(
    TestCovariance_structure(X, structure = "FO-autoregressive",
                             method = "MC")$Teststatistic,
    1.63857996457449
  )
  expect_equal(
    TestCovariance_structure(X, structure = "FO-ar",
                             method = "BT")$Teststatistic,
    1.63857996457449
  )
  expect_equal(
    TestCovariance_structure(X, structure = "diagonal",
                             method = "MC")$Teststatistic,
    4.88780263620412
  )
  expect_equal(
    TestCovariance_structure(X, structure = "diag",
                             method = "BT")$Teststatistic,
    4.88780263620412
  )
  expect_equal(
    TestCovariance_structure(X, structure = "sphericity",
                             method = "MC")$Teststatistic,
    3.10441188918666
  )
  expect_equal(
    TestCovariance_structure(X, structure = "spher",
                             method = "BT")$Teststatistic,
    3.10441188918666
  )
  expect_equal(
    TestCovariance_structure(X, structure = "compoundsymmetry",
                             method = "MC")$Teststatistic,
    1.65692260204835
  )
  expect_equal(
    TestCovariance_structure(X, structure = "cs", method = "BT")$Teststatistic,
    1.65692260204835
  )
  expect_equal(
    TestCovariance_structure(X, structure = "toeplitz",
                             method = "MC")$Teststatistic,
    1.63921322978212
  )
  expect_equal(
    TestCovariance_structure(X, structure = "toep",
                             method = "BT")$Teststatistic,
    1.63921322978212
  )
})

test_that("TestCovariance_structure pvalue", {
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "autoregressive",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.129
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "ar",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.142
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "FO-autoregressive",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.197
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "FO-ar",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.214
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "diagonal",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.017
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "diag",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.053
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "sphericity",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.033
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "spher",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.077
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "compoundsymmetry",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.177
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "cs",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.227
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "toeplitz",
      method = "MC",
      seed = 31415
    )$pvalue,
    0.18
  )
  expect_equal(
    TestCovariance_structure(
      X,
      structure = "toep",
      method = "BT",
      seed = 31415
    )$pvalue,
    0.232
  )
})

test_that("TestCovariance_structure wrong method", {
  expect_error(
    TestCovariance_structure(
      X = X,
      structure = "cs",
      method = "abc",
      repetitions = 1000,
      seed = 31415
    )
  )
  expect_equal(
    TestCovariance_structure(
      X = X,
      structure = "cs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$Teststatistic,
    1.65692260204835
  )
})

test_that("TestCovariance_structure input list", {
  expect_warning(expect_equal(
    TestCovariance_structure(
      X = X_list,
      structure = "cs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.177
  ))
  expect_equal(
    TestCovariance_structure(
      X = list(X),
      structure = "cs",
      method = "mc",
      repetitions = 1000,
      seed = 31415
    )$pvalue,
    0.177
  )
})


test_that("TestCovariance_structure wrong dimension / hypothesis", {
  expect_error(TestCovariance_structure(
    X = X_list[[1]][1, 1:12, drop = FALSE],
    structure = "cs",
    method = "BT"
  ))
  expect_error(TestCovariance_structure(
    X = X_list[[1]],
    structure = "a",
    method = "BT"
  ))
})


test_that("Jacobian wrong function", {
  X <- matrix(rnorm(25),5,5)
  d <- 4
  p <- d * (d + 1) / 2
  pu <- d * (d - 1) / 2
  a <- cumsum(c(1, (d):2))
  expect_no_error(Jacobian(X = X, a = a, d = d, p = p,
                           fun = "subdiagonal_mean_ratio_fct"))
  expect_error(Jacobian(X = X, a = a, d = d, p = p,
                        fun = "wrong_name"))
})

test_that("dvech square matrix", {
  d <- 4
  p <- d * (d + 1) / 2
  pu <- d * (d - 1) / 2
  a <- cumsum(c(1, (d):2))
  expect_no_error(dvech(X = matrix(rnorm(25),5,5), a = a, d = d, p = pu,
                        inc_diag = FALSE))
  expect_error(dvech(X = matrix(rnorm(50, 10,5)), a = a, d = d, p = pu,
                     inc_diag =  FALSE))
})

test_that("vechp no square matrix", {
  expect_error(vechp(X = matrix(rnorm(10,2,5))))
})

test_that("weighted directed sum one group", {
  expect_no_error(WDirect.sumL(matrix(5,1,1), 2))
})

test_that("print covtest", {
  X <- CovTest()
  expect_null(print(X))
})
