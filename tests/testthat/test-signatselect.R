context("test-tsinfer")

test_that("tsinfer works", {
  
  library(signatselect)
  
  tvec = c(0, 10, 20)
  bvec = c(2000, 4000, 6000)
  nvec = c(10000, 10000, 10000)
  
  tsinfer_test_output <- 
  tsinfer(
    tvec = tvec,
    bvec = bvec,
    nvec = nvec,
    verbose = FALSE
  )
  
  expect_equal(tsinfer_test_output$s, 0.08925489, tolerance = .002)
  expect_equal(tsinfer_test_output$LL.0, 18.04008, tolerance = .002)
  expect_equal(tsinfer_test_output$alpha, 16436.478, tolerance = .002)
})

test_that("fit works", {
  
  library(signatselect)
  
  # data slightly modded from Feder et al. Table2 so we get a p<0.05
  time <- c(415 , 505 , 585 , 665 , 745 , 825 , 910)
  freq <- c(0.06956522, 0.23125000, 0.62352941, 0.78494624, 0.93333333, 0.97979798, 0.98979592)
  
  fit_test_output <- 
    fit(
      time = time,
      v = freq
    )
  
  expect_equal(fit_test_output$fit_stat, 3.262457, tolerance = .002)
  expect_equal(fit_test_output$fit_p, 0.02238466, tolerance = .002)
})
