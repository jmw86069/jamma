

testthat::test_that("centerGeneData with NA values", {
   # define test data
   testmatrix <- matrix(1:15, ncol=3)
   colnames(testmatrix) <- head(LETTERS[1:3])
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- NA;
   testmatrix[3, 2] <- NA;
   testmatrix[5, ] <- NA;
   testmatrix

   # answer matrix
   answermatrix_1 <- testmatrix - c(6, 9.5, 8, 9, NA)
   answermatrix_2 <- testmatrix - c(1, NA, 3, 4, NA);
   answermatrix_3 <- testmatrix - c(1, 1, 3, 4, 1);
   answermatrix_4 <- testmatrix - c(1, 2, 3, 4, 2);
   answermatrix_5 <- testmatrix - c(1, 9.5, 3, 4, NA);

   # tests
   testthat::expect_equal(
      centerGeneData(testmatrix),
      answermatrix_1)

   testthat::expect_equal(
      centerGeneData(testmatrix,
         controlSamples="A"),
      answermatrix_2)

   testthat::expect_equal(
      centerGeneData(testmatrix,
         naControlAction="na",
         controlSamples="A"),
      answermatrix_2)

   testthat::expect_equal(
      centerGeneData(testmatrix,
         naControlAction="min",
         controlSamples="A"),
      answermatrix_3)

   testthat::expect_equal(
      centerGeneData(testmatrix,
         naControlAction="floor",
         naControlFloor=2,
         controlSamples="A"),
      answermatrix_4)

   testthat::expect_equal(
      centerGeneData(testmatrix,
         naControlAction="row",
         controlSamples="A"),
      answermatrix_5)
})

testthat::test_that("centerGeneData median or mean", {
   # define test data
   testmatrix <- matrix(1:15, ncol=3)
   colnames(testmatrix) <- head(LETTERS[1:3])
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- NA;
   testmatrix[3, 2] <- NA;
   testmatrix[5, ] <- NA;
   #
   testmatrix2 <- testmatrix;
   testmatrix2[,3] <- testmatrix[,3] * 2;

   # answer matrix
   answermatrix_6 <- testmatrix2 - c(6, 15.5, 14.5, 9, NA)
   answermatrix_7 <- testmatrix2 - c(9+2/3, 15.5, 14.5, 13+2/3, NA)

   # centerGeneData(testmatrix2)
   # centerGeneData(testmatrix2, useMedian=FALSE)
   # centerGeneData(testmatrix2, useMedian=TRUE)

   testthat::expect_equal(
      centerGeneData(testmatrix2),
      answermatrix_6)

   testthat::expect_equal(
      centerGeneData(testmatrix2,
         useMedian=FALSE),
      answermatrix_7)

})


testthat::test_that("centerGeneData naControlAction", {
   # define test data
   testmatrix <- matrix(1:25 + 20, ncol=5)
   colnames(testmatrix) <- paste0(rep(c("A", "B"), c(3, 2)), c(1, 2, 3, 1, 2))
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- NA;
   testmatrix[3, 2] <- NA;
   testmatrix[5, 1:3] <- NA;
   testmatrix
   controlSamples <- colnames(testmatrix)[1:3];

   expectmatrix1 <- testmatrix - c(26, 29.5, 28, 29, 20);
   expectmatrix1[5, 4:5] <- NA;
   testthat::expect_equal(
      centerGeneData(testmatrix, controlSamples=controlSamples),
      expectmatrix1)

   # naControlAction="floor"
   testthat::expect_equal(
      centerGeneData(testmatrix, controlSamples=controlSamples,
         naControlAction="row", naControlFloor=20),
      testmatrix - c(26, 29.5, 28, 29, 42.5))

   # naControlAction="min"
   testthat::expect_equal(
      centerGeneData(testmatrix, controlSamples=controlSamples,
         naControlAction="min", naControlFloor=20),
      testmatrix - c(26, 29.5, 28, 29, 26))

   # naControlAction="floor"
   testthat::expect_equal(
      centerGeneData(testmatrix, controlSamples=controlSamples,
         naControlAction="floor", naControlFloor=20),
      testmatrix - c(26, 29.5, 28, 29, 20))
})

testthat::test_that("centerGeneData SparseMatrix", {
   if (requireNamespace("Matrix", quietly=TRUE)) {
      # define test data
      testmatrix <- matrix(1:25 + 20, ncol=5)
      colnames(testmatrix) <- paste0(rep(c("A", "B"), c(3, 2)),
         c(1, 2, 3, 1, 2))
      rownames(testmatrix) <- head(letters[1:5])
      testmatrix[2, 1] <- NA;
      testmatrix[3, 2] <- NA;
      testmatrix[5, 1:3] <- NA;
      stestmatrix <- as(testmatrix, "CsparseMatrix")
      controlSamples <- colnames(testmatrix)[1:3];

      centerGeneData(stestmatrix, controlSamples=controlSamples)

      expectmatrix1 <- testmatrix - c(26, 29.5, 28, 29, 20);
      expectmatrix1[5, 4:5] <- NA;
      testthat::expect_equal(
         centerGeneData(testmatrix, controlSamples=controlSamples),
         expectmatrix1)

      # naControlAction="floor"
      testthat::expect_equal(
         centerGeneData(testmatrix, controlSamples=controlSamples,
            naControlAction="row", naControlFloor=20),
         testmatrix - c(26, 29.5, 28, 29, 42.5))

      # naControlAction="min"
      testthat::expect_equal(
         centerGeneData(testmatrix, controlSamples=controlSamples,
            naControlAction="min", naControlFloor=20),
         testmatrix - c(26, 29.5, 28, 29, 26))

      # naControlAction="floor"
      testthat::expect_equal(
         centerGeneData(testmatrix, controlSamples=controlSamples,
            naControlAction="floor", naControlFloor=20),
         testmatrix - c(26, 29.5, 28, 29, 20))
   }
})
