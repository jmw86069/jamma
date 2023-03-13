
context("centerGeneData")


testthat::test_that("centerGeneData with NA values", {
   # define test data
   testmatrix <- matrix(1:15, ncol=3)
   colnames(testmatrix) <- head(LETTERS[1:3])
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- NA;
   testmatrix[3, 2] <- NA;
   testmatrix[5, ] <- NA;

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

   centerGeneData(testmatrix2)
   centerGeneData(testmatrix2, useMedian=FALSE)
   centerGeneData(testmatrix2, useMedian=TRUE)

   testthat::expect_equal(
      centerGeneData(testmatrix2),
      answermatrix_6)

   testthat::expect_equal(
      centerGeneData(testmatrix2,
         useMedian=FALSE),
      answermatrix_7)

})
