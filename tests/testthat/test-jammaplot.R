# tests for jammaplot,ggjammaplot and input types

testthat::test_that("jammaplot DGEList input", {
   testthat::skip_if_not_installed("edgeR")

   # define test data
   set.seed(123)
   testmatrix <- matrix(1:15, ncol=3)
   colnames(testmatrix) <- head(LETTERS[1:3])
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- 0;
   testmatrix[3, 2] <- 0;
   testmatrix[5, ] <- 0;
   testmatrixDGE <- edgeR::DGEList(testmatrix)

   testjc <- jammaplot(testmatrix, doPlot=FALSE)
   testjcDGE <- jammaplot(testmatrixDGE, doPlot=FALSE)
   # tests
   testthat::expect_identical(testjc, testjcDGE)

   
   testggjc <- ggjammaplot(testmatrix, doPlot=FALSE)
   testggjcDGE <- ggjammaplot(testmatrixDGE, doPlot=FALSE)
   testthat::expect_identical(testggjc@data, testggjcDGE@data)

})

testthat::test_that("jammaplot SummarizedExperiment input", {
   testthat::skip_if_not_installed("SummarizedExperiment")

   # define test data
   set.seed(123)
   testmatrix <- matrix(1:15, ncol=3)
   colnames(testmatrix) <- head(LETTERS[1:3])
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- 0;
   testmatrix[3, 2] <- 0;
   testmatrix[5, ] <- 0;
   testmatrixSE <- SummarizedExperiment::SummarizedExperiment(
      assays=list(exprs=testmatrix))

   testjc <- jammaplot(testmatrix, doPlot=FALSE)
   testjcSE <- jammaplot(testmatrixSE, doPlot=FALSE)
   # tests
   testthat::expect_identical(testjc, testjcSE)

   
   testggjc <- ggjammaplot(testmatrix, doPlot=FALSE)
   testggjcSE <- ggjammaplot(testmatrixSE, doPlot=FALSE)
   testthat::expect_identical(testggjc@data, testggjcSE@data)

})

testthat::test_that("jammaplot DESeqDataSet input", {
   testthat::skip_if_not_installed("DESeq2")

   # define test data
   set.seed(123)
   testmatrix <- matrix(1:15, ncol=3)
   colnames(testmatrix) <- head(LETTERS[1:3])
   rownames(testmatrix) <- head(letters[1:5])
   testmatrix[2, 1] <- 0;
   testmatrix[3, 2] <- 0;
   testmatrix[5, ] <- 0;
   testmatrixDS <- DESeq2::DESeqDataSetFromMatrix(
      countData=testmatrix,
      design=~group,
      colData=data.frame(group=colnames(testmatrix)))

   testjc <- jammaplot(testmatrix, doPlot=FALSE)
   testjcDS <- jammaplot(testmatrixDS, doPlot=FALSE)
   # tests
   testthat::expect_identical(testjc, testjcDS)

   
   testggjc <- ggjammaplot(testmatrix, doPlot=FALSE)
   testggjcDS <- ggjammaplot(testmatrixDS, doPlot=FALSE)
   testthat::expect_identical(testggjc@data, testggjcDS@data)

   # Use the example DESeq2 data
   dds <- DESeq2::makeExampleDESeqDataSet(m=6)
   testjc1 <- jammaplot(DESeq2::counts(dds), doPlot=FALSE)
   testjc1DS <- jammaplot(dds, doPlot=FALSE)
   testthat::expect_identical(testjc1, testjc1DS)

})
