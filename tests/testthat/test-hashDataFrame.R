library(amsteRdam)

dfInp<-data.frame(col1=c(1,2,3),col2=c(4,5,6),col3=c(5,6,7),
                  col4=c(-1,-1,-1),col5=c(2,2,2),col6=c(3,3,3))

dfOut<-data.frame(H1=c(-2,-3,-4),H2=c(1,1,1),H3=c(-5,-7,-9))

test_that("Basic sanity check of data frame hashing on simple example", {
  expect_equal(hashDataFrame(dfInp,3,0),dfOut)
})
