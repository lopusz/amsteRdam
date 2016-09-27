library(amsteRdam)

expect_not_equal<-function(x,y) {
  expect_false(isTRUE(all.equal(x, y)))
}

dfInp<-data.frame(col1=c(1,2,3),col2=c(4,5,6),col3=c(5,6,7),
                  col4=c(-1,-1,-1),col5=c(2,2,2),col6=c(3,3,3))

dfOut<-data.frame(H1=c(-2,-3,-4),H2=c(1,1,1),H3=c(-5,-7,-9))

test_that("Basic sanity check of data frame hashing on simple example", {
  expect_equal(hashDataFrame(dfInp,3,0),dfOut)
})

dfInpNADouble<-dfInp
dfInpNADouble$col7<-c(NA,NA,NA)
dfInpNADouble$col8<-c(NA,NA,NA)
dfInpNADouble$col9<-c(NA,NA,NA)

dfInpNAInt<-dfInp
dfInpNAInt$col7<-as.integer(c(NA,NA,NA))
dfInpNAInt$col8<-as.integer(c(NA,NA,NA))
dfInpNAInt$col9<-as.integer(c(NA,NA,NA))

test_that("Sanity check for missing values (NA)", {
  expect_equal(hashDataFrame(dfInpNADouble,3,0),dfOut)
  expect_equal(hashDataFrame(dfInpNAInt,3,0),dfOut)
})

dfInpFactor<-cbind(dfInp,col7=factor(1:(nrow(dfInp))))

test_that("Basic sanity check of data frame hashing with factors", {
  expect_equal(hashDataFrame(dfInpFactor,3,0,factorMode="skip"),dfOut)
  # Complicated way of expressing expect_not_equal
  expect_false(isTRUE(all.equal(hashDataFrame(dfInpFactor,3,0,factorMode="integer"),
                   hashDataFrame(dfInpFactor,3,0,factorMode="oneHot"))))
})
