#include <string>
#include <iostream>

#include <stdlib.h>
#include <Rcpp.h>

#include "MurmurHash3.h"

using namespace std;
using namespace Rcpp;

int calcHash(const char *s, int max , uint32_t seed)  {
    const void *p=(const void *)s;
    int hashingRes;
    MurmurHash3_x86_32(p,strlen(s),seed,&hashingRes);
    int res=(abs(hashingRes) % max)+1;
    return hashingRes>0 ? res : -res;
}

inline bool isNotDoubleNA(double x) {
  return x==x;
}

inline bool isNotIntNA(int x) {
  return x!=NA_INTEGER;
}

//' Hashes the data frame \code{df}, so the output contains \code{numCols}
//'
//' @param df data frame to be hashed
//' @param numCols integer, number of columns for the output data frame
//' @param hashSeed integer, seed for the employed hash function (MurmurHash3)
//' @export
// [[Rcpp::export]]
DataFrame hashDataFrame(DataFrame df, int numCols, int hashSeed) {
    // Create empty numCols dataFrame and fill it with zeros

    int nrows=df.nrows();
    List resList(numCols);
    for(int i=0; i<numCols; i++) {
        resList[i]=NumericVector(nrows);
    }
    CharacterVector colNames(numCols);

    for(int i=0; i<numCols; i++) {
        stringstream ss;
        ss << "H" <<  (i+1);
        colNames[i]=ss.str();
    }

    // Loop over columns, add substract each column depending on hash

    CharacterVector cnames=df.attr("names");
    for(CharacterVector::iterator iter=cnames.begin(); iter!=cnames.end(); iter++) {
        int type=0;
        int h=calcHash(*iter, numCols, hashSeed);
        int targetCol=abs(h);
        double sign=h<0 ? -1.0 : 1.0;
        string s=as<string>(*iter);

      if (TYPEOF(df[s])==REALSXP) {
            NumericVector x=df[s];
        for(int i=0; i<nrows; i++) {
          if (isNotDoubleNA(x[i])) {
            as<NumericVector>(resList[targetCol-1])[i]+=sign*x[i];
            }
        }
      }
      else if (TYPEOF(df[s])==INTSXP) {
        IntegerVector x=df[s];
        for(int i=0; i<nrows; i++) {
          if (isNotIntNA(x[i])) {
             as<NumericVector>(resList[targetCol-1])[i]+=sign*x[i];
            }
        }
    }
    }

  DataFrame result(resList);
  result.attr("names") = colNames;
  return result;
}
