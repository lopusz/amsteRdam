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
    // Below, there has to be +1 because sign transports the info about substraction
    // or addition of hashed feature. This info is needed also for 0th column.
    int res=(abs(hashingRes) % max)+1;
    return hashingRes>0 ? res : -res;
}

inline bool isNotDoubleNA(double x) {
  return x==x;
}

inline bool isNotIntNA(int x) {
  return x!=NA_INTEGER;
}

void processNumericVector(const char *colName, int numCols, int hashSeed, const NumericVector &x, List &resList) {
  int h=calcHash(colName, numCols, hashSeed);
  int targetCol=abs(h);
  double sign=h<0 ? -1.0 : 1.0;

  for(int i=0; i<x.length(); i++) {
    if (isNotDoubleNA(x[i])) {
      as<NumericVector>(resList[targetCol-1])[i]+=sign*x[i];
    }
  }
}

void processIntegerVector(const char *colName, int numCols, int hashSeed, const IntegerVector &x, List &resList) {
  int h=calcHash(colName, numCols, hashSeed);
  int targetCol=abs(h);
  double sign=h<0 ? -1.0 : 1.0;

  for(int i=0; i<x.length(); i++) {
    if (isNotIntNA(x[i])) {
      as<NumericVector>(resList[targetCol-1])[i]+=sign*x[i];
    }
  }
}

bool isFactor(const IntegerVector &x) {
  bool res=false;
  if (x.hasAttribute("class")) {
    string className=as<string>(x.attr("class"));
    res=(className=="factor");
  }
  return res;
}

void processFactorVector(const char *colName, int numCols, int hashSeed, const IntegerVector &x, List &resList) {

  // Create levels to columns mapping
  CharacterVector levels=as<CharacterVector>(x.attr("levels"));
  // for(CharacterVector::iterator it=levels.begin(); it!=levels.end(); it++) {
  //   Rcout << (*it) << endl;
  // }
  vector<int> levelToColumn(levels.length()+1);  // levels are numbered from 1
  vector<double> levelToSign(levels.length()+1);  // levels are numbered from 1

  for(int i=0; i<levels.length(); i++) {
    string totalName=string(colName)+string(levels[i]);
    int h=calcHash(totalName.c_str(), numCols, hashSeed);
    int targetCol=abs(h);
    int sign=h<0 ? -1.0 : 1.0;

    levelToColumn[i+1]=targetCol;
    levelToSign[i+1]=sign;
  }

  for(int i=0; i<x.length(); i++) {
    int level=x[i];
    if (isNotIntNA(level)) {
      int targetCol=levelToColumn[level];
      double sign=levelToSign[level];
      as<NumericVector>(resList[targetCol-1])[i]+=sign;
    }
  }
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

    // Loop over columns, add/substract each source column to target column depending on hash

    CharacterVector cnames=df.attr("names");
    for(CharacterVector::iterator colNameIter=cnames.begin(); colNameIter!=cnames.end(); colNameIter++) {
        string colName=as<string>(*colNameIter);

        if (TYPEOF(df[colName])==REALSXP) {
              processNumericVector(*colNameIter,numCols,hashSeed,df[colName],resList);
        }
        else if (TYPEOF(df[colName])==INTSXP) {
          IntegerVector x=df[colName];
          if (isFactor(x)) {
            processFactorVector(*colNameIter,numCols,hashSeed,x,resList);
          } else {
            processIntegerVector(*colNameIter,numCols,hashSeed,x,resList);
          }
      }
    }

  DataFrame result(resList);
  result.attr("names") = colNames;
  return result;
}
