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

  vector<int> levelToColumn(levels.length()+1);  // levels are numbered from 1
  vector<double> levelToSign(levels.length()+1);  // levels are numbered from 1

  for(int i=0; i<levels.length(); i++) {
    string totalName=string(colName)+":"+string(levels[i]);
    int h=calcHash(totalName.c_str(), numCols, hashSeed);
    int targetCol=abs(h);
    int sign=h<0 ? -1.0 : 1.0;

    levelToColumn[i+1]=targetCol;
    levelToSign[i+1]=sign;
  }

  // Map each value to precalculated column

  for(int i=0; i<x.length(); i++) {
    int level=x[i];
    if (isNotIntNA(level)) {
      int targetCol=levelToColumn[level];
      double sign=levelToSign[level];
      as<NumericVector>(resList[targetCol-1])[i]+=sign;
    }
  }
}

//' Hashes the data frame, so the output contains much smaller number of columns.
//'
//' Calculates the reduced representation of the given data frame based on the
//' hashing trick.
//'
//' @param df Data frame to be hashed.
//' @param numCols Integer, number of columns for the output data frame.
//' @param hashSeed Integer, seed for the employed hash function (MurmurHash3).
//'
//' @return Data frame having \code{numCols} columns, containing combination of
//' columns from \code{df} according to hashing function.
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

    CharacterVector inpColNames=df.attr("names");
    for(CharacterVector::iterator inpColNameIter=inpColNames.begin(); inpColNameIter!=inpColNames.end(); inpColNameIter++) {
        string inpColName=as<string>(*inpColNameIter);

        if (TYPEOF(df[inpColName])==REALSXP) {
              processNumericVector(*inpColNameIter,numCols,hashSeed,df[inpColName],resList);
        }
        else if (TYPEOF(df[inpColName])==INTSXP) {
          IntegerVector x=df[inpColName];
          if (isFactor(x)) {
            processFactorVector(*inpColNameIter,numCols,hashSeed,x,resList);
          } else {
            processIntegerVector(*inpColNameIter,numCols,hashSeed,x,resList);
          }
      }
    }

  DataFrame result(resList);
  result.attr("names") = colNames;
  return result;
}

void updateCountsForIntegerOrNumericVector
  (const char *colName, int numCols, int hashSeed, vector<int> &bucketCount) {
  int h=calcHash(colName, numCols, hashSeed);
  int targetCol=abs(h);
  bucketCount[targetCol-1]++;
}

void updateCountsForFactorVector
  (const char *colName, int numCols, int hashSeed, const IntegerVector &x, vector<int> &bucketCount) {
  CharacterVector levels=as<CharacterVector>(x.attr("levels"));

  for(int i=0; i<levels.length(); i++) {
    string totalName=string(colName)+":"+string(levels[i]);
    int h=calcHash(totalName.c_str(), numCols, hashSeed);
    int targetCol=abs(h);
    bucketCount[targetCol-1]++;
  }
}

void updateExplanationForIntegerOrNumericVector
  (const char *colName, int numCols, int hashSeed, List &resList, vector<int> &bucketIndex) {
  int h=calcHash(colName, numCols, hashSeed);
  int targetCol=abs(h);
  int sign=h<0 ? -1.0 : 1.0;
  int ind=targetCol-1;
  as<IntegerVector>(resList[ind])[bucketIndex[ind]]=sign;
  CharacterVector n=as<IntegerVector>(resList[ind]).attr("names");
  as<CharacterVector>(n)[bucketIndex[ind]]=colName;
  bucketIndex[ind]++;
}

void updateExplanationForFactorVector
  (const char *colName, int numCols, int hashSeed, IntegerVector &x, List &resList, vector<int> &bucketIndex) {
  CharacterVector levels=as<CharacterVector>(x.attr("levels"));

  for(int i=0; i<levels.length(); i++) {
    string totalName=string(colName)+":"+string(levels[i]);
    int h=calcHash(totalName.c_str(), numCols, hashSeed);
    int targetCol=abs(h);
    int sign=h<0 ? -1.0 : 1.0;
    int ind=targetCol-1;
    as<IntegerVector>(resList[ind])[bucketIndex[ind]]=sign;
    CharacterVector n=as<IntegerVector>(resList[ind]).attr("names");
    as<CharacterVector>(n)[bucketIndex[ind]]=totalName.c_str();
    bucketIndex[ind]++;
  }
}

//' Provides explanation, how given data frame will be hashed.
//'
//' Creates a data structure describing the hashing process.
//'
//' @param df Data frame to be hashed.
//' @param numCols Integer, number of columns for the output data frame.
//' @param hashSeed Integer, seed for the employed hash function (MurmurHash3).
//'
//' @return List of integer vectors. Each entry in list corresponds to one column
//' in the output of hashing. Entries in integer vectors give signs and names of integer
//' vector entries provide names of columns in original data frame.
//'
//' @export
// [[Rcpp::export]]

List explainHashDataFrame(DataFrame df, int numCols, int hashSeed) {

  // Prepare names for resulting list
  CharacterVector colNames(numCols);

  for(int i=0; i<numCols; i++) {
    stringstream ss;
    ss << "H" <<  (i+1);
    colNames[i]=ss.str();
  }

  List resList(numCols);
  resList.attr("names") = colNames;

  // Calculate number of entries in each bucket
  vector<int> bucketCount(numCols);
  fill(bucketCount.begin(), bucketCount.end(), 0);

  CharacterVector inpColNames=df.attr("names");
  for(CharacterVector::iterator inpColNameIter=inpColNames.begin(); inpColNameIter!=inpColNames.end(); inpColNameIter++) {
    string inpColName=as<string>(*inpColNameIter);

    if (TYPEOF(df[inpColName])==REALSXP) {
      updateCountsForIntegerOrNumericVector(*inpColNameIter,numCols,hashSeed,bucketCount);
    }
    else if (TYPEOF(df[inpColName])==INTSXP) {
      IntegerVector x=df[inpColName];
      if (isFactor(x)) {
        updateCountsForFactorVector(*inpColNameIter,numCols,hashSeed,x,bucketCount);
      } else {
        updateCountsForIntegerOrNumericVector(*inpColNameIter,numCols,hashSeed,bucketCount);
      }
    }
  }

  // Prepare buckets

  for(int i=0; i<numCols; i++) {
    IntegerVector v(bucketCount[i]);
    v.attr("names")=CharacterVector(bucketCount[i]);
    resList[i]=v;
  }

  // Fill buckets
  vector<int> bucketIndex(numCols);
  fill(bucketIndex.begin(), bucketIndex.end(), 0);

  for(CharacterVector::iterator inpColNameIter=inpColNames.begin(); inpColNameIter!=inpColNames.end(); inpColNameIter++) {
    string inpColName=as<string>(*inpColNameIter);

    if (TYPEOF(df[inpColName])==REALSXP) {
      updateExplanationForIntegerOrNumericVector(*inpColNameIter,numCols,hashSeed,resList,bucketIndex);
    }
    else if (TYPEOF(df[inpColName])==INTSXP) {
      IntegerVector x=df[inpColName];
      if (isFactor(x)) {
        updateExplanationForFactorVector(*inpColNameIter,numCols,hashSeed,x,resList,bucketIndex);
      } else {
        updateExplanationForIntegerOrNumericVector(*inpColNameIter,numCols,hashSeed,resList,bucketIndex);
      }
    }
  }
  for(CharacterVector::iterator inpColNameIter=inpColNames.begin(); inpColNameIter!=inpColNames.end(); inpColNameIter++) {
    string inpColName=as<string>(*inpColNameIter);

    if (TYPEOF(df[inpColName])==REALSXP) {
      updateCountsForIntegerOrNumericVector(*inpColNameIter,numCols,hashSeed,bucketCount);
    }
    else if (TYPEOF(df[inpColName])==INTSXP) {
      IntegerVector x=df[inpColName];
      if (isFactor(x)) {
        updateCountsForFactorVector(*inpColNameIter,numCols,hashSeed,x,bucketCount);
      } else {
        updateCountsForIntegerOrNumericVector(*inpColNameIter,numCols,hashSeed,bucketCount);
      }
    }
  }
  return resList;
}
