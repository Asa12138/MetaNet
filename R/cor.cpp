#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace std;

typedef vector<double> Vector;

double sum(Vector a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += a[i];
  }
  return s;
}

double mean(Vector a)
{
  return sum(a) / a.size();
}

double sqsum(Vector a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(Vector nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

Vector operator-(Vector a, double b)
{
  Vector retvect;
  for (int i = 0; i < a.size(); i++)
  {
    retvect.push_back(a[i] - b);
  }
  return retvect;
}

Vector operator*(Vector a, Vector b)
{
  Vector retvect;
  for (int i = 0; i < a.size() ; i++)
  {
    retvect.push_back(a[i] * b[i]);
  }
  return retvect;
}

// [[Rcpp::export]]
double pearsoncoeff(Vector X, Vector Y)
{
  return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

Vector rankify(Vector & X) {

  int N = X.size();

  // Rank Vector
  Vector Rank_X(N);

  for(int i = 0; i < N; i++)
  {
    int r = 1, s = 1;

    // Count no of smaller elements
    // in 0 to i-1
    for(int j = 0; j < i; j++) {
      if (X[j] < X[i] ) r++;
      if (X[j] == X[i] ) s++;
    }

    // Count no of smaller elements
    // in i+1 to N-1
    for (int j = i+1; j < N; j++) {
      if (X[j] < X[i] ) r++;
      if (X[j] == X[i] ) s++;
    }

    // Use Fractional Rank formula
    // fractional_rank = r + (n-1)/2
    Rank_X[i] = r + (s-1) * 0.5;
  }

  // Return Rank Vector
  return Rank_X;
}

// [[Rcpp::export]]
double spearmancoeff(Vector X, Vector Y)
{
  // Get ranks of vector X
  Vector rank_x = rankify(X);

  // Get ranks of vector y
  Vector rank_y = rankify(Y);
  return pearsoncoeff(rank_x, rank_y);
}




