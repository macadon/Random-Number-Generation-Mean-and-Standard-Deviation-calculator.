#pragma once

#include "stdafx.h"
#include <vector>
const double PI = 3.14159265358979;

/**
 *  Computes the cumulative
 *  distribution function of the
 *  normal distribution
 */
double normcdf( double x );

/**
 *  Computes the inverse of normcdf
 */
double norminv( double x );

void solvequadratic(double a, double b, double c, double &z1, double &z2);
/**
 *  Test function
 */
double mean(vector<double> v);
double sdev(vector <double> v);
void testMatlib();
void swap(double a, double b);
void sort(vector<double> &v);
void randuniform(int n, vector <double> &v);
double unif();
double BM(int N, vector <double> &t);
