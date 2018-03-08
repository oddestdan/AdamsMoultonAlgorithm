#pragma once
double solution(int neq, double x);
void function(int neq, double x, double *y, double *yDerivative);
void Output(double **yApproximate, double **yAnalytical, double *xk, int neq, int nOutIter);
double RungeKutta(int neq, double Xk, double Xk1, double h, double *y, double *yDerivative);