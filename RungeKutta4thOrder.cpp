#include "Header.h"

double RungeKutta(int neq, double Xk, double Xk1, double h, double *y, double *yDerivative)
{
	// coefficients for Runge Kutta 4th order
	double *k1 = new double[neq];
	double *k2 = new double[neq];
	double *k3 = new double[neq];
	double *k0 = new double[neq];

	double x;

	double finalValue;

	int numOfIterations = (int)((Xk1 - Xk) / h); // Xk1 - Xk >> h --> nt is quite big
	if (numOfIterations == 0) numOfIterations = 1; // check to avoid division by zero
	h = (Xk1 - Xk) / numOfIterations; // make h more precise

	for (int it = 1; it <= numOfIterations; it++)
	{
		double a = Xk + (it - 1) * h; // left side (reserved, unlike x)

		for (int i = 0; i < neq; i++)
		{
			k0[i] = y[i]; // make copies of an array of initial values
		}

		function(neq, a, y, yDerivative); // calculate function Y' = F(Xk, Yk)

		for (int i = 0; i < neq; i++)
		{
			k1[i] = h * yDerivative[i];
			y[i] = k0[i] + k1[i] / 2.0; // set next Yi for next Ki
		}

		x = a + h / 2.0; // set next Xi for next Ki
		function(neq, x, y, yDerivative);

		for (int i = 0; i < neq; i++)
		{
			k2[i] = h * yDerivative[i];
			y[i] = k0[i] + k2[i] / 2.0;
		}

		x = a + h / 2.0;
		function(neq, x, y, yDerivative);

		for (int i = 0; i < neq; i++)
		{
			k3[i] = h * yDerivative[i];
			y[i] = k0[i] + k3[i];
		}

		x = a + h;
		function(neq, x, y, yDerivative);

		for (int i = 0; i < neq; i++)
		{
			y[i] = k0[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + yDerivative[i] * h); // k4 = yp * h
		}
	}
	return y[0];
}

