#include "Header.h"

#include <iostream>

#define RK_ITERATIONS 3

double AdamsMoulton(int neq, double Xk, double Xk1, double h, double *y, double *yDerivative)
{

	double x;
	double finalValue;

	// coefficients
	double *k0py = new double[neq];
	double *k1 = new double[neq];
	double *k2 = new double[neq];
	double *k3 = new double[neq];
	double *k4 = new double[neq];

	int numOfIterations = (int)((Xk1 - Xk) / h); // Xk1 - Xk >> h --> nt is quite big
	if (numOfIterations < 1) numOfIterations = 1; // check to avoid division by zero
	h = (Xk1 - Xk) / numOfIterations; // make h more precise

									  // temporary array of values inside algorithm function
	double **Yarray = new double*[numOfIterations + 1];
	for (int i = 0; i < numOfIterations + 1; ++i)
	{
		Yarray[i] = new double[neq];
	}


	/* RUNGE KUTTA 4th ORDER */
	int it = 1;
	for (; it <= RK_ITERATIONS && it <= numOfIterations; it++)
	{
		double a = Xk + (it - 1) * h; // left side (reserved, unlike x)

		if (it == 1)
		{
			for (int ieq = 0; ieq < neq; ieq++)
			{
				k0py[ieq] = y[ieq]; // make copies of an array of initial values

				Yarray[0][ieq] = -0.5;
				Yarray[it][ieq] = y[ieq];

				/// NOTE: Tanya
				//Yarray[0][ieq] = 0.0;
			}
			function(neq, a, y, yDerivative); // calculate function Y' = F(Xk, Yk)
		}
		else
		{
			for (int ieq = 0; ieq < neq; ieq++)
			{
				k0py[ieq] = Yarray[it - 1][ieq]; // make copies of an array of previous values
				Yarray[it][ieq] = k0py[ieq];
			}
			function(neq, a, Yarray[it - 1], yDerivative); // calculate function Y' = F(Xk, Yk)
		}

		for (int ieq = 0; ieq < neq; ieq++)
		{
			k1[ieq] = h * yDerivative[ieq];
			Yarray[it][ieq] = k0py[ieq] + k1[ieq] / 2.0; // set next Yi for next Ki
		}

		x = a + h / 2.0; // set next Xi for next Ki
		function(neq, x, Yarray[it], yDerivative);

		for (int ieq = 0; ieq < neq; ieq++)
		{
			k2[ieq] = h * yDerivative[ieq];
			Yarray[it][ieq] = k0py[ieq] + k2[ieq] / 2.0;
		}

		x = a + h / 2.0;
		function(neq, x, Yarray[it], yDerivative);

		for (int ieq = 0; ieq < neq; ieq++)
		{
			k3[ieq] = h * yDerivative[ieq];
			Yarray[it][ieq] = k0py[ieq] + k3[ieq];
		}

		x = a + h;
		function(neq, x, Yarray[it], yDerivative);

		for (int ieq = 0; ieq < neq; ieq++)
		{
			k4[ieq] = h * yDerivative[ieq];
		}

		for (int ieq = 0; ieq < neq; ieq++)
		{
			Yarray[it][ieq] = k0py[ieq] + 1.0 / 6.0 * (k1[ieq] + 2.0 * k2[ieq] + 2.0 * k3[ieq] + k4[ieq]);
			finalValue = Yarray[it][ieq];
		}
	}

	double *Ypredicted = new double[neq];

	/* ADAMS MOULTON 4th ORDER */
	for (; it <= numOfIterations; it++)
	{
		function(neq, Xk, Yarray[it - 1], yDerivative); //fk
		for (int ieq = 0; ieq < neq; ieq++)
			k2[ieq] = 19.0 * yDerivative[ieq];


		// ~yk+1 = yk + h*f(xk, yk)
		for (int ieq = 0; ieq < neq; ieq++)
			Ypredicted[ieq] = Yarray[it - 1][ieq] + h * yDerivative[ieq]; // make prediction ~Y from previous derivative Fk
		function(neq, Xk + h * 1.0, Ypredicted, yDerivative); // make predictionary derivative from ~Y -- Fk+1
		for (int ieq = 0; ieq < neq; ieq++)
			k1[ieq] = 9.0 * yDerivative[ieq];


		function(neq, Xk - h * 1.0, Yarray[it - 2], yDerivative); //fk-1
		for (int ieq = 0; ieq < neq; ieq++)
			k3[ieq] = -5.0 * yDerivative[ieq];

		function(neq, Xk - h * 2.0, Yarray[it - 3], yDerivative); //fk-2
		for (int ieq = 0; ieq < neq; ieq++)
			k4[ieq] = 1.0 * yDerivative[ieq];


		for (int ieq = 0; ieq < neq; ieq++)
		{
			Yarray[it][ieq] = Yarray[it - 1][ieq] + h / 24.0 * (k1[ieq] + k2[ieq] + k3[ieq] + k4[ieq]);
			finalValue = Yarray[it][ieq];
		}
	}


	//for (int i = 0; i < numOfIterations + 1; i++)
	//	delete[] Yarray[i];
	//delete[] Yarray;
	//delete[] Ypredicted;
	//delete[] k0py;
	//delete[] k1;
	//delete[] k2;
	//delete[] k3;
	//delete[] k4;


	// std::cout << "Yfin = " << finalValue << std::endl;
	return finalValue;
}
