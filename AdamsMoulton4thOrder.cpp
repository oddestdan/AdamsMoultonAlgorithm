#include "Header.h"

#include <iostream>

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
	for (;/* it <= 3 && */it <= numOfIterations; it++)
	{
		double a = Xk + (it - 1) * h; // left side (reserved, unlike x)

		if (it == 1)
		{
			for (int ieq = 0; ieq < neq; ieq++)
			{
				k0py[ieq] = y[ieq]; // make copies of an array of initial values
				Yarray[it][ieq] = y[ieq];
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

	/* ADAMS MOULTON 4th ORDER */
	//for (; it <= numOfIterations; it++)
	//{
	//	function(neq, x, Yarray[it], yDerivative);
	//	
	//	for (int ieq = 0; ieq < neq; ieq++)
	//		k1[ieq] = 55.0 * yDerivative[ieq];

	//	function(neq, x - h, Yarray[it], yDerivative);
	//	
	//	for (int ieq = 0; ieq < neq; ieq++)
	//		k2[ieq] = -59.0 * yDerivative[ieq];

	//	function(neq, x - h * 2.0, Yarray[it - 1], yDerivative);
	//	
	//	for (int ieq = 0; ieq < neq; ieq++)
	//		k3[ieq] = 37.0 * yDerivative[ieq];

	//	function(neq, x - h * 3.0, Yarray[it - 2], yDerivative);
	//	
	//	for (int ieq = 0; ieq < neq; ieq++)
	//		k4[ieq] = -9.0 * yDerivative[ieq];

	//	for (int ieq = 0; ieq < neq; ieq++)
	//		Yarray[it][neq] = Yarray[it - 1][neq] * h / 24.0 * (k1[0] + k2[0] + k3[0] + k4[0]);
	//}



	for (int i = 0; i < numOfIterations + 1; i++)
	{
		delete[] Yarray[i];
	}
	delete[] Yarray;
	delete[] k0py;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;

	std::cout << "Yfin = " << finalValue << std::endl;
	return finalValue;
}
