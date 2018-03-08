/**
* Adams-Moulton 4th order
* Runge-Kutta 4th order for first 3 values
* Equation -- Y + sqrt(pow(X,2) + pow(X,2)) - XY' = 0
* Cauchy -- Y(Xo) = -0.5
* Limits -- [0.0 ; 1.0]
* Solution -- Y(X) = (pow(x,2) - 1) / 2
*/

/**
* Possible testing equation system:
* y1 = sinx;
* y2 = cosx;
* ----------
* y1' = cosx = y2
* y2' = -sinx = -y1;
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>

#include "Header.h"

void Output(double **yApproximate, double **yAnalytical, double *xk, int neq, int nOutIter);



int main(void)
{
	// limits of an interval
	double a = pow(10, -10);
	double b = 1.0;

	constexpr int neq = 1; // number of equations

						   // integrational steps for different accuracy of calculations
	constexpr double h1 = 0.2;			// h1 = 1/5
	constexpr double h2 = h1 / 5.0;		// h2 = 1/25
	constexpr double h3 = h1 / 25.0;	// h3 = 1/125
	constexpr double h = h3;

	// output preparation
	double outputStep = 0.1; // bigger step for output purpose

							 // correction (not a must-have)
	int nOutIter = (int)((b - a) / outputStep); // number of output iterations
	outputStep = (b - a) / nOutIter;

	// arrays to store OUTPUT values
	// TODO remake the whole idea to use these arrays w/o using array y (and maybe yp)
	double **yApproximate = new double*[nOutIter + 1];
	double **yAnalytical = new double*[nOutIter + 1];
	for (int i = 0; i < nOutIter + 1; ++i)
	{
		yApproximate[i] = new double[neq];
		yAnalytical[i] = new double[neq];
	}
	// array to store arguments
	double *xk = new double[outputStep + 1];

	// array to store CALCULATED values
	double *y = new double[neq];

	// Cauchy condition
	y[0] = -0.5; // TODO get rid of when remaking the whole idea
	yApproximate[0][0] = -0.5;

	// function derivative
	double *yDerivative = new double[neq];



	for (int it = 1; it <= nOutIter; it++)
	{
		// [Xki; Xki+1] -- intervals we iterate over [from a to b]
		double Xk = a + outputStep * (it - 1);
		double Xk1 = a + outputStep * (it);
		xk[it - 1] = Xk;

		for (int ieq = 0; ieq < neq; ieq++)
		{
			yApproximate[it][ieq] = AdamsMoulton(neq, Xk, Xk1, h, y, yDerivative);
			y[ieq] = yApproximate[it][ieq];

			yAnalytical[it - 1][ieq] = solution(neq, Xk);
		}
	}




	Output(yApproximate, yAnalytical, xk, neq, nOutIter);

	std::cin.get();
	for (int i = 0; i < nOutIter + 1; i++)
	{
		delete[] yApproximate[i];
		delete[] yAnalytical[i];
	}
	delete[] yApproximate;
	delete[] yAnalytical;
	delete[] xk;
	delete[] y;
	delete[] yDerivative;

	return 0;
}


void Output(double **yApproximate, double **yAnalytical, double *xk, int neq, int nOutIter)
{
	std::cout << std::setprecision(12);
	// header
	std::cout << " Xk              ";
	for (int ieq = 0; ieq < neq; ieq++)
	{
		std::cout << " | Yappr(X)        "
			<< " | Y(X)            " << " | E(x)            "
			<< " | 100*E(x)/Yk" << std::left << std::endl;
	}
	for (int i = 0; i < 19 + 76 * neq; i++) std::cout << "-"; // header line
	std::cout << std::endl;

	// body
	for (int it = 0; it < nOutIter; it++)
	{
		std::cout << " " << std::setw(16) << xk[it];
		for (int ieq = 0; ieq < neq; ieq++)
		{
			std::cout << " | " << std::setw(16) << yApproximate[it][ieq]
				<< " | " << std::setw(16) << yAnalytical[it][ieq] << " | " << std::setw(16) << abs(yApproximate[it][ieq] - yAnalytical[it][ieq])
				<< " | " << std::setw(16) << abs((yApproximate[it][ieq] - yAnalytical[it][ieq]) * 100.0 / yAnalytical[it][ieq]) << std::fixed << std::endl;
		}
	}
}

double solution(int neq, double x)
{
	return (pow(x, 2) - 1) / 2;
}

void function(int neq, double x, double *y, double *yDerivative)
{
	yDerivative[0] = (y[0] + sqrt(pow(x, 2) + pow(y[0], 2))) / x;
}