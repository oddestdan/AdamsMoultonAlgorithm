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
	double **yApproximate = new double*[neq];
	double **yAnalytical = new double*[neq];
	for (int i = 0; i < neq; ++i)
	{
		yApproximate[i] = new double[nOutIter + 1];
		yAnalytical[i] = new double[nOutIter + 1];
	}
	// array to store arguments
	double *xk = new double[32];

	// array to store CALCULATED values
	double *y = new double[neq];

	// Cauchy condition
	y[0] = -0.5; // TODO get rid of when remaking the whole idea
	yApproximate[0][0] = -0.5;

	// function derivative
	double *yDerivative = new double[neq];

	for (int it = 0; it < nOutIter; it++)
	{
		// [Xki; Xki+1] -- intervals we iterate over [from a to b]
		double Xk = a + outputStep * (it);
		double Xk1 = a + outputStep * (it + 1);
		xk[it + 1] = Xk;

		for (int i = 0; i < neq; i++)
		{
			yApproximate[i][it + 1] = RungeKutta(neq, Xk, Xk1, h, y, yDerivative);
			y[i] = yApproximate[0][it + 1];

			yAnalytical[i][it] = solution(neq, Xk);
		}

	}

	
	// BASHFORTH TEST



	//


	Output(yApproximate, yAnalytical, xk, neq, nOutIter);

	std::cin.get();
	for (int i = 0; i < neq; i++)
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
	for (int i = 0; i < neq; i++)
	{
		std::cout << " | Y(Xk)           "
			<< " | Y(X)            " << " | E(x)            "
			<< " | 100*E(x)/Yk" << std::left << std::endl;
	}
	for (int i = 0; i < 19 + 76 * neq; i++) std::cout << "-"; // header line
	std::cout << std::endl;

	// body
	for (int j = 0; j < nOutIter; j++)
	{
		std::cout << " " << std::setw(16) << xk[j];
		for (int i = 0; i < neq; i++)
		{
			std::cout << " | " << std::setw(16) << yApproximate[i][j]
				<< " | " << std::setw(16) << yAnalytical[i][j] << " | " << std::setw(16) << yAnalytical[i][j] - yApproximate[i][j]
				<< " | " << std::setw(16) << (yApproximate[i][j] - yAnalytical[i][j]) * 100.0 / yAnalytical[i][j] << std::fixed << std::endl;
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