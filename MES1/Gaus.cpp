#include <iostream>
#include <cmath>
#include "Gaus.h"
#include <vector>

using namespace std;

double GausPoint::x2[2] = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
double GausPoint::w2[2] = { 1.0, 1.0 };

double GausPoint::x3[3] = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
double GausPoint::w3[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

double GausPoint::x4[4] = {
	-sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0),
	-sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0),
	sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0),
	sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0)
};
double GausPoint::w4[4] = {
	(18.0 - sqrt(30.0)) / 36.0,
	(18.0 + sqrt(30.0)) / 36.0,
	(18.0 + sqrt(30.0)) / 36.0,
	(18.0 - sqrt(30.0)) / 36.0
};

GausPoint gausPoints;

double f(double x)
{
	return 5 * x * x + 3 * x + 6;
}

double g(double x, double y)
{
	return 5 * x * x * y * y + 5 * x * y + 6;
}

double skalowanie1d(double a, double b, double xi)
{
	return ((b - a) / 2.0) * xi + (a + b) / 2.0;
}

double kw1D(double (*f)(double), double a, double b, int n, double x[], double w[])
{
	double wynik = 0.0;
	for (int i = 0; i < n; i++)
	{
		wynik += w[i] * f(skalowanie1d(a, b, x[i]));
	}

	return ((b - a) / 2.0) * wynik;
}

double kwadratura(double (*f)(double), double a, double b, int n)
{
	if (n == 2)
	{
		return kw1D(f, a, b, n, gausPoints.x2, gausPoints.w2);
	}
	else if (n == 3)
	{
		return kw1D(f, a, b, n, gausPoints.x3, gausPoints.w3);
	}
	else if (n == 4)
	{
		return kw1D(f, a, b, n, gausPoints.x4, gausPoints.w4);
	}
	else
	{
		cout << "ERROR ";
		return 0;
	}
}

double skalowanie2d(double a, double b, double xi) {
	return ((b - a) / 2.0) * xi + (a + b) / 2.0;
}

double kw2D(double (*g)(double, double), double ax, double bx, double ay, double by, int n, double x[], double w[])
{
	double wynik = 0.0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) 
		{
			double x_scaled = skalowanie2d(ax, bx, x[i]);
			double y_scaled = skalowanie2d(ay, by, x[j]);
			wynik += w[i] * w[j] * g(x_scaled, y_scaled);
		}
	}

	return ((bx - ax) / 2.0) * ((by - ay) / 2.0) * wynik;
}

double kwadratura2D(double (*g)(double, double), double ax, double bx, double ay, double by, int n)
{
	if (n == 2)
	{
		return kw2D(g, ax, bx, ay, by, n, gausPoints.x2, gausPoints.w2);
	}
	else if (n == 3)
	{
		return kw2D(g, ax, bx, ay, by, n, gausPoints.x3, gausPoints.w3);
	}
	else if (n == 4)
	{
		return kw2D(g, ax, bx, ay, by, n, gausPoints.x4, gausPoints.w4);
	}
	else
	{
		cout << "ERROR ";
		return 0;
	}
}




void LUDecomposition(int n, vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) {
	for (int i = 0; i < n; i++) 
	{
		
		L[i][i] = 1;

		
		for (int j = i; j < n; j++) 
		{
			U[i][j] = A[i][j];

			for (int k = 0; k < i; k++) 
			{
				U[i][j] -= L[i][k] * U[k][j];
			}
		}

		
		for (int j = i + 1; j < n; j++) 
		{
			L[j][i] = A[j][i];

			for (int k = 0; k < i; k++) 
			{
				L[j][i] -= L[j][k] * U[k][i];
			}
			L[j][i] /= U[i][i];
		}
	}
}


void forwardSubstitution(int n, vector<vector<double>>& L, vector<double>& Y, vector<double>& B) 
{
	for (int i = 0; i < n; i++) 
	{
		Y[i] = B[i];

		for (int j = 0; j < i; j++) 
		{
			Y[i] -= L[i][j] * Y[j];
		}
	}
}


void backwardSubstitution(int n, vector<vector<double>>& U, vector<double>& X, vector<double>& Y) 
{
	for (int i = n - 1; i >= 0; i--) 
	{
		X[i] = Y[i];
		for (int j = i + 1; j < n; j++) 
		{
			X[i] -= U[i][j] * X[j];
		}
		X[i] /= U[i][i];
	}
}

// A * X = B LU
vector<double> solveLU(vector<vector<double>> A, vector<double> B) 
{
	int n = A.size();

	
	vector<vector<double>> L(n, vector<double>(n, 0));
	vector<vector<double>> U(n, vector<double>(n, 0));

	
	LUDecomposition(n, A, L, U);

	
	vector<double> Y(n);
	vector<double> X(n);

	
	forwardSubstitution(n, L, Y, B);

	
	backwardSubstitution(n, U, X, Y);

	return X;
}