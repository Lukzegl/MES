#pragma once
#include <cmath>
#include <vector>

struct GausPoint
{
	static double x2[2];
	static double w2[2];
	static double x3[3];
	static double w3[3];
	static double x4[4];
	static double w4[4];
};

double f(double x);
double g(double x, double y);

double kwadratura(double (*f)(double), double a, double b, int n);
double kwadratura2D(double (*g)(double, double), double ax, double bx, double ay, double by, int n);

std::vector<double> solveLU(std::vector<std::vector<double>> A, std::vector<double> B);