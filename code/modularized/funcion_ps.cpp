#include<vector>
#include<fstream>
#include<cmath>
#include"constants.h"

double ps(const std::vector<double>& y)
{
	double hs = (y[1] - std::sqrt(std::pow(y[1], 2) - k3 * y[3] * (2 * y[1] - y[3]))) / k3;
	double cs = (y[3] - hs) / 2;

	return k4 * (std::pow(hs, 2) / cs);
}