#include<vector>
#include<fstream>
#include<cmath>
#include"constants.h"

double fuel(double t)
{
	double result = 0;

	for (int i = 0; i < 3; i++) {
		result += c[i] * std::exp(-std::pow((t - t1[i]), 2) / std::pow(s[i], 2));
	}

	return result;
}