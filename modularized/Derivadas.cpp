#include<vector>
#include<fstream>
#include<cmath>
#include"constants.h"

double fuel(double t);
double ps(const std::vector<double>& y);

double f(double t, const std::vector<double>& y, int i)
{
	if (i == 0) {
		return (ps(y) - y[0]) / d + fuel(t) / u1;
	}
	else if (i == 1) {
		return (w * (y[2] - y[1]) - k1 - u2 * (ps(y) - y[0]) / d) / vs;
	}
	else if (i == 2) {
		return (k1 - w * (y[2] - y[1])) / vd;
	}
	else if (i == 3) {
		return (w * (y[4] - y[3]) - k2) / vs;
	}
	else if (i == 4) {
		return (k2 - w * (y[4] - y[3])) / vd;
	}
}