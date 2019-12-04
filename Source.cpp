/*
Using RK4 for solving the differential equation for concentration on CO2 in the atmosphere

  1. Reproduce the graphics for emission of CO_2 --> f(t)
  2. Which metod is the most apropiate ?
  3. When does the atmospheric concentration of CO_2 reach its maximum ?

*/

//-----LIBRARIES-----//

#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

//------GLOBAL CONSTANTS-----//

const double d = 8.64;
const double u1 = 4.95;
const double u2 = 0.049;
const double vs = 0.12;
const double vd = 1.23;
const double w = 0.001;
const double k1 = 0.000219;
const double k2 = 0.0000612;
const double k3 = 0.997148;
const double k4 = 0.0679;

//------DATA VECTORS FOR f(t)-----//

std::vector<double> c{ 2.0, 10.5, 2.7 };
std::vector<double> t{ 1988, 2100, 2265 };
std::vector<double> s{ 21, 96, 57 };

//------DECLARATION OF ADITIONAL FUNCTIONS-----//

double function_f(double T);  //Emission of CO_2 between years 1000 and 5000
double function_ps(std::vector<double>& y); // Shallow-ocean equilibrium equatios
double f(double T, std::vector<double>& y, int id); //internal function for RK4 (derivates)
void RK4(double TA, double TB, double H, std::vector<double>& y);

int main(void)
{
	//------INITIAL CONDITIONS-----//

	double p = 1.00;
	double Ss = 2.01;
	double Sd = 2.23;
	double As = 2.2;
	double Ad = 2.26;
	double T = 1000;

	/*
	For nomenclature we are taking:
	- Ss : Sigma sub s -->  Concentration of carbonates dissolved C in shallow wathers
	- Sd : Sigma sub d -->  Concentration of carbonates dissolved C in oceans
	- As : Alhpa sub s -->  Alkalinities in shallow wathers
	- Ad : Alpha sub d -->  Alkalinities in oceans
	- p  :			   -->  Partial pressure of CO_2 in atmosphere
	- T  :			   -->  For calculating year
	*/

	//------OPERATING DATA-----//

	const double TA = 1000;
	const double TB = 5000;
	const double H = 1;
	std::vector<double> y = { p, Ss, Sd, As, Ad };  //make y vector constant 

	/*
	For nomenclature we are taking:
	- N  -->
	- TA -->  Initial calculation point
	- TB -->  Final calculation point
	- H  -->  Presition (derivate function step for RK4
	- y[0] --> p
	- y[1] --> Ss
	- y[2] --> Sd
	- y[3] --> As
	- y[4] --> Ad
	*/
	std::cout.precision(5); std::cout.setf(std::ios::scientific);

	RK4(TA, TB, H, y);
}

double function_f(double T) //T is the dynamic parameter (changes in time)
{
	double suma = 0.0;
	for (int id = 0; id < 3; id++)
	{
		suma += c[id] * std::exp((-(T - t[id]) * (T - t[id])) / (s[id] * s[id]));
	}
	return suma;
}
double function_ps(std::vector<double>& y)
{
	double hs = (y[1] - std::sqrt(y[1] * y[1] - k3 * y[3] * (2 * y[1] - y[3]))) / k3;
	double cs = (y[3] - hs) / 2;
	return (k4 * hs * hs) / cs;
}
double f(double T, std::vector<double>& y, int id)
{
	if (id == 0)
	{
		return ((function_ps(y) - y[0]) / d) + (function_f(T) / u1);
	}
	else if (id == 1)
	{
		return (w * (y[2] - y[1]) - k1 - u2 * ((function_ps(y) - y[0]) / d)) / vs;
	}
	else if (id == 2)
	{
		return (k1 - w * (y[1] - y[2])) / vd;
	}
	else if (id == 3)
	{
		return (k1 - w * (y[4] - y[3])) / vd;
	}
	else if (id == 4)
	{
		return (w * (y[4] - y[3]) - k2) / vs;
	}
	else
	{
		std::cerr << "ERROR!!!!! it does not exist -> " << id << std::endl;
		exit(1);
	}
}
void RK4(double TA, double TB, double H, std::vector<double>& y)
{
	std::vector<double>  C1, C2, C3, C4, aux;
	C1.resize(y.size());
	C2.resize(y.size());
	C3.resize(y.size());
	C4.resize(y.size());
	aux.resize(y.size());

	const int N = (TB - TA) / H;
	for (int id = 0; id < N; id++) {
		double T = TA + (H * id);
		for (int ii = 0; ii < y.size(); ii++)
		{
			C1[ii] = H * f(T, y, ii); //una función interna
		}

		//C2 aux
		for (int ii = 0; ii < y.size(); ii++)
		{
			aux[ii] = y[ii] + C1[ii] / 2;
		}

		//C2
		for (int ii = 0; ii < y.size(); ii++)
		{
			C2[ii] = H * f(T + H / 2, aux, ii);
		}

		//C3 aux
		for (int ii = 0; ii < y.size(); ii++)
		{
			aux[ii] = y[ii] + C2[ii] / 2;
		}

		//C3
		for (int ii = 0; ii < y.size(); ii++)
		{
			C3[ii] = H * f(T + H / 2, aux, ii);
		}

		//C4 aux
		for (int ii = 0; ii < y.size(); ii++)
		{
			aux[ii] = y[ii] + C3[ii];
		}

		//C4
		for (int ii = 0; ii < y.size(); ii++)
		{
			C4[ii] = H * f(T + H, aux, ii);
		}

		//write new yy
		for (int ii = 0; ii < y.size(); ii++)
		{
			y[ii] = y[ii] + (C1[ii] + 2 * C2[ii] + 2 * C3[ii] + C4[ii]) / 6.0;
		}
		std::cout << T << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t" << y[4] << std::endl;
		//replaced fout for std::cout
	}
}
