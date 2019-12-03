/*
Using RK4 for solving the differential equation for concentration on CO2 in the atmosphere

  1. Reproduce the graphics for emission of CO_2 --> f(t)
  2. Which metod is the most apropiate ?
  3. When does the atmospheric concentration of CO_2 reach its maximum ?

*/

//------DECLARATION OF INCLUDED LIBRARIES-----//

#include<iostream>
#include<vector>
#include<fstream>
//#include"constants.h"

//------GLOBAL CONSTANTS-----//

const double d = 8.64;
const double u1 = 4.95;
const double u2 = 0.049;
const double vs = 0.12;
const double vd = 1.23;
const double w = 0.001;
const double k1 = 0.00219;
const double k2 = 0.0000612;
const double k3 = 0.997148;
const double k4 = 0.0679;

//------DATA STRUCTURE FOR f(t)-----//

std::vector<double> c{ 2.0, 10.5, 2.7 };
std::vector<double> t{ 1988, 2100, 2265 };
std::vector<double> s{ 21, 96, 57 };


//------DECLARATION OF ADITIONAL FUNCTIONS-----//

double function_f(double T);
double funciton_ps(std::vector<double>& y);
void f_int(double T, std::vector<double>& y, double id);
double f(double T, std::vector<double>& y, double id);
void RK4(double TA, double TB, double H, std::vector<double>& y);

//------DECLARATION OF ADITIONAL  %% TESTING %% FUNCTIONS-----//

void print_vector(std::vector<double>& y); //		--> for testing purpose, its not necesary
void print_function_f(double T);

void print_vector(std::vector<double>& y)
{
	for (int id = 0; id < y.size(); ++id)
	{
		std::cout << y[id] << "\t";
	}
	std::cout << "\n";
}

double function_f(double T)
{
	double suma = 0.0;
	for (int id = 0; id < 3; id++)
	{
		suma += c[id] * std::exp((-(T - t[id]) * (T - t[id])) / (s[id] * s[id]));
	}
	return suma;

}

void print_function_f(double T)
{
	double suma = 0.0;
	for (int id = 0; id < 3; id++)
	{
		suma += c[id] * std::exp((-(T - t[id]) * (T - t[id])) / (s[id] * s[id]));
	}
	std::cout << suma << std::endl; //Print the result of the sumatory function
}

double function_ps(std::vector<double>& y)
{
	const double k3 = 0.997148;
	const double k4 = 0.0679;
	double hs = 0.0;
	double cs = 0.0;
	double ps = 0.0;
	hs = (y[1] - std::sqrt(y[1] * y[1] - k3 * y[3] * (2 * y[1] - y[3]))) / k3;
	cs = (y[3] - hs) / 2;
	ps = (k4 * hs * hs) / cs;
	return ps;
}

double f(double T, std::vector<double>& y, int id) //rk4 internal function
{
	if (id == 0)
	{
		double id_0 = 0.0;
		id_0 = ((function_ps(y) - y[0]) / d) + (function_f(T) / u1);
		return id_0;
	}
	else if (id == 1)
	{
		double id_1 = 0.0;
		id_1 = (w * (y[2] - y[1]) - k1 - u2 * ((function_ps(y) - y[0]) / d)) / vd;
		return id_1;
	}
	else if (id == 2)
	{
		double id_2 = 0.0;
		id_2 = (k1 - w * (y[1] - y[2])) / vd;
		return id_2;
	}
	else if (id == 3)
	{
		double id_3 = 0.0;
		id_3 = (k1 - w * (y[4] - y[3])) / vd;
		return id_3;
	}
	else if (id == 4)
	{
		double id_4 = 0.0;
		id_4 = (w * (y[4] - y[3] - k2)) / vs;
		return id_4;;
	}
	else
	{
		std::cerr << "ERROR!!!!! it does not exist -> " << id << std::endl;
		exit(1);
	}
}

void f_int(double T, std::vector<double>& y, int id) //rk4 internal function
{
	if (id == 0)
	{
		double id_0 = 0.0;
		id_0 = ((function_ps(y) - y[0]) / d) + (function_f(T) / u1);
		//return id_0;
		std::cout << id_0 << std::endl;
	}
	else if (id == 1)
	{
		double id_1 = 0.0;
		id_1 = (w * (y[2] - y[1]) - k1 - u2 * ((function_ps(y) - y[0]) / d)) / vd;
		//return id_1;
		std::cout << id_1 << std::endl;
	}
	else if (id == 2)
	{
		double id_2 = 0.0;
		id_2 = (k1 - w * (y[1] - y[2])) / vd;
		//return id_2;
		std::cout << id_2 << std::endl;
	}
	else if (id == 3)
	{
		double id_3 = 0.0;
		id_3 = (k1 - w * (y[4] - y[3])) / vd;
		//return id_3;
		std::cout << id_3 << std::endl;
	}
	else if (id == 4)
	{
		double id_4 = 0.0;
		id_4 = (w * (y[4] - y[3] - k2)) / vs;
		//return id_4;
		std::cout << id_4 << std::endl;
	}
	else
	{
		std::cerr << "ERROR!!!!! it does not exist -> " << id << std::endl;
		exit(1);
	}
}

void rk4(double TA, double TB, double H, std::vector<double>& y) //Thi RK4 is for all points on the curve
{
	std::vector<double>  C1, C2, C3, C4, aux;
	C1.resize(y.size());
	C2.resize(y.size());
	C3.resize(y.size());
	C4.resize(y.size());
	aux.resize(y.size());

	const int N = (TB - TA) / H;
	for (int id = 0; id < N; ++id) {
		double t = TA + H * id;
		for (int ii = 0; ii < y.size(); ii++) {
			C1[ii] = H * f(t, y, ii); //una función interna
		}

		//C2 aux
		for (int ii = 0; ii < y.size(); ++ii) {
			aux[ii] = y[ii] + C1[ii] / 2;
		}

		//C2
		for (int ii = 0; ii < y.size(); ++ii) {
			C2[ii] = H * f(t + H / 2, aux, ii);
		}

		//C3 aux
		for (int ii = 0; ii < y.size(); ++ii) {
			aux[ii] = y[ii] + C2[ii] / 2;
		}

		//C3
		for (int ii = 0; ii < y.size(); ++ii) {
			C3[ii] = H * f(t + H / 2, aux, ii);
		}

		//C4 aux
		for (int ii = 0; ii < y.size(); ++ii) {
			aux[ii] = y[ii] + C3[ii];
		}

		//C4
		for (int ii = 0; ii < y.size(); ++ii) {
			C4[ii] = H * f(t + H, aux, ii);
		}

		//write new y
		for (int ii = 0; ii < y.size(); ++ii) {
			y[ii] = y[ii] + (C1[ii] + 2 * C2[ii] + 2 * C3[ii] + C4[ii]) / 6.0;
		}
		fout << t << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << /*"\t" << y[3] << "\t" << y[4] <<*/ std::endl;
		//replaced fout for std::cout
	}
}

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

	const double N = 2;
	const double TA = 1000;
	const double TB = 5000;
	const double H = 100;
	std::vector<double> y = { p, Ss, Sd, As, Ad };

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

	//------PRESITION FOR CALULATIONS CODE-----//

	std::cout.precision(5); std::cout.setf(std::ios::scientific);
	
	//------MAIN CODE-----//

	rk4(TA, TB, H, y);

	return 0;
	std::cin.get();
}