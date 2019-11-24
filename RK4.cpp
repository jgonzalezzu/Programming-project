/*
Using RK4 for solving the differential equation for concentration on CO2 in the atmosphere
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>

const double W = 2.98768765;

double function_f(); // function of emition of CO2 in to the atmosphere
double d_p();
double d_Ss();
double d_Sd();
double d_As();
double d_Ad();


double f(double t, const std::vector<double> &y, int id);
void rk4 (double ta, double tb, double h, std::vector<double> & y);
std::ofstream fout ("emissions.txt");
//std::ofstream fout ("concentration.txt");


int main (void)
{
  // NUMERICAL CONSTANTS //
  
  const double d = 8.64;
  const double u1 = 495.0;
  const double u2 = 0.0495;
  const double vs = 0.12;
  const double vd = 1.23;
  const double w = 0.001;
  const double k1 = 0.00219;
  const double k2 = 0.0000612;
  const double k3 = 0.997148;
  const double k4 = 0.0679;
  
  // INICIAL CONDITIONS //

  
  double p = 1.00;
  double Ss = 2.01; //Ss equals to sigma sub s
  double Sd = 2.23; 
  double As = 2.2; //As equals to alpha sub s
  double Ad = 2.26;
  double t = 1000;   //calculating year
  const double N = 2; //
  const double TA = 1000; //initial point //modify for user input for assigned year
  const double TB = 5000; // End point
  const double H = 0.01; //step
  std::vector<double> y = {0.12, 0}; //x0, y0

  rk4(TA, TB, H, y);

  return 0;
}



void rk4 (double ta, double tb, double h, std::vector<double> & y)
{
  std::vector<double>  k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  const int N = (tb-ta)/h;
  for (int nt = 0; nt < N; ++nt){
    double t = ta + h*nt;
    for (int ii = 0; ii < y.size(); ++ii){
      k1[ii] = h*f(t, y, ii);
    }
    
    //k2 aux
    for (int ii = 0; ii < y.size(); ++ii){
      aux[ii] = y[ii] + k1[ii]/2;
    }
    
    //k2
    for (int ii = 0; ii < y.size(); ++ii){
      k2[ii] = h*f(t + h/2, aux, ii);
    }
    
    //k3 aux
    for (int ii = 0; ii < y.size(); ++ii){
      aux[ii] = y[ii] + k2[ii]/2;
    }
    
    //k3
    for (int ii = 0; ii < y.size(); ++ii){
      k3[ii] = h*f(t + h/2, aux, ii);
    }
    
    //k4 aux
    for (int ii = 0; ii < y.size(); ++ii){
      aux[ii] = y[ii] + k3[ii];
    }
    
    //k4
    for (int ii = 0; ii < y.size(); ++ii){
      k4[ii] = h*f(t + h, aux, ii);
    }
    
    //write new y
    for (int ii = 0; ii < y.size(); ++ii){
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }
    fout << t << "\t" << y[0] << "\t" << y[1] << std::endl;
  }
}

double f(double t, const std::vector<double> &y, int id)
{
  if(0 == id){
    return y[1];
  }
  else if (1 == id){
    return -W*W*y[0];
  }
  else{
    std::cerr << "ERROR!!!!! it does not exist -> " << id << std::endl;
    exit(1);
  }
}
