/*
Using RK4 for solving the differential equation for concentration on CO2 in the atmosphere

  // REPRODUCIR GRAFICAS
  // QUE METODO ES MEJOR(COMPARAR CONCENTRACIONES INICIALES Y FINALES)
  // CUANDO LA CONCENTRACION DE CO2 ES MAXIMA?

*/

#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>

double function_f(std::vector<double> &c, std::vector<double> &t, std::vector<double> &s, double t); // function of emition of CO2 in to the atmosphere
double h_s(double Ss);
double c_s(double & h_s);
double p_s(double & h_s, double % c_s);
double d_p();
double d_Ss();
double d_Sd();
double d_As();
double d_Ad();

double f(double t, const std::vector<double> &y, int id);
void rk4 (double ta, double tb, double h, std::vector<double> & y);
std::ofstream fout ("dato.txt");

int main (void)
{
  // NUMERICAL CONSTANTS //
  
  const double d = 8.64 , u1 = 495.0 , u2 = 0.049 , vs = 0.12 , vd = 1.23 , w = 0.001 , k1 = 0.00219 , k2 = 0.0000612 , k3 = 0.997148 , k4 = 0.0679;
  
  // INICIAL CONDITIONS //

  
  double p = 1.00;
  double Ss = 2.01; //Ss equals to sigma sub s
  double Sd = 2.23; 
  double As = 2.2; //As equals to alpha sub s
  double Ad = 2.26;
  double t = 1000;   //calculating year

  // DATA STRUCTURE //

  std::vector<double> c{2.0, 10.5, 2.7}
  std::vector<double> t{1988, 2100, 2265}
  std::vector<double> s{21, 96, 57}
  
  // OPERATING DATA //
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

double function_f(std::vector<double> &c, std::vector<double> &t, std::vector<double> &s, double t)
{
  double suma = 0
  
  for(int ii=0; ii<3; ++ii)
  {
    suma+= c[ii]*std::exp((-(T-t[ii])*(T-t[ii]))/(s[ii]*s[ii]));
  }
  std::cout<<suma<<std::endl;
}

double p_s(double & h_s, double % c_s)
{
  double p_s = 0.0;
  p_s = k4 * ((h_s*h_s)/c_s);
  return p_s;
}

double c_s(double & h_s)
{
  double c_s = 0.0;
  c_s = (As-h_s)/2;
  return c_s;
}

double h_s(double Ss)
{
  double h_s = 0.0;
  h_s = (Ss-std::sqrt(Ss*Ss-k3*As*(2*Ss-As))/k3;
  return h_s;
}
  
