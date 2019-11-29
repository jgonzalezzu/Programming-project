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

using  vec = double(std::vector<double>) ;
using  dou = double(vec, std::vector <double>));
using  raro = double(vec, dou, double));


double function_f(std::vector<double> &c, std::vector<double> &t, std::vector<double> &s); // function of emition of CO2 in to the atmosphere
double h_s(const double & y );
double c_s(double & h_s);
double p_s(double & h_s, double & c_s);
double f(double t, const std::vector<double> &y, int id);
void rk4 (double ta, double tb, double h, std::vector<double> & y);
std::ofstream fout ("dato.txt");


// NUMERICAL CONSTANTS //
  
  const double d = 8.64 , u1 = 495.0 , u2 = 0.049 , vs = 0.12 , vd = 1.23 , w = 0.001 , k1 = 0.00219 , k2 = 0.0000612 , k3 = 0.997148 , k4 = 0.0679;
  
  // INICIAL CONDITIONS //
  double p = 1.00;
  double Ss = 2.01; //Ss equals to sigma sub s
  double Sd = 2.23; 
  double As = 2.2; //As equals to alpha sub s
  double Ad = 2.26;
  double T = 1000;   //calculating year

  // DATA STRUCTURE //
  std::vector<double> c{2.0, 10.5, 2.7};
  std::vector<double> t{1988, 2100, 2265};
  std::vector<double> s{21, 96, 57};
  
  // OPERATING DATA //
  const double N = 2; //
  const double TA = 1000; //initial point //modify for user input for assigned year
  const double TB = 5000; // End point
  const double H = 0.01; //step
  std::vector<double> y = {0.12, 0}; //x0, y0


//INT MAIIIIIMMMM//
int main (void)
{


  // LLAMAR A P_S RAROOO
  rk4(TA, TB, H, y, p_s);

  return 0;
}






void rk4 (double ta, double tb, double h, std::vector<double> & y, raro ps)
{
  std::vector<double>  k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  const int N = (tb-ta)/h;
  for (int id = 0; nt < N; ++id){
    double t = ta + h*id;//MODIFICA EL TIEMPO EN EL QUE ESTA MODIFICAR!!!!!!!
    for (int ii = 0; ii < y.size(); ++ii){
      k1[ii] = h*f(t, y, ii, ps);
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
    std::cout << t << "\t" << y[0] << "\t" << y[1] << "\t"<< y[2]<< "\t" <<y[3] <<  "\t" <<y[4] <<   std::endl;
  }
}
//FUNCIONES DENTRO DEL ROUGE KUTTA
double f(double t, const std::vector<double> &y, raro ps, vec hs, dou cs)
{
  double PS=ps(hs(y),cs(hs(y),y),k4)-y[0]);
  
  if(0 == id){
    
    return (PS/d + function_f(c,t,s)/u1 ;
  }
  else if (1 == id){
    return (w(Y[2 ;
  }
   else if (1 == id){
    return ;
  }
   else if (1 == id){
    return ;
  }
   else if (1 == id){
    return ;
  }
  else{
    std::cerr << "ERROR!!!!! it does not exist -> " << id << std::endl;
    exit(1);
  }
}

double function_f(std::vector<double> &c, std::vector<double> &t, std::vector<double> &s)
{
  double suma = 0;
  
  for(int ii=0; ii<3; ++ii)
  {
    suma+= c[ii]*std::exp((-(T-t[ii])*(T-t[ii]))/(s[ii]*s[ii]));
  }
  std::cout<<suma<<std::endl;
}

double p_s(vec h_s, dou  c_s, double k4)
{
  double rp_s = 0.0;
  rp_s = k4 * ((h_s(y)*h_s(y))/c_s(h_s ,y));
  return rp_s;
}

 double c_s(vec h_s, const sdt::vector <double> & y)
{
  double rc_s = 0.0;
  rc_s = (y[3]-h_s(y)/2;
  return rc_s;
}

    double h_s(const sdt::vector <double> & y )
{
  double rh_s = 0.0;
  rh_s = (y[1]-std::sqrt(y[1]*y[1]-k3*y[3]*(2*y[1]-y[3]))/k3;
  return rh_s;
}
  
