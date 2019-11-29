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
usigin fptr = double(double, double, std:: vector<double>)
//using  vec = double( std::vector<double>) ;
//using  dou = double(vec, std::vector <double>));
//using  raro = double(vec, dou, double));


double function_f(std::vector<double> &c, std::vector<double> &t, std::vector<double> &s); // function of emition of CO2 in to the atmosphere
double p_s(double & h_s, double & c_s);
double f(double t, const std::vector<double> &y, int id);
void rk4 (double ta, double tb, double h, std::vector<double> & y);
std::ofstream fout ("dato.txt");





//INT MAIIIIIMMMM//
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
  std::vector<double> y = {1.0, 2.01, 2,23, 2.2, 2,26}; //x0, y0


  // LLAMAR A P_S RAROOO
  rk4(TA, TB, H, y, p_s);

  return 0;
}






void rk4 (double ta, double tb, double h, std::vector<double> & y, fptr p_s)
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
double f(double t, const std::vector<double> &y, double k1, double k2, double k3, double k4, double u1, double u2, fptr p_s)
{
  
  
  if(0 == id){
    
    return (p_s(k4,k3,y)/d + function_f(c,t,s)/u1 ;
  }
  else if (1 == id){
    return (w*(Y[2]-y[1])-k1-u2*((rp_s-y[0])/d))/vd ; //Reevisar sintaxis, pero teóricamente diana está bien :v
  }
   else if (2 == id){
     return (k1-w*(Y[1]-y[2]))/vd ; //y[1]=sigma_d y[2]=sigma_s
  }
   else if (3 == id){
     return (w*(y[4]-y[3]-k2)/vs ; //y[4]= alpha_d , y[3]= alpha_s
  }
   else if (4 == id){     return (k2-w*(y[4]-y[3])/vd;
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


	    double p_s(double k4, double k3, const sdt::vector <double> & y)
      {
	double rp_s = 0.0;
	double rc_s = 0.0;
	double rh_s = 0.0;
	
	rh_s = (y[1]-std::sqrt(y[1]*y[1]-k3*y[3]*(2*y[1]-y[3]))/k3;     
		
	rc_s = (y[3]-rh_s/2;
			
	rp_s = k4 * ((rh_s*rh_s)/c_s);
			return rp_s;		
}

  
  
