#include<fstream>
#include<vector>
#include<cmath>

void rkta(double t0, double t1, double dt, std::vector<double> &y);
double f(double t, const std::vector<double> &y, int i);
double fuel(double t);
double ps(const std::vector<double> &y);
std::ofstream fout("Datos.txt");

//constants
const double d=8.64, u1=495.0 , u2=0.0495, vs=0.12, vd=1.23, w=0.001, k1=0.000219, k2=0.0000612, k3=0.997148, k4 = 0.0679;

int main(void)
{
	const int N=5;
	double t0=1000.0;
	double t1=5000;
	double dt=1;
	std::vector<double> y(N);
	y={1,2.01,2.23,2.20,2.26};
	
	rkta(t0,t1,dt,y);
	
	return 0;
}

void rkta(double t0, double t1, double dt, std::vector<double> &y)
{
	int N=(t1-t0)/dt;
	double max=0;
	double tmax=0;
	double t=0;
	std::vector<double> k1, k2, k3, k4, aux, kaux;
	k1.resize(y.size());
	k2.resize(y.size());
	k3.resize(y.size());
	k4.resize(y.size());
	aux.resize(y.size());
	kaux.resize(y.size());

	
	fout<< t0 <<"\t";
	for(int i=0; i<3; i++){
	  fout<< y[i] <<"\t";
	}
	fout<<fuel(t0)<< "\n";
	
	for(int n=1; n<=N; n++){
	  t=t0+dt*n;
	  std::copy(y.begin(),y.end(),aux.begin());
	  for(int i=0; i<y.size(); i++){
	    //k1
	    k1[i]=dt*f(t,aux,i);
	    //k2
	    for(int j=0; j<y.size(); j++){
	      kaux[j]=aux[j]+k1[j]/2;
	    }
	    k2[i]=dt*f(t+dt/2,kaux,i);
	    //k3
	    for(int j=0; j<y.size(); j++){
	      kaux[j]=aux[j]+k2[j]/2;
	    }
	    k3[i]=dt*f(t+dt/2,kaux,i);
	    //k4
	    for(int j=0; j<y.size(); j++){
	      kaux[j]=aux[j]+k3[j];
	    }
	    k4[i]=dt*f(t+dt,kaux,i);
	    //t evolution
	    y[i]=y[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
	    
	  }
	  if(y[0]-max > 1e-6){
	    max=y[0];
	    tmax=t;
	  }
	  
	  
	  
	  //impresion de t,p, Ss,Sd,f(t)
	  fout<< t <<"\t";
	  
	  
	  for(int i=0; i<3; i++){
	    fout<< y[i] <<"\t";
	    
	  }		
	  fout<< fuel(t)<< "\n";	
	}	
	fout<< tmax << "\t"<< max << "\n";
fout.close();
}

double f(double t, const std::vector<double> &y, int i)
{
	if(i==0){
		return (ps(y)-y[0])/d+fuel(t)/u1;
	}
	else if(i==1){
		return (w*(y[2]-y[1])-k1-u2*(ps(y)-y[0])/d)/vs;
	}
	else if(i==2){
		return (k1-w*(y[2]-y[1]))/vd;
	}
	else if(i==3){
		return (w*(y[4]-y[3])-k2)/vs;
	}
	else if(i==4){
		return (k2-w*(y[4]-y[3]))/vd;
	}
}

double fuel(double t)
{
	double result=0;
	const std::vector<double> c{2,10.5,2.7}, t1{1988,2100,2265}, s{21,96,57};
	
	for(int i=0; i<3; i++){
		result+=c[i]*std::exp(-std::pow((t-t1[i]),2)/std::pow(s[i],2));
	}
	
	return result;
}

double ps(const std::vector<double> &y)
{
	double hs=(y[1]-std::sqrt(std::pow(y[1],2)-k3*y[3]*(2*y[1]-y[3])))/k3;
	double cs=(y[3]-hs)/2;
	
	return k4*(std::pow(hs,2)/cs);
}
//

