#include<iostream>
#include<vector>
#include<cmath>

int main(void)
{

 double T=1000;
 double suma=0;
 std::vector<double> c{2.0 ,10.5 , 2.70};
 std::vector<double> t{1988 , 2100 , 2265 };
 std::vector<double> s{21.0 , 96.0 , 57.0};
  
  for(int ii=0; ii<3; ++ii)
  {
    suma+= c[ii]*std::exp((-(T-t[ii])*(T-t[ii]))/(s[ii]*s[ii]));
  }
  std::cout<<suma<<std::endl;
}
