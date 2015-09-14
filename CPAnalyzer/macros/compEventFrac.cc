#include <iostream>
#include <math.h>

void computeFrac(double mean1, double sigma1, double mean2, double sigma2, double binwidth)
{
  double N1 = mean1 * fabs(sigma1) * sqrt(2*M_PI) / binwidth;
  double N2 = mean2 * fabs(sigma2) * sqrt(2*M_PI) / binwidth;

  std::cout<<"Number of events = "<<N1<<" , "<<N2<<std::endl;
  std::cout<<"Fraction of events = "<<N1/(N1+N2)<<" , "<<N2/(N1+N2)<<std::endl;

}

int main ()
{
  //Without using BS 
  /*
  std::cout<<"Events for DeltaZ"<<std::endl;
  computeFrac(0.0223012, 0.00233353, 0.00414710, 0.00611880, 0.0002);

  std::cout<<"Events for DeltaX"<<std::endl;
  computeFrac(0.0296489, 0.00164184, 0.00740740, 0.00394297, 0.0002);

  std::cout<<"Events for DeltaY"<<std::endl;
  computeFrac(0.0284937, 0.00158672, 0.00841945, 0.00386476, 0.0002);
  */

  //std::cout<<"Events for DeltaPCA_tdp"<<std::endl;
  //computeFrac(0.0183726, 0.00157436, 0.00296559, 0.00347566, 0.0001);
  
  std::cout<<"Events for DeltaPCA_t"<<std::endl;
  computeFrac(0.0190395, 0.00144168, 0.00375728, 0.00313530, 0.0001);

}
