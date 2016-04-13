//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * darksusy_int.cpp *                                              4-13-2016 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/* First must compile dsinit.f (do this in working directory) using 
gfortran -c ~/research/darksusy-5.1.2/src/ini/dsinit.f

*/


/*
in order to compile and run, use
g++ -o <filename> <filename>.cpp dsinit.o constants.o particle.o cluster_params.o astro.o synchrotron.o
 -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -ldarksusy -lFH -lHB -lgfortran

*/



#include <iostream>

using namespace std;




extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	


extern"C" {
double __synchrotron_MOD_ssyn(double *r);
}










int main()
{

   double r = .01;
   dsinit_();							//initialize darksusy
   std::cout << __synchrotron_MOD_ssyn(&r)<< std::endl;
  
   return 0;
}