//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * darksusy_int.cpp *                                              4-13-2016 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/* First must compile dsinit.f (do this in working directory) using 
gfortran -c ~/research/darksusy-5.1.2/src/ini/dsinit.f  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -ldarksusy -lFH -lHB -lgfortran

*/


/*
in order to compile and run, use
g++ -o darksusy_int darksusy_int.cpp dsinit.o constants.o particle.o cluster_params.o astro.o synchrotron.o -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -ldarksusy -lFH -lHB -lgfortran

*/


#include <iostream>

using namespace std;




extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {
double dshayield_(double *mwimp, double *emuthr,int *ch,  int *yieldk, int *istat);
}



int main()
{

   double r = .01;
   double mx = 50;
   double emuthr1 = 0.1;
   int ch = 25;
   int yieldk = 151;
   int istat;




   dsinit_();	//initialize darksusy

   for(int i = 0 ; i  < 10 ; ++i){
   	double emuthr1 = i * 0.1; 
   dshayield_(&mx, &emuthr1, &ch,  &yieldk, &istat);			
   std::cout << dshayield_(&mx, &emuthr1, &ch,  &yieldk, &istat)<< std::endl;
  };
   return 0;
}