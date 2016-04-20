//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * darksusy_int.cpp *                                              4-13-2016 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/* First must compile dsinit.f (do this in working directory) using 
gfortran -c ~/research/darksusy-5.1.2/src/ini/dsinit.f  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -ldarksusy -lFH -lHB -lgfortran

*/


/*
in order to compile and run, use
g++ -o darksusy_int darksusy_int.cpp dsinit.o constants.o particle.o cluster_params.o astro.o synchrotron.o -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran

*/


#include <iostream>
#include <stdio.h> //only for printf
#include <gsl/gsl_integration.h>

using namespace std;




extern"C" {	 								//interface with fortran code to initialize darksusy

void dsinit_();

}	

extern"C" {
double dshayield_(double *mwimp, double *emuthr,int *ch,  int *yieldk, int *istat);
}



int Set_ch(){
   double ch = 25;
   return ch;
}



double ds (double Ep, void * params){

   int ch = 13;
   int yieldk = 151;
   int istat;
   double mx = *(double *)params; 
   
   double ds = dshayield_(&mx, &Ep, &ch, &yieldk, &istat);
   //std::cout << ds << std::endl;

   return ds;
}











void gslInt_ds( double i, double f, double mx){
     gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (10000);

  double result, error;
  double expected = -4.0;
  

  gsl_function F;
  F.function = &ds;
  F.params = &mx;

  gsl_integration_qags (&F, i, f, 0, 1e-4, 10000,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);

  printf ("estimated error = % .18f\n", error);

  printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);

}



















int main()
{
    double mx_max = 1000;
    double mx = 5; 


    dsinit_();	//initialize darksusy

    gslInt_ds(0,mx_max, mx);
    return 0;
}