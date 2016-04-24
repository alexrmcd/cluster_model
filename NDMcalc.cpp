//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * darksusy_int.cpp *                                              4-13-2016 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

/* First must compile dsinit.f (do this in working directory) using 
gfortran -c ~/research/darksusy-5.1.2/src/ini/dsinit.f  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -ldarksusy -lFH -lHB -lgfortran

*/


/*
in order to compile and run, use
g++ -o NDMcalc NDMcalc.cpp  -I/home/alex/research/darksusy-5.1.2/include -L/home/alex/research/darksusy-5.1.2/lib  -lgsl -lgslcblas -ldarksusy -lFH -lHB -lgfortran

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

class Cluster{  
  
  public:
  std::string name;
  double z ;
  double rh; //halo radius
  int ch ;
  double B0 ;         //microGauss
  double rcore ;    //0.291; //core radius in Mpc for coma
  double root_dv;

  Cluster() : name(""),
        z(0),             //redshift
        rh(0),            //halo radius Mpc
        ch(0),            //darksusy channel
        B0(0), rcore(0) ,     //Bfield Params
        root_dv(0.035)          //DIffusion parameter, not really a cluster thing but easy access is good

  { //everything labelled or Coma, sgould work in user options
    std::cout << "creating cluster... " << std::endl;
  }

};

Cluster c;


double ddarksusy (double Ep, void * params){

   int yieldk = 151;
   int istat;
   double mx = *(double *)params; 
   
   double ds = dshayield_(&mx, &Ep, &c.ch, &yieldk, &istat);


   return ds;
}



double inject( double Ep, double mx){
    
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (5000);

  double result, error;

  gsl_function F;
  F.function = &ddarksusy;
  F.params = &mx;

  gsl_integration_qags (&F, Ep, mx, 0, 1e-3, 5000, 
                      w, &result, &error); 

  gsl_integration_workspace_free (w);
  return result;

}



double DM_profile(double r){

    double rhos = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
    double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)


    double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100

    return rho;

};


double bloss(double E , double r){

  double bloss = 0.0253*pow(bfield_model(r), 2)*E*E           //bsyn bfield_model(r)
          + 0.265 * pow(1 + c.z, 4 )*E*E            //bIC
          + 1.51*n*(0.36 + log(E/me/n) )      //+ 1.51*n*(0.36 + log(E/me) )*E            //bbrem Note this is most likely incorrect, no factor of E, and should be E/me/nin log
          + 6.13*n*( 1 + log(E/me/n)/75);           //bcoul

  bloss = bloss *1e-16; 

  return bloss;
};




double spectrum(){  // dn/dE = <sv> rho^2 / (2 * mx^2 * bloss(){})inject()



}












main()
{
    double mx_max = 1000;
    double mx = 5; 


    dsinit_();	//initialize darksusy

    //gslInt_ds(0,mx_max, mx);
}