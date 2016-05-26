/////////////////////////////////////////////////////////
///				GreensLUT.cpp					  5-25-16 ///
/////////////////////////////////////////////////////////

/*
g++ -o Greens Greens.cpp  -lgsl -lgslcblas

*/



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <iomanip> 


#include "Constants.h"




double DM_profile(double r){

    double rhos = 0.039974 ;//DM char density in Gev /cm^3 Storm13 (Mathematica)
    double rs = 0.404*mpc2cm;  //DM scale radius in Mpc Storm13  (Mathematica)
    double rc = 0.241 *mpc2cm;
    
    // NFW
    double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100
	

    /*
    // Berkert
    double rho = rhos / ((1+ r/rc) * (1 + pow(r/rc , 2 ))); //shouldnt use rhos here, different rho
	*/
/*
	// N04
	double rho = rhos * exp(-2.0/0.17 * ( pow(r/rc, 0.17) - 1 ));
 */




    return rho;

};



double dGreens(double rp, void * params ){ 

	std::vector<double> greenParam = *(std::vector<double> *)params;
	double ri = greenParam[0];
	double r = greenParam[1];
	double root_dv = greenParam[2]; // 


	double dGreens = pow(root_dv , -1) * rp/ri * (exp( - pow( (rp-ri)/(2*root_dv) , 2)) 
		- exp( - pow( ( rp + ri)/(2*root_dv) , 2)) ) * pow( DM_profile(rp),2)/pow( DM_profile(r),2);

	return dGreens;

}



double gslInt_Greens(double ri , double r, double root_dv, double rh){


	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);

	double result, error;

	std::vector<double> greenParam (3);


	greenParam[0] = ri;
	greenParam[1] = r;
	greenParam[2] = root_dv;

	gsl_function F;
	F.function = &dGreens;
	F.params = &greenParam;
	gsl_set_error_handler_off();
	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-5, 1000, //x?
	                w, &result, &error); 

	gsl_integration_workspace_free (w);

	return result;

}


double Greens (double r, double root_dv) {  //called by ddsyn

	double rh = 0.415 * mpc2cm ;
	int imNum = 7; //number of image pairs
	double Gsum = 0 ;
	for (int i = - imNum; i < imNum + 1; ++i ){


		double ri;
		
		if (i == 0)
			ri = r;
		else
			ri = (pow(-1 , i)*r + 2*i*rh);

		Gsum += pow(-1, i) * gslInt_Greens(ri, r, root_dv, rh);

	}

	double Greens = pow(4*pi , -1.0/2.0)*Gsum ;


	return Greens;

};



std::vector<double> createLUT(double n_r, double n_rootdv){
	// iteration timer start
	std::clock_t Gstart;
	double Greensduration;
	Gstart = std::clock();
	int Ga ; 
	///////before algorithm
	std::cout << "creating LUT..." <<std::endl; 

	std::vector<double> Greens_lookup( n_r * n_rootdv );

	double rh = 415*mpc2cm/1000;
	double r_scale = rh/n_r;

	double rootdv_max = 40*mpc2cm/1000;
	double rootdv_scale = rootdv_max/n_rootdv;
	

	std::cout << r_scale/mpc2cm*1000 << " , rdv_max = "<<rootdv_max<<" n_rootdv = " <<n_rootdv <<" rdv_scale = " << rootdv_scale/mpc2cm*1000 << std::endl; 
	

	for (int j = 1 ; j < n_r ; ++j ){
		
		double dr = rh/n_r;
		
		for(int i = 1; i < n_rootdv ; ++i){
			double drootdv = rootdv_max/n_rootdv;

			Greens_lookup[j + n_r*i] = Greens(j*dr, i*rootdv_scale);

			
		}

		std::cout << j << "/" << n_r <<std::endl;

	////////after algorithm
	Greensduration = (std::clock()  -  Gstart)/(double) CLOCKS_PER_SEC;
	std::cout << "Glookup time = " << Greensduration <<std::endl;

	return Greens_lookup;
}


void runGreens(double r){ //runs through values of root_dv
	int n = 10000;
	int n_r = 415;
	int n_rootdv = 1000;
	double root_dv = 0.04*mpc2cm;

	double rk = mpc2cm/1000*r;
	double h = root_dv/n;

	double rh = 415*mpc2cm/1000;
	double rootdv_max = 40*mpc2cm/1000;
	double r_scale = rh/n_r;
	double rootdv_scale = rootdv_max/n_rootdv;


	std::vector<double> glookup = createLUT(n_r, n_rootdv);
	std::ostringstream makefilename;
	makefilename << "Greens_NFW_r" << r << "kpc" <<"GLUTtest.txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int ix = 0 ; ix < n ; ++ix  ){
		double Gi  = (int)(rk/r_scale) +(int)( h*ix/rootdv_scale)*n_r;
		
		double data = glookup[Gi];

		file << h*ix/mpc2cm*1000 << "\t" <<  data <<std::endl;
		std::cout << h*ix/mpc2cm*1000 << "\t" << data << std::endl;


		if(ix % 1000 == 0)
			std::cout << ix/1000 << "/" << n /100000<< std::endl;;		


	};
	
}

main(){





//createLUT(415, 1000);
/*
	runGreens(1);
	std::cout << "1/2" <<std ::endl;
	runGreens(5);
	std::cout << "2/7" <<std ::endl;
	runGreens(10);
	std::cout << "3/7" <<std ::endl;
	runGreens(50);
	std::cout << "4/7" <<std ::endl;
	runGreens(100);
	std::cout << "2/2" << std::endl;*/
	//std::cout << "5/7" <<std ::endl;
	runGreens(200);
	std::cout << "6/7" <<std ::endl;
	//runGreens(500);
	//std::cout << "7/7" <<std ::endl;


}