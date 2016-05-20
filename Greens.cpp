/////////////////////////////////////////////////////////
///				Greens.cpp					  5-19-16 ///
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
    /* 
    // NFW
    double rho = rhos / ( r/rs * pow(1 + r/rs , 2 )); //x? + 1e-100
	*/

    /*
    // Berkert
    double rho = rhos / ((1+ r/rc) * (1 + pow(r/rc , 2 ))); //shouldnt use rhos here, different rho
	*/

	// N04
	double rho = rhos * exp(-2.0/0.17 * ( pow(r/rc, 0.17) - 1 ));
 




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
	gsl_integration_qags (&F, 1e-16, rh, 0, 1e-3, 1000, //x?
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
		//std::cout << i << std::endl;
		Gsum += pow(-1, i) * gslInt_Greens(ri, r, root_dv, rh);

	}

	double Greens = pow(4*pi , -1.0/2.0)*Gsum ;
	//std::cout << Greens << std::endl;

	return Greens;

};

void runGreens(double r){ //runs through values of root_dv
	
	double root_dv = 0.04*mpc2cm;
	int n = 1000;
	double rk = mpc2cm/1000*r;
	double h = root_dv/n;

	std::ostringstream makefilename;
	makefilename << "Greens_N04_r" << r << "kpc" <<".txt" ;
	std::string filename = makefilename.str();
	std::ofstream file(filename.c_str());

	for (int ix = 1 ; ix < n + 1; ++ix  ){

		double data = Greens(rk, h*ix );
		file << h*ix/mpc2cm*1000 << "\t" <<  data <<std::endl;
		//std::cout << h*ix/mpc2cm*1000 << "\t" << data << std::endl;

	};
	
}

main(){

	runGreens(1);
	std::cout << "1/7" <<std ::endl;
	runGreens(5);
	std::cout << "2/7" <<std ::endl;
	runGreens(10);
	std::cout << "3/7" <<std ::endl;
	runGreens(50);
	std::cout << "4/7" <<std ::endl;
	runGreens(100);
	std::cout << "5/7" <<std ::endl;
	runGreens(200);
	std::cout << "6/7" <<std ::endl;
	runGreens(500);
	std::cout << "7/7" <<std ::endl;
}