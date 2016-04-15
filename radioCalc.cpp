////  adioCalc.cpp
/////c++ version of storm code 4-15-16
////


#include <iostream>
#include <math.h>
#include <cmath>



double distint(double z){
	double mpc2cm = 3.085678e24; 	//Mpc to cm
	double clight = 299792.458;		// km/s
	double H0 = 70.4;				// km/s/Mpc
	double OmegaM = 0.272 ; 		
	double OmegaL = 0.728 ;

	double distint = mpc2cm * clight / ( H0 * sqrt( OmegaM * pow(1 + z , 3)  + OmegaL ) ); 

	return distint;

}


//distance as a function of redshift
double Dist(double z){ 

	double mpc2cm = 3.085678e24; 	//Mpc to cm
	double clight = 299792.458;		// km/s
	double H0 = 70.4;				// km/s/Mpc
	double OmegaM = 0.272 ; 		
	double OmegaL = 0.728 ;

	int ns = 1000000;

	double h = z/ns;
	double simp1 = distint(h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += distint(h * i + h/2);
		simp2 += distint( h * i);

	}

	double result =  h/6 * (distint(0) + distint(z) + 4 * simp1 + 2 * simp2);


	



	return result;
 		
}




//SYnchrotron emmission spectral function form Cola2006.
double fff(double x){

	double fff = 1.25 * pow( x , 1/3) * exp( -x ) * pow((648 + x*x) , 1/12);


	return fff;
}






double psynIntegrate(double theta){


	double bmu = 1; 	//muGauss
	double E = 50;
	double z = 0.0232;


	double nu = 1.4;	//Ghz

	double psyn0 = 1.463232e-25 ; // Gev/s/Hz
	double x0 = 62.1881 ;			// dimensionless constant
	double nu_em = ( 1 + z )* nu; // (observing freq)*(1+z)



	double x = x0 *nu_em / (bmu * E*E) ;
	double psynIntegrate = psyn0 * bmu * 0.5 * pow( sin(theta) , 2) * fff( x  / std::abs(sin(theta)) ); // fff(x) doesnt like negative arguments, sin(pi) ~ -2e-13, take abs for now(requires cmath)

	return psynIntegrate;


}








double psyn(){

	double pi = 3.14159265359;
	double dndeeq = 1;
	int ns = 1000;

	double h = (pi)/ns;
	double simp1 = psynIntegrate( h/2);

	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += psynIntegrate( h * i + h/2);
		simp2 += psynIntegrate(  h * i);

	}

	double result = h/6 * (psynIntegrate(0) + psynIntegrate(pi) + 4 * simp1 + 2 * simp2);

	std::cout << result << std::endl;
	



	double psyn = 2* result * dndeeq;

	return psyn;

}


































main () {

	//std ::cout << Dist(0.0232) << std::endl;
	double pi = 3.14159265359;
	
	//std::cout << psyn()<<std::endl;	


}