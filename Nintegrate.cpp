#include <iostream>
#include <math.h>


double f(double x){
	
	double val = 1/x;
	
	return val;


}


double NIntegrate(double rf, double r0, int ns){

	double h = (rf - r0)/ns;
	double simp1 = f( r0 + h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += f( r0 + h * i + h/2);
		simp2 += f( r0 + h * i);

	}

	double result = h/6 * (f(r0) + f(rf) + 4 * simp1 + 2 * simp2);

	std::cout << result << std::endl;
	return result;




}

double whileNIntegrate(double rf, double r0, int ns){

	double h = (rf - r0)/ns;
	double simp1 = f( r0 + h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += f( r0 + h * i + h/2);
		simp2 += f( r0 + h * i);

	}

	double result = h/6 * (f(r0) + f(rf) + 4 * simp1 + 2 * simp2);

	std::cout << result << std::endl;
	return result;




}





main(){
	int ns;
	std::cout << "Number of steps?" <<std::endl;
	std::cin >> ns;
	NIntegrate(100,1, ns);
	//std::cout << f(2) <<std::endl;
}