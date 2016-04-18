///////////// gsl_test.cc ///////////


////////// need to compile with g++ -o <file> <file>.cc -lgsl -lcblas 

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <ctime>

double f (double x, void * params) {
  double alpha = *(double *) params;
  double f = log(alpha*x) / sqrt(x);
  return f;
}



double b(double x, void * params){
	
	double val = x*x*log(x + 1);
	
	return val;


}

double g(double x){
	
	double val = x*x* log(x+1);
	
	return val;


}


double NIntegrate(double rf, double r0, int ns){

	double h = (rf - r0)/ns;
	double simp1 = g( r0 + h/2);
	double simp2 = 0;

	for (int i = 0 ; i < ns ; ++i){
		simp1 += g( r0 + h * i + h/2);
		simp2 += g( r0 + h * i);

	}

	double result = h/6 * (g(r0) + g(rf) + 4 * simp1 + 2 * simp2);

	std::cout << result << std::endl;
	return result;




}


int
main (void)
{

	std::clock_t gsl_start;

	double gsl_duration;

	gsl_start = std::clock();

	///// start gsl
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (100);
  
  double result, error;
  double expected = -4.0;
  double alpha = 1.0;

  gsl_function F;
  F.function = &b;
  F.params = &alpha;

  gsl_integration_qags (&F, 0, 1, 0, 1e-7, 100,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);

 	////////after algorithm
	gsl_duration = (std::clock()  -  gsl_start)/(double) CLOCKS_PER_SEC;

	std::cout << "gsl time:  " << gsl_duration <<std::endl;



//// time Nintegrate

  	std::clock_t start;

	double duration;

	start = std::clock();


	int ns = 100;

	NIntegrate(1,1e-16, ns);
	//std::cout << f(2) <<std::endl;

	////////after algorithm
	duration = (std::clock()  -  start)/(double) CLOCKS_PER_SEC;

	std::cout << "simp time:  " << duration <<std::endl;


	 return 0;
}
