#include <iostream>
#include <math.h>
#include <fstream>



main(){

		//outout path for files ~/research/cluster_model/output/diffuisionEq
		std::ofstream test("test.txt");	
		//write test file
		for (int i =0 ; i < 11 ; ++i){

			test << i << std::endl;
		}

		test.close();


}



//////Scrap greens work


//inegrand


