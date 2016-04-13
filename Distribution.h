#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include<vector>


class Distribution {
	private:
		
		std::vector <double> data;
	public:
		int nx , ny ;
		Distribution(int nx, int ny) : nx(10), ny(10), data(nx*ny, 0) {}  // access data elements by (x, y) = (nx*y + x)

};




#endif