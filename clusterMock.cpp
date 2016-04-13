





#include<vector>
#include<iostream>



class Distribution {
    public:


		int nx , ny ;
		std::vector <double> data;

		
		Distribution(int nx, int ny) : nx(10), ny(10), data(nx*ny, 0) {}  // access data elements by (x, y) = (nx*y + x)

};



class Cluster {
    private:
        Distribution bfield_;
        Distribution energies_;

    public:
    	Cluster();
    	Distribution bfield;
    	
};



int main() {
    Cluster c;
    c.bfield.data[0] = 1;
};