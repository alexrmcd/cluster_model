#ifndef CLUSTER
#define CLUSTER

#include<vector>

class Grid {
    private:
        int nx, int ny;
        std::vector<double> data;
    public:
        Grid(int nx, int ny) : nx(nx), ny(ny), data(nx*ny, 0) {}
        double at(int x, int y) {
            int index = y*nx + x;
            if (index > data.size()) {
                throw std::runtime_error("Hey thats not a valid index!");
            }
            return data[index];
        }
        int nx() { return nx; }
        int ny() { return ny; }
};

class Cluster {
    private:
        Grid bfield_;
        Grid energies_;

        std::vector<std::vector<double>> grid;
    public:
        Cluster(const Grid& bfield) :
            bfield_(bfield),
            energies_(bfield.nx(), bfield.ny());
            {}
        Cluster() : Cluster(Grid(10, 10)) {}
        void simulate(); // run the simulation
        double energy_at(int x, int y) { return energies_.at(x, y) };
        Grid energies() { return energies_; }
};

#endif
