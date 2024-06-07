#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <functional>
#include <string>

class Jacobi{
    private:
        // I didn't put max_it and tolr into the data.json file just for less communication and easier code
        int n;
        int it = 1;
        int max_it = 1000;
        double tolr = 1e-6;
        double error = 1.;
        double h;
        std::function<double(double&, double&)> f;

        std::vector<double> M;      // result matrix
        std::vector<std::vector<double>> out_M; // matrix only for output

    public:
        Jacobi(int n_, const std::function<double(double&, double&)> &func): n(n_), f(func){
            h=1/double(n-1);
            M.resize(n*n);
        };
        
        double get_error(){
            return error;
        }

        int get_it(){
            return it;
        }

        // MPI parallel solver
        void solve_parallel();

        // compute the error between my result M and the input function
        double test(std::function<double(double&, double&)>);

        // print M
        void print_res();

        // write the VTK file
        void output_VTK(std::string out_name);

};

#endif