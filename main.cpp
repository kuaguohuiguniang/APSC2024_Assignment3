#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <functional>

#include "json.hpp"
#include "Jacobi.hpp"
#include "muparser_fun.hpp"


// the RHS function
double f(double& x, double& y){
    return 8*std::pow(M_PI, 2)*std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
}

// the exact solution
double ex_f(double& x, double& y){
    return std::sin(2*M_PI*x)*std::sin(2*M_PI*y);
}

using json = nlohmann::json;

int main(int argc, char** argv){
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n;                      // nodes number in one row
    int fun_str_len;
    std::string fun_str;        // string for muparser to read f
    std::string ex_fun_str;     // string for muparser to read ex_f

    // rank 0 read the data.json file which contains n and the expression of f and ex_f
    if (rank == 0){
        std::ifstream file("data.json");
        json data = json::parse(file);

        n = data.value("n", 11);
        fun_str = data.value("func", "");
        ex_fun_str = data.value("ex_func", "");
        fun_str_len = fun_str.length();
        file.close();
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&fun_str_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // rank 0 boradcast the function string to other rank
    char* fun_str_buf = new char[fun_str_len + 1];
    if (rank == 0) {
        std::strcpy(fun_str_buf, fun_str.c_str());
    }
    MPI_Bcast(fun_str_buf, fun_str_len + 1, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank != 0){
        fun_str = std::string(fun_str_buf);
    }
    delete[] fun_str_buf;

    MuparserFun f(fun_str); // read f using muparser, but I get an error when compile

    std::function<double(double&, double&)> fun=f;  // instead, I use std::function
    Jacobi Jacobi_solver(n, fun);                   
    Jacobi_solver.solve_parallel();                 // mpi parallel function
    
    if (rank == 0){
        MuparserFun ex_f(ex_fun_str);
        std::function<double(double&, double&)> ex_fun=ex_f;
        std::cout<<"test error: "<<Jacobi_solver.test(ex_fun)<<std::endl; // error between my result and the exact solution
        std::cout<<"error: "<<Jacobi_solver.get_error()<<std::endl;       // error between my result and the result in last iteration
        std::cout<<"it: "<<Jacobi_solver.get_it()<<std::endl;             // iteration number 
        Jacobi_solver.output_VTK("output.VTK");                           // write the result in 
    }
    
    MPI_Finalize();
    return 0;
}
