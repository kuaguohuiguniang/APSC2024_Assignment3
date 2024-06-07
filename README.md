# APSC2024_Assignment3
To compile, just type "make";
To run, "mpirun -np 3 -x OMP_NUM_THREADS=3 ./main". One can freely choose the number of ranks and threads.

data.json contains the nodes number in a row, expression of f and expression of the exact solution. I didn't put the tolrance and maximum iteration number in it just to simply the main file. If one want to change them, you can go to Jacobi.hpp. 
Jacobi.hpp and Jacobi.cpp are the core of my code. Jacobi.solve_parallel() implement the Jacobi iteration mathod with MPI and Openmp. The idea for communication is that add extra rows to the local matrix and they are update by communicating with neighbor rank. 
Moreover, I modify the muparser_fun.hpp wich is given in Lab3 to make it can read a function with two varibles. 

Problem I encounter: I can't use Muparser to read the string which contains the function expression. Therefore, I comment the line "//MuparserFun fun(fun_str);" in main.cpp and use std::function instead. If I use Muparser to read the function, when I compile my files with "make", a long and unreadble error will raise. I believe the problem is either in my Makefile or in the Muparser folder. But I don't know how to fix. If you can give me some hint, I would be realy appreciate. 
