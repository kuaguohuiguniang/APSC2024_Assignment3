#include "Jacobi.hpp"
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <omp.h>
#include <fstream>

void Jacobi::solve_parallel() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int nrows = n;
    int local_ncols = n;
    int local_nrows;
    double local_error_square_sum = 0.0;

    // each rank deal with a balanced rows number
    // for example, n=11, size=4, rank 0 update 2 rows, rank 1 update 3 rows, rank 2 update 3 rows, rank 3 update 1 row
    if (rank == 0 || rank == size - 1) {
        local_nrows = (nrows%size > rank) ? nrows/size + 2 : nrows/size + 1;
    } else {
        local_nrows = (nrows%size > rank) ? nrows/size + 3 : nrows/size + 2;
    }

    // local_M has 1 or 2 rows more than the rows need to be updated, for example, local_M for rank 1 has 5 rows
    // the information from its neighbor rank are stored in these extra rows
    std::vector<double> local_M(local_ncols*local_nrows, 0);
    std::vector<double> local_M_new(local_ncols*local_nrows, 0);

    while (error > tolr && it < max_it) {
        // update local_M_new
#pragma omp parallel for shared(local_M, local_M_new), reduction(+:local_error_square_sum)
        for (int j = 1; j < local_ncols - 1; j++) {
            for (int i = 1; i < local_nrows - 1; i++) {
                double x = j * h;
                double y = (nrows % size > rank) ? ((n/size+1)*rank+i - 1*(rank > 0))*h :
                                                       1-((n/size)*(size-rank)-i)*h;
                int index = i * local_ncols + j;
                local_M_new[index] = 0.25 * (local_M[index-1]+local_M[index+1] + local_M[index-local_ncols]
                                             + local_M[index+local_ncols] + std::pow(h, 2)*f(x, y));
                // compute the local_error
                local_error_square_sum += h * std::pow(local_M_new[index] - local_M[index], 2);
            }
        }

        // pass local_M_new to local_M
        local_M = local_M_new;

        // send and receive the extra rows to its neighbor rank
        MPI_Request requests[4];
        int request_count = 0;

        if (rank == 0) {
            MPI_Isend(local_M.data() + (local_nrows - 2) * local_ncols, local_ncols, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(local_M.data() + (local_nrows - 1) * local_ncols, local_ncols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        } else if (rank == size - 1) {
            MPI_Isend(local_M.data() + local_ncols, local_ncols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(local_M.data(), local_ncols, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &requests[request_count++]);
        } else {
            MPI_Isend(local_M.data() + local_ncols, local_ncols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Isend(local_M.data() + (local_nrows - 2) * local_ncols, local_ncols, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(local_M.data(), local_ncols, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &requests[request_count++]);
            MPI_Irecv(local_M.data() + (local_nrows - 1) * local_ncols, local_ncols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &requests[request_count++]);
        }

        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
        
        // compute the totall error
        MPI_Allreduce(&local_error_square_sum, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        error = std::sqrt(error);
        it++;
        local_error_square_sum = 0.;
    }

    // gather the local_M to M
    int send_count;
    if (rank == 0 || rank == size - 1) {
        send_count = (local_nrows - 1) * local_ncols;
    } else {
        send_count = (local_nrows - 2) * local_ncols;
    }
    std::vector<int> recv_counts(size);
    MPI_Gather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> displs(size, 0);
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; ++i) {
            displs[i] = displs[i-1] + recv_counts[i-1];
        }
    }

    MPI_Gatherv(local_M.data() + (rank == 0 ? 0 : local_ncols),
                send_count, MPI_DOUBLE,
                M.data(),
                recv_counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

double Jacobi::test(std::function<double(double&, double&)> func) {
    double test_error = 0.;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            double x = j*h;
            double y = i*h;
            test_error += h*std::pow(M[i*n+j] - func(x,y), 2);
        }
    }
    test_error = std::sqrt(test_error);
    return test_error;
}

void Jacobi::print_res(){
    for (int i=0; i < n; i++){
        for (int j=0; j < n; j++){
            int index = i*n + j;
            std::cout<<M[index]<<" ";
        }
        std::cout<<";"<<std::endl;
    }
}

// this function is almost the same with the one in lab11
void Jacobi::output_VTK(std::string out_name){
    out_M.clear();
    for (int i = 0; i < n; i++){
        std::vector<double> row;
        for (int j = 0; j < n; j++){
            row.push_back(M[i*n+j]);
        }
        out_M.push_back(row);
    }
    //generateVTKFile(out_name, out_M, n, n, h, h);
    // opens the file
    std::ofstream vtkFile(out_name);

    // check if the file was opened
    if (!vtkFile.is_open()) {
        std::cerr << "Error: could not open file " << out_name << std::endl;
        return;
    }

    // Write VTK header
    vtkFile <<  "# vtk DataFile Version 3.0\n";
    vtkFile << "Scalar Field Data\n";
    vtkFile << "ASCII\n";                                // file format
    

    // Write grid data
    vtkFile << "DATASET STRUCTURED_POINTS\n";                             // format of the dataset
    vtkFile << "DIMENSIONS " << n << " " << n << " " << 1 << "\n";  // number of points in each direction
    vtkFile << "ORIGIN 0 0 0\n";                                          // lower-left corner of the structured grid
    vtkFile << "SPACING" << " " << h << " " << h << " " << 1 << "\n";   // spacing between points in each direction
    vtkFile << "POINT_DATA " << (n) * (n) << "\n";                  // number of points
                                                                
    
    // Write scalar field data
    vtkFile << "SCALARS scalars double\n";               // description of the scalar field
    vtkFile << "LOOKUP_TABLE default\n";                 // color table

    // Write vector field data
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            vtkFile <<  out_M[i][j] << "\n";
        }
    }
}
