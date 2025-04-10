//------------------------------------- Code for m-sized System of 1D Wave Equations using Finite Volume Method ------------------------------------
//------------------------------------------------------------- Scheme used - Upwind ---------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <iomanip>
#include <vector>

using namespace std;
using namespace Eigen;

/*  A guide to reading this code:
    -> The entire code has been written in a modular way, with very few operations done in the main() portion of the code.
    -> All the operations are done using the group of user-defined functions defined in the successive sections.
    -> The best way to understand this code is to first go the main() section. Constant values are defined in the initial part.
    -> Thereafter, as you come across a user defined function in main(), refer to the segment where the function is defined.
    -> The description and working of the user defined function is mentioned inside the respective function.
    -> This code makes extensive usage of the Eigen Library for C++. It has quite a few different datatypes like VectorXd and MatrixXd
       which are different from the conventional vector and array datatypes. */

// The following function evaluates maximum of two values
double find_max(double var1, double var2)
{
    if((var1>var2)||(var1==var2))
    {
        return var1;
    }
    else
    {
        return var2;
    }
}

// The following function evaluates minimum of two values
double find_min(double var1, double var2)
{
    if((var1<var2)||(var1==var2))
    {
        return var1;
    }
    else
    {
        return var2;
    }
}

// The following function is used for finding the size of the A matrix for the m x m system of PDEs
int count_matrix_size(const std::string& filename)
{
    /*  There is a special instruction for storing A.
        Firstly, store it in a .txt file
        Secondly, at the end of the last row, there should not be \n
        The size is determined as number of \n + 1
        Hence, if there is a \n at the end of A, it will count 1 extra  */

    std::ifstream file(filename);
    char ch;
    int count_var = 0;
    while(file.get(ch))
    {
        if (ch == '\n')
        {
            count_var = count_var + 1;
        }
    }
    file.close();
    return count_var;
}

// The following function stores the A matrix from the A.txt file and stores it into a vector
void matrix_A_storing(const std::string& filename, int matrix_size, vector<vector<double>>& A)
{
    std::ifstream file(filename);
    A.resize(matrix_size, std::vector<double>(matrix_size));
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            file >> A[i][j];
        }
    }
    file.close();
}

// The following function stores the single valued data inputs like L, N, t from their respective files and storing them into respective variables
void read_single_data(const std::string& filename, double& dataname)
{
    std::ifstream file(filename);
    file >> dataname;
    file.close();
}

// The following function stores the UL and UR values by reading them and storing them into VectorXd datatype - this is special to the Eigen library
void read_UL_UR(const std::string& filename, int matrix_size, VectorXd& UL_or_UR)
{
    std::ifstream file(filename);
    for(int i=0;i<matrix_size;i++)
    {
        file >> UL_or_UR(i);
    }
    file.close();
}

// The following function stores the A matrix from a vector datatype into a MatrixXd datatype which is needed when using the Eigen library
void create_vectorized_A(int matrix_size, vector<vector<double>>& A, MatrixXd& A_vectorized)
{
    for(int i=0;i<matrix_size;i++)
    {
        for(int j=0;j<matrix_size;j++)
        {
            A_vectorized(i,j) = A[i][j];
        }
    }
}

/* The following function sorts the eigenvalues, and arranges the corresponding eigenvectors also.
Sorting is necessary, because when the eigenvalues and eigenvectors are evaluated using the library functions, they are not in an ascending order
The eigenvalues and hence the eigenvectors are stored in a random order.
For the K matrix it is necessary that the first column should contain the eigenvector for the least eigenvalue, and so on thereafter.
Hence when the eigenvalues are swapped while comparison during sorting, the eigenvectors are also sorted.
The easiest sorting algorithm is Bubble Sort and although computationally inefficient, it is easy to implement and hence was done in this code. */

void bubble_sort(int matrix_size, VectorXd& eigenvalues, MatrixXd& eigenvectors)
{
    for(int i=0;i<matrix_size-1;i++)
    {
        for(int j=0;j<matrix_size-i-1;j++)
        {
            if (eigenvalues(j) > eigenvalues(j+1))   // This compares two consecutive eigenvalues
            {
                double temp_eigen_val = eigenvalues(j+1); // Temporary variable to store one eigenvalue before swapping
                eigenvalues(j+1) = eigenvalues(j);
                eigenvalues(j) = temp_eigen_val;

                double temp_eigen_vec_array[matrix_size]; // temporary 1D array to store the corresponding eigenvector of the eigenvalue stored in temp
                for(int m=0;m<matrix_size;m++)
                {
                    temp_eigen_vec_array[m] = eigenvectors(m,j+1);
                }
                for(int m=0;m<matrix_size;m++)
                {
                    eigenvectors(m,j+1) = eigenvectors(m,j);
                }
                for(int m=0;m<matrix_size;m++)
                {
                    eigenvectors(m,j) = temp_eigen_vec_array[m];
                }
            }
        }
    }
}

// The following function finds the maximum of the absolute value of the eigenvalues
double max_mod_eigenvalue(int matrix_size, VectorXd& eigenvalues)
{
    double max_eigenvalue = fabs(eigenvalues(0));
    for(int i=1;i<matrix_size;i++)
    {
        if (fabs(eigenvalues(i)) > max_eigenvalue)
        {
            max_eigenvalue = fabs(eigenvalues(i));
        }
    }
    return max_eigenvalue;
}

// The following function evaluates the delta_t value based on the CFL criteria
double evaluate_delta_t(double CFL, double delta_x, double max_eigenvalue)
{
    return CFL*delta_x/max_eigenvalue;
}

// The following function initializes the solution in W-space with WL and WR values
void initialize_domain(int matrix_size, int N, double L, double delta_x, VectorXd& WL, VectorXd& WR, MatrixXd& W_total)
{
    for(int i=0;i<matrix_size;i++)
    {
        double x_counter = -L/2 + delta_x/2; // Counter that traverses from one cell centre to another. Initialized with Centre of first cell
        for(int j=0;j<N;j++)
        {
            if(x_counter < 0)
            {
                W_total(i,j) = WL(i); // As long as cell centre is located left of origin, it is initialized with WL
            }
            else
            {
                W_total(i,j) = WR(i); // As long as cell centre is located right of origin, it is initialized with WR
            }
            x_counter = x_counter + delta_x; // Updating counter to move to next cell centre
        }
    }
}

// The following function initializes the solution in W-space with WL and WR values
void initialize_U(int matrix_size, int N, double L, double delta_x, VectorXd& UL, VectorXd& UR, MatrixXd& U_initial)
{
    for(int i=0;i<matrix_size;i++)
    {
        double x_counter = -L/2 + delta_x/2; // Counter that traverses from one cell centre to another. Initialized with Centre of first cell
        for(int j=0;j<N;j++)
        {
            if(x_counter < 0)
            {
                U_initial(i,j) = UL(i); // As long as cell centre is located left of origin, it is initialized with WL
            }
            else
            {
                U_initial(i,j) = UR(i); // As long as cell centre is located right of origin, it is initialized with WR
            }
            x_counter = x_counter + delta_x; // Updating counter to move to next cell centre
        }
    }
}

// The following function executes the Upwind scheme
void Upwind(int matrix_size, int N, double delta_x, double delta_t, int integral_number_of_time_iterations, double extra_fractional_time_to_cover, VectorXd& eigenvalues, MatrixXd& W_total, MatrixXd& W_total_new)
{
    for(int i=0;i<matrix_size;i++)
    {
        double wave_speed = eigenvalues(i); // current wave speed (lambda) is current eigenvalue
        double wave_speed_plus = find_max(wave_speed,0.0); // lambda+ evaluation
        double wave_speed_minus = find_min(wave_speed,0.0); // lambda- evaluation

        // Time marching this solution for Integral number of time iterations
        for(int t=0;t<integral_number_of_time_iterations;t++)
        {
            for(int j=1;j<N-1;j++)
            {
                W_total_new(i,j) = W_total(i,j) - (wave_speed_minus*delta_t/delta_x)*(W_total(i,j+1) - W_total(i,j)) - (wave_speed_plus*delta_t/delta_x)*(W_total(i,j) - W_total(i,j-1));
            }
            for(int j=1;j<N-1;j++)

            {
                W_total(i,j) = W_total_new(i,j); // Copying to W_total for next timestep
            }
        }

        /* This part evaluates the extra time marching needed for the time, which is less than a timestep value, and could not be completed will the last
        integral time iteration*/
        for(int j=1;j<N-1;j++)
        {
            W_total_new(i,j) = W_total(i,j) - (wave_speed_minus*extra_fractional_time_to_cover/delta_x)*(W_total(i,j+1) - W_total(i,j)) - (wave_speed_plus*extra_fractional_time_to_cover/delta_x)*(W_total(i,j) - W_total(i,j-1));
        }
    }
}

// The following function stores the results in a CSV file
void put_in_csv(const std::string& filename, int matrix_size, int N, MatrixXd& U_final)
{
    std::ofstream file(filename);
    for(int i=0;i<matrix_size;i++)
    {
        for(int j=0;j<N;j++)
        {
            file << U_final(i,j);
            if(j<N-1)
            {
                file << ",";
            }
        }
        if (i<matrix_size-1)
        {
            file << "\n";
        }
    }
    file << std::endl;
    file.close();
}

int main()
{
    // Reading Matrix A from A.txt and storing it into a vector
    std::string input_A = "input_A.txt";
    int matrix_size = count_matrix_size(input_A) + 1;
    vector<vector<double>> A;
    matrix_A_storing(input_A, matrix_size, A);

    // Creating a MatrixXd datatype and storing A into it so that Eigen library functions can be used
    MatrixXd A_vectorized(matrix_size, matrix_size);
    create_vectorized_A(matrix_size,A,A_vectorized);

    // Reading inputs from corresponding input files
    std::string input_L = "input_L.txt";
    std::string input_t = "input_t.txt";
    std::string input_N = "input_N.txt";
    double L,total_time,N_mid;
    read_single_data(input_L, L);
    read_single_data(input_t, total_time);
    read_single_data(input_N, N_mid);
    int N = (int)N_mid;

    // Reading UL and UR from the respective input files
    std::string UL_file_p = "input_UL.txt";
    std::string UR_file_p = "input_UR.txt";
    VectorXd UL(matrix_size), UR(matrix_size);
    read_UL_UR(UL_file_p,matrix_size,UL);
    read_UL_UR(UR_file_p,matrix_size,UR);

    // The following evaluates the eigenvalues and eigenvectors
    EigenSolver<MatrixXd> solver(A_vectorized, true);
    VectorXd eigenvalues = solver.eigenvalues().real();         // This is a type of a 1D Vector
    MatrixXd eigenvectors = solver.eigenvectors().real();       // This is a type of a 2D Vector

    /*  Sorting is necessary, because when the eigenvalues and eigenvectors are evaluated using the library functions, they are not in an ascending order
        The eigenvalues and hence the eigenvectors are stored in a random order.
        For the K matrix it is necessary that the first column should contain the eigenvector for the least eigenvalue, and so on thereafter.
        The easiest sorting algorithm is Bubble Sort and although computationally inefficient, it is easy to implement and hence was done in this code. */

    bubble_sort(matrix_size, eigenvalues, eigenvectors);
    // After sorting the eigenvectors matrix is the same as the K matrix

    MatrixXd K_inverse = eigenvectors.inverse();        // Computing K inverse matrix
    VectorXd WL(matrix_size), WR(matrix_size);

    // Transforming Initial Riemann Conditions from U space to W Space
    WL = K_inverse*UL;
    WR = K_inverse*UR;

    double delta_x = L/N;                               // Grid size
    double CFL = 0.8;                                   // CFL number
    double endtime = 0.2;                               // Simulation endtime

    MatrixXd W_total(matrix_size,N);                                // 2D Matrix to store the W1, W2, ...
    MatrixXd W_total_new(matrix_size,N);                            // Copy matrix

    initialize_domain(matrix_size,N,L,delta_x,WL,WR,W_total);       // Initializing W domain for W1, W2, ...
    W_total_new = W_total;                                          // Initializing copy of W domain for W1, W2, ...

    double max_eigenvalue = max_mod_eigenvalue(matrix_size, eigenvalues);       // Evaluating maximum modulus of eigenvalues
    double delta_t = evaluate_delta_t(CFL,delta_x,max_eigenvalue);              // Evaluating timestep

    double total_number_of_time_iterations = endtime/delta_t;
    int integral_number_of_time_iterations = (int)total_number_of_time_iterations;  // Time intervals that are with complete timesteps


    // Time value remaining that was not covered till the last complete timestep and is lesser than a timestep
    double extra_fractional_time_to_cover = endtime - delta_t*integral_number_of_time_iterations;

    Upwind(matrix_size,N,delta_x,delta_t,integral_number_of_time_iterations,extra_fractional_time_to_cover,eigenvalues,W_total,W_total_new);

    MatrixXd U_final(matrix_size,N);        // Matrix to store final U values
    U_final = eigenvectors*W_total_new;     // U = K*W

    //If there are no errors
    cout << "Execution Successful\n";

    // Initial state in U space
    MatrixXd U_initial(matrix_size,N);
    initialize_U(matrix_size,N,L,delta_x,UL,UR,U_initial);
    put_in_csv("u_initial.csv",matrix_size,N,U_initial);

    // Storing final results in CSV file which are to be read and plotted in MATLAB
    put_in_csv("output.csv",matrix_size,N,U_final);

    return 0;
}

