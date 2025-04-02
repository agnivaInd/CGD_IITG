//----------------------------------------- Code for 1D Wave Equation using Finite Volume Method -----------------------------------------
//--------------------------------- Schemes implemented - Upwind, FTFS, FTBS, Lax-Friedrichs, Lax-Wendroff -------------------------------

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

/*  A guide to reading this code:
    -> The entire code has been written in a modular way, with no operations done in the main() portion of the code.
    -> All the operations are done using the group of user-defined functions defined in the successive sections.
    -> The best way to understand this code is to first go the main() section. Constant values are defined in the initial part.
    -> Thereafter, as you come across a user defined function in main(), refer to the segment where the function is defined.
    -> The description and working of the user defined function is mentioned inside the respective function. */

void initialize_U(int total_points, double delta_x, double sigma, double U[])
{
    double x_start = -1;                    // Starting of the domain i.e. x = -1
    double xi = x_start + delta_x/2;        // xi is the cell centre of the first cell
    double e = exp(1.0);

    /* The next few lines of code evaluate the averaged u value across each cell.
       A brief explanation about the following lines:
       -> The integral is done numerically, using the trapezoidal rule, and each individual cell,
          has been divided into 100 parts to carry out this integral.
       -> Outer for loop moves the entire domain for the current cell to the next cell.
       -> Inner for loop is for evaluating trapezoidal rule integral, which is assigned to the cell centre at end of its completion.
       Please contact me if there is any further confusion related to this part of the code. */

    for(int i=0;i<total_points;i++)
    {
        double cell_lower_limit = xi - delta_x/2;           // left edge for the current cell which is located at cell centre minus grid size/2
        double cell_upper_limit = xi + delta_x/2;           // right edge for the current cell which is located at cell centre plus grid size/2

        double intra_cell_spacing = (cell_upper_limit - cell_lower_limit)/100;      // Dividing current cell into 100 divisions for finding averaged u
        double intra_x = cell_lower_limit;          // Initialization of counter that helps to traverse through inside the current cell
        double sum = 0;

        for(int j=0;j<100;j++)                                          //
        {                                                               //
            if((j==0)||(j==99))                                         //
            {                                                           //
                sum = sum + pow(e,-pow((intra_x/(2*sigma)),2));         //  This block is for the trapezoidal rule integral for inside each cell
            }                                                           //
            else                                                        //
            {                                                           //
                sum = sum + 2*pow(e,-pow((intra_x/(2*sigma)),2));       //
            }
            intra_x = intra_x + intra_cell_spacing;     // Counter that helps to traverse through inside the current cell
        }

        U[i] = 0.5*sum*intra_cell_spacing/(cell_upper_limit - cell_lower_limit); // Assigning averaged cell value to each cell centre
        xi = xi + delta_x;      // Moves onto next cell centre, and the process of getting averaged u for that cell is done similarly
    }
}

double L2_norm(int total_points, double sigma, double delta_x, double U[])
{
    double x_start = -1;
    double xi = x_start + delta_x/2;
    double e = exp(1.0);
    double sum = 0;
    double analytical;

    for(int i=0;i<total_points;i++)
    {
        analytical = pow(e,-pow((xi/(2*sigma)),2));     // Evaluating the analytical u value at the current cell centre
        xi = xi + delta_x;                              // Moving to next cell centre
        sum = sum + pow((analytical - U[i]),2);         // Summing up the square of the errors at each grid point
    }

    return sqrt(sum/total_points);      // Returns the L2 norm value for the particular case
}

void put_in_csv(const std::string& filename, int total_points, double U[])
{
    std::ofstream file(filename);
    for (int i=0;i<total_points;i++)
    {
        file << U[i];
        if (i < total_points-1)
        {
            file << ",";
        }
    }
    file << std::endl;
    file.close();
}

void FTBS(int total_points, double wave_speed, double delta_t, double delta_x, double total_number_of_timesteps, double U[], double U_new[])
{
    for(int t=0;t<total_number_of_timesteps;t++)
    {
        for(int j=0;j<total_points;j++)
        {
            if(j==0)
            {
                U_new[j] = U[j] - (wave_speed*delta_t/delta_x)*(U[j] - U[total_points-1]);
                // For j-1 cell at j = 0 boundary value, last cell taken as j-1 cell (Periodic Boundary condition)
            }
            else
            {
                U_new[j] = U[j] - (wave_speed*delta_t/delta_x)*(U[j] - U[j-1]);
            }
        }

        for(int j=0;j<total_points;j++)
        {
            U[j] = U_new[j];        // Copying new U values into old array before next time iteration
        }
    }
}

void FTFS(int total_points, double wave_speed, double delta_t, double delta_x, double total_number_of_timesteps, double U[], double U_new[])
{
    for(int t=0;t<total_number_of_timesteps;t++)
    {
        for(int j=0;j<total_points;j++)
        {
            if(j==total_points-1)
            {
                U_new[j] = U[j] - (wave_speed*delta_t/delta_x)*(U[0] - U[j]);
                // For j+1 cell at j = total_points-1 boundary value, first cell taken as j+1 cell (Periodic Boundary condition)
            }
            else
            {
                U_new[j] = U[j] - (wave_speed*delta_t/delta_x)*(U[j+1] - U[j]);
            }
        }
        for(int j=0;j<total_points;j++)
        {
            U[j] = U_new[j];        // Copying new U values into old array before next time iteration
        }
    }
}

void Lax_Friedrichs(int total_points, double wave_speed, double delta_t, double delta_x, double total_number_of_timesteps, double U[], double U_new[])
{
    for(int t=0;t<total_number_of_timesteps;t++)
    {
        for(int j=0;j<total_points;j++)
        {
            // Periodic boundary condition is implemented for inaccessible cells as in FTBS and FTFS
            if(j==0)
            {
                U_new[j] = 0.5*(U[j+1] + U[total_points-1]) - 0.5*(wave_speed*delta_t/delta_x)*(U[j+1] - U[total_points-1]);
            }
            else if (j==total_points-1)
            {
                U_new[j] = 0.5*(U[0] + U[j-1]) - 0.5*(wave_speed*delta_t/delta_x)*(U[0] - U[j-1]);
            }
            else
            {
                U_new[j] = 0.5*(U[j+1] + U[j-1]) - 0.5*(wave_speed*delta_t/delta_x)*(U[j+1] - U[j-1]);
            }
        }
        for(int j=0;j<total_points;j++)
        {
            U[j] = U_new[j];        // Copying new U values into old array before next time iteration
        }
    }
}

void Lax_Wendroff(int total_points, double wave_speed, double delta_t, double delta_x, double total_number_of_timesteps, double U[], double U_new[])
{
    for(int t=0;t<total_number_of_timesteps;t++)
    {
        for(int j=0;j<total_points;j++)
        {
            // Periodic boundary condition is implemented for inaccessible cells as in FTBS and FTFS
            if(j==0)
            {
                U_new[j] = U[j] - 0.5*(wave_speed*delta_t/delta_x)*(U[j+1] - U[total_points-1]) + 0.5*pow((wave_speed*delta_t/delta_x),2)*(U[j+1] - 2*U[j] + U[total_points-1]);
            }
            else if (j==total_points-1)
            {
                U_new[j] = U[j] - 0.5*(wave_speed*delta_t/delta_x)*(U[0] - U[j-1]) + 0.5*pow((wave_speed*delta_t/delta_x),2)*(U[0] - 2*U[j] + U[j-1]);
            }
            else
            {
                U_new[j] = U[j] - 0.5*(wave_speed*delta_t/delta_x)*(U[j+1] - U[j-1]) + 0.5*pow((wave_speed*delta_t/delta_x),2)*(U[j+1] - 2*U[j] + U[j-1]);
            }
        }
        for(int j=0;j<total_points;j++)
        {
            U[j] = U_new[j];        // Copying new U values into old array before next time iteration
        }
    }
}

void Upwind(int total_points, double wave_speed, double delta_t, double delta_x, double total_number_of_timesteps, double U[], double U_new[])
{
    double wave_speed_plus;                     // positive wave speed i.e. a+
    double wave_speed_minus;                    // negative wave speed i.e. a-

    /* Next if-else condition is for assigning a+ and a- values
       a- = minimum(a,0) and a+ = maximum(a,0) */

    if ((wave_speed > 0)||(wave_speed = 0))
    {
        wave_speed_minus = 0;
        wave_speed_plus = wave_speed;
    }
    else
    {
        wave_speed_minus = wave_speed;
        wave_speed_plus = 0;
    }

    for(int t=0;t<total_number_of_timesteps;t++)
    {
        for(int j=0;j<total_points;j++)
        {
            // Periodic boundary condition is implemented for inaccessible cells as in FTBS and FTFS
            if(j==0)
            {
                U_new[j] = U[j] - (wave_speed_minus*delta_t/delta_x)*(U[j+1] - U[j]) - (wave_speed_plus*delta_t/delta_x)*(U[j] - U[total_points-1]);
            }
            else if (j==total_points-1)
            {
                U_new[j] = U[j] - (wave_speed_minus*delta_t/delta_x)*(U[0] - U[j]) - (wave_speed_plus*delta_t/delta_x)*(U[j] - U[j-1]);
            }
            else
            {
                U_new[j] = U[j] - (wave_speed_minus*delta_t/delta_x)*(U[j+1] - U[j]) - (wave_speed_plus*delta_t/delta_x)*(U[j] - U[j-1]);
            }
        }

        for(int j=0;j<total_points;j++)
        {
            U[j] = U_new[j];        // Copying new U values into old array before next time iteration
        }
    }
}

int main()
{
    // Assignment of variable values
    double wave_speed = 2;                              // Wave Speed i.e. "a"
    double domain_length = 2;                           // Length of domain from x = -1 to x = 1
    double delta_x = 0.01;                              // Grid size i.e. delta x
    int total_points = domain_length/delta_x;           // Total no of cells (Also, a count of the number of the cell centre points)
    double delta_t = 0.0025;                            // Timestep value
    double endtime = 1;                                 // Final time upto which marching is done
    int total_number_of_timesteps = endtime/delta_t;    // Intuitive
    double sigma = 1.0/8.0;                             // Given data

    // 1D arrays to store current and previous timestep values
    double U[total_points]; // For storing previous timestep values
    double U_new[total_points]; // For storing current timestep value after each iteration

    // Initializing the value of U by assigning averaged "u" value across the cell and assigning it to the cell centre.
    initialize_U(total_points,delta_x,sigma,U);

    /* Before running the code, chose the desired scheme by not commenting out the same.
       Here Upwind is not commented out, hence is executed in the code. */

    // FTBS(total_points,wave_speed,delta_t,delta_x,total_number_of_timesteps,U,U_new);
    // FTFS(total_points,wave_speed,delta_t,delta_x,total_number_of_timesteps,U,U_new);
    // Lax_Friedrichs(total_points,wave_speed,delta_t,delta_x,total_number_of_timesteps,U,U_new);
    // Lax_Wendroff(total_points,wave_speed,delta_t,delta_x,total_number_of_timesteps,U,U_new);
    Upwind(total_points,wave_speed,delta_t,delta_x,total_number_of_timesteps,U,U_new);

    // Displaying the result for immediate inference
    for(int i=0;i<total_points;i++)
    {
        cout << setprecision(6) << U[i] << " ";
    }

    // Displaying the L2 norm
    cout << "\n"
         << "L2 norm of this scheme is: " << L2_norm(total_points,sigma,delta_x,U);

    // Storing the final results in a csv file, which is then accessed through MATLAB for further post processing
    put_in_csv("output.csv",total_points,U);
}
