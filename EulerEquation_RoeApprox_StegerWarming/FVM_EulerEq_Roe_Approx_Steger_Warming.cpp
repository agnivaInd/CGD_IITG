//----------------------------------------------------------- Code for 1D Sod Shock Tube Problem ----------------------------------------------------------
//--------------------- Schemes implemented - Roe Approximate Solver WITH and WITHOUT Entropy Fix, Steger Warming Flux Splitting Scheme -------------------

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

// This function returns lambda+ values
double return_eigenvalue_plus(double eigenvalue)
{
    return 0.5*(eigenvalue + fabs(eigenvalue));
}

// This function returns lambda- values
double return_eigenvalue_minus(double eigenvalue)
{
    return 0.5*(eigenvalue - fabs(eigenvalue));
}

/* This function assigns the initial values (conserved quantities - rho, rho_u and rho_e to the domain
   The first row stores rho, The second row stores rho_u, The third row stores rho_e */
void initialize_domain(int gridpts, double x_start, double x_end, double delta_x, double domain_length, MatrixXd& Domain, double rho_l, double rho_r, double p_l, double p_r, double u_l, double u_r, double gamma)
{
    double rho_u_l = rho_l*u_l;
    double rho_u_r = rho_r*u_r;

    double rho_e_l = 0.5*rho_l*pow(u_l,2) + p_l/(gamma - 1);
    double rho_e_r = 0.5*rho_r*pow(u_r,2) + p_r/(gamma - 1);

    double x_counter = x_start + delta_x/2; // Initializing counter that travels cell centres to the cell centre of the leftmost cell

    for(int j=0;j<gridpts;j++)
    {
        if (x_counter<(domain_length/2))
        {
            Domain(0,j) = rho_l;
            Domain(1,j) = rho_u_l;
            Domain(2,j) = rho_e_l;
        }
        else
        {
            Domain(0,j) = rho_r;
            Domain(1,j) = rho_u_r;
            Domain(2,j) = rho_e_r;
        }
        x_counter = x_counter + delta_x;    // Updating cell centre traversing counter
    }
}

// The following user defined function is for implementing the Roe Approximate Solver without the entropy fix
void Roe_Solver_WITHOUT_entropy_fix(int gridpts, MatrixXd& Domain, MatrixXd& Copy_Domain, double initial_time, double final_time, double CFL, double delta_x, double gamma)
{
    double t = initial_time;        // Initializing current time with final time
    while(t < final_time)
    {
        VectorXd P(gridpts);            // This 1D vector stores Pressure values for the all the gridpoints
        for(int j=0;j<gridpts;j++)
        {
            P(j) = (Copy_Domain(2,j) - pow(Copy_Domain(1,j),2)/(2*Copy_Domain(0,j)))*(gamma - 1);
        }

        // The following 2D Vector stores the tilda values for each pair of continuous cells, and there are gridpts-1 of such boundaries where tilda values are to be evaluated
        MatrixXd Store_tilda_values(3,gridpts-1);
        for(int j=0;j<gridpts-1;j++)
        {
            double local_ul = Copy_Domain(1,j)/Copy_Domain(0,j);
            double local_ur = Copy_Domain(1,j+1)/Copy_Domain(0,j+1);

            double local_Hl = (Copy_Domain(2,j) + P(j))/Copy_Domain(0,j);
            double local_Hr = (Copy_Domain(2,j+1) + P(j+1))/Copy_Domain(0,j+1);

            double local_root_rho_l = pow(Copy_Domain(0,j),0.5);
            double local_root_rho_r = pow(Copy_Domain(0,j+1),0.5);

            Store_tilda_values(0,j) = (local_root_rho_r*local_ur + local_root_rho_l*local_ul)/(local_root_rho_l + local_root_rho_r);    // u tilda
            Store_tilda_values(1,j) = (local_root_rho_r*local_Hr + local_root_rho_l*local_Hl)/(local_root_rho_l + local_root_rho_r);    // H tilda
            Store_tilda_values(2,j) = pow((gamma - 1)*(Store_tilda_values(1,j) - pow(Store_tilda_values(0,j),2)/2),0.5);                // a tilda
        }

        // The following vector evaluates and stores the eigenvalues using the tilda values evaluated previously
        MatrixXd Store_eigenvalues(3,gridpts-1);
        for(int j=0;j<gridpts-1;j++)
        {
            Store_eigenvalues(0,j) = Store_tilda_values(0,j) - Store_tilda_values(2,j);
            Store_eigenvalues(1,j) = Store_tilda_values(0,j);
            Store_eigenvalues(2,j) = Store_tilda_values(0,j) + Store_tilda_values(2,j);
        }

        double max_mod_eigenvalue = fabs(Store_eigenvalues(0,0));
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<gridpts-1;j++)
            {
                if (fabs(Store_eigenvalues(i,j)) > max_mod_eigenvalue)
                {
                    max_mod_eigenvalue = fabs(Store_eigenvalues(i,j));
                }
            }
        }
        // Maximum of the modulus of the eigenvalues is evaluated in the previous for loops

        double delta_t = CFL*delta_x/max_mod_eigenvalue;       // timestep or delta_t is evaluated using the CFL criteria
        if((t+delta_t)>final_time)
        {
            delta_t = final_time - t;
            // This is done to ensure that in the final timestep the total timestep does not exceed the final time
        }

        MatrixXd Flux(3,gridpts-1);     // This is a 2D vector that stores the Flux value at the cell interfaces
        for(int j=0;j<gridpts-1;j++)
        {
            double del_u1 = Copy_Domain(0,j+1) - Copy_Domain(0,j);  // Difference between the first row entities for the right and left cell at each cell boundary
            double del_u2 = Copy_Domain(1,j+1) - Copy_Domain(1,j);  // Difference between the second row entities for the right and left cell at each cell boundary
            double del_u3 = Copy_Domain(2,j+1) - Copy_Domain(2,j);  // Difference between the third row entities for the right and left cell at each cell boundary

            double delta_tilda_values[3];
            // This 1D array stores the delta_tilda_values as given in the literature for Roe Approximate Solver
            delta_tilda_values[1] = (gamma-1)*(del_u1*(Store_tilda_values(1,j) - pow(Store_tilda_values(0,j),2)) + Store_tilda_values(0,j)*del_u2 - del_u3)/pow(Store_tilda_values(2,j),2);
            delta_tilda_values[0] = (del_u1*(Store_tilda_values(0,j)+Store_tilda_values(2,j)) - del_u2 - Store_tilda_values(2,j)*delta_tilda_values[1])/(2*Store_tilda_values(2,j));
            delta_tilda_values[2] = del_u1 - delta_tilda_values[0] - delta_tilda_values[1];

            VectorXd FL(3); // Between the current pair of cells, stores the flux value of the left cell
            FL(0) = Copy_Domain(1,j);
            FL(1) = pow(Copy_Domain(1,j),2)/Copy_Domain(0,j) + P(j);
            FL(2) = (Copy_Domain(1,j)/Copy_Domain(0,j))*(Copy_Domain(2,j)+P(j));

            MatrixXd Eigenvectors(3,3); // This 2D vector stores the eigenvalues for the current pair of cells
            Eigenvectors(0,0) = 1; Eigenvectors(0,1) = 1; Eigenvectors(0,2) = 1;
            Eigenvectors(1,0) = Store_tilda_values(0,j) - Store_tilda_values(2,j); Eigenvectors(1,1) = Store_tilda_values(0,j); Eigenvectors(1,2) = Store_tilda_values(0,j) + Store_tilda_values(2,j);
            Eigenvectors(2,0) = Store_tilda_values(1,j) - Store_tilda_values(0,j)*Store_tilda_values(2,j); Eigenvectors(2,1) = pow(Store_tilda_values(0,j),2)/2; Eigenvectors(2,2) = Store_tilda_values(1,j) + Store_tilda_values(0,j)*Store_tilda_values(2,j);

            Flux.col(j) = FL;   // Initializing the flux column of interface between two cells with the flux of left cell
            for(int i=0;i<3;i++)
            {
                if((Store_eigenvalues(i,j)<0)||(Store_eigenvalues(i,j)==0))     // If eigenvalues are less than or equal to zero
                {
                    Flux.col(j) = Flux.col(j) + Store_eigenvalues(i,j)*delta_tilda_values[i]*Eigenvectors.col(i);
                }
            }
        }

        for(int m=1;m<gridpts-1;m++)
        {
            // Updating the contents of the Domain, using the flux values obtained for each cell interface between two cells
            Domain.col(m) = Copy_Domain.col(m) - (delta_t/delta_x)*(Flux.col(m) - Flux.col(m-1));
        }

        Domain.col(0) = Domain.col(1);                      // Flux boundary condition at left boundary
        Domain.col(gridpts-1) = Domain.col(gridpts-2);      // Flux boundary condition at right boundary

        for(int j=0;j<gridpts;j++)
        {
            Copy_Domain.col(j) = Domain.col(j);     // Storing updated values to Copy_Domain for usage in next timestep
        }

        t = t + delta_t;    // Updating current time values
    }
}

// The following user defined function is for implementing the Roe Approximate Solver with the entropy fix
void Roe_Solver_WITH_entropy_fix(int gridpts, MatrixXd& Domain, MatrixXd& Copy_Domain, double initial_time, double final_time, double CFL, double delta_x, double gamma, double epsilon)
{
    double t = initial_time;        // Initializing current time with final time
    while(t < final_time)
    {
        VectorXd P(gridpts);            // This 1D vector stores Pressure values for the all the gridpoints
        for(int j=0;j<gridpts;j++)
        {
            P(j) = (Copy_Domain(2,j) - pow(Copy_Domain(1,j),2)/(2*Copy_Domain(0,j)))*(gamma - 1);
        }

        // The following 2D Vector stores the tilda values for each pair of continuous cells, and there are gridpts-1 of such boundaries where tilda values are to be evaluated
        MatrixXd Store_tilda_values(3,gridpts-1);
        for(int j=0;j<gridpts-1;j++)
        {
            double local_ul = Copy_Domain(1,j)/Copy_Domain(0,j);
            double local_ur = Copy_Domain(1,j+1)/Copy_Domain(0,j+1);

            double local_Hl = (Copy_Domain(2,j) + P(j))/Copy_Domain(0,j);
            double local_Hr = (Copy_Domain(2,j+1) + P(j+1))/Copy_Domain(0,j+1);

            double local_root_rho_l = pow(Copy_Domain(0,j),0.5);
            double local_root_rho_r = pow(Copy_Domain(0,j+1),0.5);

            Store_tilda_values(0,j) = (local_root_rho_r*local_ur + local_root_rho_l*local_ul)/(local_root_rho_l + local_root_rho_r);    // u tilda
            Store_tilda_values(1,j) = (local_root_rho_r*local_Hr + local_root_rho_l*local_Hl)/(local_root_rho_l + local_root_rho_r);    // H tilda
            Store_tilda_values(2,j) = pow((gamma - 1)*(Store_tilda_values(1,j) - pow(Store_tilda_values(0,j),2)/2),0.5);                // a tilda
        }

        // The following vector evaluates and stores the eigenvalues using the tilda values evaluated previously
        MatrixXd Store_eigenvalues(3,gridpts-1);
        for(int j=0;j<gridpts-1;j++)
        {
            Store_eigenvalues(0,j) = Store_tilda_values(0,j) - Store_tilda_values(2,j);
            Store_eigenvalues(1,j) = Store_tilda_values(0,j);
            Store_eigenvalues(2,j) = Store_tilda_values(0,j) + Store_tilda_values(2,j);
        }

        double max_mod_eigenvalue = fabs(Store_eigenvalues(0,0));
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<gridpts-1;j++)
            {
                if (fabs(Store_eigenvalues(i,j)) > max_mod_eigenvalue)
                {
                    max_mod_eigenvalue = fabs(Store_eigenvalues(i,j));
                }
            }
        }
        // Maximum of the modulus of the eigenvalues is evaluated in the previous for loops

        double delta_t = CFL*delta_x/max_mod_eigenvalue;        // timestep or delta_t is evaluated using the CFL criteria
        if((t+delta_t)>final_time)
        {
            delta_t = final_time - t;
            // This is done to ensure that in the final timestep the total timestep does not exceed the final time
        }

        MatrixXd Flux(3,gridpts-1);     // This is a 2D vector that stores the Flux value at the cell interfaces
        for(int j=0;j<gridpts-1;j++)
        {
            double del_u1 = Copy_Domain(0,j+1) - Copy_Domain(0,j);      // Difference between the first row entities for the right and left cell at each cell boundary
            double del_u2 = Copy_Domain(1,j+1) - Copy_Domain(1,j);      // Difference between the second row entities for the right and left cell at each cell boundary
            double del_u3 = Copy_Domain(2,j+1) - Copy_Domain(2,j);      // Difference between the third row entities for the right and left cell at each cell boundary

            double delta_tilda_values[3];
            // This 1D array stores the delta_tilda_values as given in the literature for Roe Approximate Solver
            delta_tilda_values[1] = (gamma-1)*(del_u1*(Store_tilda_values(1,j) - pow(Store_tilda_values(0,j),2)) + Store_tilda_values(0,j)*del_u2 - del_u3)/pow(Store_tilda_values(2,j),2);
            delta_tilda_values[0] = (del_u1*(Store_tilda_values(0,j)+Store_tilda_values(2,j)) - del_u2 - Store_tilda_values(2,j)*delta_tilda_values[1])/(2*Store_tilda_values(2,j));
            delta_tilda_values[2] = del_u1 - delta_tilda_values[0] - delta_tilda_values[1];

            VectorXd FL(3);     // Between the current pair of cells, stores the flux value of the left cell
            FL(0) = Copy_Domain(1,j);
            FL(1) = pow(Copy_Domain(1,j),2)/Copy_Domain(0,j) + P(j);
            FL(2) = (Copy_Domain(1,j)/Copy_Domain(0,j))*(Copy_Domain(2,j)+P(j));

            MatrixXd Eigenvectors(3,3);     // This 2D vector stores the eigenvalues for the current pair of cells
            Eigenvectors(0,0) = 1; Eigenvectors(0,1) = 1; Eigenvectors(0,2) = 1;
            Eigenvectors(1,0) = Store_tilda_values(0,j) - Store_tilda_values(2,j); Eigenvectors(1,1) = Store_tilda_values(0,j); Eigenvectors(1,2) = Store_tilda_values(0,j) + Store_tilda_values(2,j);
            Eigenvectors(2,0) = Store_tilda_values(1,j) - Store_tilda_values(0,j)*Store_tilda_values(2,j); Eigenvectors(2,1) = pow(Store_tilda_values(0,j),2)/2; Eigenvectors(2,2) = Store_tilda_values(1,j) + Store_tilda_values(0,j)*Store_tilda_values(2,j);

            Flux.col(j) = FL;       // Initializing the flux column of interface between two cells with the flux of left cell
            for(int i=0;i<3;i++)
            {
                if((Store_eigenvalues(i,j)<0)||(Store_eigenvalues(i,j)==0))     // If eigenvalues are less than or equal to zero
                {
                    if (fabs(Store_eigenvalues(i,j))<epsilon)
                    {
                        Store_eigenvalues(i,j) = 0.5*(pow(Store_eigenvalues(i,j),2)/epsilon + epsilon);     // Entropy fix is applied in this step
                    }
                    Flux.col(j) = Flux.col(j) + Store_eigenvalues(i,j)*delta_tilda_values[i]*Eigenvectors.col(i);
                }
            }
        }

        for(int m=1;m<gridpts-1;m++)
        {
            // Updating the contents of the Domain, using the flux values obtained for each cell interface between two cells
            Domain.col(m) = Copy_Domain.col(m) - (delta_t/delta_x)*(Flux.col(m) - Flux.col(m-1));
        }

        Domain.col(0) = Domain.col(1);                          // Flux boundary condition at left boundary
        Domain.col(gridpts-1) = Domain.col(gridpts-2);          // Flux boundary condition at right boundary

        for(int j=0;j<gridpts;j++)
        {
            Copy_Domain.col(j) = Domain.col(j);     // Storing updated values to Copy_Domain for usage in next timestep
        }

        t = t + delta_t;    // Updating current time value
    }
}

void Steger_Warming(int gridpts, MatrixXd& Domain, MatrixXd& Copy_Domain, double initial_time, double final_time, double CFL, double delta_x, double gamma)
{
    double t = initial_time;        // Initializing current time with final time
    while(t < final_time)
    {
        VectorXd P(gridpts);        // This 1D vector stores Pressure values for the all the gridpoints
        for(int j=0;j<gridpts;j++)
        {
            P(j) = (Copy_Domain(2,j) - pow(Copy_Domain(1,j),2)/(2*Copy_Domain(0,j)))*(gamma - 1);
        }

        MatrixXd Store_eigenvalues(3,gridpts);      // This 2D vector stores all the 3 eigenvalues of each cell
        VectorXd H(gridpts);                        // This 1D vector stores the H (specific enthalpy) values of each cell
        VectorXd a(gridpts);                        // This 1D vector stores the a (speed of sound) values of each cell
        for(int j=0;j<gridpts;j++)
        {
            Store_eigenvalues(0,j) = Copy_Domain(1,j)/Copy_Domain(0,j) - sqrt(gamma*P(j)/Copy_Domain(0,j));     // Row 1 stores u-a eigenvalue
            Store_eigenvalues(1,j) = Copy_Domain(1,j)/Copy_Domain(0,j);                                         // Row 2 stores u eigenvalue
            Store_eigenvalues(2,j) = Copy_Domain(1,j)/Copy_Domain(0,j) + sqrt(gamma*P(j)/Copy_Domain(0,j));     // Row 3 stores u+a eigenvalues
            H(j) = (Copy_Domain(2,j) + P(j))/Copy_Domain(0,j);          // Storing H value for the cell
            a(j) = sqrt(gamma*P(j)/Copy_Domain(0,j));                   // Storing a value for the cell
        }

        double max_mod_eigenvalue = fabs(Store_eigenvalues(0,0));
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<gridpts;j++)
            {
                if (fabs(Store_eigenvalues(i,j)) > max_mod_eigenvalue)
                {
                    max_mod_eigenvalue = fabs(Store_eigenvalues(i,j));
                }
            }
        }
        // Maximum of the modulus of the eigenvalues is evaluated in the previous for loops

        double delta_t = CFL*delta_x/max_mod_eigenvalue;        // timestep or delta_t is evaluated using the CFL criteria
        if((t+delta_t)>final_time)
        {
            delta_t = final_time - t;
            // This is done to ensure that in the final timestep the total timestep does not exceed the final time
        }

        MatrixXd Flux(3,gridpts+1); // This is a 2D vector that stores the Flux value at the cell interfaces
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<gridpts+1;j++)
            {
                Flux(i,j) = 0.0;    // Initializing all values of the Flux vector with 0
            }
        }

        for(int j=0;j<gridpts;j++)
        {
            // For each cell the lambda+ and lambda- values are evaluated for all the 3 eigenvalues associated with it
            VectorXd eigenvalues_plus(3), eigenvalues_minus(3);     // Two 2D vectors to store the lambda+ and lambda- values
            for(int i=0;i<3;i++)
            {
                eigenvalues_plus(i) = return_eigenvalue_plus(Store_eigenvalues(i,j));
                eigenvalues_minus(i) = return_eigenvalue_minus(Store_eigenvalues(i,j));
            }

            /* In this section, flux is split into two parts and is stored into 2 1D vectors,
               Flux_minus is the part of the cell flux that goes towards the left cell boundary
               Flux_plus is the part of the cell flux that goes towards the right cell boundary
               The formulas/expressions of corresponding components of the entries of the Flux_plus and Flux_minus are taken from the book of EL Toro */

            VectorXd Flux_plus(3), Flux_minus(3);
            Flux_plus(0) = (Copy_Domain(0,j)/(2*gamma))*(eigenvalues_plus(0) + 2*(gamma-1)*eigenvalues_plus(1) + eigenvalues_plus(2));
            Flux_minus(0) = (Copy_Domain(0,j)/(2*gamma))*(eigenvalues_minus(0) + 2*(gamma-1)*eigenvalues_minus(1) + eigenvalues_minus(2));

            Flux_plus(1) = (Copy_Domain(0,j)/(2*gamma))*((Copy_Domain(1,j)/Copy_Domain(0,j) - a(j))*eigenvalues_plus(0) + 2*(gamma-1)*(Copy_Domain(1,j)/Copy_Domain(0,j))*eigenvalues_plus(1) + (Copy_Domain(1,j)/Copy_Domain(0,j) + a(j))*eigenvalues_plus(2));
            Flux_minus(1) = (Copy_Domain(0,j)/(2*gamma))*((Copy_Domain(1,j)/Copy_Domain(0,j) - a(j))*eigenvalues_minus(0) + 2*(gamma-1)*(Copy_Domain(1,j)/Copy_Domain(0,j))*eigenvalues_minus(1) + (Copy_Domain(1,j)/Copy_Domain(0,j) + a(j))*eigenvalues_minus(2));

            Flux_plus(2) = (Copy_Domain(0,j)/(2*gamma))*((H(j) - (Copy_Domain(1,j)/Copy_Domain(0,j))*a(j))*eigenvalues_plus(0) + (gamma-1)*pow(Copy_Domain(1,j)/Copy_Domain(0,j),2)*eigenvalues_plus(1) + (H(j) + (Copy_Domain(1,j)/Copy_Domain(0,j))*a(j))*eigenvalues_plus(2));
            Flux_minus(2) = (Copy_Domain(0,j)/(2*gamma))*((H(j) - (Copy_Domain(1,j)/Copy_Domain(0,j))*a(j))*eigenvalues_minus(0) + (gamma-1)*pow(Copy_Domain(1,j)/Copy_Domain(0,j),2)*eigenvalues_minus(1) + (H(j) + (Copy_Domain(1,j)/Copy_Domain(0,j))*a(j))*eigenvalues_minus(2));

            Flux.col(j) = Flux.col(j) + Flux_minus;
            Flux.col(j+1) = Flux.col(j+1) + Flux_plus;
        }

        for(int m=1;m<gridpts-1;m++)
        {
            // Updating the contents of the Domain, using the flux values obtained for each cell interface between two cells
            Domain.col(m) = Copy_Domain.col(m) - (delta_t/delta_x)*(Flux.col(m+1) - Flux.col(m));
        }

        Domain.col(0) = Domain.col(1);                      // Flux boundary condition at left boundary
        Domain.col(gridpts-1) = Domain.col(gridpts-2);      // Right boundary condition at left boundary

        for(int j=0;j<gridpts;j++)
        {
            Copy_Domain.col(j) = Domain.col(j);         // Storing updated values to Copy_Domain for usage in next timestep
        }

        t = t + delta_t;        // Updating current time value
    }
}

// The following function reads input values from the input files storing UL and UR values
void read_rho_p_u(const std::string& filename_l, const std::string& filename_r, double& rho_l, double& rho_r, double& p_l, double& p_r, double& u_l, double& u_r)
{
    std::ifstream file_l(filename_l);
    for(int i=0;i<3;i++)
    {
        if(i==0)
        {
            file_l >> rho_l;
        }
        if(i==1)
        {
            file_l >> p_l;
        }
        if(i==2)
        {
            file_l >> u_l;
        }
    }
    file_l.close();

    std::ifstream file_r(filename_r);
    for(int i=0;i<3;i++)
    {
        if(i==0)
        {
            file_r >> rho_r;
        }
        if(i==1)
        {
            file_r >> p_r;
        }
        if(i==2)
        {
            file_r >> u_r;
        }
    }
    file_r.close();
}

// The following function stores the result, the primitive variable results into a CSV file
void put_in_csv(const std::string& filename, int gridpts, MatrixXd& Domain)
{
    std::ofstream file(filename);
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<gridpts;j++)
        {
            file << Domain(i,j);
            if(j<gridpts-1)
            {
                file << ",";
            }
        }
        if (i<4)
        {
            file << "\n";
        }
    }
    file << std::endl;
    file.close();
}

int main()
{
    int gridpts = 400;                              // No of grid points
    double x_start = 0;                             // Start of domain
    double x_end = 1;                               // End of domain
    double domain_length = x_end - x_start;
    double delta_x = (x_end - x_start)/gridpts;     // Grid size

    double gamma = 1.4;

    double initial_time = 0;
    double final_time = 0.15;
    double CFL = 0.8;

    /* Domain and Copy_Domain are two 2D vectors of MatrixXd datatype that store the values of the conserved variables
    u is stored in the first row, rho_u in the second row and rho_e in the third row*/
    MatrixXd Domain(3,gridpts);
    MatrixXd Copy_Domain(3,gridpts);

    double rho_l, p_l, u_l;
    double rho_r, p_r, u_r;

    std::string UL = "input_UL.txt";                    // Input text file for rho_l, p_l, u_l
    std::string UR = "input_UR.txt";                    // Input text file for rho_r, p_r, u_r
    read_rho_p_u(UL,UR,rho_l,rho_r,p_l,p_r,u_l,u_r);

    double epsilon = pow(10,-6);            // Epsilon value for entropy fix in Roe Approximate solver

    initialize_domain(gridpts,x_start,x_end,delta_x,domain_length,Domain,rho_l,rho_r,p_l,p_r,u_l,u_r,gamma);        // Initializing Domain
    initialize_domain(gridpts,x_start,x_end,delta_x,domain_length,Copy_Domain,rho_l,rho_r,p_l,p_r,u_l,u_r,gamma);   // Initializing Copy_Domain

    int choice;             // Choice to user to choose desired scheme
    cout << "Enter 1: For Roe Approx Solver WITHOUT Entropy fix\n"
         << "Enter 2: For Roe Approx Solver WITH Entropy fix\n"
         << "Enter 3: For Steger Warming Scheme\n";
    cout << "Choice: "; cin >> choice;

    MatrixXd Final_results(4,gridpts);      // 2D array that stores the primitive variable values, after extracting them from conserved values - from Domain 2D Matrix

    switch(choice)
    {
        case 1:
        {
            Roe_Solver_WITHOUT_entropy_fix(gridpts,Domain,Copy_Domain,initial_time,final_time,CFL,delta_x,gamma);
            cout << "Execution successful for Roe Solver WITHOUT entropy fix\n";
            for(int j=0;j<gridpts;j++)
            {
                Final_results(0,j) = Domain(0,j);                                                           // density
                Final_results(1,j) = (Domain(2,j) - pow(Domain(1,j),2)/(2*Domain(0,j)))*(gamma - 1);        // pressure
                Final_results(2,j) = Domain(1,j)/Domain(0,j);                                               // velocity
                Final_results(3,j) = Domain(2,j)/Domain(0,j) - 0.5*pow(Domain(1,j)/Domain(0,j),2);          // internal energy
            }
            // Results stored in CSV files to be extracted for post processing in MATLAB
            put_in_csv("Roe_WITHOUT_Entropy_Fix_Output.csv", gridpts, Final_results);
            break;
        }
        case 2:
        {
            Roe_Solver_WITH_entropy_fix(gridpts,Domain,Copy_Domain,initial_time,final_time,CFL,delta_x,gamma,epsilon);
            cout << "Execution successful for Roe Solver WITH entropy fix\n";
            for(int j=0;j<gridpts;j++)
            {
                Final_results(0,j) = Domain(0,j);                                                           // density
                Final_results(1,j) = (Domain(2,j) - pow(Domain(1,j),2)/(2*Domain(0,j)))*(gamma - 1);        // pressure
                Final_results(2,j) = Domain(1,j)/Domain(0,j);                                               // velocity
                Final_results(3,j) = Domain(2,j)/Domain(0,j) - 0.5*pow(Domain(1,j)/Domain(0,j),2);          // internal energy
            }
            // Results stored in CSV files to be extracted for post processing in MATLAB
            put_in_csv("Roe_WITH_Entropy_Fix_Output.csv", gridpts, Final_results);
            break;
        }
        case 3:
        {
            Steger_Warming(gridpts,Domain,Copy_Domain,initial_time,final_time,CFL,delta_x,gamma);
            cout << "Execution successful for Steger Warming Scheme\n";
            for(int j=0;j<gridpts;j++)
            {
                Final_results(0,j) = Domain(0,j);                                                           // density
                Final_results(1,j) = (Domain(2,j) - pow(Domain(1,j),2)/(2*Domain(0,j)))*(gamma - 1);        // pressure
                Final_results(2,j) = Domain(1,j)/Domain(0,j);                                               // velocity
                Final_results(3,j) = Domain(2,j)/Domain(0,j) - 0.5*pow(Domain(1,j)/Domain(0,j),2);          // internal energy
            }
            // Results stored in CSV files to be extracted for post processing in MATLAB
            put_in_csv("Steger_Warming_output.csv", gridpts, Final_results);
            break;
        }
    }
}
