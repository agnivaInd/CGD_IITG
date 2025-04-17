#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

void initialize_U(int total_points, double delta_x, double U[])
{
    double x_start = -10;
    double xi = x_start + delta_x/2;
    double e = exp(1.0);
    double pi = 3.1415;

    for(int i=0;i<total_points;i++)
    {
        double cell_lower_limit = xi - delta_x/2;
        double cell_upper_limit = xi + delta_x/2;

        double intra_cell_spacing = (cell_upper_limit - cell_lower_limit)/100;
        double intra_x = cell_lower_limit;
        double sum = 0;

        for(int j=0;j<100;j++)
        {
            if((j==0)||(j==99))
            {
                sum = sum + pow(e,-(pow(intra_x-2,2))/2)/pow(2*pi,0.5) - pow(e,-(pow(intra_x+2,2))/2)/pow(2*pi,0.5);
            }
            else
            {
                sum = sum + 2*(pow(e,-(pow(intra_x-2,2))/2)/pow(2*pi,0.5) - pow(e,-(pow(intra_x+2,2))/2)/pow(2*pi,0.5));
            }
            intra_x = intra_x + intra_cell_spacing;
        }

        U[i] = 0.5*sum*intra_cell_spacing/(cell_upper_limit - cell_lower_limit);
        xi = xi + delta_x;
    }
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

void Godunov(int total_points, double delta_t, double delta_x, double total_number_of_timesteps, double U[], double U_new[])
{
    for(int t=0;t<total_number_of_timesteps;t++)
    {
        for(int j=0;j<total_points;j++)
        {
            if(j==0)
            {
                U_new[j] = U[j+1];
            }
            else if(j==total_points-1)
            {
                U_new[j] = U[j-1];
            }
            else
            {
                double flux_right;
                double flux_left;

                if(U[j]>U[j+1])
                {
                    double shock_speed = (U[j]+U[j+1])/2;
                    if(shock_speed<0)
                    {
                        flux_right = pow(U[j+1],2)/2;
                    }
                    else
                    {
                        flux_right = pow(U[j],2)/2;
                    }
                }
                else
                {
                    if((U[j]>0)&&(U[j+1]>0))
                    {
                        flux_right = pow(U[j],2)/2;
                    }
                    else if((U[j]<0)&&(U[j+1]<0))
                    {
                        flux_right = pow(U[j+1],2)/2;
                    }
                    else if ((U[j]<0)&&(U[j+1]>0))
                    {
                        flux_right = pow((U[j]+U[j+1])/2,2)/2;
                    }
                }

                if(U[j-1]>U[j])
                {
                    double shock_speed = (U[j-1]+U[j])/2;
                    if(shock_speed<0)
                    {
                        flux_left = pow(U[j],2)/2;
                    }
                    else
                    {
                        flux_left = pow(U[j-1],2)/2;
                    }
                }
                else
                {
                    if((U[j-1]>0)&&(U[j]>0))
                    {
                        flux_left = pow(U[j-1],2)/2;
                    }
                    else if((U[j-1]<0)&&(U[j]<0))
                    {
                        flux_left = pow(U[j],2)/2;
                    }
                    else if ((U[j-1]<0)&&(U[j]>0))
                    {
                        flux_left = pow((U[j-1]+U[j])/2,2)/2;
                    }
                }
                U_new[j] = U[j] - (delta_t/delta_x)*(flux_right - flux_left);
            }
        }

        for(int j=0;j<total_points;j++)
        {
            U[j] = U_new[j];
        }
    }
}

int main()
{
    double domain_length = 20;

    int total_points = 1000;
    double delta_x = domain_length/total_points;
    double delta_t = 0.001;
    double endtime = 20;
    int total_number_of_timesteps = endtime/delta_t;

    double U[total_points];
    double U_new[total_points];

    initialize_U(total_points,delta_x,U);
    initialize_U(total_points,delta_x,U_new);

    Godunov(total_points,delta_t,delta_x,total_number_of_timesteps,U,U_new);

    put_in_csv("output.csv",total_points,U);
}

