#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include "constants.h"
#include "butcher.h"
#include "functions.h"
#include "particle.h"
#include <iostream>

#include "opacities.h"

using namespace std;

void Opacities::read_data(const char *s_temp_table, const char *s_size_table,
                     const char *s_stellar_abs_table, const char *s_stellar_scat_table,
                     const char *s_particle_abs_table, const char *s_particle_scat_table,
                     bool log_tables_n)
{
    // open files
    
    FILE *f_temp_table = fopen(s_temp_table, "r");
    FILE *f_size_table = fopen(s_size_table, "r");
    FILE *f_stellar_abs_table = fopen(s_stellar_abs_table, "r");
    FILE *f_stellar_scat_table = fopen(s_stellar_scat_table, "r");
    FILE *f_particle_abs_table = fopen(s_particle_abs_table, "r");
    FILE *f_particle_scat_table = fopen(s_particle_scat_table, "r");

    
    // declare variables
    double x;

    // read temperature and size tables from files
    temp_table = {};
    while (fscanf(f_temp_table, "%lf ", &x) != EOF)
        temp_table.push_back(x);
        
    size_table = {};
    while (fscanf(f_size_table, "%lf ", &x) != EOF)
        size_table.push_back(x);
       
   
    // keep table dimensions
    temp_table_length = temp_table.size();
    size_table_length = size_table.size();
    
    // read absorption and scattering coefficients to stellar radiation
    stellar_abs_table = {};
    for (int i=0; i<size_table_length; i++)
    {   
       
        fscanf(f_stellar_abs_table, "%lf", &x);
        stellar_abs_table.push_back(x);
        
    }
    
    stellar_scat_table = {};
    for (int i=0; i<size_table_length; i++)
    {
        fscanf(f_stellar_scat_table, "%lf", &x);
        stellar_scat_table.push_back(x);
    }
    
    // read absorption and scattering coefficients to particle own radiation
    particle_abs_table = {};
    for (int i=0; i<size_table_length; i++)
    {
        vector<double> opac_line = {};
        for (int j=0; j<temp_table_length; j++)
        {
            fscanf(f_particle_abs_table, "%lf", &x);
            opac_line.push_back(x);
        }
        particle_abs_table.push_back(opac_line);
    }

    particle_scat_table = {};
    for (int i=0; i<size_table_length; i++)
    {
        vector<double> opac_line = {};
        for (int j=0; j<temp_table_length; j++)
        {
            fscanf(f_particle_scat_table, "%lf", &x);
            opac_line.push_back(x);
        }
        particle_scat_table.push_back(opac_line);
    }
   
    // set whether the temperature and size tables are log-spaced
    log_tables = log_tables_n;

    // close files
    fclose(f_temp_table);
    fclose(f_size_table);
    fclose(f_stellar_abs_table);
    fclose(f_stellar_scat_table);
    fclose(f_particle_abs_table);
    fclose(f_particle_scat_table);
}

double Opacities::stellar_abs(double A) const
{
    return interpolate_1d(stellar_abs_table, A);
}

double Opacities::stellar_scat(double A) const
{
    return interpolate_1d(stellar_scat_table, A);
}

double Opacities::particle_abs(double A, double T) const
{
    return interpolate_2d(particle_abs_table, A, T);
}

double Opacities::particle_scat(double A, double T) const
{
    return interpolate_2d(particle_scat_table, A, T);
}

double Opacities::interpolate_1d(const vector<double> &opac_table, double A) const
{
    // find the size interval in which to interpolate
    int x = size_table_length/2; // approx. index
    if (log_tables)
        x = (int)floor((log10(A)-log10(size_table[0])) * (size_table_length-1) / (log10(size_table[size_table_length-1])-log10(size_table[0])));
    else
        x = (int)floor((A-size_table[0]) * (size_table_length-1) / (size_table[size_table_length-1]-size_table[0]));
    if (x < 0)
        x = 0;

    // if out of bounds, return extrapolated values, otherwise check the interval
    if (A <= size_table[0])
    {
        return opac_table[0];
    }
    else if (A >= size_table[size_table_length-1]-1e-11)
    {
        return opac_table[size_table_length-1];
    }
    else
    {
        while (A > size_table[x+1])
            x++;
        while (A < size_table[x])
            x--;
    }

    // linear interpolation
    return opac_table[x] + (A-size_table[x])*(opac_table[x+1]-opac_table[x])/(size_table[x+1]-size_table[x]);
}

double Opacities::interpolate_2d(const vector<vector<double>> &opac_table, double A, double T) const
{
    // find the size interval in which to interpolate
    int x = size_table_length/2; // approx. index
    if (log_tables)
        x = (int)floor((log10(A)-log10(size_table[0])) * (size_table_length-1) / (log10(size_table[size_table_length-1])-log10(size_table[0])));
    else
        x = (int)floor((A-size_table[0]) * (size_table_length-1) / (size_table[size_table_length-1]-size_table[0]));
    if (x < 0)
        x = 0;
    if (x >= size_table_length)
        x = size_table_length-2;

    if (A <= size_table[0])
    {
        x = 0;
    }
    else if (A >= size_table[size_table_length-1])
    {
        x = size_table_length-1;
    }
    else
    {
        while (A > size_table[x+1])
            x++;
        while (A < size_table[x])
            x--;
    }

    // find the temperature interval in which to interpolate
    int y = temp_table_length/2; // approx. index
    if (log_tables)
        y = (int)floor((log10(T)-log10(temp_table[0])) * (temp_table_length-1) / (log10(temp_table[temp_table_length-1])-log10(temp_table[0])));
    else
        y = (int)floor((T-temp_table[0]) * (temp_table_length-1) / (temp_table[temp_table_length-1]-temp_table[0]));
    if (y < 0)
        y = 0;
    if (y >= temp_table_length)
        y = temp_table_length-2;

    if (T <= temp_table[0])
    {
        y = 0;
    }
    else if (T >= temp_table[temp_table_length-1])
    {
        y = temp_table_length-1;
    }
    else
    {
        while (T > temp_table[y+1])
            y++;
        while (T < temp_table[y])
            y--;
    }

    // handling cases where input is outside the limits of temperature, size tables
    if ((A <= size_table[0] || A >= size_table[size_table_length-1]) && !(T <= temp_table[0] || T >= temp_table[temp_table_length-1]))
    {
        return opac_table[x][y] + (T-temp_table[y]) * (opac_table[x][y+1]-opac_table[x][y]) / (temp_table[y+1]-temp_table[y]);
    }
    else if (!(A <= size_table[0] || A >= size_table[size_table_length-1]) && (T <= temp_table[0] || T >= temp_table[temp_table_length-1]))
    {
        return opac_table[x][y] + (A-size_table[x]) * (opac_table[x+1][y]-opac_table[x][y]) / (size_table[x+1]-size_table[x]);
    }
    else if ((A <= size_table[0] || A >= size_table[size_table_length-1]) && (T <= temp_table[0] || T >= temp_table[temp_table_length-1]))
    {
        return opac_table[x][y];
    }

    // bilinear interpolation
    double fQ11 = opac_table[x][y];
    double fQ12 = opac_table[x+1][y];
    double fQ21 = opac_table[x][y+1];
    double fQ22 = opac_table[x+1][y+1];

    double fxy1 = (temp_table[y+1] - T)/(temp_table[y+1] - temp_table[y])*fQ11 + (T - temp_table[y])/(temp_table[y+1] - temp_table[y])*fQ21;
    double fxy2 = (temp_table[y+1] - T)/(temp_table[y+1] - temp_table[y])*fQ12 + (T - temp_table[y])/(temp_table[y+1] - temp_table[y])*fQ22;

    double kappa = (size_table[x+1] - A)/(size_table[x+1] - size_table[x])*fxy1 + (A - size_table[x])/(size_table[x+1] - size_table[x])*fxy2;

    // a check kappa is positive, just in case something went wrong
    if (kappa <= 0.)
    {
        printf("Opacity error: opacity has non-positive value \n");
        printf("T = %e, A = %e \n", T, A);
        exit(-1);
    }

    return kappa;
}
