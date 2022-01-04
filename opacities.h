#ifndef OPACITY_H
#define OPACITY_H
/*
    Class Opacity that reads in, keeps and interpolates over tabulated opacities
*/

#include <vector>

class Opacities
{
    public:
    void read_data(const char *s_temp_table, const char *s_size_table,
              const char *s_stellar_abs_table, const char *s_stellar_scat_table,
              const char *s_particle_abs_table, const char *s_particle_scat_table,
              bool log_tables_n);
    double stellar_abs(double A) const;
    double stellar_scat(double A) const;
    double particle_abs(double A, double T) const;
    double particle_scat(double A, double T) const;

    private:
    int temp_table_length, size_table_length;
    std::vector<double> temp_table, size_table;
    std::vector<double> stellar_abs_table, stellar_scat_table;
    std::vector<std::vector<double>> particle_abs_table, particle_scat_table;
    bool log_tables;

    double interpolate_1d(const std::vector<double> &opac_table, double A) const;
    double interpolate_2d(const std::vector<std::vector<double>> &opac_table, double A, double T) const;
};

#endif

extern Opacities opac;

