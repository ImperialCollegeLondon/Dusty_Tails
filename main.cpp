#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<vector>
#include<string>
#include<iostream>

#include "opacities.h"

using namespace std;

double kappa_abs()
{
    // read the planck-mean absorption, scattering
    string opacity_dir("./opacs_jankovic/calc_dust_opac/fayalite_F01");
    opacity_dir = opacity_dir + "/opac_";

    Opacities opac = Opacities((opacity_dir+"temp.dat").c_str(), (opacity_dir+"sizes.dat").c_str(),
        (opacity_dir+"planck_abs_stellar_3830.dat").c_str(), (opacity_dir+"planck_sca_stellar_3830.dat").c_str(),
        (opacity_dir+"planck_abs.dat").c_str(), (opacity_dir+"planck_sca.dat").c_str(),
        true);

    // interpolate over a grid and save to files
    FILE *fouta = fopen("../test_1a.dat", "w");
    FILE *foutb = fopen("../test_1b.dat", "w");
    for (double logA=log10(1e-5); logA<log10(1e-3); logA+=2./400.)
    {
        double A = pow(10., logA);
        for (double logT=log10(100.); logT<log10(2000.); logT+=(log10(2000.)-log10(100.))/410.)
        {
            double T = pow(10., logT);
            fprintf(fouta, "%e ", opac.particle_abs(A, T));
            fprintf(foutb, "%e ", opac.particle_scat(A, T));
        }
        fprintf(fouta, "\n");
        fprintf(foutb, "\n");
    }
    fclose(fouta);
    fclose(foutb);
}

int main()
{
    test_1();

    return 0;
}
