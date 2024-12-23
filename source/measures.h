//
// Created by ilaria on 2019-11-16.
//

#ifndef MEASURES_H
#define MEASURES_H

#include "main.h"
#include <fstream>
#include "constants.h"
#include "montecarlo.h"
#include "initialization.h"
#include "o2.h"
#include <h5pp/h5pp.h>


struct Measures{
/*******THE ORDER SHOULD BE THE SAME OF THE H5T INSERT REGISTRATION**************/
    double E=0.; //Energy
    double E_josephson=0.;
    double E_lambda=0.;
    double E_kin=0.;
    double E_B=0.;
    double z2_m=0.; //magnetization (for the phase chirality of the two components)
    //Binder cumulant U=<m⁴>/(3*<m²>²)
    double mx_phasesum=0.; //x component of the magnetization of the single component phase
    double my_phasesum=0.; //y component of the magnetization of the single component phase
    double d_rhoz=0; //Dual stiffness along z
    double density_psi[NC] = {0.};
    double density_diff=0.;
    double Mx_nem=0.;
    double My_nem=0.;
    double Mz_nem=0.;
    double Mx_gamma=0.;
    double My_gamma=0.;
    double Mx_theta12=0.;
    double My_theta12=0.;
    double vortex_density[NC]={0.};
    double antivortex_density[NC]={0.};
    double v1v2_density[NC]={0.};
    double v1av2_density[NC]={0.};
    double DH_Ddi[NC]={0.}; //1st derivative in the twisted phase of the i component
    double D2H_Dd2i[NC]={0.}; //2nd derivative in the twisted phase of the i component
    // in this case the derivative of H with respect both the two current (mixed term) is zero
    //double D2H_Dd2ij[NC]={0}; //2nd mixed derivative in the twisted phases of the component i and j
    // see https://github.com/DavidAce/h5pp/blob/master/examples/example-04b-compound-datatype-fixed-arrays.cpp
    int my_rank = 0;
    void reset(){
        *this = Measures();
    }
};

struct Vdensity {
    double v1[2]={0.}; //first component +, second -
    double v2[2]={0.};
    double v1v2[2]; //first component ++, second --
    double v1av2[2]={0.}; //first component +-, second -+
};

void helicity_modulus(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site);
void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site);
void Z2_magnetization(struct Measures &mis, const std::vector<Node> &Site);
void magnetization_phasesum(struct Measures &mis, const std::vector<Node> &Site);
void energy(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site);
void nematic_order(struct Measures &mis, const std::vector<Node> &Site);
void save_lattice(const std::vector<Node> &Site, const fs::path & directory_write, const std::string &configuration);
void save_lattice_chargezero(const std::vector<Node> &Site, const fs::path & directory_write, const std::string  &configuration);
void new_vorticity(struct Measures &mis, std::vector<Vdensity> &local_vort_density ,struct H_parameters &Hp, const std::vector<Node> &Site );
void save_vortexlattice(const std::vector<Vdensity> &v_local, const fs::path & directory_write, const std::string & configuration);
#endif //MEASURES_H
