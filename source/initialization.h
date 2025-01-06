//
// Created by ilaria on 2019-11-13.
//

#pragma once

#include "constants.h"
#include "o2.h"
#include "rnd.h"
#include "robust_filesystem.h"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

struct Node{
    /*Fluctuating vector potential, with DIM spatial dimensions*/
    mutable std::array<double, DIM> A;
    /*Fixed guage field determined by an external rotation at fixed velocity (the equivalent of an external magnetic field when e!=0)*/
    mutable std::array<double, DIM> R_ext;
    /*complex order parameter of the NC SC components*/
    mutable std::array<O2, NC> Psi;

};

struct H_parameters{
    /*These values are externally given by an input file*/
    double K;
    double lambda;
    double e;
    double h;
    double b_low; //lowest beta for the parallel tempering procedure
    double b_high; // highest beta for the parallel tempering procedure
    /* A_fixed = (2\pi fx y_i, 2 \pi fy x_i)*/
    double fx;
    double fy;
    int init; //Initial conditions
    int london_app; //If 0 (false): GL model; if 1 (true): London model
    int phase_update; //If 0 (false): only density updates, fixed phases

};



struct MC_parameters{
    /*These values are externally given by an input file*/
    int nmisu; //total number of independent measures
    int tau; // estimate of the auto-correlation time
    int transient; // estimate of the thermalization time
    int freq_autosave; //frequency at which intermediate configuration are saved
    double lbox_l; //length of the box for the uniform distribution of l (cartesian transformation of Psi)
    double lbox_theta; //length of the box for the uniform distribution of theta (polar transformation of Psi --> phase)
    double lbox_A; //length of the box for the uniform distribution of dA (transformation of the vector potential)
//    The resulting files were too large. It seems better to save the spin configurations.
//    bool measure_corr; //if that is true than measure the correlation functions.
};

void initialize_lattice(const std::vector<Node> &Site, const fs::path & directory_read, struct H_parameters &Hp);
void initialize_Hparameters(struct H_parameters &Hp);
void initialize_MCparameters(struct MC_parameters &MCp);

