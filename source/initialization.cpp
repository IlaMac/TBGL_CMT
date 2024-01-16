//
// Created by ilaria on 2019-11-13.
//
#include "initialization.h"
#include "robust_filesystem.h"
#include "tid/tid.h"
#include "cfg/cfg.h"
#include "fmt/core.h"
void initialize_Hparameters(struct H_parameters &Hp){
    auto t_init = tid::tic_scope(__FUNCTION__);

    fs::path hp_init_file = fmt::format("{}/HP_init.txt", cfg::paths_dir::directory_parameters);


    if(fs::exists(hp_init_file)){
        FILE *fin= nullptr;
        /* clang-format off */
        if((fin=fopen(hp_init_file.c_str(), "r"))) {
            if(fscanf(fin, "%lf" , &Hp.K) == 0 ) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.lambda) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.e) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.h) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.b_low) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.b_high) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.fx) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf" , &Hp.fy) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%d", &Hp.init) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%d", &Hp.london_app) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            fclose(fin);
            /* clang-format on */
            //With this modification Hp.beta in not anymore part of the Hamiltonian parameters list
        }
    }else{

        Hp.K=-1;
        Hp.lambda=0.1;
        Hp.e=0;
        Hp.h= 1;
        Hp.b_low=0.1;
        Hp.b_high=10.0;
        /* A_fixed = (2\pi fx y_i, 2 \pi fy x_i)*/
        Hp.fx=0;
        Hp.fy=0.;
        Hp.init=4;
        Hp.london_app=0;

    }

}

void initialize_MCparameters(struct MC_parameters &MCp){
    auto t_init = tid::tic_scope(__FUNCTION__);
    fs::path mc_init_file = fmt::format("{}/MC_init.txt", cfg::paths_dir::directory_parameters);
    if(fs::exists(mc_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(mc_init_file.c_str(), "r"))) {
            /* clang-format off */
            if(fscanf(fin, "%d", &MCp.nmisu) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%d", &MCp.tau) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%d", &MCp.transient) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%d", &MCp.freq_autosave) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf", &MCp.lbox_l) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
	        if(fscanf(fin, "%lf", &MCp.lbox_theta) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            if(fscanf(fin, "%lf", &MCp.lbox_A) == 0) {fclose(fin); throw std::runtime_error("fscanf failed");}
            /* clang-format on */
            fclose(fin);
        }
    }else{
        MCp.nmisu=5000;
        MCp.tau=32;
        MCp.transient=1000;
        MCp.freq_autosave=6000;
        MCp.lbox_l=1.0;
        MCp.lbox_theta=C_PI*0.25;
        MCp.lbox_A=0.1;
    }

}

void initialize_lattice(const std::vector<Node> &Site, const fs::path & directory_read, struct H_parameters &Hp){
    using namespace cfg;
    auto t_init = tid::tic_scope(__FUNCTION__);
    fs::path psi_init_file = directory_read / "Psi_restart.bin";
    fs::path a_init_file = directory_read / "A_restart.bin";

    if(RESTART==1){
        psi_init_file = directory_read / "Psi_restart.bin";
        a_init_file = directory_read / "A_restart.bin";
    }
    else if (RESTART==2){
        psi_init_file = directory_read / "Psi_final.bin";
        a_init_file = directory_read / "A_final.bin";
    }

    /*Initialize external field*/
    for (size_t ix = 0; ix < Lx; ix++) {
        for (size_t iy = 0; iy < Ly; iy++) {
            size_t i= ix +iy*Lx;
            Site[i].R_ext[0]= (double)C_TWO_PI*(Hp.fx)* static_cast<double>(iy);
            Site[i].R_ext[1]= (double)C_TWO_PI*(Hp.fy)* static_cast<double>(ix);
            }
        }


        if((fs::exists(psi_init_file)) and (fs::exists(a_init_file)) and (RESTART!=0)){

        FILE *fPsi= nullptr;
        FILE *fA= nullptr;

        if((fPsi=fopen(psi_init_file.c_str(), "r")) and (fA=fopen(a_init_file.c_str(), "r")) ) {
            for(auto & s : Site){
                auto res_psi = fread(s.Psi.data(), sizeof(struct O2), NC, fPsi);
                auto res_A = fread(s.A.data(), sizeof(double), DIM, fA);
                if(res_psi < NC) {
                    fclose(fPsi);
                    throw std::runtime_error("res_psi < NC"); }
                if(res_A < NC) {
                    fclose(fA);
                    throw std::runtime_error("res_A < DIM"); }
            }
            fclose(fA);
            fclose(fPsi);
        }
    }else{
        if(Hp.init==0) { //ground state K>0
            double gamma=C_PI/2;
            for (size_t i = 0; i < N; i++) {
                Site[i].Psi[0].r = cos(gamma/2);
                Site[i].Psi[0].t = 0.;
                polar_to_cartesian(Site[i].Psi[0]);
                Site[i].Psi[1].r = sin(gamma/2);
                Site[i].Psi[1].t =  -C_PI/2;;
                polar_to_cartesian(Site[i].Psi[1]);
            }
        }
        else if(Hp.init==1) { //ground state K<0; lambda>0; gamma=C_PI/3
            double gamma=C_PI/3;
            for (size_t i = 0; i < N; i++) {
                Site[i].Psi[0].t = 0.;
                Site[i].Psi[0].r = sqrt(cos(gamma/2)*cos(gamma/2));
                polar_to_cartesian(Site[i].Psi[0]);
                Site[i].Psi[1].t = 0;
                Site[i].Psi[1].r = sqrt(sin(gamma/2)*sin(gamma/2));
                polar_to_cartesian(Site[i].Psi[1]);
            }
        }
        else if(Hp.init==2) { //ground state K<0; lambda>0; gamma=C_PI
            double gamma=C_PI;
            for (size_t i = 0; i < N; i++) {
                Site[i].Psi[0].t =  sqrt(rnd::uniform_double_box(0, 1));
                Site[i].Psi[0].r = sqrt(cos(gamma/2)*cos(gamma/2));
                polar_to_cartesian(Site[i].Psi[0]);
                Site[i].Psi[1].t = sqrt(rnd::uniform_double_box(0, 1));
                Site[i].Psi[1].r = sqrt(sin(gamma/2)*sin(gamma/2));
                polar_to_cartesian(Site[i].Psi[1]);
            }
        }
        else if(Hp.init==3) { //ground state for the anisotropic London model
            for (size_t i = 0; i < N; i++) {
                Site[i].Psi[0].t = 0.;
                Site[i].Psi[0].r = 1;
                polar_to_cartesian(Site[i].Psi[0]);
                Site[i].Psi[1].t = 0;
                Site[i].Psi[1].r = 0.5;
                polar_to_cartesian(Site[i].Psi[1]);
            }
        }
        else if((Hp.init<0) or (Hp.init>3)) { //Random
            for (size_t i = 0; i < N; i++) {
                Site[i].Psi[0].r = sqrt(rnd::uniform_double_box(0, 1));
                Site[i].Psi[1].r= sqrt(1. - Site[i].Psi[0].r*Site[i].Psi[0].r);
                for (size_t alpha = 0; alpha < NC; alpha++) {
                    Site[i].Psi[alpha].t = rnd::uniform_double_box(0, C_TWO_PI);
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }
    }

}



