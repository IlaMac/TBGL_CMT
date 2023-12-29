//
// Created by ilaria on 2019-11-13.
//
#include "initialization.h"
#include "robust_filesystem.h"

void initialize_Hparameters(struct H_parameters &Hp, const fs::path & directory_parameters){

    fs::path hp_init_file = directory_parameters / "HP_init.txt";


    if(fs::exists(hp_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(hp_init_file.c_str(), "r"))) {
            fscanf(fin, "%lf" , &Hp.K);
            fscanf(fin, "%lf" , &Hp.lambda);
            fscanf(fin, "%lf" , &Hp.e);
            fscanf(fin, "%lf" , &Hp.h);
            fscanf(fin, "%lf" , &Hp.b_low);
            fscanf(fin, "%lf" , &Hp.b_high);
            fscanf(fin, "%lf" , &Hp.fx);
            fscanf(fin, "%lf" , &Hp.fy);
            fscanf(fin, "%d", &Hp.init);
            fscanf(fin, "%d", &Hp.london_app);
            fclose(fin);
            //With this modification Hp.beta in not anymore part of the Hamiltonian parameters list
        }
    }else{

        Hp.K=-1;
        Hp.lambda=1;
        Hp.e=0;
        Hp.h= 1;
        Hp.b_low=1.5;
        Hp.b_high=3.0;
        /* A_fixed = (2\pi fx y_i, 2 \pi fy x_i)*/
        Hp.fx=0.;
        Hp.fy=0.;
        Hp.init=4;
        Hp.london_app=0;

    }

}

void initialize_MCparameters(struct MC_parameters &MCp, const fs::path & directory_parameters){
    fs::path mc_init_file = directory_parameters / "MC_init.txt";
    if(fs::exists(mc_init_file)){
        FILE *fin= nullptr;
        if((fin=fopen(mc_init_file.c_str(), "r"))) {
            fscanf(fin, "%d", &MCp.nmisu);
            fscanf(fin, "%d", &MCp.tau);
            fscanf(fin, "%d", &MCp.transient);
            fscanf(fin, "%d", &MCp.freq_autosave);
            fscanf(fin, "%lf", &MCp.lbox_l);
	        fscanf(fin, "%lf", &MCp.lbox_theta);
            fscanf(fin, "%lf", &MCp.lbox_A);

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

void initialize_lattice(const std::vector<Node> &Site, const fs::path & directory_read, int RESTART, struct H_parameters &Hp){

    fs::path psi_init_file = directory_read / "Psi_restart.bin";
    fs::path a_init_file = directory_read / "A_restart.bin";

    if(RESTART==1){
        fs::path psi_init_file = directory_read / "Psi_restart.bin";
        fs::path a_init_file = directory_read / "A_restart.bin";
    }
    else if (RESTART==2){
        fs::path psi_init_file = directory_read / "Psi_final.bin";
        fs::path a_init_file = directory_read / "A_final.bin";
    }

    /*Initialize external field*/
    for (size_t ix = 0; ix < Lx; ix++) {
        for (size_t iy = 0; iy < Ly; iy++) {
            size_t i= ix +iy*Lx;
            Site[i].R_ext[0]= C_TWO_PI*Hp.fx*iy;
            Site[i].R_ext[1]= C_TWO_PI*Hp.fy*ix;
            }
        }


        if((fs::exists(psi_init_file)) and (fs::exists(a_init_file)) and (RESTART!=0)){

        FILE *fPsi= nullptr;
        FILE *fA= nullptr;

        if((fPsi=fopen(psi_init_file.c_str(), "r")) and (fA=fopen(a_init_file.c_str(), "r")) ) {
            for(auto & s : Site){
                fread(s.Psi.data(), sizeof(struct O2), NC, fPsi);
                fread(s.A.data(), sizeof(double), DIM, fA);
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
                for (int alpha = 0; alpha < NC; alpha++) {
                    Site[i].Psi[alpha].t = rnd::uniform_double_box(0, C_TWO_PI);
                    polar_to_cartesian(Site[i].Psi[alpha]);
                }
            }
        }
    }

}



