//
// Created by ilaria on 2022-01-12.
//

#include "wolff_mc.h"


//void growCluster_density(int i, int* clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta){
//
//    std::array<O2, NC> NewPsi;
//    std::array<O2, NC> OldPsi;
//    int ix, iy;
//    int ip, im, vec;
//    //addProbability for the link i,j has the form: 1- exp(-beta \Delta E_ij) with \Delta E_ij
//    double rand=0., addProbability;
//    double dE;
//    double inv_h2=1./(Hp.h*Hp.h);
//
//    ix= i%Lx;
//    iy= (i/Lx)%Ly;
//
//    //phase in i added to the cluster and flipped
//    clusterSpin[i]=1;
//    //int phase_diff = Site[i].Psi[0].t - Site[i].Psi[1].t;
//
//    OldPsi[0]= Site[i].Psi[0];
//    OldPsi[1]= Site[i].Psi[1];
//
//    NewPsi[0] = Site[i].Psi[0];
//    NewPsi[1] = Site[i].Psi[1];
//    /*switch of the densities*/
//    NewPsi[0].t = Site[i].Psi[1].r;
//    NewPsi[1].t = Site[i].Psi[0].r;
//    polar_to_cartesian(NewPsi[0]);
//    polar_to_cartesian(NewPsi[1]);
//
//    Site[i].Psi[0]= NewPsi[0];
//    Site[i].Psi[1]= NewPsi[1];
//
//    //Check of the neighbours of i
//    // if the neighbor spin does not belong to the cluster, but it has the preconditions to be added, then try to add it to the cluster
//    for (vec = 0; vec < DIM; vec++) {
//        if (vec == 0) {
//            ip=(ix == Lx-1 ? 0: ix+1) + Lx * iy;
//            im = (ix == 0 ? Lx-1: ix-1)+ Lx * iy;
//        }
//        if (vec == 1) {
//            ip = ix + Lx * ((iy == Ly-1 ? 0: iy+1)) ;
//            im = ix + Lx * ((iy == 0 ? Ly-1: iy-1));
//        }
//        if (clusterSpin[im]==0){
//            dE=0.;
//            for(int alpha=0; alpha<NC; alpha++){
//                //the gauge phase remains the same
//                double gauge_phase = Site[i].Psi[alpha].t - Site[im].Psi[alpha].t + Hp.h * Hp.e * Site[i].A[vec];
//                double old_E = ;

//                dE-=inv_h2 * (Site[i].Psi[alpha].r* Site[im].Psi[alpha].r)*(cos(new_gauge_phase) - cos(old_gauge_phase));
//            }
//            addProbability=1-exp(-my_beta*dE);
//            rand= rn::uniform_real_box(0, 1);;
//            if (rand < addProbability){
//                growCluster(im, clusterSpin, Site, MCp, Hp, my_beta);}
//        }
//        if (clusterSpin[ip]==0){
//            dE=0.;
//            for(int alpha=0; alpha<NC; alpha++){
//                double new_gauge_phase = Site[ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.h * Hp.e * Site[ip].A[vec];
//                double old_gauge_phase = Site[ip].Psi[alpha].t - OldPsi[alpha].t + Hp.h * Hp.e * Site[ip].A[vec];
//                //the modulus remains the same
//                dE-=inv_h2 * (Site[i].Psi[alpha].r* Site[ip].Psi[alpha].r)*(cos(new_gauge_phase) - cos(old_gauge_phase));
//            }
//            addProbability=1-exp(-my_beta*dE);
//            rand= rn::uniform_real_box(0, 1);;
//            if (rand < addProbability){
//                growCluster(ip, clusterSpin, Site, MCp, Hp, my_beta);}
//        }
//    }
//
//    return;
//}
//
//void wolff_density(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta){
//
//    int i, iseed, count=0;
//    int clusterSpin[N]={0};
//
//    /*choose randomly a site of the lattice*/
//    iseed = rn::uniform_integer_box(0, N-1);
//
//    growCluster(iseed, clusterSpin, Site, MCp, Hp, my_beta);
//
//    for(int iy=0; iy<Ly; iy++){
//        for(int ix=0; ix<Lx; ix++){
//            i=ix+ Lx*(iy);
//            if( clusterSpin[i]==1){count++;}
//        }
//    }
//    return;
//}

void growCluster_BTRS(size_t i, std::vector<size_t> & clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta){

    std::array<O2, NC> NewPsi;
    std::array<O2, NC> OldPsi;
    size_t ix, iy;
    size_t ip, im, vec;
    //addProbability for the link i,j has the form: 1- exp(-beta \Delta E_ij) with \Delta E_ij
    double rand=0., addProbability;
    double dE;
    double inv_h2=1./(Hp.h*Hp.h);

    ix= i%Lx;
    iy= (i/Lx)%Ly;

    //phase in i added to the cluster and flipped
    clusterSpin[i]=1;
    //int phase_diff = Site[i].Psi[0].t - Site[i].Psi[1].t;

    OldPsi[0]= Site[i].Psi[0];
    OldPsi[1]= Site[i].Psi[1];

    NewPsi[0] = Site[i].Psi[0];
    NewPsi[1] = Site[i].Psi[1];
    /*switch of the phase*/
    NewPsi[0].t = Site[i].Psi[1].t;
    NewPsi[1].t = Site[i].Psi[0].t;
    polar_to_cartesian(NewPsi[0]);
    polar_to_cartesian(NewPsi[1]);

    Site[i].Psi[0]= NewPsi[0];
    Site[i].Psi[1]= NewPsi[1];

    //Check of the neighbours of i
    // if the neighbor spin does not belong to the cluster, but it has the preconditions to be added, then try to add it to the cluster
    for (vec = 0; vec < DIM; vec++) {
        if (vec == 0) {
            ip=(ix == Lx-1 ? 0: ix+1) + Lx * iy;
            im = (ix == 0 ? Lx-1: ix-1)+ Lx * iy;
        }
        if (vec == 1) {
            ip = ix + Lx * ((iy == Ly-1 ? 0: iy+1)) ;
            im = ix + Lx * ((iy == 0 ? Ly-1: iy-1));
        }
        if (clusterSpin[im]==0){
           dE=0.;
            for(int alpha=0; alpha<NC; alpha++){
                double new_gauge_phase = Site[i].Psi[alpha].t - Site[im].Psi[alpha].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
                double old_gauge_phase = OldPsi[alpha].t - Site[im].Psi[alpha].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
                //the modulus remains the same and so the value of cos(2 \phi_12)
                dE-= (Site[i].Psi[alpha].r* Site[im].Psi[alpha].r)*(cos(new_gauge_phase) - cos(old_gauge_phase));
            }
            addProbability=1-exp(-my_beta*dE);
            rand= rn::uniform_real_box(0, 1);;
            if (rand < addProbability){
                growCluster_BTRS(im, clusterSpin, Site, MCp, Hp, my_beta);}
        }
        if (clusterSpin[ip]==0){
            dE=0.;
            for(int alpha=0; alpha<NC; alpha++){
                double new_gauge_phase = Site[ip].Psi[alpha].t - Site[i].Psi[alpha].t + Site[ip].R_ext[vec] + Hp.h * Hp.e * Site[ip].A[vec];
                double old_gauge_phase = Site[ip].Psi[alpha].t - OldPsi[alpha].t + Site[ip].R_ext[vec] + Hp.h * Hp.e * Site[ip].A[vec];
                //the modulus remains the same
                dE-= (Site[i].Psi[alpha].r* Site[ip].Psi[alpha].r)*(cos(new_gauge_phase) - cos(old_gauge_phase));
            }
            addProbability=1-exp(-my_beta*dE);
            rand= rn::uniform_real_box(0, 1);;
            if (rand < addProbability){
                growCluster_BTRS(ip, clusterSpin, Site, MCp, Hp, my_beta);}
        }
    }

}

void wolff_BTRS(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta){

    size_t  count=0;
    std::vector<size_t> clusterSpin(N,0); //vector of int of size N filled with 0

    /*choose randomly a site of the lattice*/
    auto iseed = static_cast<size_t>(rn::uniform_integer_box(0, N - 1));

    growCluster_BTRS(iseed, clusterSpin, Site, MCp, Hp, my_beta);

    for(size_t iy=0; iy<Ly; iy++){
        for(size_t ix=0; ix<Lx; ix++){
            size_t i=ix+ Lx*(iy);
            if( clusterSpin[i]==1){count++;}
        }
    }
}


void growCluster_nemK(size_t i, size_t alpha_up, std::vector<size_t> & clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta){

    //here I want to go from a phase difference equal to 0 to a phase difference equal to \pi.

    std::array<O2, NC> NewPsi{};
    std::array<O2, NC> OldPsi{};

    size_t ip, im, vec;
    //addProbability for the link i,j has the form: 1- exp(-beta \Delta E_ij) with \Delta E_ij
    double rand=0., addProbability;
    double inv_h2=1./(Hp.h*Hp.h);

    size_t ix= i%Lx;
    size_t iy= (i/Lx)%Ly;

    //phase in i added to the cluster and flipped
    clusterSpin[i]=1;


    OldPsi[0]= Site[i].Psi[0];
    OldPsi[1]= Site[i].Psi[1];

    NewPsi[0] = Site[i].Psi[0];
    NewPsi[1] = Site[i].Psi[1];

    /*Rotate the phase of alpha_up by \pi*/
    NewPsi[alpha_up].t = static_cast<_FPTYPE>(Site[i].Psi[alpha_up].t + C_PI);

    polar_to_cartesian(NewPsi[0]);
    polar_to_cartesian(NewPsi[1]);

    Site[i].Psi[0]= NewPsi[0];
    Site[i].Psi[1]= NewPsi[1];

    //Check of the neighbours of i
    // if the neighbor spin does not belong to the cluster, but it has the preconditions to be added, then try to add it to the cluster
    for (vec = 0; vec < DIM; vec++) {
        if (vec == 0) {
            ip=(ix == Lx-1 ? 0: ix+1) + Lx * iy;
            im = (ix == 0 ? Lx-1: ix-1)+ Lx * iy;
        }
        if (vec == 1) {
            ip = ix + Lx * ((iy == Ly-1 ? 0: iy+1)) ;
            im = ix + Lx * ((iy == 0 ? Ly-1: iy-1));
        }
        if (clusterSpin[im]==0){
            double new_gauge_phase = Site[i].Psi[alpha_up].t - Site[im].Psi[alpha_up].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
            double old_gauge_phase = OldPsi[alpha_up].t - Site[im].Psi[alpha_up].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
            //the modulus remains the same
            double dE= -(Site[i].Psi[alpha_up].r* Site[im].Psi[alpha_up].r)*(cos(new_gauge_phase) - cos(old_gauge_phase));
            addProbability=1-exp(-my_beta*dE);
            rand= rn::uniform_real_box(0, 1);;
            if (rand < addProbability){
                growCluster_nemK(im, alpha_up, clusterSpin, Site, MCp, Hp, my_beta);}
        }
        if (clusterSpin[ip]==0){
            double new_gauge_phase = Site[ip].Psi[alpha_up].t - Site[i].Psi[alpha_up].t + Site[ip].R_ext[vec] + Hp.h * Hp.e * Site[ip].A[vec];
            double old_gauge_phase = Site[ip].Psi[alpha_up].t - OldPsi[alpha_up].t + Site[ip].R_ext[vec] + Hp.h * Hp.e * Site[ip].A[vec];
            //the modulus remains the same
            double dE= -(Site[i].Psi[alpha_up].r* Site[ip].Psi[alpha_up].r)*(cos(new_gauge_phase) - cos(old_gauge_phase));
            addProbability=1-exp(-my_beta*dE);
            rand= rn::uniform_real_box(0, 1);;
            if (rand < addProbability){
                growCluster_nemK(ip, alpha_up, clusterSpin, Site, MCp, Hp, my_beta);}
        }
    }
}


void wolff_nemK(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta){

    size_t count=0;
    std::vector<size_t> clusterSpin(N,0); //vector of int of size N filled with 0

    /*choose randomly a site of the lattice*/
    auto iseed = static_cast<size_t>(rn::uniform_integer_box(0, N - 1));

    /*choose randomly which of the two phases will be updated*/
    auto alpha_up = static_cast<size_t>(rn::uniform_integer_box(0, 1));
    std::cout<< alpha_up<< std::endl;

    growCluster_nemK(iseed, alpha_up, clusterSpin, Site, MCp, Hp, my_beta);

    for(size_t iy=0; iy<Ly; iy++){
        for(size_t ix=0; ix<Lx; ix++){
            size_t i=ix+ Lx*(iy);
            if( clusterSpin[i]==1){count++;}
        }
    }
}

