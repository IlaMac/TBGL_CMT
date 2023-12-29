#include "montecarlo.h"
#include "class_tic_toc.h"
#include "main.h"
#include "rnd.h"

void metropolis(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){

    double l, d_theta, d_A, rand;
    double acc_rate=0.5, acc_A=0., acc_theta=0.;
    std::array<O2, NC> NewPsi;
    std::array<O2, NC> OldPsi;
    double NewA, OldA;
    double newE, oldE, minus_deltaE;
    double h2=(Hp.h*Hp.h);

    for (int iy = 0; iy < Ly; iy++) {
        for (int ix = 0; ix < Lx; ix++) {
            int i = ix + Lx * (iy);
            /*choose randomly a site of the lattice*/
            //i = rnd::uniform_integer_box(0, N-1);
            //ix=i%Lx;
            //iy=i/Ly;
            if(Hp.london_app == 0) {
                /*************PSI UPDATE: density update with total density contraint **********/
                OldPsi[0] = Site[i].Psi[0];
                OldPsi[1] = Site[i].Psi[1];
                NewPsi[0] = Site[i].Psi[0];
                NewPsi[1] = Site[i].Psi[1];
                l = rnd::uniform_double_box(0, 1);
                NewPsi[0].r = sqrt(l);
                NewPsi[1].r = sqrt(1 - l);
                polar_to_cartesian(NewPsi[0]);
                polar_to_cartesian(NewPsi[1]);
                oldE = local_HPsi(OldPsi, ix, iy, Hp, Site);
                newE = local_HPsi(NewPsi, ix, iy, Hp, Site);
                minus_deltaE = h2 * (oldE - newE);
                if (minus_deltaE > 0) {
                    Site[i].Psi[0] = NewPsi[0];
                    Site[i].Psi[1] = NewPsi[1];
                    //oldE = newE; //I don't have to recompute oldE
                } else {
                    rand = rnd::uniform_double_box(0, 1);
                    //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
                    if (rand < exp(my_beta * minus_deltaE)) {
                        Site[i].Psi[0] = NewPsi[0];
                        Site[i].Psi[1] = NewPsi[1];
                        //oldE = newE; //I don't have to recompute oldE
                    }
                }
            }
//            if(Hp.london_app == 2) {
//                /************* PSI UPDATE: density update with NO total density constraint **********/
//                /************* In this case I have to add a potential term in the free energy  **********/
//                for (int alpha = 0; alpha < NC; alpha++) {
//                    OldPsi[0] = Site[i].Psi[0];
//                    OldPsi[1] = Site[i].Psi[1];
//                    NewPsi[0] = Site[i].Psi[0];
//                    NewPsi[1] = Site[i].Psi[1];
//                    l = rnd::uniform_double_box(0, 1);
//                    NewPsi[alpha].r = sqrt(l);
//                    polar_to_cartesian(NewPsi[alpha]);
//                    oldE = local_HPsi(OldPsi, ix, iy, Hp, Site);
//                    newE = local_HPsi(NewPsi, ix, iy, Hp, Site);
//                    minus_deltaE = h2 * (oldE - newE);
//                    if (minus_deltaE > 0) {
//                        Site[i].Psi[alpha] = NewPsi[alpha];
//                        //oldE = newE; //I don't have to recompute oldE
//                    } else {
//                        rand = rnd::uniform_double_box(0, 1);
//                        //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
//                        if (rand < exp(my_beta * minus_deltaE)) {
//                            Site[i].Psi[alpha] = NewPsi[alpha];
//                            //oldE = newE; //I don't have to recompute oldE
//                        }
//                    }
//                }
//            }
            /*************PSI UPDATE: phase update **********/

            for (int alpha = 0; alpha < NC; alpha++) {
                OldPsi[0] = Site[i].Psi[0];
                OldPsi[1] = Site[i].Psi[1];
                NewPsi[0] = Site[i].Psi[0];
                NewPsi[1] = Site[i].Psi[1];
                d_theta = rnd::uniform_double_box(-MCp.lbox_theta, MCp.lbox_theta);
                NewPsi[alpha].t = fmod(OldPsi[alpha].t + d_theta, C_TWO_PI);
                NewPsi[alpha].r = OldPsi[alpha].r;
                polar_to_cartesian(NewPsi[alpha]);
                oldE = local_HPsi(OldPsi, ix, iy, Hp, Site);
                newE = local_HPsi(NewPsi, ix, iy, Hp, Site);
                minus_deltaE = h2 * (oldE - newE);
                if (minus_deltaE > 0) {
                    Site[i].Psi[alpha] = NewPsi[alpha];
                    acc_theta++;
                } else {
                    rand = rnd::uniform_double_box(0, 1);
                    //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
                    if (rand < exp(my_beta * minus_deltaE)) {
                        Site[i].Psi[alpha] = NewPsi[alpha];
                        acc_theta++;
                    }
                }
            }
        }
    }

    if (Hp.e != 0) {
    /**********VECTOR POTENTIAL UPDATE********/
        for (int iy = 0; iy < Ly; iy++) {
            for (int ix = 0; ix < Lx; ix++) {
                int i = ix + Lx * (iy );
                for (int vec = 0; vec < DIM; vec++) {
                    //Update of A
                    OldA = Site[i].A[vec];
                    oldE = local_HA(OldA, ix, iy, vec, Hp, Site);
                    d_A = rnd::uniform_double_box(-MCp.lbox_A, MCp.lbox_A);
                    NewA = OldA + d_A;
                    newE = local_HA(NewA, ix, iy, vec, Hp, Site);
                    minus_deltaE = h2 * (oldE - newE);
                    if (minus_deltaE > 0.) {
                        Site[i].A[vec] = NewA;
                        acc_A++;
                    } else {
                        rand = rnd::uniform_double_box(0, 1);
                        //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
                        if (rand < exp(my_beta * minus_deltaE)) {
                            Site[i].A[vec] = NewA;
                            acc_A++;
                        }
                    }
                }
            }
        }
    }

    acc_theta=(double) acc_theta/(NC*N);
    acc_A=(double) acc_A/(DIM*N);
    MCp.lbox_theta= MCp.lbox_theta*((0.5*acc_theta/acc_rate)+0.5);
    MCp.lbox_A= MCp.lbox_A*((0.5*acc_A/acc_rate)+0.5);

}



double local_HPsi(std::array<O2, NC> &Psi, int ix, int iy, H_parameters &Hp, const std::vector<Node> &Site) {

    double h_Kinetic=0., h_Josephson=0., h_lambda=0., h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double gauge_phase1, gauge_phase2;
    double gamma=0;
    int nn_ip, nn_im;
    int i=ix +Lx*(iy);

    int ip=(ix == Lx-1 ? 0: ix+1);
    int ipx= ip+Lx*(iy);
    int jp=(iy == Ly-1 ? 0: iy+1);
    int ipy= ix+(Lx*jp);
    int imx= (ix == 0 ? Lx-1: ix-1)+Lx*(iy);
    int imy= ix+Lx*((iy == 0 ? Ly-1: iy-1));

    //Compute the local Energy respect to a given component (alpha) of the matter field Psi and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving Psi

    for(int alpha=0; alpha<NC; alpha ++) {
        //Kinetic= -(1/h²)*\sum_k=1,2 (|Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) + (|Psi_{alpha}(r-k)||Psi_{alpha}(r)|* cos(theta_{alpha}(r) - theta_{alpha}(r-k) +h*e*A_k(r-k)))
        for (int vec = 0; vec < DIM; vec++) {
            if (vec == 0) {
                nn_ip = ipx;
                nn_im = imx;
            } else if (vec == 1) {
                nn_ip = ipy;
                nn_im = imy;
            }
            gauge_phase1 = Site[nn_ip].Psi[alpha].t - Psi[alpha].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
            gauge_phase2 = Psi[alpha].t - Site[nn_im].Psi[alpha].t  + Site[nn_im].R_ext[vec]  + Hp.h * Hp.e * Site[nn_im].A[vec];
            h_Kinetic -= inv_h2 * (Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * cos(gauge_phase1);
            h_Kinetic -= inv_h2 * (Psi[alpha].r * Site[nn_im].Psi[alpha].r) * cos(gauge_phase2);
        }
    }
    //Biquadratic Josephson= 2*\beta_2 * |Psi_{1}(r)|^2|Psi_{2}(r)|^2* cos(2*(theta_{1}(r) - theta_{2}(r)))
    //h_Josephson += 2* Hp.b2*(Psi[0].r*Psi[1].r)*(Psi[0].r*Psi[1].r)*cos(2*(Psi[0].t -Psi[1].t));
    h_Josephson +=  Hp.K*(Psi[0].r*Psi[1].r)*(Psi[0].r*Psi[1].r)*(cos(2*(Psi[0].t -Psi[1].t)) -1.);

    if(Hp.lambda !=0 ) {
        gamma=2* atan2( Psi[1].r, Psi[0].r);
        h_lambda+=Hp.lambda*(cos(3*gamma) + 3*cos(gamma)*sin(gamma)*sin(gamma)*sin((Psi[0].t -Psi[1].t))*sin((Psi[0].t -Psi[1].t)));
    }

    h_tot=  h_Kinetic + h_Josephson + h_lambda;
    return h_tot;
}

double local_HA(double A, int ix, int iy,  int vec,  struct H_parameters &Hp, const std::vector<Node> &Site){

    double h_Kinetic=0., h_B, h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double F2_A=0., F_A;
    double gauge_phase1;
    int alpha, i, l;

    int nn_ip, nn_ipl, nn_iml;
    i=ix +static_cast<int>(Lx)*(iy);

    int ip=(ix == static_cast<int>(Lx-1) ? 0: ix+1);
    int ipx= ip+static_cast<int>(Lx)*(iy);
    int jp=(iy == static_cast<int>(Ly-1) ? 0: iy+1);
    int ipy= ix+(static_cast<int>(Lx)*jp);


    int imx= (ix == 0 ? static_cast<int>(Lx-1): ix-1)+static_cast<int>(Lx)*(iy);
    int imy= ix+static_cast<int>(Lx)*((iy == 0 ? static_cast<int>(Ly-1): iy-1));

    if(vec==0){
        nn_ip=ipx;
    }else if(vec==1){
        nn_ip=ipy;
    }

    //Compute the local Energy respect to a given component (alpha) of the vector potential A and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving A

    //Kinetic= -(1/h²)*\sum_k=1,2 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(alpha=0; alpha<NC; alpha++) {
        gauge_phase1=Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t +  Site[i].R_ext[vec] +  Hp.h*Hp.e*A;
        h_Kinetic -= inv_h2* Site[nn_ip].Psi[alpha].r*Site[i].Psi[alpha].r*cos(gauge_phase1);
	 }

    //All the plaquettes involving A_vec(i)
    for(l=0; l<DIM; l++){
        if(l==0){
            nn_ipl=ipx;
        }else if(l==1){
            nn_ipl=ipy;
        }
        if(l!= vec){
            F_A=(A + (Site[nn_ip].A[l] ) - Site[nn_ipl].A[vec] - Site[i].A[l]);
            F2_A+=(F_A*F_A);
        }
    }
    for(l=0; l<DIM; l++){
        if(l==0){
            nn_iml=imx;
        }else if(l==1){
            nn_iml=imy;
        }
        if(l!= vec){
            F_A=(Site[nn_iml].A[vec]+ Site[nn(nn_iml, vec, 1)].A[l] - A - Site[nn_iml].A[l]);
            F2_A+=(F_A*F_A);
        }
    }

    h_B=(0.5*inv_h2)*F2_A;
    /******/
    h_tot= h_Kinetic + h_B;
    return h_tot;
}


