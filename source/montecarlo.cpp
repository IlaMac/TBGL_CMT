#include "montecarlo.h"
#include "main.h"
#include "rnd.h"
#include "tid/tid.h"
#include "cfg/cfg.h"

void metropolis(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){
    using namespace cfg;
    auto t_metropolis = tid::tic_scope(__FUNCTION__);
    double l, d_theta, d_A, rand;
    double acc_rate=0.5, acc_A=0., acc_theta=0.;
    std::array<O2, NC> NewPsi{};
    std::array<O2, NC> OldPsi{};
    double NewA, OldA;
    double newE, oldE, minus_deltaE;
    double h2=(Hp.h*Hp.h);

    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {
            size_t i = ix + Lx * (iy);
            /*choose randomly a site of the lattice*/
            //i = rnd::uniform_integer_box(0, N-1);
            //ix=i%Lx;
            //iy=i/Ly;
            if (Hp.london_app == 0) {
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
            /*************PSI UPDATE: phase update **********/
            if (Hp.phase_update != 0) {

                for (size_t alpha = 0; alpha < NC; alpha++) {
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
    }

    if (Hp.e != 0) {
    /**********VECTOR POTENTIAL UPDATE********/
        for (size_t iy = 0; iy < Ly; iy++) {
            for (size_t ix = 0; ix < Lx; ix++) {
                size_t i = ix + Lx * (iy );
                for (size_t vec = 0; vec < DIM; vec++) {
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

    acc_theta= acc_theta/static_cast<double>(NC*N);
    acc_A= acc_A/static_cast<double>(DIM*N);
    MCp.lbox_theta= MCp.lbox_theta*((0.5*acc_theta/acc_rate)+0.5);
    MCp.lbox_A= MCp.lbox_A*((0.5*acc_A/acc_rate)+0.5);

}

//slightly different phase update procedure where first I update N times the phase of the component 1 and afterward N times that of componennt 2
void metropolis2(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){
    using namespace cfg;
    using namespace rnd;
    auto t_metropolis = tid::tic_scope(__FUNCTION__);
    double d_theta, rand;
    double acc_rate=0.5, acc_theta=0.;
    std::array<O2, NC> NewPsi{};
    std::array<O2, NC> OldPsi{};
    double newE, oldE, minus_deltaE;
    double h2=(Hp.h*Hp.h);

    for (size_t alpha=0; alpha<NC; alpha++) {
        for (size_t ni = 0; ni < N; ni++) {
            /*choose randomly a site of the lattice*/
            size_t i = uniform_integer_box(0, (int)N-1);
            size_t ix=i%Lx;
            size_t iy=i/Ly;
            /*************PSI UPDATE: phase update **********/
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
                rand = uniform_double_box(0, 1);
                //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
                if (rand < exp(my_beta * minus_deltaE)) {
                    Site[i].Psi[alpha] = NewPsi[alpha];
                    acc_theta++;
                }
            }
        }
    }

    acc_theta= acc_theta/static_cast<double>(NC*N);
    MCp.lbox_theta= MCp.lbox_theta*((0.5*acc_theta/acc_rate)+0.5);

}


double local_HPsi(std::array<O2, NC> &Psi, size_t ix, size_t iy, H_parameters &Hp, const std::vector<Node> &Site) {
    using namespace cfg;
    double h_Kinetic=0., h_Josephson=0., h_lambda=0., h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double gauge_phase1, gauge_phase2;
    double gamma=0;
    size_t nn_ip, nn_im;
    size_t i=ix +Lx*(iy);

    size_t ip=(ix == Lx-1 ? 0: ix+1);
    size_t ipx= ip+Lx*(iy);
    size_t jp=(iy == Ly-1 ? 0: iy+1);
    size_t ipy= ix+(Lx*jp);
    size_t imx= (ix == 0 ? Lx-1: ix-1)+Lx*(iy);
    size_t imy= ix+Lx*((iy == 0 ? Ly-1: iy-1));

    //Compute the local Energy respect to a given component (alpha) of the matter field Psi and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving Psi
    for(size_t alpha=0; alpha<NC; alpha ++) {
        //Kinetic= -(1/h²)*\sum_k=1,2 (|Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))) + (|Psi_{alpha}(r-k)||Psi_{alpha}(r)|* cos(theta_{alpha}(r) - theta_{alpha}(r-k) +h*e*A_k(r-k)))
        for (size_t vec = 0; vec < DIM; vec++) {
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

double local_HA(double A, size_t ix, size_t iy,  size_t vec,  struct H_parameters &Hp, const std::vector<Node> &Site){
    using namespace cfg;
    auto t_local = tid::tic_scope(__FUNCTION__);
    double h_Kinetic=0., h_B, h_tot;
    double inv_h2=1./(Hp.h*Hp.h);
    double F2_A=0., F_A;

    size_t nn_ip, nn_ipl, nn_iml;
    size_t i=ix +Lx*(iy);

    size_t ip=(ix == Lx-1 ? 0: ix+1);
    size_t ipx= ip+ Lx*(iy);
    size_t jp=(iy == Ly-1 ? 0: iy+1);
    size_t ipy= ix+(Lx*jp);

    size_t imx= (ix == 0 ? Lx-1: ix-1)+Lx*iy;
    size_t imy= ix+Lx*(iy == 0 ? Ly-1: iy-1);

    if(vec==0){
        nn_ip=ipx;
    }else if(vec==1){
        nn_ip=ipy;
    } else throw std::logic_error("unexpected value for vec");

    //Compute the local Energy respect to a given component (alpha) of the vector potential A and a given spatial position (r=(ix, iy, iz))
    //We need to compute just the part of the Hamiltonian involving A

    //Kinetic= -(1/h²)*\sum_k=1,2 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos(theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
    for(size_t alpha=0; alpha<NC; alpha++) {
        double gauge_phase1=Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t +  Site[i].R_ext[vec] +  Hp.h*Hp.e*A;
        h_Kinetic -= inv_h2* Site[nn_ip].Psi[alpha].r*Site[i].Psi[alpha].r*cos(gauge_phase1);
	 }

    //All the plaquettes involving A_vec(i)
    for(size_t l=0; l<DIM; l++){
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
    for(size_t l=0; l<DIM; l++){
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


/* For the moment implemented for lambda=0. DOES NOT WORK */
void overrelaxation(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp,  double my_beta){
    using namespace cfg;
    auto t_metropolis = tid::tic_scope(__FUNCTION__);
    std::array<O2, NC> pxPsi{};
    std::array<O2, NC> pyPsi{};
    std::array<O2, NC> mxPsi{};
    std::array<O2, NC> myPsi{};
    struct O2 localField{};
    struct O2 newPsi{};
    double rnorm, coeff;

    for(size_t i=0; i<N; i++) {
        for (size_t alpha = 0; alpha < NC; alpha++) {
            polar_to_cartesian(Site[i].Psi[alpha]);
        }
    }

    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {
            size_t i = ix + Lx * (iy);

            size_t ipx=(ix == Lx-1 ? 0: ix+1);
            size_t imx=(ix == 0 ? Lx-1: ix-1);

            size_t ipy=(iy == Ly-1 ? 0: iy+1);
            size_t imy=(iy == 0 ? Ly-1: iy-1);

            for (size_t alpha = 0; alpha < NC; alpha++) {
                pxPsi[alpha].t= Site[ipx + Lx*iy].Psi[alpha].t + Site[i].R_ext[0] + Hp.h * Hp.e * Site[i].A[0];
                pxPsi[alpha].r= Site[ipx + Lx*iy].Psi[alpha].r;

                mxPsi[alpha].t= Site[imx + Lx*iy].Psi[alpha].t - Site[imx + Lx*iy].R_ext[0] - Hp.h * Hp.e * Site[imx + Lx*iy].A[0];
                mxPsi[alpha].r= Site[imx + Lx*iy].Psi[alpha].r;

                pyPsi[alpha].t= Site[ix + Lx*ipy].Psi[alpha].t + Site[i].R_ext[1] + Hp.h * Hp.e * Site[i].A[1];
                pyPsi[alpha].r= Site[ix + Lx*ipy].Psi[alpha].r;

                myPsi[alpha].t= Site[ix + Lx*imy].Psi[alpha].t - Site[ix + Lx*imy].R_ext[1] - Hp.h * Hp.e * Site[ix + Lx*imy].A[1];
                myPsi[alpha].r= Site[ix + Lx*imy].Psi[alpha].r;

                polar_to_cartesian(pxPsi[alpha]);
                polar_to_cartesian(mxPsi[alpha]);
                polar_to_cartesian(pyPsi[alpha]);
                polar_to_cartesian(myPsi[alpha]);

                localField.y = -(pxPsi[alpha].y + mxPsi[alpha].y + pyPsi[alpha].y + myPsi[alpha].y );
                localField.x = -(pxPsi[alpha].x + mxPsi[alpha].x + pyPsi[alpha].x + myPsi[alpha].x );

                for(size_t beta=0; beta<NC; beta++){
                    if(beta!= alpha){
                        localField.x += 2*Hp.K*Site[i].Psi[beta].x;
                        localField.y += 2*Hp.K*Site[i].Psi[beta].y;
                    }
                }

                if( O2norm2(localField)!= 0){

                    rnorm= 1./O2norm2(localField);
                    coeff= 2.0*rnorm*O2prod(Site[i].Psi[alpha],localField);
                    O2comb(newPsi, -1.0, Site[i].Psi[alpha], coeff, localField);

                    Site[i].Psi[alpha]=newPsi;
                    cartesian_to_polar(Site[i].Psi[alpha]);
                }
            }
        }
    }
}

