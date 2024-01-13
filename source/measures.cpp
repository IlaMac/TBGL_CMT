//
// Created by ilaria on 2019-11-16.
//

#include "measures.h"
#include "cfg/cfg.h"
#include "main.h"
#include "rnd.h"
#include "tid/tid.h"
#include <fmt/format.h>


void energy(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site){
    using namespace cfg;
    auto t_energy = tid::tic_scope(__FUNCTION__);
    double h_Kinetic=0., h_Josephson=0., h_B=0., h_lambda=0.;
    double F_A=0;
    double gauge_phase;
    double gamma=0.;
    double inv_h2=1./(Hp.h*Hp.h);
    double h2=(Hp.h*Hp.h);
    size_t nn_ip, ip;

    for(size_t iy=0; iy<Ly; iy++){
        for(size_t ix=0; ix<Lx; ix++){
            size_t i=ix + Lx * (iy );
            for(size_t alpha=0; alpha<NC; alpha++) {
                //Kinetic= -(1/h²)*\sum_k=1,2 |Psi_{alpha}(r)||Psi_{alpha}(r+k)|* cos( theta_{alpha}(r+k) - theta_{alpha}(r) +h*e*A_k(r))
                for (size_t vec = 0; vec < DIM; vec++) {
                    if (vec == 0) {
                        ip = (ix == Lx - 1 ? 0 : ix + 1);
                        nn_ip = ip + Lx * (iy);
                    } else if (vec == 1) {
                        ip = (iy == Ly - 1 ? 0 : iy + 1);
                        nn_ip = ix + Lx * ip;
                    }
                    gauge_phase = Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
                    h_Kinetic -= inv_h2 * (Site[i].Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * cos(gauge_phase);
                }
            }
            //Biquadratic Josephson= K * |Psi_{1}(r)|^2|Psi_{2}(r)|^2* (cos(2*(theta_{1}(r) - theta_{2}(r))) -1 )
            h_Josephson +=  Hp.K*(Site[i].Psi[0].r*Site[i].Psi[1].r)*(Site[i].Psi[0].r*Site[i].Psi[1].r)*(cos(2*(Site[i].Psi[0].t -Site[i].Psi[1].t)) -1.);
            // with h_lambda= \lambda*( |Psi_{1}(r)|^6 - |Psi_{2}(r)|^6 - 3*|Psi_{1}(r)|^2|Psi_2}(r)|^2*(|Psi_{1}(r)|^2 - |Psi_{2}(r)|^2)*(3 + 2*cos(2*(theta_{1}(r) - theta_{2}(r))))
            if(Hp.lambda !=0 ) {
                gamma=2* atan2( Site[i].Psi[1].r, Site[i].Psi[0].r);
                h_lambda+=Hp.lambda*(cos(3*gamma) + 3*cos(gamma)*sin(gamma)*sin(gamma)*sin((Site[i].Psi[0].t -Site[i].Psi[1].t))*sin((Site[i].Psi[0].t -Site[i].Psi[1].t)));
            }

            if (Hp.e != 0) {
                for(size_t vec1=0; vec1<DIM; vec1++){
                    for (size_t vec2 = vec1+1; vec2 < DIM; vec2++) {
                        //F_{alpha,vec}= A_alpha(r_i) + A_vec(ri+alpha) - A_alpha(r_i+vec) - A_vec(ri)
                        F_A = (Site[i].A[vec1] + Site[nn(i, vec1, 1)].A[vec2] - Site[nn(i, vec2, 1)].A[vec1] -
                                Site[i].A[vec2]);
                        h_B += ((0.5*inv_h2) * (F_A * F_A));
                    }
                }
            }
        }
    }

    //to compute the heat capacity it is important to consider the total physical energy which is h_tot*h³
    mis.E_kin=(double)h2*h_Kinetic;
    mis.E_josephson=(double)h2*h_Josephson;
    mis.E_B= (double)h2*h_B;
    mis.E_lambda= (double)h2*h_lambda;
    mis.E=(mis.E_kin  + mis.E_josephson + mis.E_B + mis.E_lambda);
}


void helicity_modulus(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site){
    using namespace cfg;
    auto t_helicity = tid::tic_scope(__FUNCTION__);
    double J_alpha=0., DJ_alpha_Dd=0.;
    size_t vec=0; //helicity modulus computed along the x direction
    double gauge_phase1;
    double inv_h=1./(Hp.h);

    for(size_t iy=0; iy<Ly;iy++){
        for(size_t ix=0; ix<Lx; ix++){
            size_t i=ix +Lx*(iy);
            size_t ip=(ix == Lx-1 ? 0: ix+1);
            size_t nn_ip=ip+Lx*(iy);
            for(size_t alpha=0; alpha<NC; alpha++){
                gauge_phase1=Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Site[i].R_ext[vec] + Hp.h*Hp.e*Site[i].A[vec];
                J_alpha=inv_h*(Site[i].Psi[alpha].r*Site[nn_ip].Psi[alpha].r)*sin(gauge_phase1);
                DJ_alpha_Dd=inv_h*(Site[i].Psi[alpha].r*Site[nn_ip].Psi[alpha].r)*cos(gauge_phase1);
                mis.DH_Ddi[alpha] += (inv_h * J_alpha);
                mis.D2H_Dd2i[alpha]+= (inv_h*DJ_alpha_Dd );
                //in this case the derivative of H with respect both the two current (mixed term) is zero
            }
        }
    }
}


void dual_stiffness(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site){
    using namespace cfg;
    auto t_dual = tid::tic_scope(__FUNCTION__);
    double qx_min=C_TWO_PI/static_cast<double>(Lx);
    double invNorm= 1./((C_TWO_PI)*(C_TWO_PI)*static_cast<double>(N));
    double Re_rhoz=0.;
    double Im_rhoz=0.;
    double Dx_Ay, Dy_Ax;
    double inv_h=1./(Hp.h);

    for(size_t ix=0; ix<Lx;ix++){
        for(size_t iy=0; iy<Ly;iy++){
            size_t i=ix +Lx*(iy);
            Dx_Ay=(Site[nn(i, 0, 1)].A[1]- Site[i].A[1])*inv_h;
            Dy_Ax=(Site[nn(i, 1, 1)].A[0]- Site[i].A[0])*inv_h;

            Re_rhoz+=(cos((qx_min*static_cast<double>(ix)))*(Dx_Ay -Dy_Ax));
            Im_rhoz+=(sin((qx_min*static_cast<double>(ix)))*(Dx_Ay -Dy_Ax));
        }
    }
    mis.d_rhoz=(double) invNorm*((Re_rhoz*Re_rhoz) +(Im_rhoz*Im_rhoz));
}

void Z2_magnetization(struct Measures &mis, const std::vector<Node> &Site){
    //The Ising parameter m(x,y)=+/-1 indicates the chirality between the two phases.
    using namespace cfg;
    auto t_Z2 = tid::tic_scope(__FUNCTION__);

    double phi_shifted=0.;
    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {
            size_t i=ix +Lx*(iy);
		    phi_shifted= Site[i].Psi[1].t - Site[i].Psi[0].t;
		    while(phi_shifted >= C_PI){
			    phi_shifted-= C_TWO_PI;
            }
            while(phi_shifted< -C_PI){
                phi_shifted+=C_TWO_PI;
            }
            if(phi_shifted>0){
                mis.z2_m+=1;
            }else if(phi_shifted<0){
                mis.z2_m+=(-1);
            }
        }
    }

    mis.z2_m= (double)mis.z2_m / static_cast<double>(N);
}

void magnetization_singlephase(struct Measures &mis, const std::vector<Node> &Site){
    using namespace cfg;
    auto t_magnetization = tid::tic_scope(__FUNCTION__);
    [[maybe_unused]] double cos_phi[NC]={0}; //TODO: REMOVE UNUSED
    [[maybe_unused]] double sin_phi[NC]={0}; //TODO: REMOVE UNUSED
    auto inv_N= 1. / static_cast<double>(N);

    for(auto & s : Site) {
        for (size_t alpha = 0; alpha < NC; alpha++) {
            mis.mx_phase[alpha] += cos(s.Psi[alpha].t);
            mis.my_phase[alpha] += sin(s.Psi[alpha].t);
        }
    }

    for(size_t alpha=0; alpha<NC; alpha++) {
        mis.mx_phase[alpha]*=inv_N;
        mis.my_phase[alpha]*=inv_N;
    }
}

void nematic_order(struct Measures &mis, const std::vector<Node> &Site){
    using namespace cfg;
    auto t_nematic = tid::tic_scope(__FUNCTION__);
    long double gamma_temp;
    long double theta_temp;
    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {
            for (size_t alpha = 0; alpha < NC; alpha++) {
                mis.density_psi[alpha]+=(Site[ix+Lx*(iy)].Psi[alpha].r*Site[ix+Lx*(iy)].Psi[alpha].r);
            }
            mis.density_diff+= ((Site[ix+Lx*(iy)].Psi[0].r * Site[ix+Lx*(iy)].Psi[0].r) - (Site[ix+Lx*(iy)].Psi[1].r * Site[ix+Lx*(iy)].Psi[1].r));
            gamma_temp=2*(atan2(Site[ix+Lx*(iy)].Psi[1].r, Site[ix+Lx*(iy)].Psi[0].r));
            while(gamma_temp > C_PI){
                gamma_temp-= C_PI;
            }
            while(gamma_temp< 0){
                gamma_temp+=C_PI;
            }
//            mis.gamma+=(double)gamma_temp;
//            mis.Mx_gamma+= (double)cos(gamma_temp);
//            mis.My_gamma+= (double)sin(gamma_temp);

            theta_temp=(Site[ix+Lx*(iy)].Psi[0].t - Site[ix+Lx*(iy)].Psi[1].t);
            while(theta_temp >= C_PI){
                theta_temp-= C_TWO_PI;
            }
            while(theta_temp< -C_PI){
                theta_temp+=C_TWO_PI;
            }
//            mis.theta12+= (double)theta_temp;
//            mis.Mx_theta12+= (double)cos(theta_temp);
//            mis.My_theta12+= (double)sin(theta_temp);

            mis.Mx_nem+= (double)sin(gamma_temp)* (double)cos(theta_temp);
            mis.My_nem+= (double)sin(gamma_temp)* (double)sin(theta_temp);
            mis.Mz_nem+= (double)cos(gamma_temp);
        }
    }

    for (double & density_psi_alpha : mis.density_psi) {
        density_psi_alpha= density_psi_alpha / static_cast<double>(N);
    }

    mis.density_diff/= static_cast<double>(N);
//    mis.gamma/= static_cast<double>(N);
//    mis.theta12/= static_cast<double>(N);

}

void vorticity(struct Measures &mis,struct H_parameters &Hp, const std::vector<Node> &Site ){
    using namespace cfg;

    auto t_vorticity = tid::tic_scope(__FUNCTION__);
    std::vector<std::pair<double, double>> coord_v1;
    std::vector<std::pair<double, double>> coord_av1;

    std::vector<std::pair<double, double>> coord_v2;
    std::vector<std::pair<double, double>> coord_av2;

    coord_v1.resize(N);
    coord_v2.resize(N);
    coord_av1.resize(N);
    coord_av2.resize(N);

    size_t nv1=0, nv2=0, nav1=0, nav2=0;

    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {
            double n=0;
            //Convention gauge phase
            //gauge_phase = Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Site[i].R_ext[vec] + Hp.h * Hp.e * Site[i].A[vec];
            //Circuitation of the phase component 0 around the i-th plaquette
            size_t ipx= (ix == Lx-1 ? 0: ix+1);
            size_t i= ix + Lx*(iy);
            size_t nn_ipx=ipx + Lx*(iy);
            double phi_1=  Site[nn_ipx].Psi[0].t -Site[i].Psi[0].t + Site[i].R_ext[0] + Hp.h * Hp.e * Site[i].A[0];
            while(phi_1< -C_PI){
                phi_1+= C_TWO_PI;
            }
            while(phi_1>= C_PI){
                phi_1-= C_TWO_PI;
            }
            size_t ipy= (iy == Ly-1 ? 0: iy+1);
            size_t nn_ipx_ipy= ipx + Lx*(ipy);
            double phi_2= Site[nn_ipx_ipy].Psi[0].t - Site[nn_ipx].Psi[0].t +  Site[nn_ipx].R_ext[1] + Hp.h * Hp.e * Site[nn_ipx].A[1];
            while(phi_2< -C_PI){
                phi_2+= C_TWO_PI;
            }
            while(phi_2>= C_PI){
                phi_2-= C_TWO_PI;
            }
            size_t nn_ipy= ix + Lx*(ipy);
            double phi_3= Site[nn_ipy].Psi[0].t - Site[nn_ipx_ipy].Psi[0].t -  Site[nn_ipy].R_ext[0] - Hp.h * Hp.e * Site[nn_ipx].A[0] ;
            while(phi_3< -C_PI){
                phi_3+= C_TWO_PI;
            }
            while(phi_3>= C_PI){
                phi_3-= C_TWO_PI;
            }
            double phi_4= Site[i].Psi[0].t - Site[nn_ipy].Psi[0].t -  Site[i].R_ext[1] - Hp.h * Hp.e * Site[i].A[1];
            while(phi_4< -C_PI){
                phi_4+= C_TWO_PI;
            }
            while(phi_4>= C_PI){
                phi_4-= C_TWO_PI;
            }
            n=(phi_1 + phi_2+ phi_3 + phi_4)/C_TWO_PI;
            if(n>0.01){
               mis.vortex_density[0]+=1;
               coord_v1[nv1]= std::pair(ix, iy);
               nv1++;
            }
            else if(n<-0.01){
                mis.antivortex_density[0]+=1;
                coord_av1[nav1]=(std::pair(ix, iy));
                nav1++;
            }

            //Circuitation of the phase component 1 around the i-th plaquette
            phi_1=  Site[nn_ipx].Psi[1].t -Site[i].Psi[1].t + Site[i].R_ext[0] + Hp.h * Hp.e * Site[i].A[0];
            while(phi_1< -C_PI){
                phi_1+= C_TWO_PI;
            }
            while(phi_1>= C_PI){
                phi_1-= C_TWO_PI;
            }
            phi_2= Site[nn_ipx_ipy].Psi[1].t - Site[nn_ipx].Psi[1].t +  Site[nn_ipx].R_ext[1] + Hp.h * Hp.e * Site[nn_ipx].A[1];
            while(phi_2< -C_PI){
                phi_2+= C_TWO_PI;
            }
            while(phi_2>= C_PI){
                phi_2-= C_TWO_PI;
            }
            phi_3= Site[nn_ipy].Psi[1].t - Site[nn_ipx_ipy].Psi[1].t -  Site[nn_ipy].R_ext[0] - Hp.h * Hp.e * Site[nn_ipx].A[0] ;
            while(phi_3< -C_PI){
                phi_3+= C_TWO_PI;
            }
            while(phi_3>= C_PI){
                phi_3-= C_TWO_PI;
            }
            phi_4= phi_4= Site[i].Psi[0].t - Site[nn_ipy].Psi[0].t -  Site[i].R_ext[1] - Hp.h * Hp.e * Site[i].A[1];
            while(phi_4< -C_PI){
                phi_4+= C_TWO_PI;
            }
            while(phi_4>= C_PI){
                phi_4-= C_TWO_PI;
            }
            n=(phi_1 + phi_2+ phi_3 + phi_4)/C_TWO_PI;
            if(n>0.01){
                mis.vortex_density[1]+=1.;
                coord_v2[nv2]=(std::pair(ix, iy));
                nv2++;
            }
            else if(n<-0.01){
                mis.antivortex_density[1]+=1.;
                coord_av2[nav2]=(std::pair(ix, iy));
                nav2++;
            }

        }
    }

    coord_v1.resize(nv1);
    coord_v2.resize(nv2);
    coord_av1.resize(nav1);
    coord_av2.resize(nav2);

    for(size_t alpha=0; alpha<NC; alpha++){
        mis.vortex_density[alpha]/=(double) N;
        mis.antivortex_density[alpha]/= (double) N;
//        std::cout << "alpha: " << alpha << "vdensity: " << mis.vortex_density[alpha]<<std::endl ;
    }

    //Minimum size of composite vortices formed by a vortex (antivortex) in the component 0 and a vortex (antivortex) in the component 1
    if(nv1==0 || nv2==0){
        mis.composite_vortex2_size=-1; //in the data analysis, I will consider only positive entries. In this way, I avoid to overcount the size "0".
        }else{
            for(size_t v1p=0; v1p<nv1; v1p++){
                double temp_dist2 ,dist2=10000;
                for(size_t v2p=0; v2p<nv2; v2p++) {
                    double dx= abs(coord_v1[v1p].first - coord_v2[v2p].first);
                    while (dx>(double)Lx/2){
                        dx=(double)Lx-dx;
                    }
                    double dy=abs(coord_v1[v1p].second - coord_v2[v2p].second);
                    while (dy>(double)Ly/2){
                        dy=(double)Ly-dy;
                    }
                    temp_dist2= dx*dx +dy*dy;
                    if (temp_dist2< dist2){
                        dist2= temp_dist2;
                    }
                    if(sqrt(dist2)>Lx){
                        std::cout<<"v1p - v2p. dx=" << dx << " dy="<< dy<< " dist= "<< sqrt(dist2)  << std::endl;
                    }
                }
                mis.composite_vortex2_size+= sqrt(dist2);
            }
            mis.composite_vortex2_size/=(double)nv1;
        }

    if (nav1 == 0 || nav2 == 0) {
        mis.composite_vortex2_size=-1; //in the data analysis, I will consider only positive entries. In this way, I avoid to overcount the size "0".
        }else{
            for(size_t v1m=0; v1m<nav1; v1m++){
                double temp_dist2 ,dist2=10000;
                for(size_t v2m=0; v2m<nav2; v2m++) {
                    double dx= abs(coord_av1[v1m].first - coord_av2[v2m].first);
                    while (dx>(double)Lx/2){
                        dx=(double)Lx-dx;
                    }
                    double dy=abs(coord_av1[v1m].second - coord_av2[v2m].second);
                    while (dy>(double)Ly/2){
                        dy=(double)Ly-dy;
                    }
                    temp_dist2= dx*dx +dy*dy;
                    if (temp_dist2< dist2){
                        dist2= temp_dist2;
                    }
                    if(sqrt(dist2)>Lx){
                        std::cout<<"v1m - v2m. dx=" << dx << " dy="<< dy<< " dist= "<< sqrt(dist2)  << std::endl;
                    }
                }
                mis.composite_vortex2_size+=sqrt(dist2);
            }
            mis.composite_vortex2_size/=(double)nav1;
        }


    //Minimum size of composite vortices formed by a vortex (antivortex) in the component 0 and an antivortex (vortex) in the component 1
        if (nv1 == 0 || nav2 == 0) {
            mis.composite_vortex1_size=-1; //in the data analysis, I will consider only positive entries. In this way, I avoid to overcount the size "0".
        }
        else{
            for(size_t v1p=0; v1p<nv1; v1p++){
                double temp_dist2, dist2=10000;
                for(size_t v2m=0; v2m<nav2; v2m++) {
                    double dx= abs(coord_v1[v1p].first - coord_av2[v2m].first);
                    while (dx>(double)Lx/2){
                        dx=(double)Lx-dx;
                    }
                    double dy= abs(coord_v1[v1p].second - coord_av2[v2m].second);
                    while (dy>(double)Ly/2){
                        dy=(double)Ly-dy;
                    }
                    temp_dist2= dx*dx +dy*dy;
                    if (temp_dist2< dist2){
                        dist2= temp_dist2;
                    }
                    if(sqrt(dist2)>Lx){
                        std::cout<<"v1p - v2m. dx=" << dx << " dy="<< dy<< " dist= "<< sqrt(dist2)  << std::endl;
                    }
                }
                mis.composite_vortex1_size+=sqrt(dist2);
            }
            mis.composite_vortex1_size/=(double) nv1;
        }

        if (nav1 == 0 || nv2 == 0) {
        mis.composite_vortex1_size = -1; //in the data analysis, I will consider only positive entries. In this way, I avoid to overcount the size "0".
        }
        else {
            for (size_t v1m = 0; v1m < nav1; v1m++) {
                double temp_dist2, dist2 = 10000;
                for (size_t v2p = 0; v2p < nv2; v2p++) {
                    double dx = abs(coord_av1[v1m].first - coord_v2[v2p].first);
                    while (dx > (double) Lx / 2) {
                        dx = (double) Lx - dx;
                    }
                    double dy = abs(coord_av1[v1m].second - coord_v2[v2p].second);
                    while (dy > (double) Ly / 2) {
                        dy = (double) Ly - dy;
                    }
                    temp_dist2 = dx * dx + dy * dy;
                    if (temp_dist2 < dist2) {
                        dist2 = temp_dist2;
                    }
                    if(sqrt(dist2)>Lx){
                        std::cout<<"v1m - v2p. dx=" << dx << " dy="<< dy<< " dist= "<< sqrt(dist2)  << std::endl;
                    }
                }
                mis.composite_vortex1_size += sqrt(dist2);
            }
            mis.composite_vortex1_size /= (double) nav1;
        }
}


    void save_lattice(const std::vector<Node> &Site, const fs::path & directory_write, const std::string & configuration){
        auto t_save = tid::tic_scope(__FUNCTION__);

        auto sPsi = fmt::format("Psi_{}.bin", configuration);
        auto sA = fmt::format("A_{}.bin",configuration);
        fs::path psi_init_file = directory_write / sPsi;
        fs::path a_init_file = directory_write / sA;

        FILE *fPsi= nullptr;
        FILE *fA= nullptr;

        if((fPsi=fopen(psi_init_file.c_str(), "w"))) {
            for (auto & s: Site) {
                fwrite(s.Psi.data(), sizeof(struct O2), NC, fPsi);
            }
            fclose(fPsi);
        }

        if((fA=fopen(a_init_file.c_str(), "w"))) {
            for (auto & s: Site) {
                fwrite(s.A.data(), sizeof(double), DIM, fA);
            }
            fclose(fA);
        }

    }

    void save_lattice_chargezero(const std::vector<Node> &Site, const fs::path & directory_write, const std::string &configuration){
        auto t_save = tid::tic_scope(__FUNCTION__);
        auto sPsi= fmt::format("Psi_{}.bin", configuration);
        fs::path psi_init_file = directory_write / sPsi;

        FILE *fPsi= nullptr;

        if((fPsi=fopen(psi_init_file.c_str(), "w"))) {
            for (auto & s: Site) {
                fwrite(s.Psi.data(), sizeof(struct O2), NC, fPsi);
            }
            fclose(fPsi);
        }

    }

