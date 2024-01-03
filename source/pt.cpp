//
// Created by ilaria on 2021-11-05.
//

#include "pt.h"

void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp){
    int p;
    double beta_low, beta_high, delta_beta;

    if(Hp.b_high>Hp.b_low){ //Paranoic check
        beta_low=Hp.b_low;
        beta_high=Hp.b_high;
    }else{
        beta_low=Hp.b_high;
        beta_high=Hp.b_low;
    }
    PTroot.beta.resize(PTp.np, 0.0);
    PTroot.All_Energies.resize(PTp.np, 0.0);
    PTroot.ind_to_rank.resize(PTp.np, 0);
    PTroot.rank_to_ind.resize(PTp.np, 0);
    delta_beta=(beta_high-beta_low)/(PTp.np-1);
    for(p=0; p<PTp.np; p++){
        PTroot.rank_to_ind[p]=p;
        PTroot.ind_to_rank[p]=p;
        PTroot.beta[p]=beta_low + p*delta_beta;
    }
}

void parallel_temp(double &my_E , double &my_beta, int &my_ind, struct PT_parameters &PTp, struct PTroot_parameters &PTroot){

    double coin;
    double n_rand, delta_E, delta_beta;
    double oldbeta_i, oldbeta_nn;
    int i=0, nn=0, ind_nn=0;
    int oldrank_i=0, oldrank_nn=0;
    int newrank_i=0, newrank_nn=0;


    MPI_Gather(&my_E, 1, MPI_DOUBLE, PTroot.All_Energies.data(), 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    if (PTp.rank == PTp.root) { //Root forms the pairs and decides (given the energies and the betas) which pairs will swap
        //Pair Formation
        coin = rnd::uniform_double_box(0,1);
        if(coin < 0.5) { //each even rank wil be paired with its right neighbour
            nn= +1;
        }else if(coin >= 0.5){ //each even rank wil be paired with its left neighbour
            nn=-1;
        }
        while (i < PTp.np) {
            n_rand=rnd::uniform_double_box(0,1);
            ind_nn=(PTp.np + i + nn) % PTp.np;
            oldrank_i=PTroot.ind_to_rank[i];
            oldrank_nn=PTroot.ind_to_rank[ind_nn];
            delta_E = PTroot.All_Energies[oldrank_i] - PTroot.All_Energies[oldrank_nn];
            delta_beta = PTroot.beta[oldrank_i] - PTroot.beta[oldrank_nn];
            //swapping condition
            //Boltzmann weight: exp(-\beta E) E= h³ \sum_i E(i)
            if (n_rand < exp(delta_beta * delta_E)) {

                //swap indices in the rank_to_ind array
                PTroot.rank_to_ind[oldrank_i] = ind_nn;
                PTroot.rank_to_ind[oldrank_nn] = i;

                //swap indices in the ind_to_rank array
                newrank_i= oldrank_nn;
                PTroot.ind_to_rank[i]= newrank_i;
                newrank_nn=oldrank_i;
                PTroot.ind_to_rank[ind_nn] =newrank_nn;
                //swap beta
                oldbeta_i= PTroot.beta[oldrank_i];
                oldbeta_nn= PTroot.beta[oldrank_nn];
                PTroot.beta[oldrank_i] = oldbeta_nn;
                PTroot.beta[oldrank_nn] = oldbeta_i;
            }
            i+= 2;
        }
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);

}
