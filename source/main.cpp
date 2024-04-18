#include "main.h"
#include "cfg/cfg.h"
#include "cli/cli.h"
#include "initialization.h"
#include "measures.h"
#include "memory_check.h"
#include "rnd.h"
#include "tid/tid.h"
#include <csignal>
#include <cstdlib>
#include <h5pp/h5pp.h>
#include <iostream>
#include <string>

void clean_up() {
    using namespace cfg;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string src_file= h5pp::format("{}/beta_{}/Output.h5",paths_dir::directory_parameters_temp,  rank);
    std::string tgt_file= h5pp::format("{}/beta_{}/Output.h5",paths_dir::directory_parameters,  rank);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    if(src_file == tgt_file) return;

    fs::copy(paths_dir::directory_parameters_temp, paths_dir::directory_parameters , fs::copy_options::overwrite_existing | fs::copy_options::recursive );
    h5pp::hdf5::moveFile(src_file, tgt_file, h5pp::FilePermission::REPLACE);
    std::cout<<"Exit"<<std::endl;
}


void signal_callback_handler(int signum) {
    switch(signum) {
        case SIGTERM: {
            std::cout << "Caught SIGTERM" << std::endl;
            break;
        }
        case SIGKILL: {
            std::cout << "Caught SIGKILL" << std::endl;
            break;
        }
        case SIGINT: {
            std::cout << "Caught SIGINT" << std::endl;
            break;
        }
        case SIGHUP: {
            std::cout << "Caught SIGHUP" << std::endl;
            break;
        }
        case SIGQUIT: {
            std::cout << "Caught SIGQUIT" << std::endl;
            break;
        }
        default: break;
    }
    std::cout << "Exiting" << std::endl << std::flush;
    #if defined(__APPLE__) //this function is not present on apple clang
        clean_up();
        std::exit(signum);
    #else
        std::quick_exit(signum);
    #endif
}


int main(int argc, char *argv[]){
    std::vector <Node> Lattice;
    struct H_parameters Hp{};
    struct MC_parameters MCp{};
    struct PT_parameters PTp;
    struct PTroot_parameters PTroot;
    double my_beta=0.244;
    int my_ind=0;
    int NSTART=0;
    using namespace cfg;

    cli::parse(argc, argv);

    //Safe exit
    // Register termination codes and what to do in those cases
    // Basically, we just pass the termination code such as SIGKILL to the callback handler which in turn gives it to quick_exit, for instance, std::quick_exit(SIGKILL)
    signal(SIGTERM, signal_callback_handler);
    signal(SIGINT, signal_callback_handler);
    signal(SIGKILL, signal_callback_handler);
    signal(SIGHUP, signal_callback_handler);
    signal(SIGQUIT, signal_callback_handler);

    // std::at_quick_exit is called by "std::quick_exit(int)".
    // Note that std::quick_exit does not by itself catch termination codes. We have to do it ourselves with signal(), which is found in #include<csignal>
    #if !defined(__APPLE__) //this function is not present on apple clang
    std::at_quick_exit(clean_up);
    // std::atexit is called when program terminates
    std::atexit(clean_up);
    #endif

    //initialization of the random number generator
    rnd::seed(seednumber);

    //Declaration of structure Lattice
    Lattice.resize(N);

    //Initialize H_parameters: file "H_init.txt"
    initialize_Hparameters(Hp);
    //Initialize MC_parameters: file "MC_init.txt"
    initialize_MCparameters(MCp);

    MPI_Init(nullptr, nullptr); /* START MPI */
/*DETERMINE RANK OF THIS PROCESSOR*/
    MPI_Comm_rank(MPI_COMM_WORLD, &PTp.rank);
/*DETERMINE TOTAL NUMBER OF PROCESSORS*/
    MPI_Comm_size(MPI_COMM_WORLD, &PTp.np);

    auto t_tot = tid::tic_scope(fmt::format("TBGL", PTp.rank));

    if(PTp.rank == PTp.root) {
        //Initialization ranks arrays
        initialize_PTarrays( PTp, PTroot, Hp);
    }
    MPI_Scatter(PTroot.beta.data(), 1, MPI_DOUBLE, &my_beta, 1, MPI_DOUBLE, PTp.root, MPI_COMM_WORLD);
    MPI_Scatter(PTroot.rank_to_ind.data(), 1, MPI_INT, &my_ind, 1, MPI_INT, PTp.root, MPI_COMM_WORLD);

    printf("I'm rank %d and this is my beta %lf\n", PTp.rank, my_beta);
    auto directory_read = fmt::format("{}/beta_{}",cfg::paths_dir::directory_parameters, my_ind);
    initialize_lattice(Lattice, directory_read, Hp);

    if(RESTART==1){
        std::fstream restart_file(directory_read+"/restart-0", std::ios::in);
        restart_file >> NSTART;
        std::cout << NSTART << std::endl;
        restart_file.close();
    }

    //Mainloop
    mainloop(Lattice, MCp, Hp, my_beta, my_ind, PTp, PTroot, NSTART);

    t_tot.toc();

    std::cout << "Proccess current resident ram usage: " << process_memory_in_mb("VmRSS") << " MB" << std::endl;
    std::cout << "Proccess maximum resident ram usage: " << process_memory_in_mb("VmHWM") << " MB" << std::endl;
    std::cout << "Proccess maximum virtual  ram usage: " << process_memory_in_mb("VmPeak") << " MB" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if(PTp.rank == PTp.root)
        for(const auto &t : tid::get_tree("TBGL", tid::level::normal)) fmt::print("{}\n", t.str());
    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}

void mainloop(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double &my_beta, int &my_ind, struct PT_parameters PTp, struct PTroot_parameters PTroot, int NSTART) {
    using namespace cfg;

    /*Measurements*/
    Measures mis;
    /*Auxiliary Measure struct for vortices*/
    std::vector <Vdensity> local_vort_density;
    local_vort_density.resize(N);


    std::cout<< "check" << local_vort_density[2].v1[0] << std::endl;

    auto directory_write_temp = fmt::format("{}/beta_{}", cfg::paths_dir::directory_parameters_temp , my_ind);
    h5pp::File file;

    // Initialize a file
    if(NSTART==0) {
        file=h5pp::File(directory_write_temp + "/Output.h5", h5pp::FilePermission::REPLACE);
    }
    // Initialize a file in append mode
    if(NSTART>0){
        std::cout <<"NSTART >0"<< std::endl;
        file=h5pp::File(directory_write_temp+"/Output.h5", h5pp::FilePermission::READWRITE);
    }

    std::cout << directory_write_temp << "\t" << NSTART << std::endl;
    // Enable compression
    file.setCompressionLevel(0);
    // Register the compound type
    std::array<hsize_t, 1> rho_dims = {NC};
    h5pp::hid::h5t HDF5_RHO_TYPE = H5Tarray_create(H5T_NATIVE_DOUBLE, rho_dims.size(), rho_dims.data());

    h5pp::hid::h5t MY_HDF5_MEASURES_TYPE = H5Tcreate(H5T_COMPOUND, sizeof(Measures));

    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E", HOFFSET(Measures, E), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_josephson", HOFFSET(Measures, E_josephson), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_lambda", HOFFSET(Measures, E_lambda), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_kin", HOFFSET(Measures, E_kin), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "E_B", HOFFSET(Measures, E_B), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "m", HOFFSET(Measures, z2_m), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "mx_phase", HOFFSET(Measures, mx_phase),  HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "my_phase", HOFFSET(Measures, my_phase),  HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "ds", HOFFSET(Measures, d_rhoz), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "rho", HOFFSET(Measures, density_psi), HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "rho_diff", HOFFSET(Measures, density_diff), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "Mx_nem", HOFFSET(Measures, Mx_nem), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "My_nem", HOFFSET(Measures, My_nem), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "Mz_nem", HOFFSET(Measures, Mz_nem), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "Mx_gamma", HOFFSET(Measures, Mx_gamma), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "My_gamma", HOFFSET(Measures, My_gamma), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "Mx_theta12", HOFFSET(Measures, Mx_theta12), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "My_theta12", HOFFSET(Measures, My_theta12), H5T_NATIVE_DOUBLE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "vortex_density", HOFFSET(Measures, vortex_density), HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "antivortex_density", HOFFSET(Measures, antivortex_density), HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "v1v2_density", HOFFSET(Measures, v1v2_density), HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "v1av2_density", HOFFSET(Measures, v1av2_density), HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "DH_Ddi", HOFFSET(Measures, DH_Ddi), HDF5_RHO_TYPE);
    H5Tinsert(MY_HDF5_MEASURES_TYPE, "D2H_Dd2i", HOFFSET(Measures, D2H_Dd2i), HDF5_RHO_TYPE);
    //H5Tinsert(MY_HDF5_MEASURES_TYPE, "D2H_Dd2ij", HOFFSET(Measures, D2H_Dd2ij), HDF5_RHO_TYPE);

    H5Tinsert(MY_HDF5_MEASURES_TYPE, "rank", HOFFSET(Measures, my_rank), H5T_NATIVE_INT);

    file.createTable(MY_HDF5_MEASURES_TYPE, "Measurements", "Measures");

    //Initial configuration
    save_lattice(Site, directory_write_temp, std::string("init"));

    if(NSTART==0){
        /**Thermalization**/
        for (int t = 0; t < MCp.transient; t++) {
            metropolis(Site, MCp, Hp, my_beta);
            metropolis2(Site, MCp, Hp, my_beta);

        }
    }
    for (int nM = NSTART; nM<MCp.nmisu; nM++) {
        for (int t = 0; t < MCp.tau; t++) {

            metropolis(Site, MCp, Hp, my_beta);
//            metropolis2(Site, MCp, Hp, my_beta);

            (Hp.K>4) ? wolff_BTRS(Site, MCp, Hp, my_beta) : void() ;
            (Hp.K<-4) ? wolff_nemK(Site, MCp, Hp, my_beta) : void();
        }
        {
            //Measure
            auto t_measure = tid::tic_scope("measure");
            mis.reset();
            energy(mis, Hp, Site);
            Z2_magnetization(mis, Site);
            magnetization_singlephase(mis,  Site);
            new_vorticity(mis, local_vort_density, Hp, Site);
            (Hp.e == 0) ? helicity_modulus(mis, Hp, Site) : dual_stiffness(mis, Hp, Site);
            (Hp.K < 0) ? nematic_order(mis, Site) : void();

            mis.my_rank=PTp.rank;
        }

        {
            auto t_h5pp = tid::tic_token("h5pp");
            file.appendTableRecords(mis, "Measurements");
        }
        MPI_Barrier(MPI_COMM_WORLD);

        std::ofstream restart_file(directory_write_temp+"/restart-0");
        restart_file << nM <<std::endl;
        restart_file.close();

        //Save a configuration for the restarting
        save_lattice(Site, directory_write_temp, std::string("restart"));
        if(((nM+1)%MCp.freq_autosave)==0){
            if(Hp.e!=0){
                save_lattice(Site, directory_write_temp, std::string("n")+std::to_string(nM));
            }else{
                save_lattice_chargezero(Site, directory_write_temp, std::string("n")+std::to_string(nM));
            }
            save_vortexlattice(local_vort_density,directory_write_temp, std::string("n")+std::to_string(nM));
            }

        //Parallel Tempering swap
        parallel_temp(mis.E, my_beta, my_ind, PTp, PTroot);
        //Files and directory
        directory_write_temp = fmt::format("{}/beta_{}", cfg::paths_dir::directory_parameters_temp, my_ind);
        file = h5pp::File(directory_write_temp+"/Output.h5", h5pp::FilePermission::READWRITE);
    }
    save_lattice(Site, directory_write_temp, std::string("final"));
    save_vortexlattice(local_vort_density,directory_write_temp, std::string("final"));
    MPI_Barrier(MPI_COMM_WORLD);
}


size_t nn(size_t i, size_t coord, int dir){
    using namespace cfg;
    size_t ix=i%Lx;
    size_t iy=(i/Lx)%Ly;

    if(coord==0){
        int ix_new= static_cast<int>(static_cast<int>(ix) + dir / sqrt(dir * dir));
        if(ix_new==static_cast<int>(Lx)){ ix_new=0;}
        if(ix_new < 0){ ix_new=static_cast<int>(Lx-1);}
        int iy_new=static_cast<int>(iy);
        return (static_cast<size_t>(ix_new) + Lx * (static_cast<size_t>(iy_new)));

    }
    if(coord==1){
        int iy_new= static_cast<int>(static_cast<int>(iy) + dir/sqrt(dir*dir));
        if(iy_new==static_cast<int>(Ly)){ iy_new=0;}
        if(iy_new<0){ iy_new=static_cast<int>(Ly-1);}
        int ix_new=static_cast<int>(ix);
        return (static_cast<size_t>(ix_new) + Lx * (static_cast<size_t>(iy_new)));
    }
    return 1;
}
