
#pragma once
#include "constants.h"
#include "initialization.h"
#include "montecarlo.h"
#include "o2.h"
#include "pt.h"
#include "rnd.h"
#include "robust_filesystem.h"
#include "wolff_mc.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <random>

namespace paths_dir{
    inline std::string TEMP_DIROUT;
    inline std::string DIROUT;
}

void mainloop(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double &my_beta, int &my_ind, struct PT_parameters PTp, struct PTroot_parameters PTroot, std::string directory_parameters, int NSTART);
size_t nn(size_t i, size_t coord, int dir);
void myhelp(int argd, char** argu);


