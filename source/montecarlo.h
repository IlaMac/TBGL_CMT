#pragma once

#include "constants.h"
#include "main.h"
#include "initialization.h"
#include "rng.h"
#include "class_tic_toc.h"
#include <iostream>
#include <cstring>
#include "o2.h"
#include "measures.h"

void metropolis(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
double local_HPsi(std::array<O2, NC> &Psi, int ix, int iy, struct H_parameters &Hp, const std::vector<Node> &Site);
double local_HA(double A, int ix, int iy, int alpha,  struct H_parameters &Hp, const std::vector<Node> &Site);
