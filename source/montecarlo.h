#pragma once

#include "constants.h"
#include "initialization.h"
#include "main.h"
#include "measures.h"
#include "o2.h"
#include "rnd.h"
#include <cstring>
#include <iostream>

void   metropolis(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
double local_HPsi(std::array<O2, NC> &Psi, int ix, int iy, struct H_parameters &Hp, const std::vector<Node> &Site);
double local_HA(double A, int ix, int iy, int alpha, struct H_parameters &Hp, const std::vector<Node> &Site);
