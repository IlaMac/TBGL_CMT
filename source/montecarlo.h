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
double local_HPsi(std::array<O2, NC> &Psi, size_t ix, size_t iy, struct H_parameters &Hp, const std::vector<Node> &Site);
double local_HA(double A, size_t ix, size_t iy, size_t alpha, struct H_parameters &Hp, const std::vector<Node> &Site);
