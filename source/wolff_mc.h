#pragma once

#include "class_tic_toc.h"
#include "constants.h"
#include "initialization.h"
#include "main.h"
#include "measures.h"
#include "o2.h"
#include "rnd.h"
#include <cstring>
#include <iostream>

void wolff_BTRS(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
void growCluster_BTRS(size_t i, size_t* clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
void wolff_nemK(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
void growCluster_nemK(size_t i, size_t alpha_up, size_t* clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
void wolff_density(const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
void growCluster_density(int i,  std::vector<int>  &clusterSpin, const std::vector<Node> &Site, struct MC_parameters &MCp, struct H_parameters &Hp, double my_beta);
