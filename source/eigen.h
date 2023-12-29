
#ifndef MAIN_H
#define MAIN_H
#include "initialization.h"
#include "montecarlo.h"
#include "o2.h"
#include "rnd.h"
#include "robust_filesystem.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <random>

#define NC 3 /*Number of components*/

//#define Lx 8
//#define Ly 8
//#define Lz 8
//
//#define N (Lx*Ly*Lz)

extern unsigned int Lx, Ly, Lz, N;

void mainloop(struct Node* Site, struct NN_Node* NN_Site,  struct MC_parameters &MCp, struct H_parameters &Hp, double &my_beta, int &my_ind, struct PT_parameters PTp, struct PTroot_parameters PTroot, std::string directory_parameters);
void parallel_temp(double &my_E , double &my_beta,  int &my_ind, struct PT_parameters &PTp, struct PTroot_parameters &PTroot);
unsigned int nn(unsigned int i, unsigned int coord, int dir);
void myhelp(int argd, char** argu);
#endif
