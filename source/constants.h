//
// Created by ilaria on 2021-01-23.
//
#pragma once
#include <cmath>

/*Number of components*/
inline constexpr int DIM = 2;
inline constexpr int NC = 2;
extern size_t Lx, Ly, N; // unsigned long (it is an integral)

//#define C_TWO_PI (6.2831853071795864769252867665590058L)
//#define C_PI (3.1415926535897932384626433832795029L)

static constexpr auto C_PI= 3.1415926535897932384626433832795029l;
static constexpr auto C_TWO_PI= 6.2831853071795864769252867665590058l;

//        std::acos(-1.0l);

namespace settings {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif

}
