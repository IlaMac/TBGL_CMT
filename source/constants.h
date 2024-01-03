//
// Created by ilaria on 2021-01-23.
//
#pragma once
#include <cmath>

/*Number of components*/
inline constexpr size_t DIM = 2;
inline constexpr size_t NC = 2;
//extern size_t Lx, Ly, N; // unsigned long (it is an integral)

//#define C_TWO_PI (6.2831853071795864769252867665590058L)
//#define C_PI (3.1415926535897932384626433832795029L)

static constexpr auto C_PI= 3.1415926535897932384626433832795029;
static constexpr auto C_TWO_PI= 6.2831853071795864769252867665590058;
[[maybe_unused]] static constexpr auto C_PI_LD= 3.1415926535897932384626433832795029l;
[[maybe_unused]] static constexpr auto C_TWO_PI_LD= 6.2831853071795864769252867665590058l;

namespace settings {

#ifdef NDEBUG
    inline constexpr bool debug = false;
#else
    inline constexpr bool debug = true;
#endif

}
