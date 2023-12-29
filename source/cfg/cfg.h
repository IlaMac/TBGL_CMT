#pragma once
#include <string>

namespace cfg {

    inline size_t   Lx, Ly, N;
    inline long int seednumber = -1; /*by default, it is a negative number which means that rng will use random_device*/
    inline int      RESTART    = 0;
    namespace paths_dir {
        inline std::string directory_parameters;
        inline std::string directory_parameters_temp;

    }
}
