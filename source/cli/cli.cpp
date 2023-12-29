//
// Created by david on 2021-10-12.
//

#include "cfg/cfg.h"
#include "cli.h"
#include <CLI/CLI.hpp>

// MWE: https://godbolt.org/z/jddxod53d
int cli::parse(int argc, char **argv) {
    CLI::App app;
    app.description("TBGL");
    app.get_formatter()->column_width(90);
    app.option_defaults()->always_capture_default();
    app.allow_extras(false);
    /* clang-format off */
    app.add_option("--outdir"             , cfg::paths_dir::directory_parameters         , "Output directory")->required();
    app.add_option("--tempdir"            , cfg::paths_dir::directory_parameters_temp    , "Temporary output directory")->required();
    app.add_option("-s,--seed"            , cfg::seednumber                              , "Positive number seeds the random number generator");
    app.add_option("-r,--restart"         , cfg::RESTART                                 , "Restarting mode")->check(CLI::Range(0,2));
    app.add_option("-L,--system-size"     , cfg::Lx                                      , "System size L (N = L * L)")->required();
    /* clang-format on */
    app.parse(argc, argv);
    cfg::Ly = cfg::Lx;
    cfg::N = cfg::Lx * cfg::Ly;

    return 0;
}