# TBGL_CMT

**TBGL_CMT** is a C++ simulation code for a two-component (two-band) lattice Ginzburg–Landau model of superconductivity. It studies coupled complex order parameters interacting with a fluctuating gauge field; the code can also use the London approximation, in which the order-parameter amplitudes are fixed. The model includes inter-component coupling, electromagnetic coupling, an external field/flux, and optional phase updates.

The simulation uses Monte Carlo updates, including Metropolis and Wolff-based updates where applicable, and MPI parallel tempering to sample across inverse temperatures. It records thermodynamic, ordering, stiffness, and vortex observables to HDF5.

## Model and observables

Each lattice site stores two complex order parameters and the fluctuating vector potential. Model parameters include the inter-component coupling `K`, penetration-depth-related parameter `lambda`, gauge coupling `e`, external-field parameter `h`, and optional external flux (`fx`, `fy`). The `london_app` setting selects between the Ginzburg–Landau and London treatments; `phase_update` controls whether phases are updated.

Measurements include energy contributions, the `Z₂` magnetization associated with the relative phase chirality, phase-sum magnetization, dual stiffness/helicity response, nematic order, and component-resolved vortex, antivortex, and paired-vortex densities. Results and saved configurations are written to HDF5/output files using the `h5pp` library.

## Requirements

- C++17 compiler
- CMake 3.24 or newer
- Ninja
- Conan 2.x for the bundled `*-conan` presets
- MPI implementation (enabled by the provided presets)

The Conan requirements are listed in `conanfile.txt`: h5pp, Eigen, spdlog, fmt, PCG, and CLI11. OpenMP is available as a CMake option; the supplied presets enable MPI and OpenMP.

## Build

```bash
pip install conan
conan profile detect

git clone https://github.com/IlaMac/TBGL_CMT.git
cd TBGL_CMT
cmake --list-presets
cmake --preset=release-conan
cmake --build --preset=release-conan
```

This builds the `CMT` executable at `build/release-conan/CMT`. Other compiler and debug presets are available; see `CMakePresets.json` and `cmake --list-presets` for the presets supported by your environment.

## Run

The command-line interface requires an output directory, a temporary directory, and the linear system size `L`:

```bash
mpirun -np <process-count> build/release-conan/CMT \
  --outdir <output-directory> \
  --tempdir <temporary-directory> \
  --system-size <L> \
  --seed <integer> \
  --restart <0|1|2>
```

`--system-size L` creates an `L × L` lattice. `--seed` is optional, and restart mode defaults to `0` (valid values: 0–2). The executable writes measurements and configuration/restart data to the specified directories. The `launch_scripts/` directory contains examples for running on a cluster scheduler.

## Tests

The supplied presets enable CTest. After building, run:

```bash
ctest --test-dir build/release-conan --output-on-failure
```

## Repository layout

- `source/` — lattice model, initialization, Monte Carlo updates, parallel tempering, measurements, and CLI
- `tests/` — CTest test programs
- `launch_scripts/` — example cluster job and run scripts
- `cmake/` — CMake support modules

## License

See [`LICENSE`](LICENSE) for the MIT License.
