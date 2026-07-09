# DRAGON

**Divergence-Regulation by A-field Grid Operations & Numerics**

DRAGON is a modern C++ Eulerian hydrodynamics and magnetohydrodynamics code implementing finite-volume Godunov methods. The code features configurable Riemann solvers, CFL timestep calculation, MUSCL reconstruction, and boundary-condition handling for 1D, 2D, and 3D grids.  Currenlty, DRAGON is written in C++23.


## Current Features

- Dimensionally split hydrodynamics and unsplit multidimensional hydrodynamic/MHD schemes, with optional CTU support for unsplit mode.
- Constrained transport MHD, implemented using edge-centered vector potentials
- Hydrodynamic Riemann solvers: exact, HLL, HLLC, HLLE, and Roe.
- MHD Riemann solvers: HLL, HLLD, and HLLE.
- Riemann safety checks and fallback behaviour for non-physical approximate fluxes.
- First-order Godunov and MUSCL-Hancock reconstruction.
- TVD limiters: minmod, monotonized central, van Leer, superbee, and van Albada.
- Boundary conditions: outflow, gated outflow, fixed, reflective, periodic. Additional options can be added by subclassing GhostFill.
- Boundary options can be combined using a BoundaryList, which can be invoked by syntax such as 
`Boundary::Reflective("X") + Boundary::Periodic("YZ")`
- A unit test suite covering fluid elements, Riemann Solvers, boundary conditions, Godunov methods, constrained transport, and more


## Project Layout

```text
DRAGON/
  Boundary/         Boundary-conditions
    BoundaryTypes/  Pre-defined options (Reflective, Periodic, Outflow, etc)
  FluidElement/     Primitive/conservative states, vectors, fluxes, arithmetic
  Hydro/            Godunov Scheme Components
    Grid.hpp        Grid structure which holds fluid elements
    CFL/            CFL Timestep Constraint
    Godunov/        Sweep logic, split/unsplit updates
    Riemann/        Riemann solvers
    TVD/            MUSCL Reconstruction + limiter options
  IO/               HDF5 output + simulation restart
  main/             Program Main
  MHD/              Constrained Transport
  Refinement/       Partitioning Grids into parallelisable subgrids (will eventually evolve into AMR)
  
  Config.h          Solver, reconstruction, CFL, MHD, and threading configuration
  Constants.h       Physical constants
  Problem.cpp       Example problem initialization and output handling

DRAGONWING/         Multithreading
DRAGONGAZE/         Plotting Tools
Examples/           Example Problems
Testing/            Unit test suite
```

## Configuration

Most numerical choices are compile-time settings in `DRAGON/Config.h`.

Important switches include:

- `MHD`: enables magnetohydrodynamics
- `RIEMANN_DEFAULT`: the default Riemann solver. The default `RIEMANN_HLLX` uses HLLC for pure hydrodynamics and HLLD for MHD.
- `RIEMANN_VERIFY_FALLBACK`: if enabled, physicality of the Riemann solution is verified before returning. On failure, will attempt to pivot to more robust (albeit more diffusive or more expensive) solvers if needed.
- `CFL_CALCULATION`: chooses CFL speed calculation method for multidimensional problems.
- `CFL_coeff`: global CFL coefficient. The current value is `0.3`.
- `MUSCL_Hancock`: enables second-order MUSCL-Hancock reconstruction. Comment it out for a first-order Godunov scheme.
- `MUSCL_DEFAULT_LIMITER`: the TVD limiter used by MUSCL reconstruction.
- `DIMENSION_UNSPLIT`: unsplit multidimensional advancement (turn off for Strang splitting).
- `core_count`: when using DistGrid, helps the root choose how many child grids to create.

Physical parameters meanwhile can be found in `DRAGON/Constants.h`. Currently, the only such value is the specific heat ratio `_gamma`, which by default is set to 5/3.

## Setting up your Problem

Along with `Config.h`, `Problem.cpp` contains the information needed to set up the problem you wish to run. Until HDF5 output is implemented, this file also includes the code to write selected simulation data to CSV files. 

To run a problem, edit `Config.h` to choose the numerical options and edit `Problem.cpp` to define the initial conditions, boundary conditions, and desired output.

## Running an Example Problem


The Examples folder includes several example problems. To run an exmaple problem, replace Config.h and Problem.cpp with the desired example files. 

Some problems may include additional files. For example, CollidingFlows includes a Jet class which demonstrates how to use the Boundary API to implement mixed boundary conditions on the same face. Copy these alongside the Config.h and Problem.cpp files.

## Building DRAGON on macOS

DRAGON uses HDF5 for file output, so HDF5 must be installed before building and running the code.
On macOS, HDF5 can be installed with Homebrew:
```bash
brew install hdf5
brew --prefix hdf5
```
The second command prints the path to the HDF5 installation on your computer.  This is usually `/opt/homebrew/opt/hdf5` on Apple Silicon Macs or `/usr/local/opt/hdf5` on Intel Macs.

If you are using Xcode, open the project's build settings and add (using the path you got from the last step in place of `<Path_To_HDF5>`)
- `<Path_To_HDF5>/include` to *Header Search Paths*
- `<Path_To_HDF5>/lib` to *Library Search Paths*
- `<Path_To_HDF5>/lib` to *Runpath Search Paths*. 
- `-lhdf5_cpp` to *Other Linker Flags*
- `-lhdf5` to *Other Linker Flags*

Once you do this, you should be able to build and run the DRAGON target, or the DRAGON_TESTS target for the unit test suite.

If the project builds but fails when running, go to *Signing & Capabilities* for each executable target and ensure that `Disable Library Validation` is enabled under *Hardened Runtime*. If this was disabled, clean and build before trying to run the code again.


If you aren't using Xcode, you can build DRAGON from command line. From the project's root directory, run the following:
```bash
HDF5_PATH=$(brew --prefix hdf5)
```
```bash
clang++ -std=c++23 -O3 \
    -I"$HDF5_PATH/include" \
    -L"$HDF5_PATH/lib" \
    -Wl,-rpath,"$HDF5_PATH/lib" \
    -o DRAGON \
    $(find DRAGON DRAGONWING -name "*.cpp") \
    -lhdf5_cpp -lhdf5
```
The test suite can be built similarly as follows:
```bash
HDF5_PATH=$(brew --prefix hdf5)
```
```bash
clang++ -std=c++23 -O3 \
    -DTESTMODE \
    -I"$HDF5_PATH/include" \
    -L"$HDF5_PATH/lib" \
    -Wl,-rpath,"$HDF5_PATH/lib" \
    -o DRAGON_TESTS \
    $(find DRAGON DRAGONWING Testing -name "*.cpp") \
    -lhdf5_cpp -lhdf5
```

## Building DRAGON on Linux

On Linux, use your system package manager to install HDF5. For example, on Debian/Ubuntu:
```bash
sudo apt install libhdf5-dev
```

You can compile with your normal C++ compiler and manually provide the HDF5 include and library paths. If HDF5 is installed in a standard system location, this may be sufficient:
```bash
g++ -std=c++23 -O3 \
    -o DRAGON \
    $(find DRAGON DRAGONWING -name "*.cpp") \
    -lhdf5_cpp -lhdf5
```

The test suite can be built similarly
```bash
g++ -std=c++23 -O3 \
    -DTESTMODE \
    -o DRAGON_TESTS \
    $(find DRAGON DRAGONWING Testing -name "*.cpp") \
    -lhdf5_cpp -lhdf5
```

If your system provides the HDF5 compiler wrapper, you can use `h5c++` instead of manually specifying the linker flags:
```bash
h5c++ -std=c++23 -O3 \
    -o DRAGON \
    $(find DRAGON DRAGONWING -name "*.cpp")
```
```bash
h5c++ -std=c++23 -O3 \
    -DTESTMODE \
    -o DRAGON_TESTS \
    $(find DRAGON DRAGONWING Testing -name "*.cpp")
```

If the compiler cannot find `H5Cpp.h`, add the appropriate include path with `-I`. 

If the linker cannot find the HDF5 libraries, add the appropriate library path with `-L`.


You can check whether HDF5 is available with:
```bash
h5c++ -show
```
or:
```bash
h5ls --version
```
## Making Plots

DRAGON writes simulation output as HDF5 files. The `DRAGONGAZE/` directory contains Python plotting scripts for turning those files into PNG images. Alternatively, you can import your hdf5 files into visualization software of your choice.

To set up the Python dependencies:
```bash
pip install numpy matplotlib h5py
```
For macOS users, you may need to use `pip3` instead of `pip`.


Before making your plots, edit `DRAGONGAZE/Config.py` with the following:
- `hdf_dir` should match `output_dir` from `Config.h`
- `h5_base_filename` should match `output_base_name` from `Config.h`
- `img_dir` controls where the plots are written
- Any other options you wish to edit. These typically the contents of the plots

After running DRAGON, navigate to the DRAGONGAZE directory and run the relevant python command (for macOS users, if you replaced `pip` with `pip3` earlier, replace `python` with `python3` below)

For 1D hydrodynamics:
```bash
python GAZE1D.py
```
For 1D MHD:
```bash
python GAZE1DMHD.py
```
For 2D problems, pass one or more field keys:
```bash
python GAZE2D rho p E By
```
Valid 2D keys are rho, vx, vy, vz, Bx, By, Bz, p, and E.

For 3D problems, pass one or more field keys, and optionally specify which pair(s) of axes you want ploted
```bash
python GAZE3D rho E-xz E-xy By-xz
```
Valid 3D keys are rho, vx, vy, vz, Bx, By, Bz, p, and E (same as 2D).  Valid axis options are xy, xz, and yz. A key may be specified more than once to specify multiple axis options, or the axis option can be omitted to plot all 3.

Currently, 3D only supports plotting midplane slices.


## Status and Roadmap

DRAGON is an actively developed research code. The current version is suitable for demonstrating the numerical framework, solver implementations, boundary-condition infrastructure, and unit-test coverage. Basic plotting tools are also included. 

The following features are in development for v1.0:
- Example/validation problems
- Documentation

Additional features planned for future versions include:
- Adaptive Mesh Refinement
- Methods for suppressing the carbuncle instability and other grid artifacts
- Source Terms
- Passive Scalars
- Resistivity and Extended MHD


## License

DRAGON is licensed under the Apache License 2.0. See `LICENSE` for details.


## Use of Generative AI

I, Bobbie Markwick, developed DRAGON as an independent C++ MHD code.  During development I have made use of Generative AI tools for taks such as learning modern C++ features, exploring established numerical algorithms, typo-hunting, and drafting unit tests. All produciton code for DRAGON itself was written manually by me.
