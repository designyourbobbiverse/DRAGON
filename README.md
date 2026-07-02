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
  main/             Program Main
  MHD/              Constrained Transport
  Refinement/       Partitioning Grids into parallelisable subgrids (will eventually evolve into AMR)
  
  Config.h          Solver, reconstruction, CFL, MHD, and threading configuration
  Constants.h       Physical constants
  Problem.cpp       Example problem initialization and output handling

DRAGONWING/       Multithreading
Testing/          Unit test suite
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

## Running a Problem

Along with `Config.h`, `Problem.cpp` contains the information needed to set up the problem you wish to run. Until HDF5 output is implemented, this file also includes the code to write selected simulation data to CSV files. 

To run a problem, edit `Config.h` to choose the numerical options and edit `Problem.cpp` to define the initial conditions, boundary conditions, and desired output.

In future versions, several example problem files will also be provided as templates for common validation problems.

## Running Tests

If you are using Xcode, build and run `DRAGON_TESTS`.
If you aren't using Xcode, you'll want to do the following
- Exclude DRAGON/main/main.cpp and DRAGON/Problem.cpp from your build 
- Add all files in Testing to your build
- #define TESTMODE

In either case, the test runner calls the verification routines and prints progress to standard output. A successful run ends with:

```text
All tests passed.
```


## License

DRAGON is licensed under the Apache License 2.0. See `LICENSE` for details.


## Status and Roadmap

DRAGON is an actively developed research code. The current version is suitable for demonstrating the numerical framework, solver implementations, boundary-condition infrastructure, and unit-test coverage. 

The following features are planned for v1.0:
- HDF5 file output and simulation restarts
- A cleaner public interface for a production problem setup.
- Example/validation problems
- Plotting tools
- Documentation

Additional features planned for future versions include:
- Adaptive Mesh Refinement
- Methods for suppressing the carbuncle instability
- Source Terms
- Passive Scalars
- Resistivity and Extended MHD

## Use of Generative AI

I, Bobbie Markwick, developed DRAGON as an independent C++ MHD code.  During development I have made use of Generative AI tools for taks such as learning modern C++ features, exploring established numerical algorithms, typo-hunting, and drafting unit tests. All produciton code for DRAGON itself was written manually by me.
