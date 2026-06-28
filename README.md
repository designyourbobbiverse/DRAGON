# DRAGON

**Divergence-Regulation by A-field Grid Operations & Numerics**

DRAGON is a modern C++ Eulerian hydrodynamics and magnetohydrodynamics code implementing finite-volume Godunov methods. The code features configurable Riemann solvers, CFL timestep calculation, MUSCL reconstruction, and boundary-condition handling for 1D, 2D, and 3D grids. 


## Current Features

- Dimensionally split (hydro only) and unsplit multidimensional schemes, with optional CTU support for unsplit mode.
- Constrained Transport MHD, implemented via vector potential along cell edges
- Hydrodynamic mode: Exact, HLL, HLLC, HLLE, and Roe Riemann solvers.
- MHD mode: HLL, HLLD, and HLLE Riemann solvers
- Riemann safety checks and fallback behaviour for non-physical approximate fluxes.
- First-order Godunov and MUSCL-Hancock reconstruction.
- TVD limiters: minmod, monotonized central, van Leer, superbee, and van Albada.
- Boundary conditions: outflow, gated outflow, fixed, reflective, periodic. Additional options can be added by subclassing GhostFill
- Boundary options can be combined using a BoundaryList, which can be invoked by syntax such as Reflective(X) + Periodic(Y)
- A unit test suite covering fluid elements, Riemann Solvers, Boundary Condtions, Godunov methods, Constrained Transport, and more



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
- `RIEMANN_VERIFY_FALLBACK`: if enabled, phsyicality of the Riemann solution is verified before returning. On failure, will attempt to pivot to more robust (albeit more diffusive or more expensive) solvers if needed.
- `CFL_CALCULATION`: chooses CFL speed calculation method for multidimensional problems.
- `CFL_coeff`: global CFL coefficient. The current value is `0.3`.
- `MUSCL_Hancock`: enables second-order MUSCL-Hancock reconstruction. Comment it out for a first-order Godunov scheme.
- `MUSCL_DEFAULT_LIMITER`: the TVD limiter used by MUSCL reconstruction.
- `DIMENSION_UNSPLIT`: unsplit multidimensional advancement (turn off for strang splitting).
- `core_count`: when using DistGrid, helps the root choose how many child grids to create.

## Running a Problem

TODO: Fill me in once you've reconfigured this


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

## Development Notes

- The code uses C++23 features
- The default gas index is `_gamma = 5/3` in `Constants.h`.
- Grid indices outside the physical range access ghost cells, which are filled by the configured boundary list.
- Boundary lists are composable with `+`; later boundary conditions override earlier ones on overlapping ghost cells.
- The demo driver is intentionally simple and is not yet a general runtime configuration interface.


## Roadmap

- Add HDF5 file output and simulation restarts (same format for v1.0)
- Replace hard-coded demo with proper Problem.cpp
- Add example/validation problems, plotting tools, and documentation
- Implement AMR

## Use of Generative AI

I, Bobbie Markwick, developed DRAGON as an independent C++ MHD code.  I have made some use of Generative AI tools (speicfically ChatGPT & Codex) during development. Example uses include algorithm discussion, typo-hunting, and unit test writing (I wrote some of the tests, AI wrote some). However, 100% of DRAGON itself was written manually by me.
