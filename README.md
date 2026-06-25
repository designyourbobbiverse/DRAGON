# DRAGON

**Divergence-Regulation by A-field Grid Operations & Numerics**

DRAGON is a C++ Eulerian hydrodynamics and magnetohydrodynamics code for finite-volume simulations on structured grids. The core solver uses Godunov methods with configurable Riemann solvers, CFL timestep calculation, TVD/MUSCL reconstruction, and boundary-condition handling for 1D, 2D, and 3D grids.

Adaptive mesh refinement is planned for future development. The current `AMRGrid` types wrap the base grid types and support child-grid synchronization, but actual AMR behavior has not yet been implmented.

## Current Features

- Primitive and conservative fluid-state types with ideal-gas equation-of-state helpers.
- Hydrodynamic and optional MHD state support controlled by `Config.h`.
- Exact, HLL, HLLC, HLLD, HLLE, and Roe Riemann solvers, depending on hydro/MHD mode.
- Riemann safety checks and fallback behaviour for non-physical approximate fluxes.
- First-order Godunov and MUSCL-Hancock reconstruction.
- TVD limiters: minmod, monotonized central, van Leer, superbee, and van Albada.
- CFL-constrained timestep calculation for 1D, 2D, and 3D grids.
- Dimensionally split and unsplit multidimensional updates, with optional CTU support for unsplit mode.
- Boundary conditions: outflow, gated outflow, fixed, reflective, periodic.
- Unit-style test executable covering fluid elements, Riemann solvers, grids, boundaries, TVD, CFL, and Godunov updates.

## Project Layout

```text
DRAGON/
  Boundary/         Boundary-condition implementations and composition helpers
    BoundaryTypes/  Pre-defined boundary conditions (Reflective, Periodic, Outflow, etc)
  FluidElement/     Primitive/conservative states, fluxes, arithmetic, physicality checks
  Godunov/          Structured grids, sweep logic, split/unsplit updates
    CFL/            CFL Timestep Constraint
    Sweeping/       Godunov Scheme (including 1D, Split, and Unsplit) 
    TVD/            MUSCL Reconstruction and limiter options
  main/             Program Main
  MHD/              Constrained Transport
  Refinement/       Experimental AMR grid wrappers
  Riemann/          Exact and approximate Riemann solvers
  
  Config.h          Solver, reconstruction, CFL, MHD, and threading configuration
  Constants.h       Physical constants
  Problem.cpp       Example problem initialization and output handling

DRAGONWING/       Parallel launch/synchronization helpers for AMR grids
Testing/          Test executable and component-level verification routines
```

## Configuration

Most numerical choices are compile-time settings in `DRAGON/Config.h`.

Important switches include:

- `MHD`: enable magnetohydrodynamic state variables and MHD solvers.
- `RIEMANN_DEFAULT`: selects the default flux solver. The current default is `RIEMANN_HLLX`, which uses HLLC for pure hydrodynamics and HLLD for MHD.
- `RIEMANN_VERIFY_FALLBACK`: verifies approximate fluxes and falls back to more robust solvers when needed.
- `CFL_CALCULATION`: chooses CFL speed aggregation. The current default is `CFL_ADD`.
- `CFL_coeff`: global CFL coefficient. The current value is `0.3`.
- `MUSCL_Hancock`: enables second-order MUSCL-Hancock reconstruction. Comment it out for first-order Godunov updates.
- `MUSCL_DEFAULT_LIMITER`: selects the TVD limiter used by MUSCL reconstruction.
- `DIMENSION_UNSPLIT`: enables unsplit multidimensional advancement. When enabled with `CTU`, multidimensional updates use corner-transport-upwind terms.
- `core_count`: helps the root AMR grid choose how many child grids to create.

When `TESTMODE` is enabled, several normally compile-time choices can be selected at runtime through the `CONFIG` namespace so tests can exercise multiple solver modes.

## Running the Example

The `DRAGON` executable currently runs the 2D blast-wave setup defined in `Problem.cpp`. Problem construction is exposed through `Problem::makeProblem()`, while per-cycle reporting and output are handled by `Problem::cycleComplete()`:

- grid size: `1000 x 1000`
- cell size: `dx = dy = 0.01`
- boundary condition: reflective
- ambient pressure: `0.1`
- blast pressure: `10.0`
- blast radius: `0.3`
- output: one CSV density field per frame

The output path is currently hard-coded in `Problem.cpp`:

```text
/Users/bobbiemarkwick/DRAGON_OUT/frame-####.csv
```

Change the path in `Problem.cpp` for your local run. To configure a different simulation, update the problem initialization and cycle-completion behavior there; `main.cpp` now contains only the simulation loop.

## Running Tests

Build and run the `DRAGON_TESTS` product from Xcode. The test runner calls the verification routines in `Testing/Testing.cpp` and prints progress to standard output. A successful run ends with:

```text
All tests passed.
```

The current test suite covers the main numerical components, but `Testing.cpp` still lists AMR grid and DRAGONWING coverage as TODOs.

## Minimal Usage Example

```cpp
#include "Grid.hpp" //Or DistGrid.hpp to run across multiple cores

int main() {
    Grid2D grid(128, 128, 0.01, 0.01); //Or DistGrid2D to run across multiple cores
    grid.boundary = Boundary::Reflective();

    for (int i = 0; i < grid.getSizeX(); ++i) {
        for (int j = 0; j < grid.getSizeY(); ++j) {
            grid[i, j].rho = 1.0;
            grid[i, j].p = 0.1;
            grid[i, j].v = {0.0, 0.0, 0.0};
        }
    }

    grid.advance(0.001);
}
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

- Finish Implementing Constrained Transport 
- Add proper file output and simulation restarts
- Complete and validate AMR behavior.
- Replace hard-coded demo parameters and output paths with runtime configuration.
- Add documented example problems and plotting/post-processing helpers.


## Use of Generative AI

DRAGON was developed by Bobbie Markwick as an independent C++ hydrodynamics/MHD code. Generative AI tools were used during development as an auxiliary aid for discussion, review, documentation, and some test-writing support, but not as an automated system for producing the main codebase.

AI assistance was used for tasks such as:
* checking reasoning about hydrodynamics, MHD, constrained transport, CTU updates, boundary conditions, and AMR design;
* discussing numerical algorithms and implementation strategies;
* reviewing selected code snippets for possible bugs, inconsistencies, or edge cases;
* drafting or helping draft portions of the unit-test suite;

The main DRAGON implementation was written and integrated manually by the author. Design decisions, algorithm selection, debugging, validation, and interpretation of results remain the author’s responsibility. AI-generated suggestions were treated as advisory rather than authoritative, and in multiple cases required correction, rejection, or reinterpretation before use.

The numerical methods implemented in DRAGON are based on standard published algorithms in computational fluid dynamics and magnetohydrodynamics, including Godunov-type finite-volume methods, approximate Riemann solvers, MUSCL/CTU-style reconstruction and time integration, and constrained transport for magnetic-field evolution. Where appropriate, these methods are cited directly in the source code and documentation.

In short: generative AI was used as an interactive reference, code-review assistant, documentation aid, and unit-test drafting assistant. DRAGON itself is an original implementation, with responsibility for correctness, testing, validation, and scientific interpretation resting with the author.
