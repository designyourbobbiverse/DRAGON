# DRAGON 1D Examples

This directory contains one-dimensional shock-tube and wave-interaction problems for
checking DRAGON's hydrodynamic and magnetohydrodynamic solvers. Each example is a
self-contained problem setup with its own `Config.h`, `Problem.cpp`, and reference
plot image.

Instructions for running example problems can be found in the project's base README.
For the Brio-Wu problem, be sure to also copy `Constants.h`, since that example uses `gamma = 2`. Be sure also to restore the original file before switching to a different problem.

## Problems

The examples initialize piecewise states directly in `Problem.cpp`.

| Example |  MHD |  Left state | Both states | Right state |
| --- | --- | --- | --- | --- |--- |
| Sod | NO | $\rho = 1.0$, $p = 1.0$| $v = 0$ | $\rho = 0.125$, $p = 0.1$ |
| Lax | NO |  $\rho = 0.445$, $v_x = 0.698$ , $p = 3.528$| $v_y = v_z = 0$ |$\rho = 0.5$, $v_x = 0$, $p = 0.571$ |
| Strong Shock | NO |  $p = 1000.0$ | $\rho = 1.0$, $v = 0$ | $p = 0.01$ |
| Shu & Osher | NO |  $\rho = 3.857143$, $v_x = 2.629369$, $p = 10.333333$ | $v_y = v_z = 0$ | $\rho = 1.0 + 0.2 \sin(5x)$, $v_x = 0$, $p = 1.0$ |
| Brio & Wu |  YES | $\rho = 1.0$, $p = 1.0$, $B_y = \sqrt{4\pi}$ | $v = 0$, $B_x = 0.75\sqrt{4\pi}$,  $B_z = 0$| $\rho = 0.125$, $p = 0.1$,  $B_y = -\sqrt{4\pi}$ |
| Dai & Woodward |  YES  | $\rho = 1.08$, $v = (1.2, 0.01, 0.5)$, $p = 0.95$, $B_y = 3.6$ | $B_x = 2.0$, $B_z = 2.0$  |$\rho = 1.0$, $v = (0, 0, 0)$, $p = 1.0$, $B_y = 4.0$ |
| Ryu & Jones |  YES |   $v_x = 10.0$, $p = 20.0$ | $\rho = 1.0$, $v_y = v_z = 0$, $B = (5, 5, 0)$ |  $v_x = -10.0$, $p = 1.0$,  |

All examples use outflow boundary conditions. 

The Riemann Solver was HLLC for pure hydrodynamics and HLLD for MHD. Hydrodynamic examples use the MC limiter, while the MHD examples instead use MINMOD since MC produced spurious oscillations near strong shocks in those.



## Plotting Results

After running an example, configure `DRAGONGAZE/Config.py` so `hdf_dir` matches
the example's output directory. Also set `plot_title` to the name of the test you are running.  Then run the appropriate plotting script from
`DRAGONGAZE/`. For hydrodynamic examples:

```bash
python GAZE1D.py
```

For MHD examples:

```bash
python GAZE1DMHD.py
```

The PNG files included beside each example show output generated from the example as provided. For MHD examples with no $B_z$ components, the plotting script was modified to omit the $v_z$ and $B_z$ plots.
