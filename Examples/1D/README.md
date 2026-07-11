# DRAGON 1D Examples

This directory contains one-dimensional shock-tube and wave-interaction problems for
checking DRAGON's hydrodynamic and magnetohydrodynamic solvers. Each example is a
self-contained problem setup with its own `Config.h`, `Problem.cpp`, and reference
plot image.

Instructions for running example problems can be found in the project's base README.
For the Brio-Wu problem, be sure to also copy `Constants.h`, since that example uses `gamma = 2`. Be sure also to restore the original file before switching to a different problem.

## Example Inventory

| Example | MHD | Final time | Resolution | Limiter |
| --- | --- | ---: | ---: | --- |
| Sod | NO | 0.2 | 8192 | MC |
| Lax | NO | 0.1 | 8192 | MC |  
| Strong Shock | NO | 0.012 | 8192 | MC |
| Shu-Osher | NO | 1.8 | 10240 | MC |
| Brio-Wu | YES | 0.1 | 8192 | MC |
| Dai-Woodward | YES | 0.2 | 8192 | minmod |
| Ryu-Jones | YES | 0.1 | 8192 | minmod |

All examples use outflow boundary conditions. The Riemann Solver was HLLC for pure hydrodynamics and HLLD for MHD. Most examples use the MC limiter, though the last two were run with MINMOD since MC produced spurious oscillations near strong shocks in those.

## Initial Conditions

The examples initialize piecewise states directly in `Problem.cpp`.

### Hydrodynamics

| Example | Left state | Both states | Right state |
| --- | --- | --- |--- |
| Sod | $\rho = 1.0$, $p = 1.0$| $v = (0, 0, 0)$ | $\rho = 0.125$, $p = 0.1$ |
| Lax | $\rho = 0.445$, $p = 3.528$, $v_x = 0.698$ | $v_y = v_z = 0$ |$\rho = 0.5$, $p = 0.571$, $v_x = 0$ |
| Strong Shock | $p = 1000.0$ | $\rho = 1.0$, $v = (0, 0, 0)$ | $p = 0.01$ |
| Shu & Osher | $\rho = 3.857143$, $v_x = 2.629369$, $p = 10.333333$ | $v_y = v_z = 0$ | $\rho = 1.0 + 0.2 \sin(5x)$, $v_x = 0$, `p = 1.0` |
| Brio & Wu | $\rho = 1.0$, $p = 1.0$, $B_y = \sqrt{4\pi}$ | $v = 0$, $B_x = 0.75\sqrt{4\pi}$,  $B_z = 0$| $\rho = 0.125$, $p = 0.1$,  $B_y = -\sqrt{4\pi}$ |
| Dai & Woodward | $\rho = 1.08$, $p = 0.95$, $v = (1.2, 0.01, 0.5)$, $B_y = 3.6$ | $B_x = 2.0$, $B_z = 2.0$  |$\rho = 1.0$, $p = 1.0$, $v = (0, 0, 0)$, $B_y = 4.0$ |
| Ryu & Jones |  $v_x = 10.0$, $p = 20.0$ | $\rho = 1.0$, $B = (5, 5, 0)$, $v_y = v_z = 0$ |  $v_x = -10.0$, $p = 1.0$,  |

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

The PNG files included beside each example show output generated from the example as provided.
