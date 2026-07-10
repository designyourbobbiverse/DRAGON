# DRAGON diagonal advection validation data

## Problem definition
- Periodic domain.
- Smooth sinusoidal density perturbation advected diagonally.
- Errors are measured in density against the analytic/return-state solution.
- Reported norms: `L1`, `L2`, and `Linf`/max error.

## Convergence-rate formula
For two successive resolutions $n$ and $2n$, the observed convergence rate is

$$
q = \log_2(E_n/E_{2n})
$$

where $E_n$ is the error observed at resolution $n$.
## Results
### 2D HLLC + Minmod
| n | L1 error | L1 rate | L2 error | L2 rate | Linf error | Linf rate |
|---:|---:|---:|---:|---:|---:|---:|
| 256 | 0.000286346 | — | 0.000449423 | — | 0.00151213 | — |
| 512 | 7.83388e-05 | 1.87 | 0.000142185 | 1.66 | 0.000610277 | 1.31 |

### 2D HLLC + MC
| n | L1 error | L1 rate | L2 error | L2 rate | Linf error | Linf rate |
|---:|---:|---:|---:|---:|---:|---:|
| 128 | 0.000224884 | — | 0.000310695 | — | 0.000942581 | — |
| 256 | 5.63051e-05 | 2.00 | 9.3146e-05 | 1.74 | 0.000387981 | 1.28 |
| 512 | 1.40658e-05 | 2.00 | 2.81185e-05 | 1.73 | 0.000159312 | 1.28 |

### 2D HLLC + Superbee
| n | L1 error | L1 rate | L2 error | L2 rate | Linf error | Linf rate |
|---:|---:|---:|---:|---:|---:|---:|
| 128 | 0.000693177 | — | 0.00090973 | — | 0.00278401 | — |
| 256 | 0.000198661 | 1.80 | 0.00029218 | 1.64 | 0.00127189 | 1.13 |

### 3D HLLC + MC
| n | L1 error | L1 rate | L2 error | L2 rate | Linf error | Linf rate |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 0.00116811 | — | 0.00139005 | — | 0.00291 | — |
| 128 | 0.000321543 | 1.86 | 0.000424485 | 1.71 | 0.00121935 | 1.25 |
| 256 | 8.19464e-05 | 1.97 | 0.000128312 | 1.73 | 0.000514205 | 1.25 |

## Notes
- The MC limiter gives the cleanest smooth-advection result among the limiter data recorded here.
- Superbee is more compressive and performs worse on this smooth sine-wave advection problem, especially in `Linf`.
- The overall pattern is consistent with a limiter-based finite-volume method: near-second-order convergence in `L1`, somewhat lower order in `L2`, and lower but still convergent `Linf` behavior.
