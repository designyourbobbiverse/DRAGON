# DRAGON diagonal advection validation data

The following is a summary of convergence data calculated from running the B-Field Loop Advection example at different resolutions

## Problem definition
- Periodic domain, n x n/2 x n/2
- Errors are measured in $|B|$ against the analytic solution.
- Reported norms: `L1`, `L2`, and `Linf`.

## Convergence-rate formula
For two successive resolutions $n$ and $2n$, the observed convergence rate is

$$
q = \log_2(E_n/E_{2n})
$$

where $E_n$ is the error observed at resolution $n$.
## Results
### 2D, Uniform $B_\phi$

| n | L1 error | L1 rate |
|---:|---:|---:|
| 64 | 5.21079e-08 | - |
| 128 | 2.71831e-08 | 0.94 |
| 256 | 1.61524e-08 | 0.74 |
| 512 | 1.0518e-08 | 0.62 |

### 2D, Gaussian $A_z$

| n | L1 error | L1 rate |
|---:|---:|---:|
| 64 | 2.08777e-07 | - |
| 128 | 8.03226e-08 | 1.38 |
| 256 | 2.7019e-08 | 1.57 |
| 512 | 9.0612e-09 | 1.57 |


### 3D, Uniform $B_\phi$

| n | L1 error | L1 rate |
|---:|---:|---:|
| 64 | 5.65724e-08 | - |
| 128 | 3.0396e-08 | 0.90 |
| 256 |  |  |

### 3D, Gaussian $A_z$

| n | L1 error | L1 rate |
|---:|---:|---:|
| 64 | 2.30943e-07 | - |
| 128 | 9.19098e-08 | 1.32 |
| 256 | 3.21374e-08 | 1.51 |
