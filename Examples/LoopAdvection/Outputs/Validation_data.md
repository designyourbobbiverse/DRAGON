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

| n | L1 error | L1 rate | L2 error | L2 rate | Linf error |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.0286387 | — |  0.152378 | — | 0.996704 | 
| 128 | 0.0172391 | 0.73 | 0.0949274 | 0.68 | 0.997576 |
| 256 | 0.00968561 | 0.83 | 0.0666294 | 0.51 | 0.996234 |
| 512 | 0.00641105 | 0.59 | 0.0512477 | 0.38 | 0.99921 |


### 2D, Gaussian $A_z$

| n | L1 error | L1 rate | L2 error | L2 rate | Linf error |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.601428 | — |  3.53229 | — | 31.4691 | 
| 128 | 0.273266 | 1.14 | 1.7344 | 1.03 | 20.7554 |
| 256 | 0.0753625 | 1.86 | 0.521772 | 1.73 | 8.74698 |
| 512 | 0.0261989 | 1.52 | 0.187809 | 1.47 | 3.77506 |


### 3D, Uniform $B_\phi$

| n | L1 error | L1 rate | L2 error | L2 rate | Linf error |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.0139403 | — |  0.107475 | — | 0.999567 | 
| 128 | 0.00866564 | 0.69 | 0.066301 | 0.70 | 0.998004 |
| 256 | 0.00485955 | 0.83 | 0.0461968 | 0.52 | 0.993108 |


### 3D, Gaussian $A_z$

| n | L1 error | L1 rate | L2 error | L2 rate | Linf error |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.323758 | — |  2.65926 | — | 32.7321 | 
| 128 | 0.161105 | 1.01 | 1.43543 | 0.90 | 22.2412 |
| 256 | 0.0459888 | 1.81 | 0.438878 | 1.71 | 9.38652 |
