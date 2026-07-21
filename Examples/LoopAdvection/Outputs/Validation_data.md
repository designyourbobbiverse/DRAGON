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
| 64 | 0.0286377 | — |  0.152362 | — | 0.99671 | 
| 128 | 0.017229 | 0.73 | 0.0949007 | 0.68 | 0.997568 |
| 256 | 0.00969727 | 0.83 | 0.0666787 | 0.51 | 0.996292 |
| 512 | 0.0064274 | 0.59 | 0.0513134 | 0.38 | 0.999162 |


### 3D, Uniform $B_\phi$

| n | L1 error | L1 rate | L2 error | L2 rate | Linf error |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.0148854 | — |  0.112979 | — | 0.99957 | 
| 128 | 0.00964739 | 0.63 | 0.0728129 | 0.63 | 0.997978 |
| 128 | 0.00487479 | 0.98 | 0.0461529 | 0.66 | 0.997978 |


### 3D, Gaussian $A_z$

| n | L1 error | L1 rate | L2 error | L2 rate | Linf error |
|---:|---:|---:|---:|---:|---:|
| 64 | 0.323645 | — |  2.65623 | — | 32.7052 | 
| 128 | 0.160526 | 1.01 | 1.42665 | 0.90 | 22.1701 |
| 256 | 0.0457757 | 1.81 | 0.43565 | 1.71 | 9.34788 |


