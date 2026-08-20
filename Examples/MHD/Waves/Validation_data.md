# DRAGON Linear Wave Convergence data

The following is a summary of convergence data calculated from running the MHD Linear Wave example at different resolutions

## Problem definition
- Periodic domain, 2n x n x n
- Errors are measured on all components against the initial state, then combined into a single error magnitude

## Convergence-rate formula
For two successive resolutions $n$ and $2n$, the observed convergence rate is

$$
q = \log_2(E_n/E_{2n})
$$

where $E_n$ is the error observed at resolution $n$.

## Results

### 3D, Fast

| n | (+) L1 error | (+) L1 rate | (-) L1 error | (-) L1 rate |
|---:|---:|---:|---:|---:|
| 16 | 5.38918e-07 | — | 5.38918e-07 | — |
| 32 | 1.17147e-07 | 2.20 | 1.17511e-07 | 2.20 |
| 64 | 3.46558e-08 | 1.76 | | 3.46258e-08 | 1.76 |

### 3D, Alfven

| n | (+) L1 error | (+) L1 rate | (-) L1 error | (-) L1 rate |
|---:|---:|---:|---:|---:|
| 16 | 4.65683e-07 | — | 4.65683e-07 | — |
| 32 | 9.57591e-08 | 2.28 | 9.57327e-08 | 2.28 |
| 64 | 2.48090e-08 | 1.95 | 2.48183e-08 | 1.95 |

### 3D, Slow

| n | (+) L1 error | (+) L1 rate | (-) L1 error | (-) L1 rate |
|---:|---:|---:|---:|---:|
| 16 | 1.72296e-07 | — | 1.72296e-07 | — |
| 32 | 4.45664e-08 | 1.95 | 4.45664e-08 | 1.95 |
| 64 | 1.32143e-08 | 1.75 | 1.32151e-08 | 1.75 |

### 3D, Entropy

| n | L1 error | L1 rate |
|---:|---:|---:|
| 16 | 1.25004e-07 | — |
| 32 | 4.51067e-08 | 1.47 |
| 64 | 1.38116e-08 | 1.71 |





