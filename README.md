# Flux Fermi Model

A Monte Carlo simulation of particle flux using Fermi-Dirac statistics to model exclusion principles in a lattice environment.

## Overview

This project implements a stochastic simulation of particles moving on a 2D lattice. It employs the Fermi-Dirac distribution to calculate transition probabilities for particle movements, effectively modeling systems where site occupancy is limited (exclusion principle). The simulation tracks two types of particles (A and B) and measures various physical quantities such as mobility, order parameters, and Gini coefficients to analyze the system's state.

## How it Works

The simulation uses a Monte Carlo method with the following key characteristics:

- **Lattice**: A 2D grid of dimensions `l_x` by `l_y` with periodic boundary conditions.
- **Particles**: Two species of particles, A and B. A total of `n_total` particles are present.
- **Dynamics**:
  - Particles attempt to move to adjacent sites.
  - Movement probabilities are governed by a Fermi-Dirac-like function: $P \propto \frac{1}{1 + e^{\alpha(\Delta E)}}$.
  - $\Delta E$ depends on the local density of A and B particles (`dif` and `deq` parameters).
  - The exclusion principle is enforced (site capacity limits).
- **Measurements**:
  - **Mobility**: The fraction of attempted moves that are successful.
  - **Order Parameter ($\phi$)**: Measures the degree of particle segregation or ordering.
  - **Current Density ($J$)**: The net flux of particles.
  - **Gini Coefficient**: A measure of inequality in the spatial distribution of particles.

## Compilation and Usage

The project is a single-file C++ application. You can compile it using a C++ compiler like `g++`.

### Prerequisites

- A C++ compiler (e.g., `g++`, `clang++`)
- Standard C++ library

### Build

```bash
g++ -O3 main.cpp -o flux_fermi
```

### Run

```bash
./flux_fermi
```

## Development Environment

This project uses [Nix](https://nixos.org/) to provide a reproducible development environment.

### Using Nix Flakes

To enter the development shell with all dependencies (`gcc`, `gnumake`):

```bash
nix develop
```

### Using direnv

This project is configured with `direnv`. If you have `direnv` and `nix-direnv` installed, the environment will be automatically loaded when you enter the directory.

To allow the environment:

```bash
direnv allow
```

## Output

The program generates two output files containing statistical data for different system densities ($\rho_{total}$):

### `data1.dat` - Averages

Contains sample averages of the measured quantities.
Columns:

1. $\rho_{total}$ (Total density)
2. $\rho_a$ (Density of A)
3. $\rho_b$ (Density of B)
4. Average Mobility
5. Order Parameter X ($\phi_x$)
6. Order Parameter Y ($\phi_y$)
7. Current Density Difference ($J_a - J_b$)
8. Current Density A ($J_a$)
9. Current Density B ($J_b$)
10. Gini Index A
11. Gini Index B
12. Probability of Lane Formation
13. Probability of Clogging
14. Probability of Other states

### `data2.dat` - Variances

Contains the sample variances for the corresponding quantities in `data1.dat`.
Columns match the order of `data1.dat` but represent variances (mobility, $\phi_x$, etc.).
