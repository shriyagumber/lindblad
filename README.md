# Lindblad Dynamics Codebase

## Overview
Small Fortran driver for two‑state nonadiabatic dynamics. It reads precomputed Hamiltonians (`avg-ham/Ham*_re`, `Ham*_im`), interpolates them to electronic timesteps, propagates density matrices with both time‑dependent Schrödinger and Lindblad equations, and estimates surface‑hopping populations.

## Layout
- `src/` – all Fortran sources (moved from project root)
  - `main_lindblad.f95` – sets constants, reads Hamiltonians, orchestrates propagation and averaging
  - `readHnL.f95` – reads real/imag parts of Hamiltonians and repeats files as needed
  - `inter.f95` – linear interpolation of Hamiltonians/Lindblad couplings across electronic steps
  - `solve_tdse.f95` – RK4 propagator for coherent TDSE (uses imaginary NACs as off‑diagonal terms)
  - `solve_me2.f95` – Lindblad ME RK4 propagator with simple Golden‑Rule rate for |g⟩→|e⟩
  - `fssh.f95` – fewest‑switches surface hopping population estimator
  - `printfssh.f95` – writers for electronic and surface‑hopping populations
  - `solve_me.f95`, `solve_me3.f95` – earlier Lindblad variants (not referenced by main)
- `avg-ham/` – required input data (real/imaginary Hamiltonian slices)
- `makefile` – builds `run_lindblad`

## Build
```
make            # produces ./run_lindblad
make clean      # remove objects, modules, and binary
```
Ensure `gfortran` is available and input files exist in `avg-ham/`.

## How it runs
1) `read_ham` loads the two targeted states (`state_low=16`, `state_high=17`) from each available `Ham` file, converts Rydberg→eV, and fills `ham` and `lindblad` tensors, repeating files cyclically to reach `n_files`.
2) `interpolate` linearly fills electronic sub‑steps (`nelect` per nuclear step) for both Hamiltonian and Lindblad coupling arrays.
3) For each initial condition and stochastic realization:
   - `me` (Lindblad RK4) and `tdse` (coherent RK4) propagate density matrices.
   - `fsshpop` estimates surface‑hopping populations from TDSE diagonals.
4) Populations are averaged and written via `print_se`/`print_sh` to `lindblad_pop.out`, `tdse_pop.out`, and `sh_pop.out`.

## Review Notes / Potential Issues
- `fssh` module exports only `func` (undefined); `calc_boltz` and `fsshpop` are not public—code relying on the module name alone still works because default visibility is public, but the `public func` line is a bug and should be corrected.
- `solve_tdse` allocates `drho`, `H`, `L` every step and discards `drho`; the commutator uses `H+L` so Lindblad dissipator is not included (may be intentional, but naming suggests otherwise).
- `solve_me2` allocates `nac` and `fermi_rate` each step; `fermi_rate` divides by `(H(2,2)-H(1,1))` without zero‑check and can go negative if the gap changes sign, leading to `sqrt` of a negative value.
- Output writers open with `status='new'`; rerunning without cleaning will fail if output files already exist—consider `status='replace'` for convenience.
- `main_lindblad` slices `ham_intrpl`/`lindblad_intrpl` with `-2` on the upper bound, so the final two electronic steps per trajectory are dropped; verify this truncation is intentional.

## Quick Start
1) Place `Ham*_re` and `Ham*_im` files in `avg-ham/` (Rydberg units).
2) Run `make` to build `run_lindblad`.
3) Execute `./run_lindblad`; results appear in the repo root as `lindblad_pop.out`, `tdse_pop.out`, and `sh_pop.out`.
