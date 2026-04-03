# Low-Field Metal-Insulator Transition in AB-Stacked Bilayer Graphene

**Companion code for:**
> A. Chakraborty, A. Rodin, S. Adam, G. Vignale,
> *Low-Field Metal-Insulator Transition in AB-Stacked Bilayer Graphene*,
> Physical Review B (2026). DOI: [to be inserted upon publication]

---

## Overview

This repository contains the MATLAB code used to generate the numerical results and figures in the paper. The code implements a full 4×4 tight-binding model for Bernal-stacked bilayer graphene, including skew interlayer couplings (γ₃, γ₄) that produce trigonal warping of the low-energy Fermi surface. An in-plane magnetic field enters via the Peierls substitution, and a transverse displacement field enters as an interlayer potential asymmetry V. Together, these drive an insulator–metal (IM) transition at critical fields approximately two orders of magnitude smaller than predicted without trigonal warping.

---

## Repository Contents

| File | Description |
|---|---|
| `BandStructure_FermiSurface.m` | Computes the low-energy band structure, estimates the Fermi energy via charge neutrality, and produces 2D Fermi surface contour plots and 3D band surface visualizations |
| `PhaseDiagram_V_Phi.m` | Maps the IM phase boundary in the (V, φ) plane by locating the critical displacement field V_c at each flux value via bisection search over the indirect band gap |
| `README.md` | This file |

---

## Requirements

- MATLAB R2020b or later
- **Parallel Computing Toolbox** (required for `parfor` loops; code can be run without it by replacing `parfor` with `for`, at significantly increased runtime)

No additional toolboxes or external dependencies are needed.

---

## Usage

### 1. Band Structure and Fermi Surface (`BandStructure_FermiSurface.m`)

Open the script and set the physical parameters as you please in the clearly marked sections at the top:

```matlab
v    =;   % Displacement field in eV (try 0, 0.0004, 0.0007)
phi1 =;   % In-plane magnetic flux (x-component)
phi2 =;   % In-plane magnetic flux (y-component)
```

Run the script. It will:
1. Compute eigenvalues of the 4×4 Hamiltonian on a 2000×2000 k-grid around the K valley
2. Estimate the Fermi energy by matching conduction and valence pocket areas
3. Display the 2D Fermi surface contour plot
4. Display the 3D band surface with Fermi plane
5. Save a high-resolution PNG of the contour plot

**Typical runtimes** (with Parallel Computing Toolbox, 8 cores):
- Eigenvalue calculation: ~5–10 minutes
- Fermi energy search: ~2–5 minutes

**Key outputs:** Fermi energy printed to console; figures displayed and saved.

**To reproduce Fig. 1 of the paper:** run with `v = 0`, `v = 0.0004`, and `v = 0.0007` eV in succession (corresponding to panels b, c, d respectively).

---

### 2. Phase Diagram (`PhaseDiagram_V_Phi.m`)

Set the flux sweep range and grid resolution at the top of the script:

```matlab
phi_vals     = linspace(0.02, 1.0, 100);   % Flux values to sweep
N_k          = 500;                         % k-grid resolution
v_min_search = 0.0  / 1000;                % Bisection lower bound (eV)
v_max_search = 1000 / 1000;                % Bisection upper bound (eV)
```

Run the script. It will:
1. For each flux value φ, precompute structure factors on the full k-grid
2. Run a bisection search to find V_c(φ) where the indirect gap closes
3. Plot the phase boundary V_c vs φ
4. Save a high-resolution PNG

**Typical runtimes** (with Parallel Computing Toolbox, 8 cores):
- Full sweep of 100 flux values at N_k = 500: ~30–60 minutes
- For a quick test, reduce to `phi_vals = linspace(0.02, 0.1, 20)` and `N_k = 200`

**To reproduce Fig. 2(c) of the paper:** use the default parameters as set in the script. The low-field region (φ < 0.01) corresponds to B < 25 T.

---

## Physical Parameters

All hopping parameters are in eV and follow the values of Kuzmenko et al. (PRB 2009) as used in Seiler et al. (Nat. Commun. 2024):

| Parameter | Value (eV) | Description |
|---|---|---|
| γ₀ | −2.610 | Nearest-neighbor intralayer hopping |
| γ₁ |  0.361 | Interlayer dimer hopping (B_top – A_bot) |
| γ₃ | −0.283 | Skew interlayer hopping (trigonal warping) |
| γ₄ |  0.138 | Skew interlayer hopping (electron–hole asymmetry) |
| Δ  |  0.015 | Dimer/non-dimer on-site energy difference |

The in-plane magnetic flux is related to the physical field by φ = B·a·ℓ/Φ₀, where a ≈ 1.42 Å is the carbon bond length, ℓ ≈ 3.35 Å is the interlayer spacing, and Φ₀ = h/e ≈ 4.14 × 10⁵ T·Å².

---

## Citation

If you use this code, please cite the associated paper:

```
@article{Chakraborty2026BLG,
  author  = {Chakraborty, Amarnath and Rodin, Aleksandr and Adam, Shaffique and Vignale, Giovanni},
  title   = {Low-Field Metal-Insulator Transition in {AB}-Stacked Bilayer Graphene},
  journal = {Physical Review B},
  year    = {2026},
  doi     = {to be inserted}
}
```

---

## License

This code is released under the MIT License. See `LICENSE` for details.

---

## Contact

Amarnath Chakraborty — amarnathchakraborty95@gmail.com
