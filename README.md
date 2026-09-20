# Topological Spectral Susceptibility Model (v8.0)
An open-system non-equilibrium thermodynamic and matrix perturbation framework mapping Edge Spectral Leverage landscapes, mechanical fragility, and multi-scale vibronic routing in 3D biomolecular networks.

## 🌟 Paradigm Shift: From Descriptive Feedback to Discovery Mechanics
This repository documents the multi-scale formalization and rigorous structural evolution of an independent computational biophysics pipeline, moving from legacy 1D graph models (v1-v3.0) and standardized 3D Elastic Network Models (v5.1R) to dynamic structural perturbation null models (v6-v7) and first-principles matrix perturbation theory (v8.0).

## 🧬 Multiscale Mathematical Architecture
* **Macroscopic Operator Chassis (3D ENM Hessian):** Parses raw Alpha Carbon coordinates (e.g., Cryptochrome-1, PDB: 1U3C) to construct a global isotropic Hessian matrix $H_{\rm ENM}$ and extracts positive eigenvalues and normalized eigenvectors after removing rigid-body zero modes.
* **First-Principles Matrix Perturbation:** Computes exact matrix differential operators to evaluate eigenvalue sensitivity and defines Global Edge Spectral Leverage ($L_e$) to extract invariants like the Spectral Gini Coefficient ($L_{\rm Gini}$) and Edge Participation Ratio ($\text{IPR}_L$).

## 📊 Causal Verification & Microscopic Closure
* **Null Models & Ablation:** Implements degree-preserving double-edge swaps and targeted edge ablation to map spectral robustness versus fragility, isolating cases where frequency profiles remain fixed while eigenvector geometry transitions.
* **Lindblad Dynamics:** Drives the master decoherence operator within a 4x4 electronic Hilbert space under the GKSL Master Equation, enforcing thermodynamic validity with a rigid numerical Positivity Guard.

## 🔬 Popperian Falsification Matrix & Repository
The framework establishes empirical criteria (F1-F4) to test multi-scale decoupling, isospectral-IPR collapse, vibronic null contrast, and additive interaction sums against laboratory data. 

For the complete implementation details, code execution, and dependency setup, please refer to the repository source files in the referenced GitHub project.
