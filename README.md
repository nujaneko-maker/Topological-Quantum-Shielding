# Topological "Decoherence Shielding" Pipeline via Inverse Design
An interdisciplinary data science project assessing structural mode density and environmental noise filtering in biological networks.

## 🌟 Executive Summary
This project implements a multiscale mathematical pipeline to investigate the **Inverse Topological Design** of proteins. By constructing a Graph Laplacian matrix from real atomic coordinates, the model evaluates how network geometry acts as a spectral filter to mitigate quantum decoherence against macro-environmental noise (37°C).

## 📊 Core Findings & Experimental Data
Using real structural configurations fetched from the RCSB Protein Data Bank, the algorithm solved the full spectrum for two verified architectures:
*   **Target (Active Photoreceptor):** Arabidopsis thaliana Cryptochrome-1 PHR domain (**PDB: 1U3C**) -> \(\gamma_{\text{eff}} = 0.1227\)
*   **Control (Non-Quantum Baseline):** E. coli Ammonium Transporter AmtB (**PDB: 1U7C**) -> \(\gamma_{\text{eff}} = 0.1769\)

**Result:** Cryptochrome exhibits a **30.6% quantitative reduction** in environmental noise penetration, mathematically validating that its spatial network topology naturally engineers a vibrational band-gap near the critical noise threshold (ω = 5.0).

## 🧬 Pipeline Architecture & Methodology
The data science pipeline is built entirely in Python and executes the following multiscale architecture:
1. **Biopython Integration:** Automatic extraction of Alpha Carbon (\(C_\alpha\)) 3D coordinates.
2. **Graph Theory Modeling:** Construction of adjacency matrix W using a strict \(7.0\text{ \AA}\) spatial threshold, followed by full symmetric Graph Laplacian solving (L = D - W).
3. **Spectral Filtering Computation:** Applying Silverman's Bandwidth Rule (1986) to compute continuous mode density \(F_{\text{topo}}(\omega)\).
4. **Numerical Integration:** Continuous integration of the overlapping spectra to determine the final effective decoherence summary statistic (\(\gamma_{\text{eff}}\)).

## 🚀 How to Run on Google Colab
The code is optimized for cloud execution. Click the file `topological_shielding_pipeline.ipynb` inside this repository and select the **"Open in Colab"** badge to replicate the charts instantly.
