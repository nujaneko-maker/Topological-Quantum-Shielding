# Topological Spectral Homeostasis Model (v2.1)
An open-system thermodynamic and computational framework evaluating hypothesized emergent spectral shielding in biological networks via a 4-state Factorial Double-Ablation Protocol.

## 🌟 Project Evolution & Epistemic Humility
This repository documents the structural evolution of a multiscale biological data science project.
*   **Phase 1 (Legacy Archive - `/archive_v1`):** Initial formulations explored naive physical shielding assumptions under discrete boundary constraints. These early models suffered from dimensional inconsistencies and attempted to claim total system-bath isolation, which violates non-equilibrium thermodynamic laws.
*   **Phase 2 (Current Model - `main.py`):** Re-engineered into a rigorous **Topological Spectral Homeostasis** framework. Rather than claiming biological systems annihilate environmental noise, the model demonstrates that open systems utilize non-linear network topology to minimize system-bath coupling via **Spectral Impedance Matching** constrained within a functional Hilbert space \(L^2(\Omega, d\omega)\).

## 🧬 Generalized Mathematical Foundation
The effective decoherence summary statistic (\(\gamma_{\rm eff}\)) of this open system is modeled as a functional leak rate that must satisfy \(L^1\) integrability to achieve flat convergence (\(\gamma_{\rm eff}(t) \to 0\) and \(d\gamma_{\rm eff}/dt \to 0\) as \(t \to \infty\)):
\[\gamma_{\rm eff}(t) = \int_0^\infty W(\omega,t) J_{\rm env}(\omega,t) \vert{}F_G(\omega,t)\vert{}^2 M(\omega,t) \vert{}F_{\rm topo}(\omega,t)\vert{}^2 d\omega\]

Where \(J_{\rm env}\) is the ambient thermal noise spectrum (Gaussian band centered at \(\omega = 5.0\), representing \(37^\circ\text{C}\)), \(\vert{}F_G\vert{}^2\) is the common-mode electrostatic boundary suppression factor, \(M(\omega)\) is the microbiome non-linear frequency modulation transfer function, and \(\vert{}F_{\rm topo}\vert{}^2\) is the continuous structural mode density derived from the symmetric graph Laplacian (\(L = D - W\)).

## 📊 Empirical Factorial Verification (Double-Ablation)
The computational pipeline validates structural spectral filtering using real atomic coordinates fetched from the RCSB Protein Data Bank, comparing the quantum-active **Arabidopsis thaliana Cryptochrome-1 PHR domain (PDB: 1U3C, 485 Cα nodes)** against an **E. coli AmtB metabolic baseline (PDB: 1U7C, 372 Cα nodes)**.

To establish strict topological causality, the pipeline executes a **Double-Ablation Factorial Experiment** by holding the grounding interface fixed (\(F_G = 0.5\)) and calculating the cross-term interaction signature:
\[\Delta_{\rm int} = \gamma_{11} - \gamma_{10} - \gamma_{01} + \gamma_{00} = \int WJ_{\rm env}\vert{}F_G\vert{}^2 (F_1-F_0)(M_1-M_0) d\omega\]

### 🛠️ Factorial Trajectory Results (See `TSH_v21_double_ablation.png`):
*   **State 00 (Intact):** Native 1U3C topology + Intact Modulator (\(M_0\)) -> \(\gamma_{\rm eff} = \mathbf{0.26355172}\) (Achieves elite flat convergence at the lowest exposure plateau).
*   **State 10 (Protein Topology Ablated):** Randomized 1U3C topology (degree-preserving edge swap) + Intact Modulator (\(M_0\)) -> \(\gamma_{\rm eff} = \mathbf{0.41409509}\)
*   **State 01 (Microbiome Ablated):** Native 1U3C topology + Ablated Modulator (\(M_1 = 1.0\)) -> \(\gamma_{\rm eff} = \mathbf{0.65471067}\)
*   **State 11 (Double Ablation):** Randomized 1U3C topology + Ablated Modulator (\(M_1 = 1.0\)) -> \(\gamma_{\rm eff} = \mathbf{0.87077734}\)

### 🔍 Causal Interaction Output:
\[\Delta_{\rm int} = 0.87077734 - 0.41409509 - 0.65471067 + 0.26355172 = \mathbf{0.06552329}\]
**Conclusion:** \(\Delta_{\rm int} > 0\) confirms a **super-additive interaction effect** (\(D_4 > D_2 + D_3\)). The failure of both defensive tiers causes a non-linear explosion in noise penetration, mathematically validating that spectral homeostasis emerges as a synergistic multiscale phenomenon.

## 🚀 Execution & Reproducibility
The core engine is optimized for Python 3.8+. To reproduce the results and the 2-panel diagnostic dashboard locally:
```bash
git clone https://github.com
cd Topological-Quantum-Shielding
pip install biopython numpy scipy networkx matplotlib pandas
python main.py
```
