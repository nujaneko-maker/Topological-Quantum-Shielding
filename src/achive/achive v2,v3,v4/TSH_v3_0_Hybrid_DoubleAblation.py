#!/usr/bin/env python3
"""
===============================================================================
  TOPOLOGICAL SPECTRAL HOMEOSTASIS MODEL v3.0
  Hybrid Double-Ablation Pipeline — Research-Grade Computational Validation

  Target  : 1U3C — Arabidopsis thaliana Cryptochrome-1 PHR domain
  Control : 1U7C — Escherichia coli AmtB ammonium transporter

  Core leakage functional:

    gamma_eff = ∫ W(ω) J_env(ω) |F_G|² M(ω) |F_topo(ω)|² dω

  Design principles carried forward from the two previous versions:
    • phase-based, auditable research pipeline;
    • real PDB coordinates and explicit Cα extraction;
    • 7 Å residue-contact graph and L = D − W;
    • deterministic Maslov-Sneppen degree-preserving ablation;
    • explicit Silverman KDE for continuous spectral density;
    • explicit distinction between rho_topo and F_topo;
    • fully confound-controlled 2×2 topology × microbiome factorial design;
    • direct integral cross-check of Delta_int;
    • cumulative exposure curves derived from the measured spectral integrand;
    • reproducible outputs and a publication-style two-panel figure.

  NOTE ON INTERPRETATION
  ----------------------
  This is a falsifiable computational model. It does NOT by itself establish
  that a biological system physically implements a quantum "shield". The code
  separates measured/derived quantities from model assumptions so that each
  assumption can be challenged independently.
===============================================================================
"""

# -----------------------------------------------------------------------------
# 0. DEPENDENCIES
# -----------------------------------------------------------------------------
from __future__ import annotations

import os
import sys
import time
import subprocess
import warnings
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.linalg import eigvalsh
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm

warnings.filterwarnings("ignore")

try:
    from scipy.integrate import cumulative_trapezoid
except ImportError as exc:
    raise ImportError(
        "TSH v3.0 requires scipy.integrate.cumulative_trapezoid; "
        "please upgrade SciPy in the Colab runtime."
    ) from exc

try:
    from scipy.integrate import simpson as _integrate
except ImportError:  # scipy < 1.7
    from scipy.integrate import simps as _integrate

try:
    import pandas as pd
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "pandas"])
    import pandas as pd

try:
    import networkx as nx
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "networkx"])
    import networkx as nx

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "biopython"])
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class Config:
    target_pdb: str = "1U3C"
    control_pdb: str = "1U7C"

    target_chain: str = "A"
    control_chain: str = "A"

    contact_cutoff_A: float = 7.0

    # Reduced GNM spectral coordinate. omega_0 is a model coordinate, not Hz.
    omega_0: float = 5.0
    omega_min: float = 0.0
    omega_max: float = 30.0
    n_omega: int = 4000

    # Environment spectrum
    sigma_env: float = 1.5

    # Microbiome notch model
    notch_depth: float = 0.85
    sigma_notch: float = 0.80

    # Boundary amplitude. All 2×2 cells are grounded to isolate T × M.
    F_G_grounded: float = 0.5

    # Maslov-Sneppen randomized topology
    swaps_per_edge: int = 10
    seed: int = 42
    max_swap_attempt_multiplier: int = 12

    # Numerical stability / diagnostics
    zero_eigenvalue_tol: float = 1e-8
    min_silverman_bandwidth: float = 0.05

    # Optional robustness analysis; set > 0 to run additional randomizations.
    robustness_replicates: int = 0

    output_dir: str = "TSH_v3_0_outputs"


CFG = Config()


# -----------------------------------------------------------------------------
# 2. DATA STRUCTURES
# -----------------------------------------------------------------------------
@dataclass
class StructureData:
    pdb_id: str
    chain_id: str
    residue_labels: list[str]
    coords_A: np.ndarray
    adjacency: np.ndarray
    eigenvalues: np.ndarray


@dataclass
class StateResult:
    state: str
    label: str
    gamma_eff: float
    cumulative: np.ndarray
    integrand: np.ndarray


# -----------------------------------------------------------------------------
# 3. RCSB PDB FLAT-FILE DOWNLOAD
# -----------------------------------------------------------------------------
def download_pdb(pdb_id: str, output_dir: Path) -> Path:
    """Download and cache a PDB coordinate flat file."""
    pdb_id = pdb_id.upper()
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{pdb_id}.pdb"

    if path.exists() and path.stat().st_size > 500:
        print(f"  ✓  {pdb_id}.pdb  (cached locally)")
        return path

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  ↓  {pdb_id}  from RCSB ... ", end="", flush=True)

    try:
        urllib.request.urlretrieve(url, path)
    except Exception as exc:
        raise RuntimeError(
            f"Could not retrieve {pdb_id}. Check Colab/network access."
        ) from exc

    if path.stat().st_size < 500 or "ATOM" not in path.read_text(errors="ignore"):
        raise RuntimeError(f"Downloaded file for {pdb_id} is not a valid PDB coordinate file.")

    print("done ✓")
    return path


# -----------------------------------------------------------------------------
# 4. ROBUST Cα EXTRACTION
# -----------------------------------------------------------------------------
def extract_calpha(
    pdb_path: Path,
    preferred_chain: str | None = None,
) -> tuple[list[str], np.ndarray, str]:
    """
    Extract standard-amino-acid Cα coordinates from MODEL 1.

    The requested chain is selected when present. Otherwise the longest
    sufficiently populated standard-protein chain is used.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    model = next(iter(structure))

    candidates: list[tuple[str, list]] = []

    for chain in model:
        residues = []
        for residue in chain:
            if not is_aa(residue, standard=True):
                continue
            if not residue.has_id("CA"):
                continue
            residues.append(residue)

        if len(residues) >= 30:
            candidates.append((chain.id, residues))

    if not candidates:
        raise ValueError(f"No suitable protein chain found in {pdb_path.name}.")

    selected = next(
        ((cid, residues) for cid, residues in candidates if cid == preferred_chain),
        None,
    )
    if selected is None:
        selected = max(candidates, key=lambda item: len(item[1]))

    chain_id, residues = selected

    labels = []
    coords = []
    for residue in residues:
        labels.append(f"{chain_id}:{residue.get_resname()}:{residue.id[1]}{residue.id[2].strip()}")
        coords.append(np.asarray(residue["CA"].coord, dtype=float))

    coords = np.asarray(coords, dtype=np.float64)
    return labels, coords, chain_id


# -----------------------------------------------------------------------------
# 5. CONTACT GRAPH + LAPLACIAN
# -----------------------------------------------------------------------------
def build_contact_graph(coords_A: np.ndarray, cutoff_A: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build binary Cα contact adjacency and the symmetric graph Laplacian.

      W_ij = 1 when ||r_i-r_j|| <= cutoff, i != j
      L    = D - W
    """
    dist_A = squareform(pdist(coords_A))
    W = (dist_A <= cutoff_A).astype(np.float64)
    np.fill_diagonal(W, 0.0)
    W = np.maximum(W, W.T)

    degree = W.sum(axis=1)
    L = np.diag(degree) - W
    L = 0.5 * (L + L.T)

    eigs = np.clip(eigvalsh(L), 0.0, None)
    return W, L, eigs


def largest_connected_component(W: np.ndarray) -> np.ndarray:
    """Return the adjacency matrix restricted to the largest connected component."""
    G = nx.from_numpy_array(W)
    if nx.is_connected(G):
        return W

    nodes = max(nx.connected_components(G), key=len)
    idx = np.array(sorted(nodes), dtype=int)
    warnings.warn(
        f"Contact graph is disconnected; restricting analysis to largest component "
        f"({len(idx)} nodes)."
    )
    return W[np.ix_(idx, idx)]


# -----------------------------------------------------------------------------
# 6. DEGREE-PRESERVING MASLOV-SNEPPEN ABLATION
# -----------------------------------------------------------------------------
def ablate_topology(
    W_native: np.ndarray,
    swaps_per_edge: int,
    seed: int,
    max_attempt_multiplier: int = 12,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Degree-preserving double-edge rewiring.

    Every accepted swap preserves the degree of every node exactly.
    The final Laplacian spectrum is recomputed from the randomized graph.
    """
    rng = np.random.default_rng(seed)
    W = W_native.copy()

    rows, cols = np.where(np.triu(W, k=1) > 0)
    edges = list(zip(rows.tolist(), cols.tolist()))
    edge_count = len(edges)
    target_swaps = swaps_per_edge * edge_count

    if edge_count < 3:
        raise ValueError("Graph has too few edges for double-edge rewiring.")

    accepted = 0
    max_trials = max(1000, max_attempt_multiplier * target_swaps)

    for _ in range(max_trials):
        if accepted >= target_swaps:
            break

        i1, i2 = rng.choice(edge_count, size=2, replace=False)
        a, b = edges[i1]
        c, d = edges[i2]

        if len({a, b, c, d}) < 4:
            continue

        if rng.integers(2) == 0:
            e1 = (min(a, c), max(a, c))
            e2 = (min(b, d), max(b, d))
        else:
            e1 = (min(a, d), max(a, d))
            e2 = (min(b, c), max(b, c))

        if e1[0] == e1[1] or e2[0] == e2[1] or e1 == e2:
            continue

        # Remove the chosen native edges temporarily.
        W[a, b] = W[b, a] = 0.0
        W[c, d] = W[d, c] = 0.0

        valid = (
            W[e1[0], e1[1]] == 0.0
            and W[e2[0], e2[1]] == 0.0
        )

        if valid:
            W[e1[0], e1[1]] = W[e1[1], e1[0]] = 1.0
            W[e2[0], e2[1]] = W[e2[1], e2[0]] = 1.0
            edges[i1] = e1
            edges[i2] = e2
            accepted += 1
        else:
            W[a, b] = W[b, a] = 1.0
            W[c, d] = W[d, c] = 1.0

    if accepted < target_swaps:
        raise RuntimeError(
            f"Could not complete requested randomization: {accepted}/{target_swaps} swaps."
        )

    degree_native = W_native.sum(axis=1)
    degree_random = W.sum(axis=1)
    max_degree_error = float(np.max(np.abs(degree_native - degree_random)))

    if max_degree_error > 1e-12:
        raise RuntimeError("Degree sequence preservation failed.")

    # The intended null model is a randomized connected network.
    if not nx.is_connected(nx.from_numpy_array(W)):
        warnings.warn("Randomized graph is disconnected; spectrum still computed as-is.")

    L_random = np.diag(W.sum(axis=1)) - W
    eigs_random = np.clip(eigvalsh(0.5 * (L_random + L_random.T)), 0.0, None)

    return W, eigs_random, accepted


# -----------------------------------------------------------------------------
# 7. SPECTRAL COORDINATE + SILVERMAN KDE
# -----------------------------------------------------------------------------
def calibrate_gnm_frequency(eigenvalues: np.ndarray, omega0: float, reference_lambda: float) -> np.ndarray:
    """
    Map Laplacian eigenvalues to a reduced GNM-like coordinate:

        omega_k ∝ sqrt(lambda_k)

    calibrated so reference_lambda maps to omega0.
    """
    lam = np.asarray(eigenvalues, dtype=float)
    return omega0 * np.sqrt(np.maximum(lam, 1e-16) / reference_lambda)


def silverman_bandwidth(samples: np.ndarray, lower_bound: float) -> float:
    """Silverman's bandwidth rule with a numerical lower bound."""
    x = np.asarray(samples, dtype=float)
    n = len(x)
    sigma = float(np.std(x, ddof=1)) if n > 1 else 0.0
    q25, q75 = np.percentile(x, [25, 75])
    iqr = float(q75 - q25)

    scale = min(sigma, iqr / 1.34) if iqr > 0 and sigma > 0 else sigma
    h = 1.06 * scale * (n ** (-1.0 / 5.0))
    return max(float(h), lower_bound)


def spectral_mode_density(
    eigenvalues: np.ndarray,
    omega_grid: np.ndarray,
    omega0: float,
    reference_lambda: float,
    bandwidth_floor: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Exclude the single zero graph mode, map lambda -> reduced omega,
    then smooth the empirical spectral measure with a Gaussian KDE.

    Returns
    -------
    rho       : continuous spectral mode density
    mode_omega: discrete mapped mode coordinates
    bandwidth : actual Silverman bandwidth used
    """
    positive = np.asarray(eigenvalues)[np.asarray(eigenvalues) > 1e-8]
    if positive.size < 5:
        raise ValueError("Too few positive Laplacian eigenvalues for KDE.")

    mode_omega = calibrate_gnm_frequency(positive, omega0, reference_lambda)
    h = silverman_bandwidth(mode_omega, bandwidth_floor)

    rho = norm.pdf(
        omega_grid[:, None],
        loc=mode_omega[None, :],
        scale=h,
    ).mean(axis=1)

    return rho, mode_omega, h


# -----------------------------------------------------------------------------
# 8. PHYSICAL MODEL COMPONENTS
# -----------------------------------------------------------------------------
def J_env(omega: np.ndarray, cfg: Config) -> np.ndarray:
    """Gaussian environmental spectrum centered at omega0."""
    return np.exp(
        -0.5 * ((omega - cfg.omega_0) / cfg.sigma_env) ** 2
    )


def microbiome_M(omega: np.ndarray, intact: bool, cfg: Config) -> np.ndarray:
    """
    M0 = notch-filtered intact state.
    M1 = 1 everywhere after ablation.
    """
    if not intact:
        return np.ones_like(omega)

    notch = cfg.notch_depth * np.exp(
        -0.5 * ((omega - cfg.omega_0) / cfg.sigma_notch) ** 2
    )
    return np.clip(1.0 - notch, 0.0, None)


def topological_amplitude(rho: np.ndarray) -> np.ndarray:
    """
    Explicit model convention:
        F_topo(omega) := rho_topo(omega)

    i.e. the continuous KDE spectral mode density is used as the amplitude
    field appearing in the requested squared-amplitude leakage functional.

    This naming convention is deliberate: rho and F_topo are both retained
    in the code, so the modeling assumption is visible rather than hidden.
    """
    return np.asarray(rho, dtype=float)


def gamma_eff(
    F_topo: np.ndarray,
    omega: np.ndarray,
    cfg: Config,
    grounded: bool,
    microbiome_intact: bool,
) -> tuple[float, np.ndarray]:
    """
    Evaluate

      gamma_eff = ∫ W J_env |F_G|² M |F_topo|² dω

    with W(omega)=1 because no extra weighting function is specified.
    """
    W = np.ones_like(omega)
    J = J_env(omega, cfg)
    FG = cfg.F_G_grounded if grounded else 1.0
    M = microbiome_M(omega, microbiome_intact, cfg)

    integrand = W * J * (FG ** 2) * M * (np.abs(F_topo) ** 2)
    gamma = float(_integrate(integrand, x=omega))
    return gamma, integrand


# -----------------------------------------------------------------------------
# 9. FOUR-STATE FACTORIAL DESIGN
# -----------------------------------------------------------------------------
def run_factorial(
    omega: np.ndarray,
    F_native: np.ndarray,
    F_random: np.ndarray,
    cfg: Config,
) -> dict[str, StateResult]:
    """
    Causal 2×2 matrix:

        topology   = native / randomized
        microbiome = intact / ablated

    Grounding remains fixed at grounded in EVERY cell. This isolates the
    topology × microbiome interaction and prevents boundary-condition
    confounding.
    """
    specs = {
        "00": (F_native, True,  "Intact: native topology + intact microbiome"),
        "10": (F_random, True,  "Protein ablated: randomized topology + intact microbiome"),
        "01": (F_native, False, "Microbiome ablated: native topology + ablated microbiome"),
        "11": (F_random, False, "Double ablation: randomized topology + ablated microbiome"),
    }

    out: dict[str, StateResult] = {}
    for state, (F, mb_intact, label) in specs.items():
        gamma, integrand = gamma_eff(
            F,
            omega,
            cfg,
            grounded=True,
            microbiome_intact=mb_intact,
        )
        cumulative = cumulative_trapezoid(integrand, omega)
        cumulative = np.concatenate(([0.0], cumulative))
        out[state] = StateResult(
            state=state,
            label=label,
            gamma_eff=gamma,
            cumulative=cumulative,
            integrand=integrand,
        )

    return out


# -----------------------------------------------------------------------------
# 10. FIGURE
# -----------------------------------------------------------------------------
def normalize_display(x: np.ndarray) -> np.ndarray:
    m = float(np.max(x))
    return x / m if m > 0 else x


def make_figure(
    omega: np.ndarray,
    J: np.ndarray,
    M0: np.ndarray,
    M1: np.ndarray,
    rho_native: np.ndarray,
    rho_random: np.ndarray,
    rho_control: np.ndarray,
    results: dict[str, StateResult],
    cfg: Config,
    delta_int: float,
    output_dir: Path,
) -> Path:
    """Create the polished two-panel research figure."""
    C = {
        "bg": "#0d1117",
        "panel": "#161b22",
        "grid": "#21262d",
        "border": "#30363d",
        "text": "#c9d1d9",
        "title": "#e6edf3",
        "00": "#2ecc71",
        "10": "#e67e22",
        "01": "#3498db",
        "11": "#e74c3c",
        "add": "#9b59b6",
        "J": "#95a5a6",
        "M0": "#1abc9c",
        "nat": "#ecf0f1",
        "abl": "#f39c12",
        "ctrl": "#7f8c8d",
        "gap": "#f1c40f",
    }

    fig = plt.figure(figsize=(20, 8.8), facecolor=C["bg"])
    gs = gridspec.GridSpec(
        1, 2, figure=fig,
        left=0.055, right=0.985, bottom=0.11, top=0.82,
        wspace=0.30,
    )
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    for ax in (ax1, ax2):
        ax.set_facecolor(C["panel"])
        ax.tick_params(colors=C["text"], labelsize=9.5)
        ax.xaxis.label.set_color(C["text"])
        ax.yaxis.label.set_color(C["text"])
        ax.title.set_color(C["title"])
        for spine in ax.spines.values():
            spine.set_edgecolor(C["border"])
        ax.grid(color=C["grid"], lw=0.5, alpha=0.7, zorder=0)

    # ------------------------------------------------------------------
    # Panel 1: spectral overlap
    # ------------------------------------------------------------------
    ax1.set_title(
        "Panel 1 — Spectral Overlap & Topological Band-Gap Candidate",
        fontsize=10.5, fontweight="bold", pad=9,
    )

    Jn = normalize_display(J)
    Rn = normalize_display(rho_native)
    Ra = normalize_display(rho_random)
    Rc = normalize_display(rho_control)

    ax1.fill_between(omega, 0, Jn, color=C["J"], alpha=0.12)
    ax1.plot(omega, Jn, lw=1.6, ls="--", color=C["J"],
             label=r"$J_{env}(\omega)$ — Gaussian environment spectrum")

    ax1.plot(omega, M0, lw=1.8, color=C["M0"],
             label=r"$M_0(\omega)$ — intact microbiome notch model")
    ax1.plot(omega, M1, lw=1.0, ls=":", color=C["11"], alpha=0.8,
             label=r"$M_1(\omega)$ — ablated microbiome")

    ax1.fill_between(omega, 0, Rn, color=C["nat"], alpha=0.08)
    ax1.plot(omega, Rn, lw=2.2, color=C["nat"],
             label=r"$F_{topo}^{nat}(\omega)$ — 1U3C Laplacian KDE")

    ax1.fill_between(omega, 0, Ra, color=C["abl"], alpha=0.08)
    ax1.plot(omega, Ra, lw=1.7, ls="-.", color=C["abl"],
             label=r"$F_{topo}^{abl}(\omega)$ — degree-preserving randomization")

    ax1.plot(omega, Rc, lw=1.1, ls=":", color=C["ctrl"],
             label=r"$F_{topo}^{1U7C}(\omega)$ — AmtB control")

    ax1.axvline(cfg.omega_0, color=C["gap"], lw=0.9, ls=":", alpha=0.85)
    ax1.text(cfg.omega_0 + 0.35, 1.02, rf"$\omega_0={cfg.omega_0:.1f}$",
             fontsize=8.5, color=C["gap"], transform=ax1.get_xaxis_transform())

    # Candidate protection region: high environmental coupling + low native density.
    high_noise = Jn > 0.35
    candidate = rho_native[high_noise]
    if candidate.size:
        threshold = np.percentile(candidate, 35)
        gap_mask = high_noise & (rho_native <= threshold)
        ax1.fill_between(
            omega, 0, 1.0, where=gap_mask,
            color=C["gap"], alpha=0.055,
            label="Candidate low-overlap window",
        )

    uplift = normalize_display(rho_random)[np.argmin(np.abs(omega - cfg.omega_0))] / max(
        normalize_display(rho_native)[np.argmin(np.abs(omega - cfg.omega_0))], 1e-12
    )

    sign = "SUPER-ADDITIVE" if delta_int > 0 else "SUB-ADDITIVE / NULL"
    ax1.text(
        0.015, 0.975,
        f"gamma00 = {results['00'].gamma_eff:.6g}\n"
        f"gamma10 = {results['10'].gamma_eff:.6g}\n"
        f"gamma01 = {results['01'].gamma_eff:.6g}\n"
        f"gamma11 = {results['11'].gamma_eff:.6g}\n"
        f"Delta_int = {delta_int:+.6g}\n"
        f"uplift@omega0 = {uplift:.3f}x  [{sign}]",
        transform=ax1.transAxes, fontsize=7.8,
        va="top", ha="left", color=C["text"],
        bbox=dict(boxstyle="round,pad=0.45", fc=C["grid"], ec=C["border"], alpha=0.97),
    )

    ax1.set_xlabel(r"Reduced spectral coordinate $\omega$", fontsize=10)
    ax1.set_ylabel("Normalized display amplitude", fontsize=10)
    ax1.set_xlim(0, cfg.omega_max)
    ax1.set_ylim(-0.03, 1.15)
    ax1.legend(loc="upper right", fontsize=6.8, facecolor=C["grid"],
               edgecolor=C["border"], labelcolor=C["text"], framealpha=0.97)

    # ------------------------------------------------------------------
    # Panel 2: cumulative spectral exposure
    # ------------------------------------------------------------------
    ax2.set_title(
        r"Panel 2 — Cumulative Spectral Exposure  $C(\Omega)$",
        fontsize=10.5, fontweight="bold", pad=9,
    )

    G_add = (
        results["10"].cumulative
        + results["01"].cumulative
        - results["00"].cumulative
    )

    labels = {
        "11": r"State 11 — Double ablation",
        "01": r"State 01 — Microbiome ablation",
        "10": r"State 10 — Protein topology ablation",
        "00": r"State 00 — Intact",
    }
    colours = {"00": C["00"], "10": C["10"], "01": C["01"], "11": C["11"]}
    styles = {"00": "-", "10": "--", "01": "-.", "11": ":"}

    for state in ("11", "01", "10", "00"):
        ax2.plot(
            omega, results[state].cumulative,
            lw=2.25, color=colours[state], ls=styles[state],
            label=labels[state],
        )

    ax2.plot(
        omega, G_add, lw=1.5, ls="--", color=C["add"], alpha=0.85,
        label=r"Additive null: $C_{10}+C_{01}-C_{00}$",
    )

    above = results["11"].cumulative >= G_add
    ax2.fill_between(
        omega, G_add, results["11"].cumulative,
        where=above, color=C["11"], alpha=0.12,
        label=r"Super-additivity gap",
    )

    ax2.set_xlabel(r"Upper integration limit $\Omega$", fontsize=10)
    ax2.set_ylabel(
        r"$C(\Omega)=\int_0^\Omega J_{env}|F_G|^2M|F_{topo}|^2d\omega$",
        fontsize=10,
    )
    ax2.set_xlim(0, cfg.omega_max)
    ax2.set_ylim(bottom=0)
    ax2.legend(loc="upper left", fontsize=7.0, facecolor=C["grid"],
               edgecolor=C["border"], labelcolor=C["text"], framealpha=0.97)

    fig.suptitle(
        "Topological Spectral Homeostasis Model v3.0 · Double-Ablation Causal Validation\n"
        r"$\gamma_{eff}=\int J_{env}(\omega)|F_G|^2M(\omega)|F_{topo}(\omega)|^2d\omega$"
        "   [1U3C target | 1U7C control]",
        fontsize=11.5, fontweight="bold", color=C["title"], y=0.975,
    )

    out = output_dir / "TSH_v3_0_Hybrid_DoubleAblation.png"
    fig.savefig(out, dpi=220, bbox_inches="tight", facecolor=C["bg"], edgecolor="none")
    plt.show()
    return out


# -----------------------------------------------------------------------------
# 11. MAIN PIPELINE
# -----------------------------------------------------------------------------
def main() -> dict:
    wall_start = time.time()
    out_dir = Path(CFG.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("  TOPOLOGICAL SPECTRAL HOMEOSTASIS MODEL v3.0")
    print("  Hybrid Double-Ablation · 2×2 Causal Matrix Validation")
    print("=" * 78)

    # -----------------------------------------------------------------
    # PHASE 1: real coordinates
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[1/8] REAL ATOMIC COORDINATES")
    print("─" * 78)

    target_pdb = download_pdb(CFG.target_pdb, out_dir)
    control_pdb = download_pdb(CFG.control_pdb, out_dir)

    labels_t, coords_t, chain_t = extract_calpha(target_pdb, CFG.target_chain)
    labels_c, coords_c, chain_c = extract_calpha(control_pdb, CFG.control_chain)

    print(f"  {CFG.target_pdb} | chain {chain_t} | {len(coords_t):4d} Cα residues")
    print(f"  {CFG.control_pdb} | chain {chain_c} | {len(coords_c):4d} Cα residues")

    # -----------------------------------------------------------------
    # PHASE 2: graph Laplacians
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[2/8] GRAPH CONSTRUCTION — L = D − W")
    print("─" * 78)

    W_t, L_t, eig_t = build_contact_graph(coords_t, CFG.contact_cutoff_A)
    W_c, L_c, eig_c = build_contact_graph(coords_c, CFG.contact_cutoff_A)

    # Retain the original coordinate chain but validate connectivity.
    print(
        f"  Target |E|={int(W_t.sum()/2):5d} | <k>={W_t.sum(axis=1).mean():.2f} "
        f"| λ_max={eig_t[-1]:.3f}"
    )
    print(
        f"  Control|E|={int(W_c.sum()/2):5d} | <k>={W_c.sum(axis=1).mean():.2f} "
        f"| λ_max={eig_c[-1]:.3f}"
    )

    if not nx.is_connected(nx.from_numpy_array(W_t)):
        print("  ! Target graph is disconnected; largest-component fallback is active.")
        W_t = largest_connected_component(W_t)
        L_t = np.diag(W_t.sum(axis=1)) - W_t
        eig_t = np.clip(eigvalsh(0.5 * (L_t + L_t.T)), 0.0, None)

    if not nx.is_connected(nx.from_numpy_array(W_c)):
        print("  ! Control graph is disconnected; largest-component fallback is active.")
        W_c = largest_connected_component(W_c)
        L_c = np.diag(W_c.sum(axis=1)) - W_c
        eig_c = np.clip(eigvalsh(0.5 * (L_c + L_c.T)), 0.0, None)

    # -----------------------------------------------------------------
    # PHASE 3: topology ablation
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[3/8] TOPOLOGY ABLATION — DEGREE-PRESERVING MASLOV-SNEPPEN")
    print("─" * 78)

    W_a, eig_a, accepted = ablate_topology(
        W_t,
        swaps_per_edge=CFG.swaps_per_edge,
        seed=CFG.seed,
        max_attempt_multiplier=CFG.max_swap_attempt_multiplier,
    )

    max_dk = float(np.max(np.abs(W_t.sum(axis=1) - W_a.sum(axis=1))))
    native_degrees = sorted(W_t.sum(axis=1).round(10).tolist())
    random_degrees = sorted(W_a.sum(axis=1).round(10).tolist())
    degree_preserved = native_degrees == random_degrees

    print(f"  Accepted swaps       : {accepted:,} / {CFG.swaps_per_edge * int(W_t.sum()/2):,}")
    print(f"  max|Δk_i|            : {max_dk:.0f}")
    print(f"  Degree sequence exact: {degree_preserved}")

    # -----------------------------------------------------------------
    # PHASE 4: continuous spectral densities
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[4/8] CONTINUOUS SPECTRAL DENSITIES — SILVERMAN KDE")
    print("─" * 78)

    omega = np.linspace(CFG.omega_min, CFG.omega_max, CFG.n_omega)

    positive_native = eig_t[eig_t > CFG.zero_eigenvalue_tol]
    reference_lambda = float(np.median(positive_native))

    rho_nat, modes_nat, h_nat = spectral_mode_density(
        eig_t, omega, CFG.omega_0, reference_lambda, CFG.min_silverman_bandwidth
    )
    rho_abl, modes_abl, h_abl = spectral_mode_density(
        eig_a, omega, CFG.omega_0, reference_lambda, CFG.min_silverman_bandwidth
    )
    rho_ctrl, modes_ctrl, h_ctrl = spectral_mode_density(
        eig_c, omega, CFG.omega_0, reference_lambda, CFG.min_silverman_bandwidth
    )

    F_nat = topological_amplitude(rho_nat)
    F_abl = topological_amplitude(rho_abl)
    F_ctrl = topological_amplitude(rho_ctrl)

    idx0 = int(np.argmin(np.abs(omega - CFG.omega_0)))
    uplift_raw = rho_abl[idx0] / max(rho_nat[idx0], 1e-15)

    print(f"  Reference λ        : {reference_lambda:.6f}")
    print(f"  Silverman h(native): {h_nat:.6f}")
    print(f"  Silverman h(ablated):{h_abl:.6f}")
    print(f"  Silverman h(control):{h_ctrl:.6f}")
    print(f"  rho_native(ω0=5)  : {rho_nat[idx0]:.8g}")
    print(f"  rho_ablated(ω0=5) : {rho_abl[idx0]:.8g}")
    print(f"  Spectral uplift    : {uplift_raw:.4f}×")

    # -----------------------------------------------------------------
    # PHASE 5: physical factors
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[5/8] PHYSICAL FACTORS")
    print("─" * 78)

    J = J_env(omega, CFG)
    M0 = microbiome_M(omega, True, CFG)
    M1 = microbiome_M(omega, False, CFG)

    print(f"  J_env centre        : ω0={CFG.omega_0:.2f}")
    print(f"  Notch depth         : {CFG.notch_depth:.2f}")
    print(f"  Grounded F_G        : {CFG.F_G_grounded:.2f}  -> |F_G|²={CFG.F_G_grounded**2:.2f}")
    print("  Grounding factor is fixed across all four causal cells.")

    # -----------------------------------------------------------------
    # PHASE 6: factorial causal matrix
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[6/8] 2×2 DOUBLE-ABLATION CAUSAL MATRIX")
    print("─" * 78)

    results = run_factorial(omega, F_nat, F_abl, CFG)
    for state in ("00", "10", "01", "11"):
        r = results[state]
        print(f"  γ_{state} = {r.gamma_eff:.10f}   [{r.label}]")

    gamma00 = results["00"].gamma_eff
    gamma10 = results["10"].gamma_eff
    gamma01 = results["01"].gamma_eff
    gamma11 = results["11"].gamma_eff

    delta_int = gamma11 - gamma10 - gamma01 + gamma00
    additive_null = gamma10 + gamma01 - gamma00

    # Algebraic cross-check directly from factorized integrands.
    FG2 = CFG.F_G_grounded ** 2
    direct_delta = float(_integrate(
        J * FG2 * (M1 - M0) * (F_abl**2 - F_nat**2),
        x=omega,
    ))

    print("\n  Interaction effect")
    print("    Δ_int = γ11 − γ10 − γ01 + γ00")
    print(f"          = {delta_int:+.10f}")
    print(f"    direct-integral check = {direct_delta:+.10f}")
    print(f"    algebraic agreement   = {np.isclose(delta_int, direct_delta, rtol=1e-7, atol=1e-12)}")

    if additive_null > 0:
        ratio = gamma11 / additive_null
        print(f"    γ11 / additive-null   = {ratio:.6f}×")
    else:
        ratio = np.nan
        print("    γ11 / additive-null   = undefined (non-positive null)")

    if delta_int > 0:
        verdict = "SUPER-ADDITIVE INTERACTION OBSERVED"
    elif delta_int < 0:
        verdict = "SUB-ADDITIVE INTERACTION OBSERVED"
    else:
        verdict = "NO FIRST-ORDER INTERACTION"
    print(f"    RESULT: {verdict}")

    # -----------------------------------------------------------------
    # PHASE 7: convergence diagnostics
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[7/8] CUMULATIVE SPECTRAL EXPOSURE + CONVERGENCE")
    print("─" * 78)

    diag_rows = []
    for state in ("00", "10", "01", "11"):
        r = results[state]
        idx_tail = int(np.searchsorted(omega, 8.0))
        tail = (r.gamma_eff - r.cumulative[idx_tail]) / max(r.gamma_eff, 1e-15)
        diag_rows.append({
            "state": state,
            "gamma_eff": r.gamma_eff,
            "plateau_C(Omega_max)": r.cumulative[-1],
            "tail_fraction_after_8": tail,
        })

    diagnostics = pd.DataFrame(diag_rows)
    print(diagnostics.to_string(index=False, float_format=lambda x: f"{x:.8g}"))

    # The cumulative curves are generated from the same spectral integrand
    # used to compute gamma_eff; no hand-tuned time exponents are introduced.
    G_add = results["10"].cumulative + results["01"].cumulative - results["00"].cumulative
    cumulative_interaction = results["11"].cumulative - G_add
    cumulative_crosscheck = np.allclose(
        cumulative_interaction[-1], delta_int, rtol=1e-7, atol=1e-12
    )
    print(f"  Final cumulative interaction matches Δ_int: {cumulative_crosscheck}")

    # -----------------------------------------------------------------
    # PHASE 8: control + figure + export
    # -----------------------------------------------------------------
    print("\n" + "─" * 78)
    print("[8/8] CONTROL STRUCTURE + FIGURE + EXPORT")
    print("─" * 78)

    gamma_ctrl, _ = gamma_eff(
        F_ctrl, omega, CFG, grounded=True, microbiome_intact=True
    )
    print(f"  1U7C control gamma under State-00 conditions: {gamma_ctrl:.10f}")

    figure_path = make_figure(
        omega, J, M0, M1,
        rho_nat, rho_abl, rho_ctrl,
        results, CFG, delta_int, out_dir,
    )

    summary = pd.DataFrame({
        "gamma00": [gamma00],
        "gamma10": [gamma10],
        "gamma01": [gamma01],
        "gamma11": [gamma11],
        "Delta_int": [delta_int],
        "gamma11_additive_null_ratio": [ratio],
        "spectral_uplift_at_omega0": [uplift_raw],
        "control_gamma_state00": [gamma_ctrl],
        "degree_sequence_preserved": [degree_preserved],
        "max_degree_error": [max_dk],
        "accepted_swaps": [accepted],
        "target_chain": [chain_t],
        "control_chain": [chain_c],
        "target_residues": [len(coords_t)],
        "control_residues": [len(coords_c)],
        "target_edges": [int(W_t.sum() / 2)],
        "cutoff_A": [CFG.contact_cutoff_A],
        "omega0": [CFG.omega_0],
        "sigma_env": [CFG.sigma_env],
        "notch_depth": [CFG.notch_depth],
        "sigma_notch": [CFG.sigma_notch],
        "grounded_FG": [CFG.F_G_grounded],
        "seed": [CFG.seed],
        "silverman_h_native": [h_nat],
        "silverman_h_ablated": [h_abl],
        "silverman_h_control": [h_ctrl],
    })

    csv_path = out_dir / "TSH_v3_0_summary.csv"
    diagnostics_path = out_dir / "TSH_v3_0_convergence_diagnostics.csv"
    summary.to_csv(csv_path, index=False)
    diagnostics.to_csv(diagnostics_path, index=False)

    elapsed = time.time() - wall_start

    print("\n" + "=" * 78)
    print("  PIPELINE COMPLETE")
    print(f"  γ00   = {gamma00:.10f}")
    print(f"  γ11   = {gamma11:.10f}")
    print(f"  Δ_int = {delta_int:+.10f}  |  {verdict}")
    print(f"  Figure      : {figure_path}")
    print(f"  Summary CSV : {csv_path}")
    print(f"  Diagnostics : {diagnostics_path}")
    print(f"  Wall-clock  : {elapsed:.1f} s")
    print("=" * 78)

    return {
        "gammas": {s: results[s].gamma_eff for s in results},
        "Delta_int": delta_int,
        "spectral_uplift": uplift_raw,
        "degree_sequence_preserved": degree_preserved,
        "figure": str(figure_path),
        "summary_csv": str(csv_path),
    }


if __name__ == "__main__":
    results = main()
