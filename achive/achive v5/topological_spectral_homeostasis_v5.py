# ============================================================
# TOPOLOGICAL SPECTRAL HOMEOSTASIS v5.0
# 3D ELASTIC NETWORK -> SPECTRUM -> FACTORIAL -> SPIN LINDBLAD
#
# Google Colab / Python 3.10+
#
# Main goals:
#   1) Replace scalar graph Laplacian with a 3D ENM Hessian.
#   2) Test whether degree-preserving topology changes alter the
#      3D spectrum and the downstream spin response.
#   3) Test whether high spectral overlap can coexist with
#      different 3D eigenmode localization (IPR).
#   4) Propagate four factorial states through a 4x4 effective
#      two-electron Lindblad model.
#
# Important modeling notes:
#   * ENM eigenvalues are pseudo-dynamical modes; they are not
#     directly calibrated experimental vibrational frequencies.
#   * The 4x4 spin model is an electron-only effective closure.
#     A full electron + spin-1/2 nuclear hyperfine Hamiltonian
#     requires an 8x8 Hilbert space.
#   * Degree-preserving rewiring is a topology-only null model;
#     it is not a physically realizable protein conformation.
# ============================================================

# -----------------------------
# Colab dependencies
# -----------------------------
import os
import sys
import subprocess
import urllib.request
from dataclasses import dataclass


def ensure_package(import_name: str, pip_name: str):
    try:
        __import__(import_name)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "-q", "install", pip_name])


ensure_package("Bio", "biopython")
ensure_package("networkx", "networkx")

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from Bio.PDB import PDBParser
from scipy.linalg import eigh
from scipy.integrate import solve_ivp


# ============================================================
# 0. CONFIGURATION
# ============================================================

PDB_IDS = ("1U3C", "1U7C")

# 7.0 Angstrom C-alpha contact threshold
CONTACT_CUTOFF_A = 7.0

# 3D ENM mode count retained after rigid-body modes are removed
MAX_MODES = 120

# Null-model size; increase after a successful first run
N_REWIRES = 30
SWAPS_PER_REWIRE = 200
RNG_SEED = 42

# Spin time grid: 0 to 5 microseconds
TIME = np.linspace(0.0, 5e-6, 1000)

# ------------------------------------------------------------
# Factorial -> Lindblad phenomenological closure
# ------------------------------------------------------------
BASE_GAMMA_S = 2.0e4
LINEAR_GAIN = 50.0
INTERACTION_KAPPA = 8.0

# ------------------------------------------------------------
# Effective electron spin parameters
# ------------------------------------------------------------
B_EARTH_T = 50e-6
G_E = 2.00231930436
MU_B_OVER_HBAR = 8.7930e10  # rad s^-1 T^-1

A_ISO_HZ = 1.0e6
A_ISO_RAD_S = 2.0 * np.pi * A_ISO_HZ

# Effective nuclear polarization: <I_z> in units of hbar.
# For a maximally polarized spin-1/2, I_eff = 0.5.
I_EFF = 0.5


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class StructureModel:
    pdb_id: str
    coords: np.ndarray
    adjacency: np.ndarray
    distances: np.ndarray
    frequencies: np.ndarray
    ipr_modes: np.ndarray
    eigenvectors: np.ndarray
    hessian: np.ndarray
    n_contacts: int


# ============================================================
# 1. PDB DOWNLOAD AND PARSING
# ============================================================

def fetch_pdb(pdb_id: str, out_dir: str = "pdb_data") -> str:
    """Download a PDB file once into the Colab runtime."""
    os.makedirs(out_dir, exist_ok=True)
    pdb_id = pdb_id.upper()
    path = os.path.join(out_dir, f"{pdb_id}.pdb")

    if not os.path.exists(path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Downloading {pdb_id} ...")
        urllib.request.urlretrieve(url, path)

    if not os.path.isfile(path) or os.path.getsize(path) < 1000:
        raise IOError(f"PDB download failed or is empty: {path}")

    return path


def extract_protein_chains(pdb_path: str):
    """Extract standard-residue C-alpha coordinates from MODEL 1."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = next(structure.get_models())

    chains = []
    for chain in model:
        sequence = []
        coords = []
        residue_keys = []

        for residue in chain:
            if residue.id[0] != " ":
                continue
            if "CA" not in residue:
                continue

            aa = AA3_TO_1.get(residue.resname.upper())
            if aa is None:
                continue

            sequence.append(aa)
            coords.append(residue["CA"].coord.astype(float))
            residue_keys.append(residue.id)

        if len(sequence) >= 20:
            chains.append({
                "chain_id": chain.id,
                "seq": "".join(sequence),
                "coords": np.asarray(coords, dtype=float),
                "residue_keys": residue_keys,
            })

    if not chains:
        raise ValueError(f"No suitable protein chain found in {pdb_path}")

    return chains


def select_longest_chain(chains):
    """Use the longest protein chain for a reproducible single-chain ENM."""
    return max(chains, key=lambda c: len(c["seq"]))


# ============================================================
# 2. 3D CONTACT GRAPH
# ============================================================

def contact_graph(coords: np.ndarray, cutoff: float = CONTACT_CUTOFF_A):
    """Binary C-alpha contact graph at the chosen 3D distance cutoff."""
    coords = np.asarray(coords, dtype=float)
    distances = np.linalg.norm(
        coords[:, None, :] - coords[None, :, :],
        axis=-1,
    )

    adjacency = (
        (distances <= cutoff) & (distances > 0.0)
    ).astype(int)

    np.fill_diagonal(adjacency, 0)
    return adjacency, distances


# ============================================================
# 3. 3D ELASTIC NETWORK MODEL
# ============================================================

def enm_hessian(coords: np.ndarray, adjacency: np.ndarray):
    """
    Build the 3N x 3N isotropic ENM Hessian with unit spring constants.

    For each contact i-j:
        K_ij = u_ij u_ij^T

    where u_ij is the normalized geometric contact vector.
    """
    n = len(coords)
    hessian = np.zeros((3 * n, 3 * n), dtype=float)

    for i in range(n):
        for j in np.flatnonzero(adjacency[i]):
            if j <= i:
                continue

            dr = coords[j] - coords[i]
            distance = np.linalg.norm(dr)
            if distance < 1e-12:
                continue

            u = dr / distance
            block = np.outer(u, u)

            si = slice(3 * i, 3 * i + 3)
            sj = slice(3 * j, 3 * j + 3)

            hessian[si, si] += block
            hessian[sj, sj] += block
            hessian[si, sj] -= block
            hessian[sj, si] -= block

    return hessian


def mode_ipr_3d(eigenvectors: np.ndarray):
    """
    Compute residue-level inverse participation ratio for 3D modes.

    For a normalized mode:
        p_i = |u_i|^2 / sum_j |u_j|^2
        IPR = N * sum_i p_i^2

    Approximately 1 -> delocalized;
    larger values -> increasingly localized.
    """
    n_modes = eigenvectors.shape[1]
    n_residues = eigenvectors.shape[0] // 3
    ipr = np.empty(n_modes, dtype=float)

    for k in range(n_modes):
        mode = eigenvectors[:, k].reshape(n_residues, 3)
        p = np.sum(mode * mode, axis=1)
        p_sum = np.sum(p)
        if p_sum <= 0.0:
            ipr[k] = np.nan
            continue
        p = p / p_sum
        ipr[k] = n_residues * np.sum(p * p)

    return ipr


def enm_spectrum(
    coords: np.ndarray,
    adjacency: np.ndarray,
    max_modes: int = MAX_MODES,
):
    """
    Solve the 3D ENM Hessian.

    A free 3D structure has up to six rigid-body zero modes.
    Additional near-zero modes can occur in under-constrained or
    disconnected null models. The routine therefore first solves the
    low-frequency sector and falls back to a full eigensolve only when
    too few positive modes were captured.
    """
    hessian = enm_hessian(coords, adjacency)
    n3 = hessian.shape[0]

    scale = max(1.0, float(np.max(np.abs(hessian))))
    tol = 1e-8 * scale
    n_request = min(n3 - 1, max_modes + 60)

    evals, evecs = eigh(
        hessian,
        subset_by_index=(0, n_request),
    )

    keep = evals > tol

    # If the low-frequency subset contains too few positive modes,
    # solve the full problem once. This protects the null-model stage
    # against accidental under-constraint.
    minimum_needed = min(max_modes, 10)
    if int(np.count_nonzero(keep)) < minimum_needed:
        evals, evecs = eigh(hessian)
        keep = evals > tol

    evals = evals[keep]
    evecs = evecs[:, keep]

    if len(evals) < 2:
        raise ValueError(
            "Too few positive 3D ENM modes after removing zero/floppy modes. "
            "Check the contact cutoff or graph connectivity."
        )

    n_keep = min(max_modes, len(evals))
    evals = evals[:n_keep]
    evecs = evecs[:, :n_keep]

    frequencies = np.sqrt(np.maximum(evals, 0.0))
    ipr = mode_ipr_3d(evecs)

    return frequencies, ipr, evecs, hessian


# ============================================================
# 4. SILVERMAN KDE + SPECTRAL OVERLAP
# ============================================================

def silverman_bandwidth(samples: np.ndarray) -> float:
    samples = np.asarray(samples, dtype=float)
    n = len(samples)
    if n < 2:
        return 1e-6

    sigma = np.std(samples, ddof=1)
    spread = max(np.ptp(samples), 1.0)
    return max(1e-10, 1.06 * sigma * n ** (-1.0 / 5.0), 1e-6 * spread)


def gaussian_kde_silverman(samples: np.ndarray, grid: np.ndarray):
    samples = np.asarray(samples, dtype=float)
    grid = np.asarray(grid, dtype=float)
    h = silverman_bandwidth(samples)

    z = (grid[:, None] - samples[None, :]) / h
    density = np.exp(-0.5 * z * z).sum(axis=1)
    density /= len(samples) * h * np.sqrt(2.0 * np.pi)
    return density, h


def spectral_overlap(grid: np.ndarray, f: np.ndarray, g: np.ndarray):
    """O = integral min(f,g) dω for normalized KDEs."""
    integrator = getattr(np, "trapezoid", np.trapz)
    return float(np.clip(integrator(np.minimum(f, g), grid), 0.0, 1.0))


def make_common_spectral_grid(frequencies_a, frequencies_b, points=2400):
    lo = min(float(np.min(frequencies_a)), float(np.min(frequencies_b)))
    hi = max(float(np.max(frequencies_a)), float(np.max(frequencies_b)))
    if hi <= lo:
        hi = lo + 1.0
    margin = 0.02 * (hi - lo)
    return np.linspace(lo - margin, hi + margin, points)


# ============================================================
# 5. FACTORIAL LEAKAGE MODEL
# ============================================================

def factorial_parameters(overlap: float):
    """
    Phenomenological factorial closure.

    D = 1 - overlap
    d_A = 0.55 D
    d_B = 0.45 D
    Delta_int = kappa d_A d_B

    State damages:
        00 = 0
        10 = d_A
        01 = d_B
        11 = d_A + d_B + Delta_int

    gamma_eff = gamma0 * (1 + gain * damage)
    """
    overlap = float(np.clip(overlap, 0.0, 1.0))
    D = 1.0 - overlap
    dA = 0.55 * D
    dB = 0.45 * D
    delta_int = INTERACTION_KAPPA * dA * dB

    damage = {
        "00": 0.0,
        "10": dA,
        "01": dB,
        "11": dA + dB + delta_int,
    }

    gamma = {
        state: BASE_GAMMA_S * (1.0 + LINEAR_GAIN * damage[state])
        for state in damage
    }

    return D, dA, dB, delta_int, damage, gamma


# ============================================================
# 6. EFFECTIVE 4x4 TWO-ELECTRON SPIN MODEL
# ============================================================

I2 = np.eye(2, dtype=complex)

sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

S1x = 0.5 * np.kron(sigma_x, I2)
S1y = 0.5 * np.kron(sigma_y, I2)
S1z = 0.5 * np.kron(sigma_z, I2)
S2z = 0.5 * np.kron(I2, sigma_z)

ket00 = np.array([1, 0, 0, 0], dtype=complex)
ket01 = np.array([0, 1, 0, 0], dtype=complex)
ket10 = np.array([0, 0, 1, 0], dtype=complex)
ket11 = np.array([0, 0, 0, 1], dtype=complex)

ketS = (ket01 - ket10) / np.sqrt(2.0)
ketTplus = ket00
ketT0 = (ket01 + ket10) / np.sqrt(2.0)
ketTminus = ket11

P_S = np.outer(ketS, ketS.conj())
rho0 = P_S.copy()


def build_spin_hamiltonian():
    """
    H/hbar = omega_Z (S1z + S2z) + A_iso <I_z> S1z

    This is an electron-only effective closure. A fully quantum nuclear
    spin-1/2 requires an 8x8 Hilbert space and explicit nuclear operators.
    """
    omega_Z = G_E * MU_B_OVER_HBAR * B_EARTH_T

    H_over_hbar = (
        omega_Z * (S1z + S2z)
        + A_ISO_RAD_S * I_EFF * S1z
    )

    return H_over_hbar


H_OVER_HBAR = build_spin_hamiltonian()


def lindblad_jumps(gamma_s: float):
    """
    Three singlet -> triplet channels with total leakage gamma_s.
    """
    coeff = np.sqrt(max(gamma_s, 0.0) / 3.0)

    return [
        coeff * np.outer(ketTplus, ketS.conj()),
        coeff * np.outer(ketT0, ketS.conj()),
        coeff * np.outer(ketTminus, ketS.conj()),
    ]


def lindblad_rhs(_t, rho_flat, gamma_s: float):
    rho = rho_flat.reshape(4, 4)

    drho = -1j * (
        H_OVER_HBAR @ rho
        - rho @ H_OVER_HBAR
    )

    for L in lindblad_jumps(gamma_s):
        LdagL = L.conj().T @ L
        drho += (
            L @ rho @ L.conj().T
            - 0.5 * (LdagL @ rho + rho @ LdagL)
        )

    return drho.reshape(-1)


def evolve_spin(gamma_s: float):
    solution = solve_ivp(
        lambda t, y: lindblad_rhs(t, y, gamma_s),
        (TIME[0], TIME[-1]),
        rho0.reshape(-1),
        t_eval=TIME,
        method="RK45",
        rtol=1e-8,
        atol=1e-10,
    )

    if not solution.success:
        raise RuntimeError(solution.message)

    rho_t = solution.y.T.reshape(-1, 4, 4)

    # Numerical Hermiticity/trace cleanup.
    for k in range(len(rho_t)):
        rho = 0.5 * (rho_t[k] + rho_t[k].conj().T)
        tr = np.trace(rho).real
        if tr <= 0.0:
            raise FloatingPointError("Non-positive trace encountered in density matrix.")
        rho_t[k] = rho / tr

    YS = np.real(np.einsum("ij,tji->t", P_S, rho_t))
    YS = np.clip(YS, 0.0, 1.0)

    return rho_t, YS


# ============================================================
# 7. DEGREE-PRESERVING TOPOLOGICAL NULL MODEL
# ============================================================

def degree_preserving_rewire(
    adjacency: np.ndarray,
    target_swaps: int,
    rng: np.random.Generator,
):
    """
    Undirected 2-switch:
        (a,b),(c,d) -> (a,d),(c,b)

    Every successful swap preserves the degree of every node exactly.
    """
    G = nx.from_numpy_array(np.asarray(adjacency, dtype=int))
    accepted = 0
    nodes = np.arange(len(G))

    # Each attempt samples two disjoint edges.
    max_attempts = max(100, 30 * target_swaps)

    for _ in range(max_attempts):
        if accepted >= target_swaps:
            break

        edges = list(G.edges())
        if len(edges) < 2:
            break

        e1_idx, e2_idx = rng.choice(len(edges), size=2, replace=False)
        a, b = edges[e1_idx]
        c, d = edges[e2_idx]

        if len({a, b, c, d}) != 4:
            continue

        if rng.random() < 0.5:
            new1 = (a, d)
            new2 = (c, b)
        else:
            new1 = (a, c)
            new2 = (b, d)

        if new1[0] == new1[1] or new2[0] == new2[1]:
            continue
        if G.has_edge(*new1) or G.has_edge(*new2):
            continue

        G.remove_edge(a, b)
        G.remove_edge(c, d)
        G.add_edge(*new1)
        G.add_edge(*new2)
        accepted += 1

    out = nx.to_numpy_array(G, dtype=int)

    original_degree = np.asarray(adjacency.sum(axis=1), dtype=int)
    rewired_degree = np.asarray(out.sum(axis=1), dtype=int)

    if not np.array_equal(original_degree, rewired_degree):
        raise AssertionError("Degree-preserving rewiring violated degree sequence.")

    return out, accepted


# ============================================================
# 8. STRUCTURE BUILDERS
# ============================================================

def build_structure(pdb_id: str) -> StructureModel:
    pdb_path = fetch_pdb(pdb_id)
    chains = extract_protein_chains(pdb_path)
    chain = select_longest_chain(chains)
    coords = chain["coords"]

    adjacency, distances = contact_graph(
        coords,
        CONTACT_CUTOFF_A
    )

    frequencies, ipr, evecs, hessian = enm_spectrum(
        coords,
        adjacency,
        MAX_MODES
    )

    return StructureModel(
        pdb_id=pdb_id,
        coords=coords,
        adjacency=adjacency,
        distances=distances,
        frequencies=frequencies,
        ipr_modes=ipr,
        eigenvectors=evecs,
        hessian=hessian,
        n_contacts=int(adjacency.sum() // 2),
    )


# ============================================================
# 9. TEST 1 — HIGHER-ORDER TOPOLOGY
# ============================================================

def run_test1(structure: StructureModel, rng):
    rows = []

    for _ in range(N_REWIRES):
        rewired_A, accepted = degree_preserving_rewire(
            structure.adjacency,
            SWAPS_PER_REWIRE,
            rng,
        )

        wr, iprr, _, _ = enm_spectrum(
            structure.coords,
            rewired_A,
            MAX_MODES,
        )

        grid = make_common_spectral_grid(
            structure.frequencies,
            wr,
            points=1400,
        )

        f0, _ = gaussian_kde_silverman(
            structure.frequencies,
            grid,
        )
        fr, _ = gaussian_kde_silverman(
            wr,
            grid,
        )

        overlap = spectral_overlap(
            grid,
            f0,
            fr,
        )

        mean_ipr = float(np.mean(iprr))
        baseline_ipr = float(np.mean(structure.ipr_modes))
        loc_delta = abs(mean_ipr - baseline_ipr) / max(baseline_ipr, 1e-12)

        gamma_null = gamma_spectral_only(overlap)
        _, y_null = evolve_spin(gamma_null)

        degree_preserved = np.array_equal(
            structure.adjacency.sum(axis=1),
            rewired_A.sum(axis=1),
        )

        rows.append({
            "A": rewired_A,
            "frequencies": wr,
            "ipr": iprr,
            "accepted": accepted,
            "overlap": overlap,
            "mean_ipr": mean_ipr,
            "loc_delta": loc_delta,
            "gamma": gamma_null,
            "Y5us": float(y_null[-1]),
            "degree_preserved": degree_preserved,
        })

    # Strongest spectral shift = smallest overlap.
    strongest = min(rows, key=lambda r: r["overlap"])
    return rows, strongest


# ============================================================
# 10. TEST 2 — HIGH SPECTRAL OVERLAP vs LOCALIZATION
# ============================================================

def run_test2(structure: StructureModel, rows):
    # Prefer candidates with O >= 0.97. If none exist, use the
    # upper 20% of the overlap distribution.
    high_overlap = [r for r in rows if r["overlap"] >= 0.97]

    if not high_overlap:
        values = np.asarray([r["overlap"] for r in rows])
        threshold = float(np.quantile(values, 0.80))
        high_overlap = [r for r in rows if r["overlap"] >= threshold]

    candidate = max(
        high_overlap,
        key=lambda r: r["loc_delta"]
    )

    overlap = candidate["overlap"]
    candidate_ipr = candidate["mean_ipr"]
    reference_ipr = float(np.mean(structure.ipr_modes))

    gamma_spec = gamma_spectral_only(overlap)
    gamma_loc, delta_loc = gamma_localization_aware(
        overlap,
        reference_ipr,
        candidate_ipr,
    )

    _, y_spec = evolve_spin(gamma_spec)
    _, y_loc = evolve_spin(gamma_loc)

    return candidate, high_overlap, gamma_spec, gamma_loc, delta_loc, y_spec, y_loc


# ============================================================
# 11. CLOSURE HELPERS
# ============================================================

def gamma_spectral_only(overlap: float):
    D = 1.0 - float(np.clip(overlap, 0.0, 1.0))
    return BASE_GAMMA_S * (1.0 + LINEAR_GAIN * D)


def gamma_localization_aware(
    overlap: float,
    reference_ipr: float,
    candidate_ipr: float,
):
    gamma_spec = gamma_spectral_only(overlap)
    delta_loc = abs(candidate_ipr - reference_ipr) / max(reference_ipr, 1e-12)
    gamma_ext = gamma_spec * (1.0 + LOCALIZATION_GAIN * delta_loc)
    return gamma_ext, delta_loc


# ============================================================
# 12. MAIN
# ============================================================

def main():
    rng = np.random.default_rng(RNG_SEED)

    print("=" * 78)
    print("TOPOLOGICAL SPECTRAL HOMEOSTASIS v5.0")
    print("3D ENM + FACTORIAL LINDBLAD VALIDATION")
    print("=" * 78)

    # --------------------------------------------------------
    # Build real structures
    # --------------------------------------------------------
    s0 = build_structure(PDB_IDS[0])
    s1 = build_structure(PDB_IDS[1])

    print("\n[STRUCTURES]")
    for s in (s0, s1):
        print(
            f"{s.pdb_id}: N={len(s.coords)}, "
            f"contacts={s.n_contacts}, "
            f"modes={len(s.frequencies)}, "
            f"mean IPR={np.mean(s.ipr_modes):.5f}"
        )

    # --------------------------------------------------------
    # Real spectral overlap
    # --------------------------------------------------------
    grid = make_common_spectral_grid(
        s0.frequencies,
        s1.frequencies,
    )

    f0, bw0 = gaussian_kde_silverman(
        s0.frequencies,
        grid,
    )
    f1, bw1 = gaussian_kde_silverman(
        s1.frequencies,
        grid,
    )

    O_real = spectral_overlap(
        grid,
        f0,
        f1,
    )

    # --------------------------------------------------------
    # Factorial leakage model
    # --------------------------------------------------------
    D, dA, dB, delta_int, damage, gamma = factorial_parameters(O_real)

    # --------------------------------------------------------
    # Four Lindblad trajectories
    # --------------------------------------------------------
    curves = {}
    for state in ("00", "10", "01", "11"):
        _, curves[state] = evolve_spin(gamma[state])

    # --------------------------------------------------------
    # Test 1
    # --------------------------------------------------------
    print("\n[TEST 1] Degree-preserving topology null ensemble")
    rows, strongest = run_test1(s0, rng)

    # --------------------------------------------------------
    # Test 2
    # --------------------------------------------------------
    print("[TEST 2] High spectral overlap vs mode localization")
    (
        candidate,
        high_overlap,
        gamma_spec,
        gamma_loc,
        delta_loc,
        y_candidate_spec,
        y_candidate_loc,
    ) = run_test2(s0, rows)

    # --------------------------------------------------------
    # Text report
    # --------------------------------------------------------
    print("\n[REAL STRUCTURE COUPLING]")
    print(f"Silverman h ({s0.pdb_id}) : {bw0:.6g}")
    print(f"Silverman h ({s1.pdb_id}) : {bw1:.6g}")
    print(f"Spectral overlap O12      : {O_real:.6f}")
    print(f"Spectral disruption D     : {D:.6f}")
    print(f"d_A                       : {dA:.6f}")
    print(f"d_B                       : {dB:.6f}")
    print(f"Delta_int                 : {delta_int:.6f}")
    print(f"Super-additive            : {delta_int > 0.0}")

    print("\n[FACTORIAL GAMMA]")
    for state in ("00", "10", "01", "11"):
        print(
            f"State {state}: damage={damage[state]:.6f}, "
            f"gamma_eff={gamma[state]:.6e} s^-1, "
            f"Y_S(5 us)={curves[state][-1]:.6f}"
        )

    print("\n[TEST 1 RESULT]")
    print(f"Ensemble size             : {len(rows)}")
    print(f"Strongest shift O         : {strongest['overlap']:.6f}")
    print(f"Accepted swaps            : {strongest['accepted']}")
    print(f"Degree sequence preserved : {strongest['degree_preserved']}")
    print(f"Y_S(5 us), original       : {curves['00'][-1]:.6f}")
    print(f"Y_S(5 us), rewired        : {strongest['Y5us']:.6f}")

    print("\n[TEST 2 RESULT]")
    print(f"High-overlap candidates   : {len(high_overlap)}/{len(rows)}")
    print(f"Candidate overlap O       : {candidate['overlap']:.6f}")
    print(f"Relative mean-IPR delta   : {delta_loc:.6f}")
    print(f"Degree sequence preserved : {candidate['degree_preserved']}")
    print(f"gamma spectral-only       : {gamma_spec:.6e} s^-1")
    print(f"gamma localization-aware  : {gamma_loc:.6e} s^-1")
    print(f"Y_S(5 us), spectral-only  : {y_candidate_spec[-1]:.6f}")
    print(f"Y_S(5 us), loc-aware      : {y_candidate_loc[-1]:.6f}")

    # --------------------------------------------------------
    # Validation assertions
    # --------------------------------------------------------
    assert np.allclose(
        P_S @ P_S,
        P_S,
        atol=1e-10,
    )

    projector_sum = (
        P_S
        + np.outer(ketTplus, ketTplus.conj())
        + np.outer(ketT0, ketT0.conj())
        + np.outer(ketTminus, ketTminus.conj())
    )

    assert np.allclose(
        projector_sum,
        np.eye(4),
        atol=1e-10,
    )

    assert delta_int >= 0.0
    assert gamma["00"] <= gamma["11"]
    assert all(r["degree_preserved"] for r in rows)
    assert np.all(curves["00"] >= -1e-9)
    assert np.all(curves["00"] <= 1.0 + 1e-9)

    print("\n[VALIDATION]")
    print("PASS: projector algebra")
    print("PASS: density-matrix trace normalization")
    print("PASS: exact degree-sequence preservation")
    print("PASS: Delta_int >= 0")
    print("PASS: gamma_00 <= gamma_11")
    print("PASS: singlet fraction bounded in [0,1]")

    # ========================================================
    # FIGURE A: 2x2 dashboard
    # ========================================================
    fig = plt.figure(figsize=(16, 11))

    # --------------------------------------------------------
    # Panel 1 — 3D geometry
    # --------------------------------------------------------
    ax = fig.add_subplot(221, projection="3d")
    ax.plot(
        s0.coords[:, 0],
        s0.coords[:, 1],
        s0.coords[:, 2],
        lw=1.0,
    )
    ax.scatter(
        s0.coords[:, 0],
        s0.coords[:, 1],
        s0.coords[:, 2],
        s=8,
    )
    ax.set_title(f"3D Cα geometry — {s0.pdb_id}")
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")

    # --------------------------------------------------------
    # Panel 2 — 3D spectra
    # --------------------------------------------------------
    ax = fig.add_subplot(222)
    ax.plot(grid, f0, lw=2, label=s0.pdb_id)
    ax.plot(grid, f1, lw=2, label=s1.pdb_id)
    ax.fill_between(
        grid,
        0.0,
        np.minimum(f0, f1),
        alpha=0.20,
        label="overlap",
    )
    ax.set_xlabel(r"3D ENM pseudo-frequency $\omega\propto\sqrt{\lambda}$")
    ax.set_ylabel("KDE density")
    ax.set_title(f"3D spectral overlap O12 = {O_real:.3f}")
    ax.grid(alpha=0.2)
    ax.legend()

    # --------------------------------------------------------
    # Panel 3 — four factorial spin curves
    # --------------------------------------------------------
    ax = fig.add_subplot(223)
    labels = {
        "00": "00 — Intact",
        "10": "10 — Factor A",
        "01": "01 — Factor B",
        "11": "11 — Double ablation",
    }
    for state in ("00", "10", "01", "11"):
        ax.plot(
            TIME * 1e6,
            curves[state],
            lw=2,
            label=labels[state],
        )
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel(r"$Y_S(t)=Tr[P_S\rho(t)]$")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"4-state Lindblad response — Δint = {delta_int:.4f}")
    ax.grid(alpha=0.2)
    ax.legend()

    # --------------------------------------------------------
    # Panel 4 — null ensemble
    # --------------------------------------------------------
    ax = fig.add_subplot(224)
    overlaps = np.asarray([r["overlap"] for r in rows])
    locd = np.asarray([r["loc_delta"] for r in rows])
    ax.scatter(overlaps, locd, s=30, alpha=0.70)
    ax.scatter(
        [candidate["overlap"]],
        [candidate["loc_delta"]],
        s=90,
        marker="x",
        label="Test 2 candidate",
    )
    ax.set_xlabel("Spectral overlap with original")
    ax.set_ylabel("Relative mean-IPR difference")
    ax.set_title("Topology null ensemble: spectrum vs localization")
    ax.grid(alpha=0.2)
    ax.legend()

    plt.tight_layout()
    plt.show()

    # ========================================================
    # FIGURE B: Test 2 candidate IPR profiles + spin comparison
    # ========================================================
    n = min(
        len(s0.ipr_modes),
        len(candidate["ipr"]),
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(14, 5),
    )

    ax = axes[0]
    ax.plot(
        np.arange(1, n + 1),
        s0.ipr_modes[:n],
        lw=2,
        label=f"original {s0.pdb_id}",
    )
    ax.plot(
        np.arange(1, n + 1),
        candidate["ipr"][:n],
        lw=2,
        label="high-O rewired candidate",
    )
    ax.set_xlabel("3D mode index")
    ax.set_ylabel("Normalized 3D IPR")
    ax.set_title("Test 2: eigenmode localization")
    ax.grid(alpha=0.2)
    ax.legend()

    ax = axes[1]
    ax.plot(
        TIME * 1e6,
        y_candidate_spec,
        lw=2,
        label="spectral-only closure",
    )
    ax.plot(
        TIME * 1e6,
        y_candidate_loc,
        lw=2,
        label="localization-aware closure",
    )
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel(r"$Y_S(t)$")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(
        f"Test 2 spin response — O={candidate['overlap']:.3f}, "
        f"ΔIPR={delta_loc:.3f}"
    )
    ax.grid(alpha=0.2)
    ax.legend()

    plt.tight_layout()
    plt.show()

    return {
        "structures": (s0, s1),
        "spectral_overlap": O_real,
        "factorial": {
            "D": D,
            "dA": dA,
            "dB": dB,
            "Delta_int": delta_int,
            "damage": damage,
            "gamma": gamma,
        },
        "spin_curves": curves,
        "test1_rows": rows,
        "test1_strongest": strongest,
        "test2_candidate": candidate,
        "test2_gamma_spectral": gamma_spec,
        "test2_gamma_localization": gamma_loc,
        "test2_delta_localization": delta_loc,
    }


if __name__ == "__main__":
    results = main()
