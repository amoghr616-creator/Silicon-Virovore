from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path.cwd()
INPUT = ROOT / (
    "prototype3/results/sequence_only_mutation_baseline/"
    "run_20261026T234823Z/single_substitution_pairs.csv"
)
OUT = ROOT / "prototype3/results/arise_experimental_position_audit"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20261027
BOOTSTRAPS = 2000
PERMUTATIONS = 5000
rng = np.random.default_rng(SEED)


def add_positions(df):
    positions = []
    lengths = []

    for a, b in zip(df["sequence_a"], df["sequence_b"]):
        a, b = str(a).strip().upper(), str(b).strip().upper()
        if len(a) != len(b):
            raise ValueError(f"Unequal peptide lengths: {a!r}, {b!r}")
        differences = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        if len(differences) != 1:
            raise ValueError(
                f"Expected one substitution, got {len(differences)}: {a!r}, {b!r}"
            )
        positions.append(differences[0] + 1)
        lengths.append(len(a))

    df["mutation_position_1based"] = positions
    df["peptide_length"] = lengths
    df["abs_delta_pKd"] = df["delta_pKd"].abs()
    df["informative"] = df["abs_delta_pKd"] >= 0.3
    return df


def residualize(frame, cluster_column="target_group_id"):
    """Center effect magnitude within target, assay/context, and peptide length."""
    keys = [cluster_column, "context_key", "peptide_length"]
    out = frame.copy()

    grouped = out.groupby(keys, dropna=False)
    n_rows = grouped["abs_delta_pKd"].transform("size")
    n_positions = grouped["mutation_position_1based"].transform("nunique")

    # Only retain contexts in which at least two mutation positions were tested.
    out = out[(n_rows >= 2) & (n_positions >= 2)].copy()
    out["centered_abs_effect"] = (
        out["abs_delta_pKd"]
        - out.groupby(keys, dropna=False)["abs_delta_pKd"].transform("mean")
    )
    return out


def window_stat(residualized_frame, width, starts_by_length):
    """Difference in centered absolute affinity effect: window minus elsewhere."""
    d = residualized_frame
    starts = d["peptide_length"].map(starts_by_length).to_numpy()
    pos = d["mutation_position_1based"].to_numpy()
    hot = (pos >= starts) & (pos < starts + width)
    y = d["centered_abs_effect"].to_numpy(dtype=float)

    if not hot.any() or (~hot).sum() == 0:
        return None
    return float(y[hot].mean() - y[~hot].mean())


def random_window_starts(lengths, width, generator):
    starts = {}
    for length in lengths:
        if length < width:
            raise ValueError(f"Window width {width} exceeds peptide length {length}")
        # One-based starts: 1 through length-width+1 inclusive.
        starts[length] = int(generator.integers(1, length - width + 2))
    return starts


def analyze_window(frame, label, start, width, informative_only=False):
    d = frame.copy()
    if informative_only:
        d = d[d["informative"]].copy()

    # The matching and centering are recalculated in each bootstrap replicate.
    base = residualize(d)
    lengths = sorted(base["peptide_length"].unique())
    starts = {length: start for length in lengths}

    if base.empty or len(lengths) == 0:
        return {"label": label, "status": "insufficient data"}

    observed = window_stat(base, width, starts)
    if observed is None:
        return {"label": label, "status": "insufficient window coverage"}

    # Descriptive unadjusted effects, retained for interpretability.
    in_window = (
        (d["mutation_position_1based"] >= start)
        & (d["mutation_position_1based"] < start + width)
    )
    raw_hot = d.loc[in_window, "abs_delta_pKd"]
    raw_other = d.loc[~in_window, "abs_delta_pKd"]

    # Permutation control: shift the same-width contiguous window independently
    # for each peptide length, keeping the measured mutation effects fixed.
    null = []
    for _ in range(PERMUTATIONS):
        random_starts = random_window_starts(lengths, width, rng)
        value = window_stat(base, width, random_starts)
        if value is not None and np.isfinite(value):
            null.append(value)

    p_greater = (
        (1 + sum(x >= observed for x in null)) / (1 + len(null))
        if null else None
    )
    p_two_sided = (
        (1 + sum(abs(x) >= abs(observed) for x in null)) / (1 + len(null))
        if null else None
    )

    # Cluster bootstrap over target sequences.
    original_targets = sorted(d["target_group_id"].astype(str).unique())
    bootstrap_values = []

    for _ in range(BOOTSTRAPS):
        selected = rng.choice(
            original_targets, size=len(original_targets), replace=True
        )
        chunks = []
        for occurrence, target in enumerate(selected):
            part = d[d["target_group_id"].astype(str) == target].copy()
            part["bootstrap_target"] = occurrence
            chunks.append(part)

        boot = pd.concat(chunks, ignore_index=True)
        boot = residualize(boot, cluster_column="bootstrap_target")
        boot_lengths = sorted(boot["peptide_length"].unique())
        boot_starts = {length: start for length in boot_lengths}

        if boot.empty:
            continue
        value = window_stat(boot, width, boot_starts)
        if value is not None and np.isfinite(value):
            bootstrap_values.append(value)

    ci = (
        [
            float(np.quantile(bootstrap_values, 0.025)),
            float(np.quantile(bootstrap_values, 0.975)),
        ]
        if len(bootstrap_values) >= 100 else None
    )

    return {
        "label": label,
        "positions_1based": list(range(start, start + width)),
        "informative_pairs_only": informative_only,
        "pairs_total": int(len(d)),
        "pairs_in_window": int(in_window.sum()),
        "pairs_outside_window": int((~in_window).sum()),
        "contexts_after_matching_filter": int(
            base.groupby(
                ["target_group_id", "context_key", "peptide_length"],
                dropna=False
            ).ngroups
        ),
        "unique_target_sequences": int(d["target_group_id"].nunique()),
        "raw_mean_abs_delta_pKd_in_window": float(raw_hot.mean()),
        "raw_mean_abs_delta_pKd_outside_window": float(raw_other.mean()),
        "raw_median_abs_delta_pKd_in_window": float(raw_hot.median()),
        "raw_median_abs_delta_pKd_outside_window": float(raw_other.median()),
        "centered_effect_difference_window_minus_elsewhere": observed,
        "target_cluster_bootstrap_95ci": ci,
        "bootstrap_valid_replicates": len(bootstrap_values),
        "random_contiguous_window_permutation_p_greater": p_greater,
        "random_contiguous_window_permutation_p_two_sided": p_two_sided,
        "permutation_valid_replicates": len(null),
        "permutation_count_requested": PERMUTATIONS,
        "bootstrap_count_requested": BOOTSTRAPS,
        "interpretation": (
            "Positive centered difference means mutations in the selected "
            "window have larger absolute measured affinity effects than the "
            "within-context expectation. This is retrospective and exploratory."
        ),
    }


if not INPUT.exists():
    raise SystemExit(f"Input not found: {INPUT}")

data = pd.read_csv(INPUT)
required = {
    "target_group_id", "context_key", "sequence_a",
    "sequence_b", "delta_pKd",
}
missing = required - set(data.columns)
if missing:
    raise SystemExit(f"Missing columns: {sorted(missing)}")

data["delta_pKd"] = pd.to_numeric(data["delta_pKd"], errors="coerce")
data = data.dropna(
    subset=["target_group_id", "context_key", "delta_pKd"]
).copy()
data["target_group_id"] = data["target_group_id"].astype(str)
data = add_positions(data)

# Primary prediction: latest logged hotspot [0,3) zero-based = positions 1–3.
primary = analyze_window(
    data,
    label="latest_logged_hotspot_positions_1_to_3",
    start=1,
    width=3,
)

# Sensitivity check: previous generation's [1,3) zero-based = positions 2–3.
previous = analyze_window(
    data,
    label="previous_generation_hotspot_positions_2_to_3",
    start=2,
    width=2,
)

# Secondary analysis restricted to pairs with |delta pKd| >= 0.3.
informative_primary = analyze_window(
    data,
    label="latest_hotspot_positions_1_to_3_informative_pairs",
    start=1,
    width=3,
    informative_only=True,
)

result = {
    "benchmark": "Retrospective PPIKB position-window validation",
    "input": str(INPUT),
    "n_single_substitution_pairs": int(len(data)),
    "n_informative_pairs_abs_delta_pKd_ge_0_3": int(data["informative"].sum()),
    "hotspot_coordinate_mapping": (
        "The logs specify zero-based start and end-exclusive indices; "
        "[0,3) maps to one-based positions 1,2,3; [1,3) maps to positions 2,3."
    ),
    "primary_window_test": primary,
    "previous_generation_window_sensitivity": previous,
    "informative_subset_sensitivity": informative_primary,
    "limitations": [
        "The experimental dataset was explored in prior analyses; this is not an untouched prospective test.",
        "PPIKB peptides have different targets and sequences from the 21-aa search seed.",
        "This tests a general peptide-position window, not HERV-K Env-specific residue importance.",
        "It does not evaluate target-aware structural contacts or prove that ARISE guided the mutations in the logged smoke run.",
        "The permutation test compares the observed window with random contiguous windows of equal width, stratified by peptide length.",
    ],
}

(OUT / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
(OUT / "REPORT.md").write_text(
    "# ARISE hotspot position-window audit\n\n"
    "This is a retrospective exploratory comparison of the logged hotspot "
    "positions with experimental mutation-affinity effects. It tests whether "
    "the predicted relative positions are associated with larger experimental "
    "effects than random equal-width contiguous windows. It does not establish "
    "HERV-K-specific hotspot importance or causal benefit from adaptive search.\n\n"
    "```json\n" + json.dumps(result, indent=2) + "\n```\n"
)

print(json.dumps(result, indent=2))
print("\nSaved:", OUT)
