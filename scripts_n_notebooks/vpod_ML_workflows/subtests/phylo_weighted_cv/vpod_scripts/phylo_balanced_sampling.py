"""
Utilities for building VPOD training sets balanced by broad phylogenetic group,
lambda_max coverage, and local phylogenetic redundancy.

The intended use is upstream of model fitting:

1. Load the metadata and sequence matrix as in ``main_vpod_wf.ipynb``.
2. Call ``make_phylo_phenotype_balanced_subset`` on the metadata.
3. Subset ``tr`` and ``y`` to the returned selected sequence IDs.
4. Train/evaluate/finalize the shallow models exactly as before.

Rows with ``Phylum == 'Chordata'`` are treated as vertebrate. Every other
phylum is treated as invertebrate.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
import pandas as pd


GroupName = Literal["vertebrate", "invertebrate"]
RangeMode = Literal["overlap", "full"]
BinMode = Literal["equal_width", "quantile"]


@dataclass(frozen=True)
class BalanceConfig:
    phenotype_col: str = "Lambda_Max"
    phylum_col: str = "Phylum"
    vertebrate_phylum: str = "Chordata"
    n_bins: int = 12
    range_mode: RangeMode = "overlap"
    bin_mode: BinMode = "equal_width"
    min_per_group_per_bin: int = 1
    max_per_group_per_bin: int | None = None
    phylo_percentile: float = 5.0
    random_state: int = 1


@dataclass
class BalancedSubset:
    selected_ids: list[str]
    selected_meta: pd.DataFrame
    annotated_meta: pd.DataFrame
    bin_summary: pd.DataFrame
    selection_summary: dict[str, object]


def _metadata_indexed_by_sequence(meta_data: pd.DataFrame) -> pd.DataFrame:
    """Return metadata indexed by sequence ID, matching deepBreaks behavior."""
    meta = meta_data.copy()
    if "Seq_Id" in meta.columns:
        meta = meta.set_index("Seq_Id", drop=False)
    meta.index = meta.index.astype(str)
    return meta


def annotate_vertebrate_status(
    meta_data: pd.DataFrame,
    phylum_col: str = "Phylum",
    vertebrate_phylum: str = "Chordata",
) -> pd.DataFrame:
    """Add a two-level vertebrate/invertebrate grouping column."""
    meta = _metadata_indexed_by_sequence(meta_data)
    if phylum_col not in meta.columns:
        raise KeyError(f"Could not find phylum column: {phylum_col!r}")

    phyla = meta[phylum_col].astype(str).str.strip()
    meta["vpod_balance_group"] = np.where(
        phyla.str.casefold() == vertebrate_phylum.casefold(),
        "vertebrate",
        "invertebrate",
    )
    return meta


def _pairwise_tree_distances(tree_file: str | Path) -> tuple[pd.DataFrame, float]:
    """Read a Newick tree and return all pairwise terminal distances."""
    try:
        from Bio import Phylo
    except ImportError as exc:
        raise ImportError("Biopython is required when tree_file is supplied.") from exc

    tree = Phylo.read(str(tree_file), "newick")
    terminals = tree.get_terminals()
    tip_names = [str(term.name) for term in terminals]
    dist = np.zeros((len(tip_names), len(tip_names)), dtype=float)

    for i, terminal_i in enumerate(terminals):
        for j in range(i + 1, len(terminals)):
            d = tree.distance(terminal_i, terminals[j])
            dist[i, j] = d
            dist[j, i] = d

    return pd.DataFrame(dist, index=tip_names, columns=tip_names), float(np.nanmax(dist))


def _connected_components_from_threshold(
    dist: pd.DataFrame,
    threshold: float,
) -> pd.Series:
    """Cluster tree tips that are connected by short phylogenetic distances."""
    ids = dist.index.to_list()
    parent = {seq_id: seq_id for seq_id in ids}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        root_a = find(a)
        root_b = find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    values = dist.to_numpy()
    for i, id_i in enumerate(ids):
        for j in range(i + 1, len(ids)):
            if values[i, j] <= threshold:
                union(id_i, ids[j])

    roots = [find(seq_id) for seq_id in ids]
    root_to_label = {root: f"phylo_cluster_{i:04d}" for i, root in enumerate(sorted(set(roots)))}
    return pd.Series([root_to_label[root] for root in roots], index=ids, name="phylo_cluster")


def add_phylo_clusters(
    meta_data: pd.DataFrame,
    tree_file: str | Path | None,
    percentile: float = 5.0,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """
    Add ``phylo_cluster`` labels.

    The cluster threshold is the requested percentile of all non-zero pairwise
    tree distances. Tips not present in the supplied metadata are ignored; rows
    not present in the tree are kept as singleton clusters.
    """
    meta = _metadata_indexed_by_sequence(meta_data)
    if tree_file is None:
        meta["phylo_cluster"] = [f"singleton::{idx}" for idx in meta.index]
        return meta, {
            "tree_file": None,
            "phylo_threshold": None,
            "tree_matched_rows": 0,
            "tree_missing_rows": int(meta.shape[0]),
        }

    dist, max_dist = _pairwise_tree_distances(tree_file)
    shared = [seq_id for seq_id in meta.index if seq_id in dist.index]
    if len(shared) < 2:
        raise ValueError(
            "Fewer than two metadata rows matched tree tips. Check that FASTA, "
            "metadata index/Seq_Id, and tree tip names use the same IDs."
        )

    shared_dist = dist.loc[shared, shared]
    upper = shared_dist.to_numpy()[np.triu_indices(len(shared), k=1)]
    upper = upper[np.isfinite(upper)]
    upper = upper[upper > 0]
    threshold = float(np.percentile(upper, percentile)) if len(upper) else max_dist
    clusters = _connected_components_from_threshold(shared_dist, threshold)

    meta["phylo_cluster"] = [f"singleton::{idx}" for idx in meta.index]
    meta.loc[clusters.index, "phylo_cluster"] = clusters
    return meta, {
        "tree_file": str(tree_file),
        "phylo_threshold": threshold,
        "tree_matched_rows": int(len(shared)),
        "tree_missing_rows": int(meta.shape[0] - len(shared)),
        "phylo_percentile": float(percentile),
    }


def _phenotype_edges(
    meta: pd.DataFrame,
    phenotype_col: str,
    n_bins: int,
    range_mode: RangeMode,
    bin_mode: BinMode,
) -> np.ndarray:
    """Build lambda_max bin edges from either shared or full phenotype range."""
    if n_bins < 1:
        raise ValueError("n_bins must be >= 1")

    values = pd.to_numeric(meta[phenotype_col], errors="coerce")
    if range_mode == "overlap":
        grouped = meta.assign(_phenotype=values).groupby("vpod_balance_group")["_phenotype"]
        min_edge = grouped.min().max()
        max_edge = grouped.max().min()
    elif range_mode == "full":
        min_edge = values.min()
        max_edge = values.max()
    else:
        raise ValueError("range_mode must be 'overlap' or 'full'")

    if not np.isfinite(min_edge) or not np.isfinite(max_edge) or min_edge >= max_edge:
        raise ValueError(
            f"Could not construct a valid phenotype range from {phenotype_col!r}: "
            f"{min_edge} to {max_edge}"
        )

    in_range = values[(values >= min_edge) & (values <= max_edge)]
    
    # equal_width = equal lambda_max intervals; better range coverage; may discard more data
    # quantile    = roughly equal data density per bin; keeps more data; less strict phenotype-range balance
    if bin_mode == "equal_width":
        edges = np.linspace(min_edge, max_edge, n_bins + 1)
    elif bin_mode == "quantile":
        edges = np.quantile(in_range.to_numpy(dtype=float), np.linspace(0, 1, n_bins + 1))
        edges = np.unique(edges)
        if len(edges) < 2:
            raise ValueError("Quantile binning collapsed to fewer than two unique edges.")
    else:
        raise ValueError("bin_mode must be 'equal_width' or 'quantile'")

    edges[0] = np.nextafter(edges[0], -np.inf)
    edges[-1] = np.nextafter(edges[-1], np.inf)
    return edges


def _add_phenotype_bins(
    meta: pd.DataFrame,
    phenotype_col: str,
    edges: np.ndarray,
) -> pd.DataFrame:
    """Add numeric phenotype and interval-bin labels."""
    out = meta.copy()
    out["vpod_balance_phenotype"] = pd.to_numeric(out[phenotype_col], errors="coerce")
    out["vpod_balance_bin"] = pd.cut(
        out["vpod_balance_phenotype"],
        bins=edges,
        include_lowest=True,
        duplicates="drop",
    )
    return out


def _inverse_cluster_weights(frame: pd.DataFrame) -> np.ndarray:
    """Weight candidates inversely to local phylogenetic redundancy."""
    cluster_sizes = frame.groupby("phylo_cluster")["phylo_cluster"].transform("size")
    weights = 1.0 / cluster_sizes.astype(float).to_numpy()
    total = weights.sum()
    if not np.isfinite(total) or total <= 0:
        return np.full(frame.shape[0], 1.0 / frame.shape[0])
    return weights / total


def _weighted_sample_ids(
    frame: pd.DataFrame,
    n: int,
    rng: np.random.Generator,
) -> list[str]:
    if n > frame.shape[0]:
        raise ValueError(f"Cannot sample {n} rows from only {frame.shape[0]} candidates.")
    if n == frame.shape[0]:
        return frame.index.astype(str).to_list()

    probs = _inverse_cluster_weights(frame)
    chosen = rng.choice(frame.index.astype(str).to_numpy(), size=n, replace=False, p=probs)
    return chosen.tolist()


def summarize_bins(annotated_meta: pd.DataFrame) -> pd.DataFrame:
    """Return group/bin counts before sampling."""
    usable = annotated_meta.dropna(subset=["vpod_balance_bin"])
    summary = (
        usable.groupby(["vpod_balance_bin", "vpod_balance_group"], observed=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for col in ["vertebrate", "invertebrate"]:
        if col not in summary.columns:
            summary[col] = 0
    summary["balanced_n_per_group"] = summary[["vertebrate", "invertebrate"]].min(axis=1)
    summary["bin_left"] = summary["vpod_balance_bin"].map(lambda x: float(x.left))
    summary["bin_right"] = summary["vpod_balance_bin"].map(lambda x: float(x.right))
    return summary[
        ["vpod_balance_bin", "bin_left", "bin_right", "vertebrate", "invertebrate", "balanced_n_per_group"]
    ]


def make_phylo_phenotype_balanced_subset(
    meta_data: pd.DataFrame,
    tree_file: str | Path | None = None,
    config: BalanceConfig | None = None,
) -> BalancedSubset:
    """
    Select a final-model training subset balanced between vertebrates and invertebrates.

    Balancing happens in three layers:

    * Broad group balance: each retained lambda_max bin contributes the same number
      of vertebrate and invertebrate rows.
    * Phenotype balance: bins span the overlapping phenotype range by default, so
      both groups cover the same lambda_max support.
    * Phylogenetic balance: sampling within each group/bin is weighted inversely
      to local tree-cluster size, reducing oversampling of dense clades.
    """
    cfg = config or BalanceConfig()
    rng = np.random.default_rng(cfg.random_state)

    meta = annotate_vertebrate_status(
        meta_data,
        phylum_col=cfg.phylum_col,
        vertebrate_phylum=cfg.vertebrate_phylum,
    )
    if cfg.phenotype_col not in meta.columns:
        raise KeyError(f"Could not find phenotype column: {cfg.phenotype_col!r}")

    meta, phylo_summary = add_phylo_clusters(meta, tree_file, percentile=cfg.phylo_percentile)
    edges = _phenotype_edges(
        meta,
        phenotype_col=cfg.phenotype_col,
        n_bins=cfg.n_bins,
        range_mode=cfg.range_mode,
        bin_mode=cfg.bin_mode,
    )
    meta = _add_phenotype_bins(meta, cfg.phenotype_col, edges)
    bin_summary = summarize_bins(meta)

    selected: list[str] = []
    skipped_bins: list[str] = []
    usable = meta.dropna(subset=["vpod_balance_bin"]).copy()
    for bin_label, bin_frame in usable.groupby("vpod_balance_bin", observed=True):
        group_frames = {
            "vertebrate": bin_frame[bin_frame["vpod_balance_group"] == "vertebrate"],
            "invertebrate": bin_frame[bin_frame["vpod_balance_group"] == "invertebrate"],
        }
        n_per_group = min(group_frames["vertebrate"].shape[0], group_frames["invertebrate"].shape[0])
        if cfg.max_per_group_per_bin is not None:
            n_per_group = min(n_per_group, cfg.max_per_group_per_bin)

        if n_per_group < cfg.min_per_group_per_bin:
            skipped_bins.append(str(bin_label))
            continue
        
        # The sampling step happens here: 
        # So if a bin has 40 vertebrates and 6 invertebrates, 
        # the balanced bin keeps 6 vertebrates + 6 invertebrates.
        # This uses rng to uniformly random sample from the bins if no phylo tree is given,
        # otherwise the sampling would be weighted by the tree distance matrix to prevent heavy sampling of one clade within a bin
        selected.extend(_weighted_sample_ids(group_frames["vertebrate"], n_per_group, rng))
        selected.extend(_weighted_sample_ids(group_frames["invertebrate"], n_per_group, rng))

    selected = list(dict.fromkeys(selected))
    selected_meta = meta.loc[selected].copy()
    selection_summary: dict[str, object] = {
        **phylo_summary,
        "n_selected": int(len(selected)),
        "n_selected_vertebrate": int((selected_meta["vpod_balance_group"] == "vertebrate").sum()),
        "n_selected_invertebrate": int((selected_meta["vpod_balance_group"] == "invertebrate").sum()),
        "n_bins_requested": int(cfg.n_bins),
        "n_bins_used": int(selected_meta["vpod_balance_bin"].nunique()),
        "skipped_bins": skipped_bins,
        "phenotype_min_selected": float(selected_meta["vpod_balance_phenotype"].min()) if len(selected_meta) else np.nan,
        "phenotype_max_selected": float(selected_meta["vpod_balance_phenotype"].max()) if len(selected_meta) else np.nan,
        "range_mode": cfg.range_mode,
        "bin_mode": cfg.bin_mode,
        "random_state": int(cfg.random_state),
    }

    return BalancedSubset(
        selected_ids=selected,
        selected_meta=selected_meta,
        annotated_meta=meta,
        bin_summary=bin_summary,
        selection_summary=selection_summary,
    )


def subset_training_matrix(
    tr: pd.DataFrame,
    meta_data: pd.DataFrame,
    selected_ids: Iterable[str],
    phenotype_col: str = "Lambda_Max",
    use_ev: bool = False,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    """
    Align a deepBreaks sequence matrix and metadata to selected IDs.

    Returns ``X_balanced, y_balanced, selected_meta``.
    """
    ids = [str(seq_id) for seq_id in selected_ids]
    meta = _metadata_indexed_by_sequence(meta_data)
    shared = [seq_id for seq_id in ids if seq_id in tr.index and seq_id in meta.index]
    missing = sorted(set(ids) - set(shared))
    if missing:
        raise KeyError(
            f"{len(missing)} selected IDs were missing from tr or metadata. "
            f"First missing IDs: {missing[:10]}"
        )

    X = tr.loc[shared].copy()
    y = pd.to_numeric(meta.loc[shared, phenotype_col], errors="raise").to_numpy(dtype=float)
    if use_ev:
        y = 1239.8 / y
    return X, y, meta.loc[shared].copy()


def write_balanced_subset_report(
    result: BalancedSubset,
    out_dir: str | Path,
    prefix: str = "phylo_phenotype_balanced",
) -> dict[str, str]:
    """Write auditable subset files and return their paths."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    selected_meta_path = out_path / f"{prefix}_metadata.tsv"
    annotated_meta_path = out_path / f"{prefix}_all_rows_annotated.tsv"
    ids_path = out_path / f"{prefix}_ids.txt"
    bin_summary_path = out_path / f"{prefix}_bin_summary.tsv"
    summary_path = out_path / f"{prefix}_summary.tsv"

    result.selected_meta.to_csv(selected_meta_path, sep="\t", index=True)
    result.annotated_meta.to_csv(annotated_meta_path, sep="\t", index=True)
    result.bin_summary.to_csv(bin_summary_path, sep="\t", index=False)
    pd.Series(result.selected_ids).to_csv(ids_path, index=False, header=False)
    pd.Series(result.selection_summary, name="value").to_csv(summary_path, sep="\t")

    return {
        "selected_meta": str(selected_meta_path),
        "annotated_meta": str(annotated_meta_path),
        "selected_ids": str(ids_path),
        "bin_summary": str(bin_summary_path),
        "summary": str(summary_path),
    }
