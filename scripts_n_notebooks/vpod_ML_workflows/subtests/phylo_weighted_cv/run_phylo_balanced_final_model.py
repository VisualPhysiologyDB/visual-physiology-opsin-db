"""
Create a vertebrate/invertebrate phylo-phenotype-balanced VPOD training set,
optionally then train/finalize shallow deepBreaks models on that subset.

Example, balance-only audit:

python run_phylo_balanced_final_model.py ^
  --metaDataFileName ..\\..\\vpod_1.3_data_splits_2025-10-06_16-50-06\\wds_meta.tsv ^
  --tree_file trees\\vpod_1.2\\vpod_1.2_wt\\wt_vpod_1.2_LG_F_R7.treefile ^
  --output_dir balanced_wds_audit ^
  --n_bins 12

Example, train final models:

python run_phylo_balanced_final_model.py ^
  --seqFileName ..\\..\\vpod_1.3_data_splits_2025-10-06_16-50-06\\wds_aligned_VPOD_1.3_het.fasta ^
  --metaDataFileName ..\\..\\vpod_1.3_data_splits_2025-10-06_16-50-06\\wds_meta.tsv ^
  --tree_file path\\to\\matching_alignment.treefile ^
  --output_dir balanced_wds_training ^
  --train_final_model
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import pandas as pd

from vpod_scripts.phylo_balanced_sampling import (
    BalanceConfig,
    make_phylo_phenotype_balanced_subset,
    subset_training_matrix,
    write_balanced_subset_report,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a VPOD final-model training subset balanced by clade, lambda_max, and tree redundancy."
    )
    parser.add_argument("--seqFileName", default=None, help="Aligned FASTA used by deepBreaks. Required for training.")
    parser.add_argument("--metaDataFileName", required=True, help="Metadata TSV with Seq_Id, Phylum, and Lambda_Max.")
    parser.add_argument("--tree_file", default=None, help="Newick/IQ-TREE .treefile matching the FASTA IDs.")
    parser.add_argument("--output_dir", default=None, help="Directory for reports and optional model outputs.")

    parser.add_argument("--mt", default="Lambda_Max", help="Phenotype column.")
    parser.add_argument("--phylum_col", default="Phylum", help="Phylum column.")
    parser.add_argument("--vertebrate_phylum", default="Chordata", help="Phylum treated as vertebrate.")
    parser.add_argument("--n_bins", type=int, default=12, help="Number of lambda_max bins.")
    parser.add_argument(
        "--range_mode",
        choices=["overlap", "full"],
        default="overlap",
        help="Use shared vertebrate/invertebrate phenotype range, or full observed range.",
    )
    parser.add_argument(
        "--bin_mode",
        choices=["equal_width", "quantile"],
        default="equal_width",
        help="Equal-width bins preserve phenotype-range coverage; quantile bins preserve row count.",
    )
    parser.add_argument("--min_per_group_per_bin", type=int, default=1)
    parser.add_argument("--max_per_group_per_bin", type=int, default=None)
    parser.add_argument("--phylo_percentile", type=float, default=5.0)
    parser.add_argument("--random_state", type=int, default=1)

    parser.add_argument("--train_final_model", action="store_true", help="Train/finalize deepBreaks models after balancing.")
    parser.add_argument("--seq_type", default="aa")
    parser.add_argument("--ana_type", default="reg")
    parser.add_argument("--gap_threshold", type=float, default=0.5)
    parser.add_argument("--encoding_method", choices=["aa_prop", "hot"], default="aa_prop")
    parser.add_argument("--dataset", default="wds", help="deepBreaks dataset key for optimized params.")
    parser.add_argument("--ignore_dataset_params", action="store_true", help="Use generic deepBreaks model params.")
    parser.add_argument("--use_ev", action="store_true", help="Train on eV instead of nm.")
    parser.add_argument(
        "--props_to_keep",
        nargs="+",
        default=["best"],
        help="Amino-acid properties, or 'best' to use deepBreaks get_best_aa_prop_combos(dataset).",
    )
    parser.add_argument("--cv", type=int, default=10, help="CV folds for model_compare_cv/finalize_top.")
    return parser.parse_args()


def _make_output_dir(args: argparse.Namespace) -> Path:
    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        stamp = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        out_dir = Path(f"phylo_phenotype_balanced_{args.mt}_{stamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def _make_preprocess_pipeline(
    encoding_method: str,
    props_to_keep: list[str],
    ana_type: str,
    dataset: str,
) :
    from deepBreaks.preprocessing import (
        AminoAcidPropertyEncoder,
        CollinearCare,
        ConstantCare,
        CustomOneHotEncoder,
        FeatureSelection,
        MisCare,
        URareCare,
    )
    from deepBreaks.utils_alt2 import get_best_aa_prop_combos, make_pipeline

    if encoding_method == "aa_prop":
        if props_to_keep == ["best"]:
            props_to_keep = get_best_aa_prop_combos(dataset)
        return make_pipeline(
            steps=[
                ("mc", MisCare(missing_threshold=0.05)),
                ("cc", ConstantCare()),
                ("aa_prop", AminoAcidPropertyEncoder(props_to_keep=props_to_keep)),
                ("feature_selection", FeatureSelection(model_type=ana_type, alpha=0.10, keep=True)),
                ("collinear_care", CollinearCare(dist_method="correlation", threshold=0.01, keep=True)),
            ]
        )

    return make_pipeline(
        steps=[
            ("mc", MisCare(missing_threshold=0.05)),
            ("cc", ConstantCare()),
            ("ur", URareCare(threshold=0.025)),
            ("cc2", ConstantCare()),
            ("one_hot", CustomOneHotEncoder()),
            ("feature_selection", FeatureSelection(model_type=ana_type, alpha=0.10, keep=True)),
            ("collinear_care", CollinearCare(dist_method="correlation", threshold=0.01, keep=True)),
        ]
    )


def train_final_models(args: argparse.Namespace, out_dir: Path, selected_ids: list[str]) -> None:
    if not args.seqFileName:
        raise ValueError("--seqFileName is required with --train_final_model")

    from deepBreaks.models import finalize_top, model_compare_cv
    from deepBreaks.preprocessing import read_data, write_fasta
    from deepBreaks.utils_alt2 import get_models, get_params, get_scores, make_pipeline

    dataset = "ignore" if args.ignore_dataset_params else args.dataset
    meta_data = read_data(args.metaDataFileName, seq_type=None, is_main=False)
    tr = read_data(
        args.seqFileName,
        seq_type=args.seq_type,
        is_main=True,
        gap_threshold=args.gap_threshold,
    )

    X_balanced, y_balanced, selected_meta = subset_training_matrix(
        tr=tr,
        meta_data=meta_data,
        selected_ids=selected_ids,
        phenotype_col=args.mt,
        use_ev=args.use_ev,
    )
    selected_meta.to_csv(out_dir / "deepbreaks_training_metadata.tsv", sep="\t", index=True)
    write_fasta(dat=X_balanced, fasta_file="deepbreaks_training_gap_dropped.fasta", report_dir=str(out_dir))

    prep_pipeline = _make_preprocess_pipeline(
        encoding_method=args.encoding_method,
        props_to_keep=args.props_to_keep,
        ana_type=args.ana_type,
        dataset=args.dataset,
    )

    report, top = model_compare_cv(
        X=X_balanced,
        y=y_balanced,
        preprocess_pipe=prep_pipeline,
        models_dict=get_models(ana_type=args.ana_type, encoding=args.encoding_method, dataset=dataset),
        scoring=get_scores(ana_type=args.ana_type),
        report_dir=str(out_dir),
        cv=args.cv,
        ana_type=args.ana_type,
        cache_dir=str(out_dir),
    )
    report.to_csv(out_dir / "balanced_model_compare_report.tsv", sep="\t", index=True)

    final_prep = _make_preprocess_pipeline(
        encoding_method=args.encoding_method,
        props_to_keep=args.props_to_keep,
        ana_type=args.ana_type,
        dataset=args.dataset,
    )
    top_models = [make_pipeline(steps=[("prep", final_prep), model.steps[-1]]) for model in top]
    finalize_top(
        X=X_balanced,
        y=y_balanced,
        top_models=top_models,
        grid_param=get_params(),
        report_dir=str(out_dir),
        cv=args.cv,
    )


def main() -> None:
    args = parse_args()
    out_dir = _make_output_dir(args)

    meta_data = pd.read_csv(args.metaDataFileName, sep="\t")
    config = BalanceConfig(
        phenotype_col=args.mt,
        phylum_col=args.phylum_col,
        vertebrate_phylum=args.vertebrate_phylum,
        n_bins=args.n_bins,
        range_mode=args.range_mode,
        bin_mode=args.bin_mode,
        min_per_group_per_bin=args.min_per_group_per_bin,
        max_per_group_per_bin=args.max_per_group_per_bin,
        phylo_percentile=args.phylo_percentile,
        random_state=args.random_state,
    )
    result = make_phylo_phenotype_balanced_subset(
        meta_data=meta_data,
        tree_file=args.tree_file,
        config=config,
    )
    paths = write_balanced_subset_report(result, out_dir=out_dir)

    print("Balanced subset created")
    print(f"Output directory: {out_dir}")
    print(f"Selected rows: {result.selection_summary['n_selected']}")
    print(f"Vertebrate rows: {result.selection_summary['n_selected_vertebrate']}")
    print(f"Invertebrate rows: {result.selection_summary['n_selected_invertebrate']}")
    print(f"Selected metadata: {paths['selected_meta']}")
    print(f"Selected IDs: {paths['selected_ids']}")

    if args.train_final_model:
        train_final_models(args, out_dir=out_dir, selected_ids=result.selected_ids)
        print("Training complete")


if __name__ == "__main__":
    main()
