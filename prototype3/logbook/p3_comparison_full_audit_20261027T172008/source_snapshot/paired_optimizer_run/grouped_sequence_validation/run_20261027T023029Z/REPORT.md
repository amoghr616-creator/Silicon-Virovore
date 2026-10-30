# Leakage-resistant sequence-only mutation-pair validation

- Input file: `/Users/centurion616/Desktop/Silicon-Virovore/prototype3/results/ppikb_kd_benchmark/filtered_kd_measurements.csv`
- Collapsed target/assay/peptide records: **2496**
- Single-substitution pairs: **693**
- Unique target sequences: **123**
- Leakage components: **111**; folds: **5**
- Split rule: connected components link all rows sharing a target sequence, an individual peptide sequence, or an exact unordered peptide pair.
- Each fold was checked for zero shared target sequences, individual peptide sequences, and exact pair signatures.
- Ridge alpha: **10.0**, fixed without tuning on held-out folds.

## Out-of-fold metrics

| Model | n | MAE (pKd units) | RMSE | Spearman rho | Informative direction coverage | Direction accuracy among covered | Direction score (ties half) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Zero-change baseline | 693 | 0.8254 | 1.3840 | NA | 0.0000 | NA | 0.5000 |
| Training-fold mean | 693 | 0.8182 | 1.3656 | 0.0459 | 1.0000 | 0.6266 | 0.6266 |
| Sequence-only Ridge | 693 | 0.8643 | 1.3957 | 0.1087 | 0.9099 | 0.5896 | 0.5815 |

## Interpretation guardrails

This is a separate sequence-only research comparator. It does not validate the native `c_score`, prove the adaptive optimizer improves, or validate the target-aware structural scorer. A positive result supports only the reported held-out statistical relationship under this split and feature definition.

Do not report target-binding, biological efficacy, or HERV-K-specific accuracy from this benchmark. A null/weak result is a valid scientific outcome.
