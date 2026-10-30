# P3 Final Defensible Evidence Package

Created: 2026-10-27T17:42:33-04:00

## Main conclusion

The existing evidence supports improved optimization against the internal score.
It does not establish improved real-world or experimentally measured binding.

## Start here

1. `logbook_DATA_RESULTS.md` — curated data/results logbook entry.
2. `logbook_ANALYSIS.md` — corresponding interpretation and limitations.
3. `poster_graphs/` — three suggested core poster figures.
4. `graphs/` — nine curated plots for poster selection or supporting material.
5. `supporting_data/` — source metrics, existing candidate records, and audit tables.
6. `source_reports/` — reports and logs from the completed run and validation scripts.
7. `source_archive/original_full_comparison_audit/` — unmodified copy of the earlier full audit.
8. `unresolved_evidence/` — source conflicts that must not be presented as confirmed results.
9. `evidence_manifest.csv` — checksums and package locations.

## Scope

This was a packaging and figure-generation step. It did not launch the optimizer,
add random seeds, rerun a seed-by-seed comparison, or rerun validation tests.

The aggregate optimizer figure uses the existing inventory. The curated figures
exclude the provisional PEPBI 125-match count and the unreconciled positive
exact-21-aa subset correlation. Original outputs are preserved in the archive.

## Important scientific distinction

Internal score improvement is not equivalent to experimentally stronger binding.
The current retrospective benchmarks do not demonstrate that the advantage
translates into better experimental affinity for the generated candidates.
