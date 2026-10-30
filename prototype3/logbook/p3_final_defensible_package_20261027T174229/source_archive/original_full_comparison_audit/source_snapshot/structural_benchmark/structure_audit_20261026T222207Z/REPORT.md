# Experimental mutation-pair structure availability audit

- Input pairs: 693
- PDB IDs requested: 139
- PDB download failures: 0
- Target/peptide/PDB mapping attempts: 265
- Mapping status counts: `{'peptide_chain_not_exact': 202, 'mapped_candidate': 51, 'target_receptor_chain_not_identified': 12}`
- Pairs with sequence-mapped structures for both variants: 18
- Informative pairs (|ΔpKd| ≥ 0.3) with both variants mapped: 14

These are sequence-level eligibility checks only. Each mapped complex must be inspected for chain identity, residues, coordinates, clashes, and matching experimental assay before its structural score can be used in a mutation benchmark.
