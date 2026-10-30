# Logbook Supplement — Structural Test Analysis

This test asks whether adaptive-versus-baseline differences persist under a target-aware geometric evaluation of existing candidate structures.

For `prototype2/structural_validation/haddock_audit_20261024_184403/{ARM}_models`, the adaptive model pass rate was 86.1% and the baseline pass rate was 92.5%; adaptive was lower on this structural criterion.
The adaptive-minus-baseline difference in mean cluster pass rate was -0.064, 95% cluster-bootstrap interval [-0.228, 0.100], exploratory two-sided cluster permutation p=0.5726.
For `prototype2/structural_validation/haddock_audit_20261024_184440/{ARM}_models`, the adaptive model pass rate was 86.1% and the baseline pass rate was 92.5%; adaptive was lower on this structural criterion.
The adaptive-minus-baseline difference in mean cluster pass rate was -0.064, 95% cluster-bootstrap interval [-0.228, 0.100], exploratory two-sided cluster permutation p=0.5726.
For `prototype2/structural_validation/haddock_audit_fixed_20261024_191309/{ARM}_models`, the adaptive model pass rate was 86.1% and the baseline pass rate was 92.5%; adaptive was lower on this structural criterion.
The adaptive-minus-baseline difference in mean cluster pass rate was -0.064, 95% cluster-bootstrap interval [-0.228, 0.100], exploratory two-sided cluster permutation p=0.5726.
For `prototype2/structural_validation/haddock_audit_retry_20261024_185024/{ARM}_models`, the adaptive model pass rate was 86.1% and the baseline pass rate was 92.5%; adaptive was lower on this structural criterion.
The adaptive-minus-baseline difference in mean cluster pass rate was -0.064, 95% cluster-bootstrap interval [-0.228, 0.100], exploratory two-sided cluster permutation p=0.5726.

These results can indicate whether an internal-score advantage is accompanied by a difference in structural plausibility. They do not establish stronger experimental binding. A direct real-world conclusion requires comparable experimental binding-affinity measurements for candidates from both arms. Candidate sequence and run lineage must also be verified before attributing the result to the current optimizer comparison.
