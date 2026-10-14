# Prototype 2 pilot runbook

This is a small, isolated preflight-gated pilot design. It has not been
executed. The machine-readable specification is
`prototype2/p2_pilot_plan.json`.

## Gate before running

Run the following from the repository root:

```bash
python prototype2/run_p2_rapidock.py \
  --preflight-only \
  --output-root prototype2/pilot_runs \
  --box 52 \
  --n-samples 20 \
  --seed 616
```

The command must exit `0`. It must report the intended `hybridock-pep`
executable, supported dock arguments, a compatible `longer_local.pt`, valid
ARG203/ARG205 mapping, and a comparable receptor coordinate frame. Exit `2`
means the pilot is blocked; do not replace the checkpoint or use the fallback.

The preflight report is written beside the requested output directory, so it
does not contaminate the directory that must remain empty for a real run.

## Fixed pilot

The six fixed 21-residue sequences and seeds are in `p2_pilot_plan.json`.
Use one subprocess at a time, 20 Stage-1 samples per sequence, and a 30-minute
per-candidate wall-time limit. The wrapper assigns seed 616 plus the candidate
index when launched with `--seed 616`.

This is a computational candidate/control comparison, not a claim that the
control arm is biologically negative. The existing sequence-only baseline
must remain a separate analysis arm.

## Evidence gates

After a successful run, independently audit every selected pose. A pose is
eligible for target-aware evidence only if its sequence and coordinates agree,
its peptide chain is unambiguous, the target mapping and coordinate frame are
valid, it has no severe clash, and it has an independently measured contact to
the intended ARG203/ARG205 site. Malformed, mismatched, ambiguous, off-pocket,
or contact-free poses remain descriptive-only. A favorable docking score is not
binding proof.
