# Expanded HADDOCK Reliability Audit
Date: 2026-11-02

## Scope

Two expanded HADDOCK3 workflows were assessed:
- Chain A receptor preparation
- Chain B receptor preparation

These are computational structural-model assessments, not experimental
binding measurements.

## Sampling and selection

| Metric | Chain A | Chain B |
|---|---:|---:|
| Rigid-body output files | 120 | 120 |
| Unique rigid-body coordinate sets | 15 | 15 |
| Selected input files | 10 | 10 |
| Unique selected input coordinate sets | 1 | 1 |
| Unique FlexRefine coordinate sets | 1 | 1 |
| Unique EMRefine coordinate sets | 1 | 1 |

Chain A selected rigidbody_44–53, ten copies of one coordinate set.
Chain B selected rigidbody_64–73, ten copies of one coordinate set.

## Scoring audit

Both expanded FlexRefine configurations used w_bsa = -0.01.
Both reported BSA = -999999.

The resulting BSA contribution is:
(-0.01) * (-999999) = +9999.99

Therefore, the expanded FlexRefine scores are invalid for pose ranking.

EMRefine used w_bsa = 0.0.

| Metric | Chain A | Chain B |
|---|---:|---:|
| HADDOCK EMRefine BSA (Å²) | 11.3948 | 9.88994 |
| HADDOCK EMRefine score | 8.95942 | 5.61418 |
| Independent FreeSASA total BSA (Å²) | approximately 0 | approximately 0 |

Independent BSA was calculated as:
SASA(receptor A) + SASA(peptide B) - SASA(complex A+B).

## Interchain geometry

| Metric | Chain A | Chain B |
|---|---:|---:|
| Minimum interchain heavy-atom distance (Å) | 22.744 | 7.212 |
| Atom pairs below 4 Å | 0 | 0 |
| Atom pairs below 5 Å | 0 | 0 |
| Atom pairs below 8 Å | 0 | 7 |

The Chain A closest pair was receptor Gln280 CB to peptide Val5 CG2.
The Chain B closest pair was receptor Val395 CG2 to peptide Arg21 NH2.

## Target-anchor geometry

Distances below are minimum heavy-atom distances from each receptor
Arg side chain to the peptide.

| Anchor | Chain A (Å) | Chain B (Å) |
|---|---:|---:|
| Arg203 | 37.217 | 18.963 |
| Arg205 | 38.848 | 17.175 |

Neither expanded structure had Arg203/Arg205 atom pairs within 8 Å.

## Conclusions

1. Expanded rigid-body output files contained substantial coordinate
   redundancy.
2. Selection retained ten copies of one coordinate set per workflow.
3. FlexRefine BSA sentinel values dominated the reported scores.
4. Independent SASA indicated approximately zero buried surface area.
5. Neither final expanded structure approached the intended Arg203/Arg205
   anchor region.

These runs do not provide positive structural evidence of HERV-K Env
binding or target-site engagement.

The findings concern these specific computational workflows and do not
establish that the peptide cannot bind experimentally or that other
docking protocols would necessarily fail.

## Evidence limitations

- Pilot pose diversity does not establish pose accuracy.
- HADDOCK scores are not binding free energies or dissociation constants.
- Independent SASA is a separate calculation and may differ from
  HADDOCK's internal implementation.
- The expanded workflows are not independent biological replicates.
- No experimental binding or inhibition validation is established.
