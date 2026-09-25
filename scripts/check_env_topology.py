from pathlib import Path

FASTA = Path("data/P61567.fasta.txt")

SEQUENCE = "".join(
    line.strip()
    for line in FASTA.read_text().splitlines()
    if line.strip() and not line.startswith(">")
).upper()

# UniProt-derived annotation for P61567.
TM_START = 522
TM_END = 542

print(f"P61567 length: {len(SEQUENCE)}")
print(f"Annotated TM: {TM_START}-{TM_END}")
print(f"TM sequence:  {SEQUENCE[TM_START - 1:TM_END]}")

assert len(SEQUENCE) == 588
assert SEQUENCE[TM_START - 1:TM_END] == "TIINLILILVCLFCLLLVCRC"

print("Topology sanity check: PASS")