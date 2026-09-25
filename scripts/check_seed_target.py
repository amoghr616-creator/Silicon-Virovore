from pathlib import Path

SEED = "MKLAVFALLVFFAGSSDLIRR"
FASTA = Path("data/P61567.fasta.txt")


def read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(
        line.strip()
        for line in lines
        if line.strip() and not line.startswith(">")
    ).upper()


env = read_fasta(FASTA)

print(f"Seed length: {len(SEED)}")
print(f"Env length:  {len(env)}")

exact = env.find(SEED)

if exact >= 0:
    print(
        f"EXACT MATCH: Env positions "
        f"{exact + 1}-{exact + len(SEED)}"
    )
else:
    print("NO EXACT 21-AA MATCH.")

# Report all exact 3+ residue substrings shared by the seed and Env.
best = []

for length in range(min(len(SEED), 8), 2, -1):
    for i in range(len(SEED) - length + 1):
        fragment = SEED[i:i + length]
        j = env.find(fragment)

        if j >= 0:
            best.append(
                (
                    length,
                    fragment,
                    i + 1,
                    j + 1,
                )
            )

    if best:
        break

if best:
    print("\nLongest exact shared fragment(s):")
    for length, fragment, seed_pos, env_pos in best:
        print(
            f"  {fragment}: "
            f"seed {seed_pos}-{seed_pos + length - 1}, "
            f"Env {env_pos}-{env_pos + length - 1}"
        )
else:
    print("No shared fragment of length >= 3.")