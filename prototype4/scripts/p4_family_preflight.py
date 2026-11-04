from pathlib import Path
import sqlite3
import pandas as pd
from collections import Counter

ROOT = Path.cwd()
SOURCE = ROOT / "prototype4/data/source/skempi_v2.csv"
AUDIT = ROOT / "prototype4/results/skempi_audit_20261029"
DB = AUDIT / "p4_evidence.sqlite"

OUT = AUDIT / "p4_family_preflight.txt"
lines = []

def emit(value=""):
    lines.append(str(value))

def show_df(name, frame, n=8):
    emit(f"\n--- {name} ---")
    emit(f"Rows: {len(frame)}")
    emit(f"Columns: {list(frame.columns)}")
    if len(frame):
        emit(frame.head(n).to_string(index=False, max_colwidth=100))

if not SOURCE.is_file():
    raise FileNotFoundError(SOURCE)
if not DB.is_file():
    raise FileNotFoundError(DB)

df = pd.read_csv(SOURCE, sep=";", dtype=str).fillna("")

emit("P4 FAMILY-MAPPING PREFLIGHT")
emit("=" * 72)
emit(f"Source: {SOURCE}")
emit(f"Source rows: {len(df)}")
emit(f"Unique PDB-chain entries: {df['#Pdb'].nunique()}")

emit("\n--- PDB CHAIN IDENTIFIERS ---")
pdb_values = df["#Pdb"].drop_duplicates().head(30)

for value in pdb_values:
    parts = value.split("_")
    emit(f"{value!r} -> {parts}")

arity = Counter(len(str(x).split("_")) for x in df["#Pdb"])
emit(f"\nNumber of underscore-separated fields: {dict(arity)}")

emit("\n--- PROTEIN PARTNER NAMES ---")
pairs = (
    df[["Protein 1", "Protein 2"]]
    .drop_duplicates()
    .head(30)
)

for row in pairs.itertuples(index=False):
    emit(f"Protein 1: {row[0]} | Protein 2: {row[1]}")

emit("\n--- MUTATION CHAIN CHECK ---")
chain_groups = []

for value in df["#Pdb"].drop_duplicates():
    parts = value.split("_")
    if len(parts) == 3:
        chain_groups.append((parts[0], parts[1], parts[2]))

emit(f"Entries with exactly three fields: {len(chain_groups)}")
emit("Examples:")
for item in chain_groups[:20]:
    emit(str(item))

emit("\n--- AUDIT OUTPUT FILES ---")
for path in sorted(AUDIT.iterdir()):
    if path.is_file():
        emit(f"{path.name}: {path.stat().st_size:,} bytes")

for filename in ["eligible_cycles.csv", "family_map_template.csv"]:
    path = AUDIT / filename
    if path.exists():
        try:
            table = pd.read_csv(path)
            show_df(filename, table)
        except Exception as exc:
            emit(f"{filename}: could not inspect: {exc}")

emit("\n--- SQLITE SCHEMA ---")

with sqlite3.connect(DB) as con:
    tables = con.execute("""
        SELECT name FROM sqlite_master
        WHERE type='table'
        ORDER BY name
    """).fetchall()

    for (table,) in tables:
        emit(f"\nTABLE: {table}")
        schema = con.execute(
            f'PRAGMA table_info("{table}")'
        ).fetchall()

        for column in schema:
            emit(f"  {column[1]} ({column[2]})")

        if table in {"skempi_raw", "cycle_nodes", "cycle_features"}:
            rows = con.execute(
                f'SELECT * FROM "{table}" LIMIT 2'
            ).fetchall()

            names = [x[1] for x in schema]

            for row in rows:
                emit("  SAMPLE:")
                for name, value in zip(names, row):
                    value = str(value)
                    emit(f"    {name} = {value[:140]}")

text = "\n".join(lines)
OUT.write_text(text)

print(text)
print(f"\nSaved preflight report: {OUT}")
