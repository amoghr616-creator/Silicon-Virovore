import requests
import json

# Read the sequence from your uploaded file (skipping the header line)
sequence = ""
with open("P61567.fasta.txt", "r") as f:
    for line in f:
        if not line.startswith(">"):
            sequence += line.strip()

print(f"Sending sequence to ESMFold (Length: {len(sequence)} AA)...")

# Call the ESMFold API
url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
response = requests.post(url, data=sequence)

if response.status_code == 200:
    with open("herv_k_env.pdb", "w") as out_file:
        out_file.write(response.text)
    print("Success! 3D structure saved as 'herv_k_env.pdb'")
else:
    print(f"Error folding protein: {response.status_code}")