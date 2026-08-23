#!/bin/bash
gmx pdb2gmx -f lead_1_pipeline_clean.pdb -o lead_1_processed.gro -water tip3p -ff charmm36-feb2026 -ignh -ter
