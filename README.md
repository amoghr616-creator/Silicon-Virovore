# Silicon Virovore

> A computational peptide engineering platform for the rational design and in silico evaluation of peptide inhibitors targeting the Human Endogenous Retrovirus K (HERV-K) envelope protein.

---

## Authors

**Amogh Ramesh**  
Lead Software Engineer • Computational Biology • Pipeline Architecture

**Paul Vu**  
Research Collaborator • Computational Biology • Scientific Development

--

## Overview

Silicon Virovore is a computational drug-discovery platform that integrates peptide generation, molecular docking, molecular dynamics, structural analysis, and hardware telemetry into a unified, reproducible workflow.

The project investigates whether computationally designed peptide inhibitors can bind conserved regions of the HERV-K envelope glycoprotein while demonstrating an end-to-end engineering pipeline suitable for rapid therapeutic candidate evaluation.

Although developed as an independent high school research project, the software is designed using modular software engineering principles commonly found in academic computational biology laboratories.

---

## Motivation

Human Endogenous Retrovirus K (HERV-K) is the most recently active endogenous retrovirus within the human genome and has been implicated in several pathological conditions, including certain cancers and neurodegenerative disorders.

Traditional peptide discovery is experimentally intensive and time consuming.

Silicon Virovore explores whether modern computational biology methods can accelerate early-stage peptide discovery through automated in silico screening before experimental validation.

---

# Features

- Modular computational biology pipeline
- Evolutionary peptide optimization
- High-performance C acceleration backend
- Automated molecular docking workflow
- Molecular dynamics simulation support
- Structural visualization
- Hardware telemetry integration
- Interactive analysis dashboard
- Reproducible pipeline architecture

---

# Computational Workflow

Protein Target
      │
      ▼
Sequence Processing
      │
      ▼
Peptide Generation
      │
      ▼
Candidate Optimization
      │
      ▼
Structure Prediction
      │
      ▼
Receptor Preparation
      │
      ▼
Molecular Docking
      │
      ▼
Pose Analysis
      │
      ▼
Molecular Dynamics
      │
      ▼
Scoring & Ranking
      │
      ▼
Visualization & Dashboard
```

---

# Repository Structure

```
Silicon-Virovore/

├── src/                # Core computational pipeline
├── c/                  # High-performance C backend
├── configs/            # Configuration files
├── dashboard/          # Streamlit visualization
├── data/               # Input datasets
├── results/            # Generated outputs
├── figures/            # Images for documentation
├── docs/               # Technical documentation
├── firmware/           # Arduino telemetry system
├── tests/              # Unit and integration tests
└── README.md
```

---

# Software Architecture

Silicon Virovore is organized as a modular workflow in which each stage performs a single computational task.

Each module produces standardized outputs that serve as inputs to the next stage, allowing the pipeline to remain reproducible, maintainable, and extensible.

The computational architecture separates

- biological data processing
- optimization algorithms
- structural prediction
- docking
- molecular dynamics
- visualization
- hardware monitoring

into independent components.

---

# Technologies

| Category | Software |
|----------|----------|
| Language | Python |
| Systems Programming | C |
| Structural Biology | PyMOL |
| Docking | AutoDock Vina |
| Molecular Dynamics | GROMACS |
| Bioinformatics | Biopython |
| Data Analysis | NumPy, Pandas |
| Visualization | Matplotlib |
| Dashboard | Streamlit |
| Embedded Systems | Arduino |

---

# Installation

Clone the repository

```bash
git clone https://github.com/yourusername/Silicon-Virovore.git
cd Silicon-Virovore
```

Install dependencies

```bash
pip install -r requirements.txt
```

Compile the C backend

```bash
make
```

Run the pipeline

```bash
python run_pipeline.py
```

Launch the dashboard

```bash
streamlit run dashboard/app.py
```

---

# Example Outputs

The pipeline produces

- ranked peptide candidates
- docking affinity tables
- docking poses
- RMSD trajectories
- RMSF analyses
- contact maps
- structural visualizations
- interactive dashboards

Example output files are located in

```
results/
```

---

# Engineering Highlights

- Modular software architecture
- Native C acceleration
- Automated computational workflow
- Reproducible analysis pipeline
- Cross-platform design
- Integrated visualization tools
- Hardware/software co-design

---

# Scientific Limitations

Silicon Virovore performs computational hypothesis generation only.

Docking scores and molecular dynamics simulations should not be interpreted as experimental evidence of therapeutic efficacy.

Experimental validation, including biochemical binding assays and cell-based studies, is required to evaluate biological activity.

---

# Future Work

Planned improvements include

- machine-learning-assisted peptide generation
- multi-objective evolutionary optimization
- GPU acceleration
- expanded structural databases
- AlphaFold integration improvements
- experimental validation
- automated statistical benchmarking

---

# Citation

If you use Silicon Virovore in academic work, please cite the repository or contact the author.

---

# License

MIT License

---

# Acknowledgments

This project was developed as an independent computational biology and software engineering research project.

The author acknowledges the developers of open-source scientific software including Biopython, AutoDock Vina, GROMACS, PyMOL, Streamlit, and the broader computational biology community.