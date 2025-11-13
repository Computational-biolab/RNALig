# RNALig: An AI-driven Scoring Function for RNA–Ligand Binding Affinity

> **Status:** Research code accompanying the JPCB submission: *“RNALig: Learning Structure-Aware Features to Predict RNA–Small Molecule Binding Free Energy”*

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#license)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](#requirements)
[![Colab Predictor](https://img.shields.io/badge/Colab-Binding%20Affinity%20Predictor-black?logo=googlecolab)](COLAB_LINK_PREDICTOR)
[![Colab CleanPDB](https://img.shields.io/badge/Colab-Clean%20PDB-black?logo=googlecolab)](COLAB_LINK_CLEANPDB)

---

## ✨ Overview

**RNALig** predicts RNA–small molecule binding free energy (ΔG, kcal·mol⁻¹) using structure-aware features derived from RNA, ligand, and complex geometries.
It performs **automated feature extraction**, **binding affinity prediction**, and **validation** with interpretable metrics.

**Supported input formats:** **PDB (`.pdb`)** and **mmCIF (`.cif`)**. Multi-model files are auto-split; hetero-ligands are detected automatically.

This repository contains:

* **Feature extraction (Linux CLI)** – batch processing of PDB/mmCIF files with automatic Results folder creation.
* **Binding-affinity predictor (Google Colab)** – ΔG prediction and visualization.
* **Clean-PDB utility (Google Colab)** – structure cleanup and ligand chain handling.

> Every feature is physically validated, unit-checked, and summarized via HTML visualization reports for transparency and reproducibility.

---

## 🧭 Repository Layout

```
RNALig/
├── Features_RNALig_Pro_Full.py     # Unified script for RNA, Ligand, and Complex features (PDB/mmCIF)
├── models/                         # Trained models and feature lists
├── notebooks/                      # Colab notebooks (Predictor, Clean-PDB)
├── data/                           # Example datasets
├── docs/                           # Schematics, flowcharts, and benchmarks
├── CITATION.cff
├── LICENSE
└── README.md
```

---

## 🚀 Quick Start (Feature Extraction)

### Step 1: Create Environment

```bash
conda create -n rnalig python=3.10 -y
conda activate rnalig
conda install -c conda-forge rdkit openbabel mdanalysis biopython freesasa -y
conda install -c bioconda viennarna -y
pip install numpy pandas scipy py3Dmol
```

Alternatively, use the provided environment file:

```bash
conda env create -f environment.yml
conda activate rnalig
```

---

### Step 2: Prepare Input Folder

Place your **PDB** or **mmCIF** files in a directory:

```bash
mkdir input_structures
cp path/to/*.pdb input_structures/
cp path/to/*.cif input_structures/
```

---

### Step 3: Run Feature Extraction

#### 🧩 Batch Mode (PDB + mmCIF)

Extract features for all structures in a directory:

```bash
python Features_RNALig_Pro_Full.py \
  --indir ./input_structures \
  --outcsv final_features.csv \
  --outdir ./viz \
  --viz_rna --viz_ligand
```

#### 🧬 Single Structure Mode

Run on one RNA–ligand complex (PDB or mmCIF):

```bash
# PDB example
python Features_RNALig_Pro_Full.py \
  --pdb ./input_structures/1f27.pdb \
  --outcsv final_features.csv \
  --outdir ./viz \
  --viz_rna --viz_ligand

# mmCIF example
python Features_RNALig_Pro_Full.py \
  --pdb ./input_structures/1f27.cif \
  --outcsv final_features.csv \
  --outdir ./viz \
  --viz_rna --viz_ligand
```

---

### Step 4: Output Folder Structure

After execution, a **viz/** directory is automatically created:

```
viz/
├── final_features.csv              # Combined features for all complexes
├── RNA/                            # RNA-specific HTML reports
│   ├── 1f27_RNA_features.html
│   └── ...
├── Ligand/                         # Ligand-specific HTML reports
│   ├── 1f27_Ligand_features.html
│   └── ...
└── logs/ (optional)                # Log files if errors occur
```

Each **HTML file** visually summarizes feature extraction (pocket detection, SASA, shape metrics, etc.).

---

## 🔧 Advanced Options

Use these flags for control over extraction and visualization:

* `--viz_rna`, `--viz_ligand` – generate per-structure HTML reports.
* `--pocket_cutoff 5.0` – distance (Å) for defining RNA pocket residues.
* `--pocket_sasa 0.05` – pocket SASA threshold.
* `--rna_label_topk 5` – top residue labels shown near the binding site.
* `--min_heavy 4` – minimum heavy atoms required for ligands.
* `--no-require_carbon` – allow ions or metal cofactors as ligands.
* `--keep_ions` – preserve ions in system.
* `--cutoff 5.0` – atomic contact cutoff for RNA–ligand interactions.

Example:

```bash
python Features_RNALig_Pro_Full.py \
  --indir Training \
  --outdir Training/Results \
  --outcsv All_Features.csv \
  --viz_rna --viz_ligand \
  --pocket_cutoff 5.0 --pocket_sasa 0.05 \
  --rna_label_topk 5 --min_heavy 4 \
  --no-require_carbon --keep_ions --cutoff 5.0
```

---

## 🧪 Features Summary

All RNA, Ligand, and Complex features are computed by **`Features_RNALig_Pro_Full.py`**.

### RNA Features

* **Shape metrics:** radius of gyration (Rg), κ² anisotropy, asphericity, acylindricity.
* **Sequence composition:** base counts (A, C, G, U), GC%.
* **Folding stability:** MFE (kcal·mol⁻¹) via ViennaRNA.
* **Solvent exposure:** SASA (Å²) using FreeSASA.
* **Pocket metrics:** residue count, pocket depth, mean and max distances.

### Ligand Features

* **Geometry/shape:** asphericity, eccentricity, inertia eigenvalues.
* **Physicochemical:** LogP, TPSA, HBD/HBA, charge stats, rotatable bonds, MMFF94 energy.
* **Surface:** SASA (polar/nonpolar), vdW volume.

### Complex Features

* **Contacts:** contact count, H-bonds, hydrophobic and vdW interactions.
* **Buried surface area (BSA)** and pocket COM distances.

---

## 🧬 Validation & Reporting

* All outputs are unit-checked and logged in **HTML reports**.
* Each computed feature is summarized for inspection.
* Reviewers can verify per-structure reports in `viz/RNA/` and `viz/Ligand/`.

---

## 🧯 Troubleshooting

* **Missing freesasa:** ensure it’s installed via conda (`conda install -c conda-forge freesasa`).
* **CIF parsing errors:** use valid mmCIF files or re-export via the PDB website.
* **RDKit/OpenBabel conflicts:** avoid installing `rdkit-pypi`; always use conda-forge RDKit.
* **User-site conflicts:** if pip warns `Defaulting to user installation`, set `export PYTHONNOUSERSITE=1`.

---

## 📬 Contact

**NextGen Computational Biology Lab**
*Maintainer:* Priyanka Sharma
*Email:* [your.email@domain](mailto:your.email@domain)
*Website:* <lab or project website>

---

> 🧡 *This version of RNALig uses the unified `Features_RNALig_Pro_Full.py` script with PDB/mmCIF support, visual validation (RNA & Ligand), and streamlined feature generation for publication-ready reproducibility.*
