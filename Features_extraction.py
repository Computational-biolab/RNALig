<<<<<<< HEAD
#Google colab link: https://colab.research.google.com/drive/1u7pWCd-Jpg1_U6xAdtR4rpJ3HI5Mr8b_#scrollTo=A0kMXXYEAJM9

# Clean uninstall
!pip uninstall -y numpy pandas scipy rdkit

# Reinstall compatible versions
!pip install numpy==1.23.5 pandas==1.5.3 scipy==1.9.3 rdkit-pypi
!pip install rdkit-pypi==2022.9.5

import numpy as np
import pandas as pd
import scipy
from rdkit import Chem

!pip install rdkit-pypi
!pip install -q condacolab
import condacolab
condacolab.install()
!mamba install -c rdkit rdkit -y
!pip install MDAnalysis biopython pubchempy scipy py3Dmol rdkit-pypi
!apt-get update
!apt-get install -y vienna-rna
!pip install biopython mdanalysis py3Dmol
!apt-get install -y openbabel

import os
import re
import csv
import subprocess
import tempfile
import argparse
import logging
import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import cdist
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
from Bio.PDB import PDBParser
import MDAnalysis as mda
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import joblib
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB import Select
import tempfile, subprocess
import MDAnalysis as mda
import py3Dmol
from IPython.display import HTML, display
from rdkit.Chem import rdmolfiles

from google.colab import files
import os

uploaded = files.upload()
for pdb_filename in uploaded.keys():
    pdb_id = os.path.splitext(pdb_filename)[0]

## RNA Specific Features
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)

def extract_pdb_id(filename):
    base = os.path.splitext(filename)[0]
    cleaned = re.sub(r'\s*\([^)]*\)$', '', base)
    return cleaned.strip()

def show_rna_structure(pdb_file):
    with open(pdb_file, 'r') as file:
        pdb_data = file.read()
    view = py3Dmol.view(width=500, height=400)
    view.addModel(pdb_data, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    return HTML(view._make_html())

def extract_rna_features(pdb_file, pdb_id, show_structure=False):
    features = {}
    structure = parser.get_structure(pdb_id, pdb_file)
    valid_residues = {"A", "U", "G", "C"}
    rna_chain_id = None
    rna_sequence = ""

    for chain in structure.get_chains():
        for residue in chain:
            if residue.resname.strip() in valid_residues:
                rna_chain_id = chain.id
                break
        if rna_chain_id:
            break

    if rna_chain_id:
        for residue in structure[0][rna_chain_id]:
            if residue.resname.strip() in valid_residues:
                rna_sequence += residue.resname.strip()

        features['nucleotide_composition_A'] = rna_sequence.count('A') / len(rna_sequence)
        features['nucleotide_composition_U'] = rna_sequence.count('U') / len(rna_sequence)
        features['nucleotide_composition_G'] = rna_sequence.count('G') / len(rna_sequence)
        features['nucleotide_composition_C'] = rna_sequence.count('C') / len(rna_sequence)
        features['gc_content'] = (rna_sequence.count('G') + rna_sequence.count('C')) / len(rna_sequence)

        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as tmp_seq_file:
            tmp_seq_file.write(rna_sequence)
            tmp_seq_filename = tmp_seq_file.name

        try:
            result = subprocess.run(['RNAfold', tmp_seq_filename], capture_output=True, text=True)
            mfe_line = result.stdout.split('\n')[1]
            mfe = float(mfe_line.split('(')[-1].strip(')'))
            features['minimum_free_energy'] = mfe
        except Exception as e:
            print("RNAfold error:", e)
            features['minimum_free_energy'] = None
        finally:
            os.remove(tmp_seq_filename)

        try:
            atom_count = len(mda.Universe(pdb_file).select_atoms("nucleic"))
            features['atom_count'] = atom_count
        except Exception as e:
            print("MDAnalysis error:", e)
            features['atom_count'] = 0
    else:
        features.update({
            'nucleotide_composition_A': 0,
            'nucleotide_composition_U': 0,
            'nucleotide_composition_G': 0,
            'nucleotide_composition_C': 0,
            'gc_content': 0,
            'minimum_free_energy': None,
            'atom_count': 0,
        })

    features['pdb_id'] = pdb_id
    features['pdb_file'] = pdb_file
    return features

def analyze_uploaded_pdbs(uploaded):
    tmp_dir = tempfile.mkdtemp()
    pdb_files = []
    pdb_ids = []

    for filename in uploaded.keys():
        filepath = os.path.join(tmp_dir, filename)
        with open(filepath, 'wb') as f:
            f.write(uploaded[filename])
        pdb_files.append(filepath)
        pdb_ids.append(extract_pdb_id(filename))

    all_features = [extract_rna_features(f, pid) for f, pid in zip(pdb_files, pdb_ids)]
    df = pd.DataFrame(all_features)

    styled = df.drop(columns=['pdb_file'], errors='ignore').style.set_table_attributes("style='display:inline'").set_caption("RNA Features Summary")
    display(HTML("<h3>RNA Features for Uploaded PDBs</h3>"))
    display(styled)

    return df

rna_df = analyze_uploaded_pdbs(uploaded)
# Extract PDB ID and path for consistent tracking
pdb_file = list(uploaded.keys())[0]
pdb_path = pdb_file
pdb_id = os.path.splitext(pdb_file)[0]

# Function to detect ligands in a PDB file
import os
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import Descriptors as Desc
from rdkit.Chem import Lipinski, Crippen, rdMolDescriptors
import pubchempy as pcp
from IPython.display import display, HTML
import pandas as pd
from rdkit.Chem.Descriptors import ExactMolWt


parser = PDBParser(QUIET=True)
common_ions = {"NA", "MG", "CL", "K", "ZN", "CA", "MN", "FE", "CU", "CO", "CD"}


def detect_ligands(pdb_file):
    structure = parser.get_structure("Structure", pdb_file)
    ligands = []

    # Loop over all residues in the PDB structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # If the residue is not a standard RNA nucleotide (A, U, G, C) and not an ion
                if (residue.get_resname() not in ['A', 'U', 'G', 'C']) and (residue.get_resname().upper() not in common_ions):
                    ligands.append(residue.get_resname())  # Add the ligand to the list

    return ligands

# Function to extract ligand features (SMILES, molecular weight, hydrogen bond donors/acceptors, etc.)
def extract_ligand_features(ligand_name):
    features = {}

    try:
        # Try to fetch ligand SMILES from PubChem
        compounds = pcp.get_compounds(ligand_name, 'name')
        if compounds:
            smiles = compounds[0].isomeric_smiles
            print(f"Found SMILES for {ligand_name}: {smiles}")

            # Generate RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                features['SMILES'] = smiles
                features['Molecular Weight (g/mol)'] = ExactMolWt(mol)
                features['Hydrogen Bond Donors'] = Lipinski.NumHDonors(mol)
                features['Hydrogen Bond Acceptors'] = Lipinski.NumHAcceptors(mol)
                features['LogP'] = Crippen.MolLogP(mol)
                features['TPSA (Å²)'] = rdMolDescriptors.CalcTPSA(mol)
                features['Atom Count'] = mol.GetNumAtoms()
                features['Chiral Atom Count'] = len([atom for atom in mol.GetAtoms() if atom.HasProp("_ChiralityPossible")])
            else:
                print(f"Failed to generate molecule for SMILES: {smiles}")
        else:
            print(f"No compounds found in PubChem for ligand: {ligand_name}")

    except Exception as e:
        print(f"Error extracting features for ligand {ligand_name}: {e}")

    return features

from IPython.display import display, HTML
import pandas as pd

def analyze_ligands_from_pdb(pdb_path):
    ligands = detect_ligands(pdb_path)
    print(f"🔍 Ligands detected: {ligands}")

    all_features = []
    for ligand in ligands:
        feats = extract_ligand_features(ligand)
        if feats:
            feats['Ligand'] = ligand
            all_features.append(feats)

    if all_features:
        df = pd.DataFrame(all_features)
        styled = df.style.set_table_attributes("style='display:inline'").set_caption("Ligand Features Summary")
        display(HTML("<h3>Ligand Features for Uploaded PDBs</h3>"))
        display(styled)
        return df
    else:
        print("⚠️ No ligand features found.")
        return pd.DataFrame()

# This assumes you uploaded a PDB file and stored it in pdb_filename
analyze_ligands_from_uploaded_pdb(pdb_filename)
ligand_df = analyze_multiple_pdbs_for_ligands(uploaded)
pdb_file = list(uploaded.keys())[0]
pdb_path = pdb_file
pdb_id = os.path.splitext(pdb_file)[0]

# Complex Features
def extract_complex_features(pdb_file, pdb_id, show_structure=True):
    features = {}
    u = mda.Universe(pdb_file)
    rna_atoms = u.select_atoms("resname A C G U")
    drug_atoms = u.select_atoms("not resname A C G U")

    if len(rna_atoms) == 0 or len(drug_atoms) == 0:
        return None

    rna_positive_atoms = rna_atoms.select_atoms("name N*")
    drug_negative_atoms = drug_atoms.select_atoms("name O*")

    electrostatic_distances = cdist(rna_positive_atoms.positions, drug_negative_atoms.positions)
    electrostatic_contacts = electrostatic_distances < 5.0
    features['num_electrostatic_contacts'] = np.sum(electrostatic_contacts)
    features['avg_electrostatic_distance'] = np.mean(electrostatic_distances[electrostatic_contacts]) if np.sum(electrostatic_contacts) > 0 else None

    distances = cdist(rna_atoms.positions, drug_atoms.positions)
    vdw_contacts = distances < 4.0
    features['num_vdw_contacts'] = np.sum(vdw_contacts)
    features['avg_vdw_distance'] = np.mean(distances[vdw_contacts]) if np.sum(vdw_contacts) > 0 else None

    return features

# Show 3D + features
def show_complex_features_and_3d(pdb_file, pdb_id="complex"):
    features = extract_complex_features(pdb_file, pdb_id)
    if features is None:
        print("No RNA or ligand atoms detected in complex.")
        return

    import pandas as pd
    from IPython.display import display
    import py3Dmol
    import MDAnalysis as mda

    df = pd.DataFrame([features])
    display(df)

    with open(pdb_file, 'r') as f:
        pdb_content = f.read()

    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_content, 'pdb')

    # Detect RNA and ligand using MDAnalysis (but do NOT change feature code)
    u = mda.Universe(pdb_file)
    rna_atoms = u.select_atoms("resname A C G U")
    drug_atoms = u.select_atoms("not resname A C G U")

    # Show RNA as cartoon with rainbow color
    view.setStyle({'resn': ['A', 'C', 'G', 'U']}, {'cartoon': {'color': 'spectrum'}})

    # Show ligand (non-RNA) as sticks in red
    ligand_resnames = list(set(atom.resname for atom in drug_atoms))
    for resn in ligand_resnames:
        view.setStyle({'resn': resn}, {'stick': {'colorscheme': 'redCarbon'}})

    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()
    view.show()

def analyze_multiple_pdbs_for_complex_features(uploaded_files, show_structure=True):
    all_complex_features = []

    for pdb_filename in uploaded_files.keys():
        print(f"\n📁 Processing Complex: {pdb_filename}")
        pdb_id = os.path.splitext(pdb_filename)[0]
        feats = extract_complex_features(pdb_filename, pdb_id, show_structure=False)

        if feats:
            feats['PDB_ID'] = pdb_id
            all_complex_features.append(feats)
        else:
            print(f"⚠️ No valid RNA-ligand complex found in {pdb_filename}.")

    if all_complex_features:
        df = pd.DataFrame(all_complex_features)
        from IPython.display import display, HTML
        styled = df.style.set_table_attributes("style='display:inline'").set_caption("Complex Interaction Features")
        display(HTML("<h3>Complex Interaction Features Summary</h3>"))
        display(styled)
        return df
    else:
        print("⚠️ No complex features extracted from any PDBs.")
        return pd.DataFrame()

complex_df = analyze_multiple_pdbs_for_complex_features(uploaded)

for pdb_filename in uploaded.keys():
    show_complex_features_and_3d(pdb_filename)
    complex_features = extract_complex_features(pdb_path, pdb_id)

# Create empty dictionaries
rna_dict = {}
ligand_dict = {}
complex_dict = {}

# Fill them using uploaded files
for pdb_filename in uploaded.keys():
    pdb_id = os.path.splitext(pdb_filename)[0]
    pdb_path = pdb_filename

    print(f"\n📦 Processing {pdb_id}")

    # Extract RNA features
    rna_df = analyze_uploaded_pdbs({pdb_filename: uploaded[pdb_filename]})
    rna_dict[pdb_id] = rna_df

    # Extract ligand features
    ligand_df = analyze_ligands_from_pdb(pdb_path)
    ligand_dict[pdb_id] = ligand_df

    # Extract complex features
    complex_feat = extract_complex_features(pdb_path, pdb_id)
    complex_dict[pdb_id] = complex_feat

from google.colab import files
import pandas as pd
import os

def merge_and_download_all_features(rna_dict, ligand_dict, complex_dict, output_file="merged_features_all_pdbs.csv"):
    all_merged = []

    for pdb_id in rna_dict.keys():
        print(f"\n🔄 Merging features for {pdb_id}")

        rna_df = rna_dict.get(pdb_id, pd.DataFrame())
        ligand_df = ligand_dict.get(pdb_id, pd.DataFrame())
        complex_feat = complex_dict.get(pdb_id, {})

        if rna_df.empty:
            print(f"⚠️ Skipping {pdb_id}: RNA features missing.")
            continue

        # Flatten ligand features
        if not ligand_df.empty:
            ligand_df.index = [f"Ligand_{i+1}" for i in range(len(ligand_df))]
            ligand_df = ligand_df.T
            ligand_df.columns = [f"{col}_{i+1}" for i, col in enumerate(ligand_df.columns)]
            ligand_flat = ligand_df.T.reset_index(drop=True).T
        else:
            ligand_flat = pd.DataFrame()

        # Complex
        complex_df = pd.DataFrame([complex_feat]) if complex_feat else pd.DataFrame()

        # Merge
        combined = pd.concat([
            rna_df.reset_index(drop=True),
            ligand_flat.reset_index(drop=True),
            complex_df.reset_index(drop=True)
        ], axis=1)

        combined.insert(0, "PDB_ID", pdb_id)
        all_merged.append(combined)

    final_df = pd.concat(all_merged, ignore_index=True)
    display(HTML("<h3>✅ Final Merged Features for All PDBs</h3>"))
    display(final_df)

    # Save and download
    final_df.to_csv(output_file, index=False)
    files.download(output_file)

    return final_df

final_df = merge_and_download_all_features(rna_dict, ligand_dict, complex_dict)
=======
# streamlit_app.py
import streamlit as st
import pandas as pd
import numpy as np
import tempfile, os, io, zipfile, traceback, joblib, difflib

# Import your existing modules
from Clean_PDB import clean_pdb
from Features_extraction import (
    extract_pdb_id, extract_rna_features,
    detect_ligands as old_detect_ligands,
    extract_ligand_features, extract_complex_features
)

# Try to import Bio.PDB for robust detection; fall back if missing
try:
    from Bio.PDB import PDBParser
    _HAS_BIOPDB = True
    _parser = PDBParser(QUIET=True)
except Exception:
    _HAS_BIOPDB = False
    _parser = None

st.set_page_config(page_title="RNALig — Option1 fixed", layout="wide")

# ---------------------------
# Load model & pipeline
# ---------------------------
@st.cache_resource
def load_pipeline():
    # choose model file; prefer retrained if present
    model_file = "retrained_rf_model_with_full_data.pkl" if os.path.exists("retrained_rf_model_with_full_data.pkl") else "best_rf_model.pkl"
    pipeline = {
        "model": joblib.load(model_file),
        "imputer": joblib.load("imputer.pkl"),
        "scaler": joblib.load("scaler.pkl"),
        "target_scaler": joblib.load("target_scaler.pkl"),
        # a fallback list — still we will attempt to infer exact expected names from fitted objects
        "features": [
            "nucleotide_composition_A", "nucleotide_composition_U",
            "nucleotide_composition_G", "nucleotide_composition_C",
            "gc_content", "minimum_free_energy", "Molecular Weight (g/mol)",
            "Hydrogen Bond Donors", "Hydrogen Bond Acceptors", "LogP",
            "TPSA (Å²)", "Atom Count", "Chiral Atom Count",
            "num_electrostatic_contacts", "avg_electrostatic_distance",
            "num_vdw_contacts", "avg_vdw_distance"
        ]
    }
    return pipeline

pipeline = load_pipeline()

# ---------------------------
# Robust ligand detection
# ---------------------------
def robust_detect_ligands_local(pdb_file, return_residue_tokens=False):
    """
    Return (ligands_list, residues_info).
    residues_info is a list of dicts with keys: resname, chain, resseq, token
    """
    ligands = []
    residues_info = []
    if _HAS_BIOPDB and _parser is not None:
        try:
            structure = _parser.get_structure("s", pdb_file)
            for model in structure:
                for chain in model:
                    ch_id = chain.id if chain.id is not None else ""
                    for residue in chain:
                        resname = residue.get_resname().strip()
                        # residue.id -> tuple (hetfield, resseq, icode)
                        hetflag = residue.id[0]
                        resseq = residue.id[1]
                        icode = residue.id[2].strip() if len(residue.id) > 2 else ""
                        token = f"{resname}_{ch_id}_{resseq}{('_'+icode) if icode else ''}"
                        residues_info.append({
                            "resname": resname,
                            "chain": ch_id,
                            "resseq": resseq,
                            "token": token
                        })
                        # Skip standard nucleotides
                        if resname in ('A','U','G','C','DA','DT','DG','DC'):
                            continue
                        # skip waters/ions simple check
                        if resname.upper() in {"HOH","WAT","NA","CL","MG","CA","K"}:
                            continue
                        # consider as ligand
                        ligands.append(token if return_residue_tokens else resname)
            # dedupe preserving order
            seen = set()
            uniq = []
            for l in ligands:
                if l not in seen:
                    seen.add(l); uniq.append(l)
            return uniq, residues_info
        except Exception as e:
            # fallback to old detector
            print("robust_detect_ligands_local error:", e)
    # fallback: use older detector (may return simple resnames)
    try:
        res = old_detect_ligands(pdb_file) or []
        residues_info = [{"resname": r, "chain": "", "resseq": "", "token": r} for r in res]
        return res, residues_info
    except Exception:
        return [], []

# ---------------------------
# Feature alignment & prediction helpers
# ---------------------------
def infer_expected_feature_names(pipeline_obj):
    """Try to extract exact feature names used for fit from fitted components."""
    for key in ("imputer", "scaler", "model"):
        comp = pipeline_obj.get(key)
        if comp is None:
            continue
        if hasattr(comp, "feature_names_in_"):
            return list(getattr(comp, "feature_names_in_"))
        try:
            if hasattr(comp, "get_feature_names_out"):
                out = list(comp.get_feature_names_out())
                if out:
                    return out
        except Exception:
            pass
    # fallback
    return list(pipeline_obj.get("features", []))

def fuzzy_map(provided_keys, expected_keys):
    """Map expected -> provided using normalized lowercase & difflib close matches."""
    prov_map = {k.lower().replace(" ", "").replace("(", "").replace(")", ""): k for k in provided_keys}
    mapping = {}
    for exp in expected_keys:
        key_norm = exp.lower().replace(" ", "").replace("(", "").replace(")", "")
        if key_norm in prov_map:
            mapping[exp] = prov_map[key_norm]
            continue
        close = difflib.get_close_matches(key_norm, prov_map.keys(), n=1, cutoff=0.80)
        mapping[exp] = prov_map[close[0]] if close else None
    return mapping

def build_input_df(feature_dict, pipeline_obj):
    expected = infer_expected_feature_names(pipeline_obj)
    if not expected:
        return pd.DataFrame([feature_dict])
    # If already has all expected keys
    if all(k in feature_dict for k in expected):
        return pd.DataFrame([{k: feature_dict.get(k, np.nan) for k in expected}])
    # Try fuzzy mapping
    mapping = fuzzy_map(list(feature_dict.keys()), expected)
    missing = [k for k,v in mapping.items() if v is None]
    if missing:
        # include debug for UI
        debug = {
            "expected_sample": expected[:30],
            "provided_sample": list(feature_dict.keys())[:60],
            "missing_expected_count": len(missing),
            "missing_expected_preview": missing[:20]
        }
        raise ValueError("Feature name mismatch (missing expected features).", debug)
    # build row
    row = {exp: feature_dict.get(mapping[exp], np.nan) for exp in expected}
    return pd.DataFrame([row], columns=expected)

def predict_with_pipeline(feature_dict, pipeline_obj):
    X = build_input_df(feature_dict, pipeline_obj)
    imputer = pipeline_obj.get("imputer")
    scaler = pipeline_obj.get("scaler")
    model = pipeline_obj.get("model")
    target_scaler = pipeline_obj.get("target_scaler", None)

    if imputer is not None:
        X_imp = imputer.transform(X)
    else:
        X_imp = X.values
    if scaler is not None:
        X_scl = scaler.transform(X_imp)
    else:
        X_scl = X_imp
    if model is None:
        raise RuntimeError("No model in pipeline.")
    y_scaled = model.predict(X_scl)
    if target_scaler is not None:
        try:
            y = target_scaler.inverse_transform(y_scaled.reshape(-1,1)).flatten()
        except Exception:
            y = target_scaler.inverse_transform(y_scaled)
    else:
        y = y_scaled
    return float(y[0])

# ---------------------------
# Helpers to process PDB -> features -> predict
# ---------------------------
def process_pdb_file(pdb_path, chosen_ligand=None):
    pdb_id = extract_pdb_id(pdb_path)
    try:
        cleaned = clean_pdb(pdb_path)
    except Exception:
        cleaned = pdb_path

    # determine ligand name to pass to extract_ligand_features (resname portion)
    ligand_name = None
    if chosen_ligand:
        ligand_name = chosen_ligand.split("_")[0] if "_" in chosen_ligand else chosen_ligand

    # extract features
    try:
        rna_feats = extract_rna_features(cleaned, pdb_id) or {}
    except Exception as e:
        rna_feats = {}
        print("RNA features error:", e)
    try:
        lig_feats = extract_ligand_features(ligand_name) if ligand_name else {}
    except Exception as e:
        lig_feats = {}
        print("Ligand features error:", e)
    try:
        complex_feats = extract_complex_features(cleaned) or {}
    except Exception as e:
        complex_feats = {}
        print("Complex features error:", e)

    merged = {}
    merged.update(rna_feats)
    merged.update(lig_feats)
    merged.update(complex_feats)
    merged["PDB_File"] = pdb_id
    merged["Used_Ligand"] = ligand_name

    # Try to predict with alignment
    try:
        pred = predict_with_pipeline(merged, pipeline)
        merged["Predicted_Binding_Affinity_kcal_per_mol"] = pred
        merged["Status"] = "SUCCESS"
    except ValueError as e:
        # feature name mismatch; provide debug info
        merged["Predicted_Binding_Affinity_kcal_per_mol"] = None
        merged["Status"] = f"FAILED: feature name mismatch"
        # attach debug info to return for UI if present
        merged["_debug"] = e.args[1] if len(e.args) > 1 else {}
    except Exception as e:
        merged["Predicted_Binding_Affinity_kcal_per_mol"] = None
        merged["Status"] = f"FAILED: {str(e)}"
        merged["_debug"] = {"exception": str(e)}
    return merged

# ---------------------------
# UI - Option 1 only (focus)
# ---------------------------
st.title("🧬 RNALig — RNA–Ligand Binding Affinity Predictor (Option 1)")

st.header("Option 1 — Upload 1–5 PDB files (quick manually-reviewed run)")
uploaded_multi = st.file_uploader("Upload 1–5 PDB files", type=["pdb"], accept_multiple_files=True, key="opt1_uploader")

# session state for choices
if "opt1_choices" not in st.session_state:
    st.session_state.opt1_choices = {}

if uploaded_multi:
    if len(uploaded_multi) > 5:
        st.error("Please upload at most 5 PDB files for Option 1.")
    else:
        files = []
        for f in uploaded_multi:
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
            tmp.write(f.read()); tmp.flush(); tmp.close()
            files.append((f.name, tmp.name))

        st.markdown("#### Per-file ligand selection (auto-select single-ligand files)")
        for idx, (display_name, local_path) in enumerate(files):
            st.markdown(f"**{idx+1}. {display_name}**")
            ligs, residues_info = robust_detect_ligands_local(local_path, return_residue_tokens=True)
            # show residue summary — now includes chain and resseq
            with st.expander("Show residue summary (debug)"):
                if residues_info:
                    # ensure columns order and presence
                    df_res = pd.DataFrame(residues_info)
                    # ensure chain column exists
                    if 'chain' not in df_res.columns:
                        df_res['chain'] = ""
                    st.dataframe(df_res[['resname','chain','resseq','token']].head(200))
                else:
                    st.info("No residue info available.")

            # Behavior:
            # - exactly 1 ligand -> auto-select
            # - >1 ligand -> dropdown and required selection
            # - 0 ligands -> manual input allowed
            if not ligs:
                st.warning("No ligands auto-detected.")
                manual_val = st.text_input(f"Manual ligand for {display_name} (resname or token)", key=f"opt1_manual_{idx}")
                use_manual = st.checkbox(f"Use manual ligand for {display_name}", key=f"opt1_use_manual_{idx}")
                if use_manual and manual_val.strip():
                    st.session_state.opt1_choices[display_name] = manual_val.strip()
                else:
                    st.session_state.opt1_choices.setdefault(display_name, None)
            elif len(ligs) == 1:
                st.success(f"Single ligand detected: {ligs[0]} (auto-selected)")
                st.session_state.opt1_choices[display_name] = ligs[0]
            else:
                st.info(f"Multiple ligands detected ({len(ligs)}). Please choose one.")
                choice = st.selectbox(f"Select ligand for {display_name}", options=ligs, key=f"opt1_select_{idx}")
                st.session_state.opt1_choices[display_name] = choice

        st.markdown("---")
        # show quick-predict when single file only and ligand chosen
        if len(files) == 1 and st.session_state.opt1_choices.get(files[0][0]) is not None:
            if st.button("Quick Predict this single PDB"):
                chosen = st.session_state.opt1_choices.get(files[0][0])
                res = process_pdb_file(files[0][1], chosen_ligand=chosen)
                if res.get("_debug"):
                    st.error("Prediction failed due to feature-name mismatch. See debug below.")
                    st.json(res.get("_debug"))
                elif res.get("Status", "").startswith("FAILED"):
                    st.error(res.get("Status"))
                else:
                    st.metric("Predicted binding affinity (kcal/mol)", f"{res.get('Predicted_Binding_Affinity_kcal_per_mol'):.4f}")
                    st.json({k: res.get(k) for k in ("PDB_File","Used_Ligand","Predicted_Binding_Affinity_kcal_per_mol","Status") if k in res})

        # Run Predictions for multiple files
        if st.button("Run Predictions (Option 1)"):
            missing = [n for n,ch in st.session_state.opt1_choices.items() if ch is None]
            if missing:
                st.error("Cannot run: some files are missing ligand selections. Fill them first.")
            else:
                results = []
                prog = st.progress(0)
                total = len(files)
                for i,(display_name, local_path) in enumerate(files, start=1):
                    chosen = st.session_state.opt1_choices.get(display_name)
                    row = process_pdb_file(local_path, chosen_ligand=chosen)
                    results.append(row)
                    prog.progress(int(i/total*100))
                df = pd.DataFrame(results)
                st.success("Option 1 predictions complete")
                # show debug column if present
                if any("_debug" in r for r in results):
                    st.warning("Some rows have debug info (feature name mismatch). See _debug column.")
                display_cols = [c for c in df.columns if c != "_debug"]
                st.dataframe(df[display_cols].fillna("").head(200))
                st.download_button("Download CSV (Option 1)", df.to_csv(index=False).encode(), file_name="rnalig_opt1_results.csv", mime="text/csv")
else:
    st.info("Upload 1–5 PDB files to run Option 1.")

# End of file — Option 2 omitted here for brevity (you can keep your existing Option2 code)
>>>>>>> 1e3dc63 (Initial commit of RNALig GUI)
