Google colab link: https://colab.research.google.com/drive/1u7pWCd-Jpg1_U6xAdtR4rpJ3HI5Mr8b_#scrollTo=A0kMXXYEAJM9

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
