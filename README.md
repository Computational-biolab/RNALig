# RNALig: AI-Driven RNA–Ligand Binding Affinity Predictor

RNALig is a structure-based machine learning model designed to predict the binding affinity between RNA molecules and small ligands using their 3D complex structures (PDB format). It enables non-computational users, especially medicinal chemists, to:

-Clean PDB files
-Extract meaningful interaction features
-Predict binding affinity (ΔG in kcal/mol)

*Powered by Random Forest Regressor, RNALig was trained on curated experimental RNA–ligand complexes and validated against real binding data and evaluated with R², RMSE, MSE, MAE, and PCC.

*Try it instantly on Google Colab – no installation required!
1. Clean PDBs: [PDB Cleaner](https://colab.research.google.com/drive/1LSxz-l2kczM9fi3W_mor72IlOP2vFq7R) | Remove unwanted molecules and prep PDBs for feature extraction having receptor and ligands.
2. Extract Features: [Feature Extractor](https://colab.research.google.com/drive/1u7pWCd-Jpg1_U6xAdtR4rpJ3HI5Mr8b_) | Extract RNA, ligand, and RNA–ligand complex features |
3. Predict Binding Affinity: [Binding Affinity Predictor](https://colab.research.google.com/drive/1ZFgGIhVFuunZtIllH1kUleCFePZYWsD8) | Predict binding affinity using a trained ML model |

## Folder Structure

RNALig/
├── notebooks/ # All Colab notebooks
├── Features.py # Feature extraction script (modular)
├── Binding Affinity Predictor.py # Prediction script using Random Forest
├── models/ # Trained ML model (joblib .pkl)
├── sample_input/ # Example PDBs and input CSVs
├── requirements.txt # Python dependencies


*Add One-Line Command for Prediction
#python Binding_Affinity_Predictor.py --input sample_input/features.csv --output prediction.csv

*Compatibility 
Tested on:
- Python 3.9+
- Biopython, RDKit, ViennaRNA, OpenBabel
- Google Colab (100% compatible)

## Installation (Optional - for local use)

```bash
git clone https://github.com/yourusername/RNALig.git
cd RNALig
pip install -r requirements.txt
sudo apt-get install -y vienna-rna openbabel
Load Model dependencies
Features Extracted
RNA-specific features
Ligand-specific features
Comples-specific features
RNA–Ligand Interaction Metrics 
Predictive Output: Binding Affinity (kcal/mol)
