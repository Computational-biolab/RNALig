# RNALig: AI-Driven RNA–Ligand Binding Affinity Predictor

RNALig is a machine learning-based pipeline designed to **predict the binding affinity** of RNA–ligand complexes. It supports structure cleaning, feature extraction, and binding affinity prediction — all through easy-to-use Google Colab notebooks!

1. Clean PDBs: [PDB Cleaner](https://colab.research.google.com/drive/1LSxz-l2kczM9fi3W_mor72IlOP2vFq7R) | Remove unwanted molecules and prep PDBs for feature extraction having receptor and ligands.
2. Extract Features: [Feature Extractor](https://colab.research.google.com/drive/1u7pWCd-Jpg1_U6xAdtR4rpJ3HI5Mr8b_) | Extract RNA, ligand, and RNA–ligand complex features |
3. Predict Binding Affinity: [Binding Affinity Predictor](https://colab.research.google.com/drive/1ZFgGIhVFuunZtIllH1kUleCFePZYWsD8) | Predict binding affinity using a trained ML model |

> Powered by **Random Forest Regressor**, trained on curated RNA–ligand data  
> Evaluated with **R², RMSE, MAE, and PCC**

## Folder Structure

RNALig/
├── notebooks/ # All Colab notebooks
├── Features.py # Feature extraction script (modular)
├── Binding Affinity Predictor.py # Prediction script using Random Forest
├── models/ # Trained ML model (joblib .pkl)
├── sample_input/ # Example PDBs and input CSVs
├── requirements.txt # Python dependencies

## Installation (Optional - for local use)

```bash
git clone https://github.com/yourusername/RNALig.git
cd RNALig
pip install -r requirements.txt
sudo apt-get install -y vienna-rna openbabel

Features Extracted
RNA Descriptors
Ligand Properties
RNA–Ligand Interaction Metrics 
Predictive Output: Binding Affinity in kcal/mol
