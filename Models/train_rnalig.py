# -*- coding: utf-8 -*-
"""RFR_Prediction.ipynb

Original file is located at
    https://colab.research.google.com/drive/1-t8j_M11koPIX4DtugbNiQiBvQqvfI7-
"""

# %%capture
!pip -q install pandas numpy scikit-learn joblib

import io, json, joblib
import numpy as np, pandas as pd
from pathlib import Path
from google.colab import files

from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.feature_selection import VarianceThreshold

SEED = 42
np.random.seed(SEED)

print("Upload Training_data.csv (must include 'Experimental_binding_affinity')...")
uploaded = files.upload()
train_fname = list(uploaded.keys())[0]

df = pd.read_csv(train_fname)
df.columns = [c.strip() for c in df.columns]
print("Loaded:", train_fname, "| shape:", df.shape)

# detect target column (case-insensitive, flexible)
target_col = next((c for c in df.columns if "experimental_binding_affinity" in c.lower()), None)
if target_col is None:
    raise ValueError("Could not find target column (needs 'Experimental_binding_affinity').")

# optional PDB column (not required for training)
pdb_col = next((c for c in df.columns if "pdb" in c.lower()), None)

# numeric features only, drop target
X_all = (
    df.drop(columns=[target_col])
      .select_dtypes(include=[np.number])
      .replace([np.inf, -np.inf], np.nan)
)
y = df[target_col].astype(float).values

# remove near-constant features
vt = VarianceThreshold(threshold=1e-4)
X_var = pd.DataFrame(vt.fit_transform(X_all), columns=X_all.columns[vt.get_support()])

# remove exact-duplicate columns (by rounded vector)
dup_cols, seen = [], {}
for c in X_var.columns:
    key = tuple(np.round(X_var[c].values, 6))
    if key in seen: dup_cols.append(c)
    else: seen[key] = c

X = X_var.drop(columns=dup_cols)
feat_names = X.columns.tolist()

print(f"Features after cleaning: {X.shape[1]} (removed {len(dup_cols)} duplicates)")

pipe = Pipeline([
    ("imp", SimpleImputer(strategy="median")),
    ("rf", RandomForestRegressor(
        n_estimators=1000,
        max_depth=30,
        max_features="sqrt",
        min_samples_leaf=1,
        random_state=SEED,
        n_jobs=-1
    ))
])

pipe.fit(X, y)
print("Trained RF pipeline on full Training_data.csv")

bundle = {
    "model": pipe,                 
    "features": feat_names,       
    "target": target_col,          
    "meta": {
        "n_train": int(len(y)),
        "random_state": SEED,
        "preprocessing": "median imputation + variance/duplicate filter",
    }
}

outpkl = "RNALig_training_model.pkl"
joblib.dump(bundle, outpkl)
print(f"Saved → {outpkl}")
files.download(outpkl)
