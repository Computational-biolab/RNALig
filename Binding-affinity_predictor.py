##https://colab.research.google.com/drive/1t__onHfzPnMIpNyn2w2Vy6e1mZMa_5bd#scrollTo=ZWrqENxFxTZl

# predict_from_pipeline.py 

import pandas as pd
import numpy as np
import joblib
import sys
from pathlib import Path

# ========= CONFIG (edit these filenames if needed) =========
PIPELINE_FILE = "binding_affinity_pipeline.pkl"   # produced by training
FEATURES_FILE = "Features.csv"                    # user-provided features with same columns as training
OUTPUT_FILE   = "RNALig_Predictions.csv"
QUIET = True  # set False to see logs

def _log(msg): 
    if not QUIET: 
        print(msg)

def main():
    # --- sanity on files ---
    for f in [PIPELINE_FILE, FEATURES_FILE]:
        if not Path(f).exists():
            raise FileNotFoundError(f"❌ File not found: {f}")

    # --- load pipeline ---
    pipe = joblib.load(PIPELINE_FILE)
    features = pipe["features"]
    calib = pipe.get("calibration", None)  # {"a":..., "b":...} if present

    # --- load features csv & normalize headers like training ---
    df = pd.read_csv(FEATURES_FILE)
    df.columns = (df.columns.str.strip()
                  .str.replace("Å", "A", regex=False)
                  .str.replace("²", "2", regex=False))

    # Keep PDB_ID if present (free-form)
    pdb_ids = df["PDB_ID"] if "PDB_ID" in df.columns else pd.Series([f"PDB_{i}" for i in range(len(df))])

    # --- verify required feature columns exist (case-sensitive, after normalization) ---
    missing = [f for f in features if f not in df.columns]
    if missing:
        # helpful debug info before failing
        print("❌ Missing required features in your CSV.\n"
              "Here are some tips:\n"
              " • Make sure your header names match exactly what the model was trained on.\n"
              " • The script normalizes 'Å'→'A' and '²'→'2'.\n"
              " • Check for trailing spaces or typos.\n")
        print("Required features (from pipeline):")
        for f in features: print("  -", f)
        print("\nYour CSV columns:")
        for c in df.columns: print("  -", c)
        raise ValueError(f"\nMissing columns: {missing}")

    # --- build feature matrix exactly in trained order & coerce to numeric ---
    X_raw = df[features].copy()
    # Convert any non-numeric entries to NaN; imputer will fill
    for col in X_raw.columns:
        X_raw[col] = pd.to_numeric(X_raw[col], errors="coerce")

    # --- transform & predict ---
    X_imp = pipe["imputer"].transform(X_raw)
    X_scl = pipe["scaler"].transform(X_imp)
    y_scaled = pipe["model"].predict(X_scl)
    y_raw = pipe["target_scaler"].inverse_transform(y_scaled.reshape(-1, 1)).ravel()

    # --- apply calibration silently if present ---
    if calib and all(k in calib for k in ("a", "b")):
        y_final = calib["a"] + calib["b"] * y_raw
        _log("(calibration applied internally)")
    else:
        y_final = y_raw
        _log("(no calibration info found; using raw predictions)")

    # --- output: ONLY PDB_ID + calibrated prediction ---
    out = pd.DataFrame({
        "PDB_ID": pdb_ids,
        "Predicted_Binding_Affinity (kcal/mol)": y_final
    })
    out.to_csv(OUTPUT_FILE, index=False)
    if not QUIET:
        print(f"✅ Saved → {OUTPUT_FILE}")
        print(out.head())

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        # make sure you see useful info even if QUIET=True
        print("\n🚨 Prediction failed.")
        print(e)
        sys.exit(1)
