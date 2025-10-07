#Google colab link: https://colab.research.google.com/drive/1erVuOR1yuk3TTIljFkcDqIPbRORpQ_8E#scrollTo=fnzg9udlQ2gJ

import re
import numpy as np
import pandas as pd
from google.colab import files

# ---- Constants ----
R = 0.0019872041    # kcal·mol^-1·K^-1
T_DEFAULT = 298.15  # K
KJ_TO_KCAL = 0.239005736
DEBUG = False

# ---- Helpers ----
def _first_number(s: str):
    m = re.search(r"[-+]?\d*\.?\d+(?:e[-+]?\d+)?", str(s))
    return float(m.group()) if m else None

def _normalize(s: str):
    return str(s).strip().lower().replace("μ","u").replace("µ","u").replace(" ","")

def _to_molar(val_str: str):
    """Parse '3.5 nM', '2.6uM', '1e-6 M' → molar (float)."""
    if val_str is None or str(val_str).strip() == "":
        return None
    s = _normalize(val_str)
    num = _first_number(s)
    if num is None:
        return None
    if "pm" in s:  return num * 1e-12
    if "nm" in s:  return num * 1e-9
    if "um" in s:  return num * 1e-6
    if "mm" in s:  return num * 1e-3
    # bare 'm' at end → molar
    if re.search(r"(?<![pnμu])m(?![a-z])", s): return num
    # no explicit unit → assume molar
    return num

def pretty_conc(M):
    """Format molar as a human-friendly string."""
    if M is None: return "NA"
    if M >= 1e-3:   return f"{M*1e3:.3g} mM"
    if M >= 1e-6:   return f"{M*1e6:.3g} µM"
    if M >= 1e-9:   return f"{M*1e9:.3g} nM"
    if M >= 1e-12:  return f"{M*1e12:.3g} pM"
    return f"{M:.3e} M"

# ---- IC50 → Ki (Cheng–Prusoff, competitive) ----
def convert_ic50_to_ki(ic50_str, substrate_str=None, km_str=None, mechanism="competitive"):
    """
    Cheng–Prusoff (competitive): Ki = IC50 / (1 + [S]/Km)
    If [S] or Km missing, falls back to Ki≈IC50 (factor=1) with a clear warning.
    Returns (Ki_M, CP_factor, provenance_text) or (None, None, reason)
    """
    mech = str(mechanism or "").strip().lower()
    if mech not in {"competitive", "comp"}:
        return (None, None, "Only competitive supported for IC50→Ki (Cheng–Prusoff).")

    IC50 = _to_molar(ic50_str)
    if IC50 is None or IC50 <= 0:
        return (None, None, "Could not parse IC50 (need a number + unit).")

    S = _to_molar(substrate_str) if substrate_str else None
    Km = _to_molar(km_str) if km_str else None

    if S is None or Km is None or Km <= 0:
        prov = ("[S] and/or Km missing → using Ki≈IC50 (Cheng–Prusoff factor=1). "
                "Provide both [S] and Km to apply Cheng–Prusoff properly.")
        return (IC50, 1.0, prov)

    cp_factor = 1.0 + (S / Km)
    Ki = IC50 / cp_factor
    prov = (f"Ki from IC50 via Cheng–Prusoff (competitive): "
            f"[S]={S:.3e} M ({pretty_conc(S)}), Km={Km:.3e} M ({pretty_conc(Km)}), "
            f"factor={cp_factor:.3f}")
    return (Ki, cp_factor, prov)

# ---- Core converter for declared type ----
def convert_by_type(value_str, value_type, T=T_DEFAULT, assume_Ki_as_Kd=True,
                    substrate=None, km=None, mechanism="competitive"):
    """
    value_type ∈ {'Kd','Ki','Ka','pKd','pKi','IC50','kJ/mol','kcal/mol'}
    Returns dict: dG, Ki_M, Used_ChengPrusoff, CP_factor, Provenance, Notes
    """
    out = dict(dG=None, Ki_M=None, Used_ChengPrusoff=False, CP_factor=None, Provenance="", Notes="")
    t = value_type.strip().lower()

    # Energies
    if t in {"kj/mol", "kjmol"}:
        num = _first_number(value_str)
        if num is None:
            out["Notes"] = "Could not parse kJ/mol value"
            return out
        out["dG"] = round(num * KJ_TO_KCAL, 4)
        out["Provenance"] = "Converted from kJ/mol"
        return out

    if t in {"kcal/mol", "kcalmol"}:
        num = _first_number(value_str)
        if num is None:
            out["Notes"] = "Could not parse kcal/mol value"
            return out
        out["dG"] = round(num, 4)
        out["Provenance"] = "Already kcal/mol"
        return out

    # Association constant
    if t == "ka":
        Ka = _first_number(value_str)
        if Ka is None or Ka <= 0:
            out["Notes"] = "Could not parse Ka or Ka <= 0"
            return out
        dG = -R * T * np.log(Ka)  # ΔG = -RT ln(Ka)
        out["dG"] = round(dG, 4)
        out["Provenance"] = "From Ka: ΔG = -RT ln(Ka)"
        return out

    # pKd
    if t == "pkd":
        pkd = _first_number(value_str)
        if pkd is None:
            out["Notes"] = "Could not parse pKd"
            return out
        Kd = 10 ** (-pkd)
        dG = R * T * np.log(Kd)   # ΔG = RT ln(Kd)
        out["dG"] = round(dG, 4)
        out["Provenance"] = "From pKd: Kd=10^-pKd; ΔG = RT ln(Kd)"
        return out

    # pKi
    if t == "pki":
        pki = _first_number(value_str)
        if pki is None:
            out["Notes"] = "Could not parse pKi"
            return out
        if not assume_Ki_as_Kd:
            out["Notes"] = "pKi given; need assumption Ki≈Kd to compute ΔG"
            return out
        Kd = 10 ** (-pki)         # treat Ki≈Kd
        dG = R * T * np.log(Kd)
        out["dG"] = round(dG, 4)
        out["Provenance"] = "From pKi; assumed Ki≈Kd (competitive); ΔG = RT ln(Kd)"
        out["Notes"] = "Assumed Ki≈Kd (competitive)."
        return out

    # IC50
    if t == "ic50":
        Ki, fac, prov = convert_ic50_to_ki(value_str, substrate, km, mechanism)
        if Ki is None:
            out["Notes"] = prov
            return out
        dG = R * T * np.log(Ki)   # ΔG = RT ln(Kd≈Ki)
        out.update(dict(
            dG=round(dG,4),
            Ki_M=Ki,
            Used_ChengPrusoff=True,
            CP_factor=round(fac, 6),
            Provenance=prov,
            Notes="ΔG from Ki (IC50 via Cheng–Prusoff, competitive)."
        ))
        return out

    # Ki (treat as Kd)
    if t == "ki":
        Kd = _to_molar(value_str)
        if Kd is None or Kd <= 0:
            out["Notes"] = "Could not parse Ki or Ki <= 0"
            return out
        dG = R * T * np.log(Kd)
        out["dG"] = round(dG, 4)
        out["Provenance"] = "From Ki treated as Kd; ΔG = RT ln(Kd≈Ki)"
        out["Notes"] = "Assumed Ki≈Kd (competitive)."
        return out

    # Kd
    if t == "kd":
        Kd = _to_molar(value_str)
        if Kd is None or Kd <= 0:
            out["Notes"] = "Could not parse Kd or Kd <= 0"
            return out
        dG = R * T * np.log(Kd)
        out["dG"] = round(dG, 4)
        out["Provenance"] = "From Kd; ΔG = RT ln(Kd)"
        return out

    out["Notes"] = f"Unsupported type: {value_type}"
    return out

# ---- Simple menus ----
def menu_select(prompt, options_dict):
    """
    Show numbered menu from dict {num: (code, label)}.
    Return selected code or None.
    """
    print(prompt)
    for k in sorted(options_dict.keys()):
        print(f"{k}. {options_dict[k][1]}")
    sel = input("Enter choice: ").strip()
    if sel.isdigit() and int(sel) in options_dict:
        return options_dict[int(sel)][0]
    return None

# ---- Single value flow with adaptive prompts ----
def single_value_flow():
    type_menu = {
        1: ("Kd", "Kd (e.g., 3.5 nM)"),
        2: ("Ki", "Ki (e.g., 5 µM)"),
        3: ("Ka", "Ka (e.g., 1e8 M^-1)"),
        4: ("pKd", "pKd (e.g., 8.2)"),
        5: ("pKi", "pKi (e.g., 9.3)"),
        6: ("IC50", "IC50 (needs [S] and Km; competitive)"),
        7: ("kJ/mol", "Energy in kJ/mol"),
        8: ("kcal/mol", "Energy in kcal/mol"),
    }

    tcode = menu_select("Select the affinity/energy format:", type_menu)
    if not tcode:
        print("❌ Invalid choice.")
        return

    val_prompt = {
        "Kd":      "Enter Kd value (with units, e.g., 3.5 nM): ",
        "Ki":      "Enter Ki value (with units, e.g., 5 µM): ",
        "Ka":      "Enter Ka value (e.g., 1e8 M^-1): ",
        "pKd":     "Enter pKd value (e.g., 8.2): ",
        "pKi":     "Enter pKi value (e.g., 9.3): ",
        "IC50":    "Enter IC50 value (with units, e.g., 2 nM): ",
        "kJ/mol":  "Enter energy in kJ/mol: ",
        "kcal/mol":"Enter energy in kcal/mol: ",
    }
    val = input(val_prompt[tcode]).strip()

    if tcode == "IC50":
        print("\n📘 IC50 conversion (Cheng–Prusoff for competitive binding)")
        print("Tip: Press Enter to skip [S]/Km and assume Ki≈IC50 (factor=1).")
        mech = input("Mechanism (competitive assumed; press Enter to accept): ").strip() or "competitive"
        substrate = input("Enter [S] with units (e.g., 50 µM): ").strip()
        km = input("Enter Km with units (e.g., 5 µM): ").strip()

        Ki, fac, prov = convert_ic50_to_ki(val, substrate or None, km or None, mech)
        if Ki is None:
            print(f"❌ {prov}")
            return

        dG = R * T_DEFAULT * np.log(Ki)  # ΔG = RT ln(Kd≈Ki)
        print(f"\n✅ Ki = {Ki:.3e} M  (factor = {fac:.3f})")
        print(f"   {prov}")
        print(f"✅ ΔG = {dG:.4f} kcal/mol   (from Ki)\n")
        return

    # All other types
    res = convert_by_type(val, tcode)
    if res["dG"] is not None:
        print(f"\n✅ ΔG = {res['dG']:.4f} kcal/mol")
        if res.get("Ki_M"):
            print(f"Ki ≈ {res['Ki_M']:.3e} M ({pretty_conc(res['Ki_M'])})")
    else:
        print("❌ ΔG not computed.")
    print("\n— Details —")
    for k in ["Provenance","Used_ChengPrusoff","CP_factor","Ki_M","Notes"]:
        if res.get(k) is not None:
            print(f"{k}: {res.get(k)}")

# ---- CSV flow ----
def csv_flow():
    print("📁 Upload your CSV")
    uploaded = files.upload()
    for filename in uploaded.keys():
        df = pd.read_csv(filename, encoding="ISO-8859-1")
        print(f"Loaded: {filename} (rows={len(df)})")

        print("\nColumns:", list(df.columns))
        val_col = input("Enter the column name with values (e.g., 'Experimental_affinity'): ").strip()
        if val_col not in df.columns:
            print(f"❌ Column '{val_col}' not found.")
            continue

        has_type = input("Does your CSV have a 'Type' column per row? (y/n): ").strip().lower() == "y"
        if not has_type:
            type_menu = {
                1: ("Kd", "Kd"),
                2: ("Ki", "Ki"),
                3: ("Ka", "Ka"),
                4: ("pKd", "pKd"),
                5: ("pKi", "pKi"),
                6: ("IC50", "IC50 (needs [S], Km, competitive)"),
                7: ("kJ/mol", "kJ/mol"),
                8: ("kcal/mol", "kcal/mol"),
            }
            chosen_type = menu_select("Select the format for ALL rows:", type_menu)
            if not chosen_type:
                print("❌ Invalid choice.")
                continue
        else:
            chosen_type = None
            if "Type" not in df.columns:
                print("❌ You said Type exists, but 'Type' column not found.")
                continue

        # IC50 parameters: columns or single values
        ic50_mode = (chosen_type == "IC50") or (has_type and (df["Type"].astype(str).str.lower() == "ic50").any())
        sub_col = km_col = mech_col = None
        sub_fallback = km_fallback = mech_fallback = None

        if ic50_mode:
            has_cols = input("Do you have [S], Km, Mechanism columns in CSV? (y/n): ").strip().lower() == "y"
            if has_cols:
                sub_col = input("Enter [S] column name (e.g., 'Substrate'): ").strip()
                km_col = input("Enter Km column name (e.g., 'Km'): ").strip()
                mech_col = input("Enter Mechanism column name (press Enter if not present): ").strip()
                if mech_col == "":
                    mech_col = None
                for need in [sub_col, km_col]:
                    if need not in df.columns:
                        print(f"❌ Required column '{need}' not found.")
                        return
                if mech_col and mech_col not in df.columns:
                    print(f"⚠️ Mechanism column '{mech_col}' not found; defaulting to 'competitive'.")
                    mech_col = None
            else:
                sub_fallback = input("Enter a single [S] for ALL IC50 rows (e.g., 50 µM): ").strip()
                km_fallback = input("Enter a single Km for ALL IC50 rows (e.g., 5 µM): ").strip()
                mech_fallback = input("Mechanism for ALL IC50 rows (default 'competitive'): ").strip() or "competitive"

        # Prepare output columns
        out_cols = ["ΔG_kcal/mol","Used_ChengPrusoff","CP_factor","Ki_M","Provenance","Notes"]
        for c in out_cols:
            if c not in df.columns:
                df[c] = np.nan

        def row_convert(row):
            tcode = chosen_type if chosen_type else str(row.get("Type","")).strip()
            if not tcode:
                return pd.Series([np.nan, False, np.nan, np.nan, "No Type", "Missing Type"], index=out_cols)

            tcode_norm = tcode.strip()
            val = row[val_col]

            if str(tcode_norm).lower() == "ic50":
                if sub_col and km_col:
                    sub = row.get(sub_col, None)
                    kmv = row.get(km_col, None)
                    mech = row.get(mech_col, "competitive") if mech_col else "competitive"
                else:
                    sub = sub_fallback
                    kmv = km_fallback
                    mech = mech_fallback or "competitive"
                res = convert_by_type(val, "IC50", substrate=sub, km=kmv, mechanism=mech)
            else:
                res = convert_by_type(val, tcode_norm)

            return pd.Series([
                res["dG"],
                res["Used_ChengPrusoff"],
                res["CP_factor"],
                res["Ki_M"],
                res["Provenance"],
                res["Notes"]
            ], index=out_cols)

        df[out_cols] = df.apply(row_convert, axis=1)

        out_name = filename.rsplit(".csv",1)[0] + "_converted.csv"
        df.to_csv(out_name, index=False)
        print(f"\n✅ Done. Downloading {out_name}")
        files.download(out_name)

# ---- Top-level menu ----
print("Select an option:")
print("1. Convert a SINGLE value")
print("2. Convert a WHOLE CSV file")
main_choice = input("Enter 1 or 2: ").strip()

if main_choice == "1":
    single_value_flow()
elif main_choice == "2":
    csv_flow()
else:
    print("❌ Invalid choice.")
