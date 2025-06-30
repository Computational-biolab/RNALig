Google Colab link: https://colab.research.google.com/drive/1n_mSo7R1eLDprwTz2SKgI8LWizaiYI23#forceEdit=true&sandboxMode=true
# Import libraries
import pandas as pd
import numpy as np
import re
from google.colab import files
from math import log

# === Constants ===
R = 0.001987  # kcal/mol·K
T = 298       # Kelvin

# === Import libraries ===
import pandas as pd
import numpy as np
import re
from google.colab import files

# === Constants ===
R = 0.001987  # kcal/mol·K
T = 298       # Kelvin

# === ΔG Conversion Function ===
def convert_affinity_to_kcal(value):
    try:
        value = str(value).strip().lower().replace(" ", "")
        value = value.replace("μ", "u").replace("µ", "u")  # Normalize to 'u'

        # Log-based units
        if 'pkd' in value or 'pki' in value:
            pk_val = float(re.findall(r"\d+\.?\d*", value)[0])
            kd_molar = 10 ** (-pk_val)

        # Energy units
        elif 'kj/mol' in value:
            kj_val = float(re.findall(r"-?\d+\.?\d*", value)[0])
            return round(kj_val * 0.239006, 4)
        elif 'kcal/mol' in value:
            return round(float(re.findall(r"-?\d+\.?\d*", value)[0]), 4)

        else:
            # Numeric extraction
            num_match = re.search(r"[-+]?\d*\.?\d+(e[-+]?\d+)?", value)
            if not num_match:
                return None
            num = float(num_match.group())

            # Concentration units
            if 'pm' in value:
                kd_molar = num * 1e-12
            elif 'nm' in value:
                kd_molar = num * 1e-9
            elif 'um' in value:
                kd_molar = num * 1e-6
            elif 'mm' in value:
                kd_molar = num * 1e-3
            elif re.fullmatch(r"[0-9.e+-]+m", value):
                kd_molar = num
            elif 'ki' in value or 'ic50' in value or 'kd' in value:
                kd_molar = num
            else:
                kd_molar = num

        # Final ΔG calculation
        if kd_molar <= 0:
            return None

        ln_kd = np.log(kd_molar)
        deltaG = R * T * ln_kd  # Forcefully ensure negative ΔG for binding
        print(f"DEBUG: kd={kd_molar}, ln(kd)={ln_kd}, ΔG={deltaG}")
        return round(deltaG, 4)

    except Exception as e:
        print(f"❌ Error: {e}")
        return None
# === Menu Interface ===
print("Select an option:")
print("1. Upload CSV file with PDB_ID and affinity values")
print("2. Convert a single affinity value manually")

choice = input("Enter 1 or 2: ")

if choice == "1":
    print("📁 Upload your CSV (must contain columns: PDB_ID, Experimental_affinity):")
    uploaded = files.upload()

    for filename in uploaded.keys():
        try:
            df = pd.read_csv(filename, encoding='ISO-8859-1')

            if 'Experimental_affinity' not in df.columns or 'PDB_ID' not in df.columns:
                print("❌ CSV must have columns: 'PDB_ID' and 'Experimental_affinity'")
            else:
                df['BA_kcal/mol'] = df['Experimental_affinity'].apply(convert_affinity_to_kcal)
                output_name = filename.replace(".csv", "_converted.csv")
                df.to_csv(output_name, index=False)
                print(f"\n✅ Conversion complete. Downloading: {output_name}")
                files.download(output_name)

        except Exception as e:
            print(f"❌ Error processing file: {e}")

elif choice == "2":
    affinity_input = input("Enter the experimental affinity (e.g., 2.6µM, 5.24 pKd, 32 kJ/mol, 3.5 Ki): ")
    converted = convert_affinity_to_kcal(affinity_input)
    if converted is not None:
        print(f"\n✅ Converted value in kcal/mol: {converted:.4f}")
    else:
        print("❌ Unable to parse or convert the input value.")

else:
    print("❌ Invalid choice. Please run the script again.")
