Google colab link for Affinity_convertor.py : https://colab.research.google.com/drive/1qyQcHqD_Qow6OKWZGbv32QrMbhtWwDU6?usp=sharing
#Import required libraries
import math
import re
# Input Constants
R = 1.987e-3  # kcal/(mol·K)
T = 298       # Kelvin

def kd_to_deltaG(Kd):
    """Calculate ΔG from Kd (in M)."""
    if Kd <= 0:
        raise ValueError("Kd must be positive")
    return R * T * math.log(Kd)
  #Conversion
def convert_to_molar(value_with_unit):
    value_with_unit = value_with_unit.strip().lower()
    value_with_unit = re.sub(r'^(kd|ki|ic50)=', '', value_with_unit)
    match = re.match(r'([0-9.]+)\s*([mun]?m)?', value_with_unit)
    if not match:
        raise ValueError(f"Invalid affinity format: {value_with_unit}")
    value = float(match.group(1))
    unit = match.group(2)
    if unit == "mm":
        return value * 1e-3
    elif unit == "um":
        return value * 1e-6
    elif unit == "nm":
        return value * 1e-9
    elif unit == "m" or unit is None:
        return value
    else:
        raise ValueError(f"Unknown unit: {unit}")

def pX_to_molar(pX):
    return 10 ** (-pX)

def molar_to_pX(molar):
    return -math.log10(molar)

def convert_and_report_single(input_value, input_type):
    try:
        if input_type in ["Kd/Ki/IC50"]:
            kd_molar = convert_to_molar(input_value)
            pX_value = molar_to_pX(kd_molar)
        elif input_type in ["pKd/pKi/pIC50"]:
            pX_value = float(input_value)
            kd_molar = pX_to_molar(pX_value)
        else:
            raise ValueError("Input type must be 'Kd/Ki/IC50' or 'pKd/pKi/pIC50'")
        
        deltaG = kd_to_deltaG(kd_molar)
        nM = kd_molar * 1e9
        uM = kd_molar * 1e6
        mM = kd_molar * 1e3
        
        print("\nConverted affinity values:")
        print(f"  Input: {input_value} ({input_type})")
        print(f"  Kd (M): {kd_molar:.3e}")
        print(f"  Kd (mM): {mM:.3e}")
        print(f"  Kd (uM): {uM:.3e}")
        print(f"  Kd (nM): {nM:.3e}")
        print(f"  pX (pKd/pKi/pIC50): {pX_value:.3f}")
        print(f"  ΔG (kcal/mol): {deltaG:.4f}\n")

    except Exception as e:
        print(f"Error: {e}\n")

def process_row(affinity_str, pX_str):
    try:
        kd_molar = convert_to_molar(affinity_str)
    except Exception as e:
        kd_molar = None
        err1 = str(e)
    else:
        err1 = None

    try:
        pX_val = float(pX_str)
    except Exception as e:
        pX_val = None
        err2 = str(e)
    else:
        err2 = None

    if kd_molar is not None:
        calc_pX = molar_to_pX(kd_molar)
        calc_dG = kd_to_deltaG(kd_molar)
    else:
        calc_pX = None
        calc_dG = None

    if pX_val is not None:
        calc_kd_from_pX = pX_to_molar(pX_val)
        calc_dG_from_pX = kd_to_deltaG(calc_kd_from_pX)
    else:
        calc_kd_from_pX = None
        calc_dG_from_pX = None

    return {
        "Affinity_input": affinity_str,
        "Kd(M)_from_affinity": kd_molar,
        "Calc_pX_from_affinity": calc_pX,
        "Calc_dG_from_affinity": calc_dG,
        "pX_input": pX_val,
        "Kd(M)_from_pX": calc_kd_from_pX,
        "Calc_dG_from_pX": calc_dG_from_pX,
        "Affinity_error": err1,
        "pX_error": err2
    }

def convert_and_report_list():
    print("\nPaste your data (tab-separated: AffinityValue <tab> pXValue), one per line.")
    print("Type 'exit' or an empty line to finish input.\n")
    lines = []
    while True:
        line = input()
        if line.strip().lower() == "exit" or line.strip() == "":
            break
        lines.append(line.strip())

    if not lines:
        print("No data entered.\n")
        return

    print("\nProcessing...\n")
    print(f"{'Affinity':15} {'Kd(M)_from_aff':12} {'pX_from_aff':10} {'ΔG_from_aff':12} {'pX_input':10} {'Kd(M)_from_pX':12} {'ΔG_from_pX':12}")
    for line in lines:
        if "\t" not in line:
            print(f"Skipping invalid line (missing tab): {line}")
            continue
        affinity_str, pX_str = line.split("\t", 1)
        res = process_row(affinity_str, pX_str)
        aff_val = f"Err:{res['Affinity_error']}" if res["Affinity_error"] else f"{res['Kd(M)_from_affinity']:.3e}"
        pX_aff = f"{res['Calc_pX_from_affinity']:.2f}" if res["Calc_pX_from_affinity"] is not None else "Err"
        dG_aff = f"{res['Calc_dG_from_affinity']:.3f}" if res["Calc_dG_from_affinity"] is not None else "Err"
        pX_inp = f"Err:{res['pX_error']}" if res["pX_error"] else f"{res['pX_input']:.2f}"
        kd_pX = f"{res['Kd(M)_from_pX']:.3e}" if res["Kd(M)_from_pX"] is not None else "Err"
        dG_pX = f"{res['Calc_dG_from_pX']:.3f}" if res["Calc_dG_from_pX"] is not None else "Err"

        print(f"{affinity_str:15} {aff_val:12} {pX_aff:10} {dG_aff:12} {pX_inp:10} {kd_pX:12} {dG_pX:12}")
    print()

def main():
    print("Binding Affinity Converter")
    print("--------------------------")
    while True:
        print("Choose an option:")
        print("1. Convert one single affinity value (Kd/Ki/IC50 or pKd/pKi/pIC50)")
        print("2. Convert a list of affinity and pX values (tab-separated per line)")
        print("Type 'exit' to quit.")
        choice = input("Enter your choice (1 or 2): ").strip().lower()
        if choice == "exit" or choice == "":
            print("Goodbye!")
            break
        elif choice == "1":
            print("\nChoose input type:")
            print("a. Kd/Ki/IC50 (with units like nM, uM, mM)")
            print("b. pKd/pKi/pIC50 (numeric values)")
            input_type_choice = input("Enter 'a' or 'b': ").strip().lower()
            if input_type_choice == "a":
                input_type = "Kd/Ki/IC50"
            elif input_type_choice == "b":
                input_type = "pKd/pKi/pIC50"
            else:
                print("Invalid choice. Try again.\n")
                continue
            value = input(f"Enter your affinity value ({input_type}): ").strip()
            if value.lower() == "exit" or value == "":
                print("Goodbye!")
                break
            convert_and_report_single(value, input_type)
        elif choice == "2":
            convert_and_report_list()
        else:
            print("Invalid choice. Please enter 1 or 2.\n")

if __name__ == "__main__":
    main()
