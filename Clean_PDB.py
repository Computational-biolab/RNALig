<<<<<<< HEAD
Google colab: https://colab.research.google.com/drive/1LSxz-l2kczM9fi3W_mor72IlOP2vFq7R


from google.colab import files
import os

def clean_pdb_colab(input_pdb_path, output_pdb_path, ligand_resnames=None):
    standard_nucleotides = {'A', 'U', 'G', 'C'}
    excluded_resnames = {'HOH', 'NA', 'K', 'CL', 'MG', 'ZN', 'CA', 'MN', 'FE', 'CU', 'BR', 'IOD',
                         'NI', 'HG', 'AG', 'CD', 'AU', 'PB', 'RB'}
=======
import sys
import os

def clean_pdb(input_pdb_path, output_pdb_path=None, ligand_resnames=None):
    standard_nucleotides = {'A', 'U', 'G', 'C'}
    excluded_resnames = {
        'HOH', 'NA', 'K', 'CL', 'MG', 'ZN', 'CA', 'MN', 'FE', 'CU', 'BR', 'IOD',
        'NI', 'HG', 'AG', 'CD', 'AU', 'PB', 'RB'
    }

    if output_pdb_path is None:
        base_name = os.path.splitext(os.path.basename(input_pdb_path))[0]
        output_pdb_path = f"{base_name}_cleaned.pdb"
>>>>>>> 1e3dc63 (Initial commit of RNALig GUI)

    cleaned_lines = []

    with open(input_pdb_path, 'r') as infile:
        for line in infile:
            record = line[:6].strip()
            resname = line[17:20].strip()
            hetatm = (record == 'HETATM')
            is_std_nt = resname in standard_nucleotides

            if record == 'ATOM' and is_std_nt:
                cleaned_lines.append(line)
            elif hetatm:
                if ligand_resnames:
                    if resname in ligand_resnames:
                        cleaned_lines.append(line)
                elif resname not in excluded_resnames:
                    cleaned_lines.append(line)
            elif record in {'TER', 'END'}:
                cleaned_lines.append(line)

    with open(output_pdb_path, 'w') as outfile:
        outfile.writelines(cleaned_lines)

<<<<<<< HEAD
    return output_pdb_path

# Upload multiple files
uploaded = files.upload()

# Process each uploaded file
cleaned_files = []
for filename in uploaded:
    base_name = os.path.splitext(filename)[0]
    cleaned_path = f"{base_name}_cleaned.pdb"
    cleaned_file = clean_pdb_colab(filename, cleaned_path)
    cleaned_files.append(cleaned_file)

# List cleaned files
print("✅ Cleaned PDBs created:")
for f in cleaned_files:
    print(f" - {f}")
for cleaned_file in cleaned_files:
    files.download(cleaned_file)
=======
    print(f"✅ Cleaned PDB saved: {output_pdb_path}")
    return output_pdb_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python Clean_PDB.py <input.pdb>")
        sys.exit(1)

    input_file = sys.argv[1]
    clean_pdb(input_file)
>>>>>>> 1e3dc63 (Initial commit of RNALig GUI)
