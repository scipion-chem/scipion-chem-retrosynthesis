import os
import argparse
from rdkit import Chem

def save_molecule_to_sdf(smiles, output_directory, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        os.makedirs(output_directory, exist_ok=True)
        output_file = os.path.join(output_directory, filename)
        
        try:
            with Chem.SDWriter(output_file) as writer:
                writer.write(mol)
            print(f"Molecule saved to: {output_file}")
        except Exception as e:
            print(f"Error writing SDF file: {e}")
    else:
        print("Error: SMILES string could not be converted to Mol.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a SMILES string to an SDF file.")
    parser.add_argument('--smiles', type=str, required=True, help="SMILES string to convert.")
    parser.add_argument('--output-directory', type=str, required=True, help="Directory to save the SDF file.")
    parser.add_argument('--filename', type=str, required=True, help="Name of the output SDF file.")
    
    args = parser.parse_args()

    save_molecule_to_sdf(args.smiles, args.output_directory, args.filename)



