from pathlib import Path
from Bio import SeqIO
import os
import torch
import pandas as pd
from rdkit import Chem


import argparse
from chai_lab.chai1 import run_inference


def canonicalize_smiles(smiles_string):
    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles
    else:
        return None

def run_chai(output_dir, filename, entry_column='Entry', seq_column='Sequence', ligand_column='Substrate'):

    df = pd.read_csv(filename)

    for label, seq, smiles in df[[entry_column, seq_column, ligand_column]].values:
        # Need to 
        seq = seq.strip().replace('*', '').replace(' ', '').upper()
        example_fasta = f">protein|{label}\n{seq}\n>ligand|{label}\n{smiles}\n"
        if not os.path.isdir(f'{output_dir}{label}'):
            os.system(f'mkdir {output_dir}{label}')
            smiles = canonicalize_smiles(smiles)
            # CHeck this is OK
            if smiles:
                fasta_path = Path(f"{output_dir}{label}/{label}.fasta")
                fasta_path.write_text(example_fasta)
                
                output_paths = run_inference(
                    fasta_file=fasta_path,
                    output_dir=Path(f"{output_dir}{label}/"),
                    # 'default' setup
                    num_trunk_recycles=3,
                    num_diffn_timesteps=200,
                    seed=42,
                    device=torch.device("cuda:0"),
                    use_esm_embeddings=True,
                )

def parse_args():
    parser = argparse.ArgumentParser(description="Run Chai on a sequence")
    parser.add_argument('-out', '--out', required=True, help='Path to the output directory')
    parser.add_argument('-df', '--df', type=str, required=True, help='Fasta of the file of interest')
    return parser.parse_args()

def main():
    args = parse_args()
    run_chai(args.out, args.df)

if __name__ == "__main__":
    main()