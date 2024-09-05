import pandas as pd
from sciutil import *
import sys
#sys.path.append('/disk1/ariane/vscode/enzymetk/')
sys.path.append('/Users/arianemora/Documents/code/docko/docko/')

from docko import *
from pathlib import Path

smiles = 'CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C=C)C)C=C)C)C(=C3CCC(=O)[O-])C)CCC(=O)[O-].[Fe]' # heme
smiles = 'O=[P@@]([H])(C)C1=CC=CC=C1' # ent1
smiles = 'O=[P@]([H])(C)C1=CC=CC=C1' # ent2

# , 270, 240
""" Test cleaning a pdb structure so that we can have """
base_dir = '/Users/arianemora/Documents/code/docko/data/hayden/'
label = 'fold_2024_08_28_hayden'

#clean_one_pdb(f'{base_dir}fold_casey_04092024_model_0.cif', f'{base_dir}{label}/{label}.pdb')

# """ Test converting a pdb to a pdbqt file. """
#pdb_to_pdbqt_protein(f'{base_dir}{label}/{label}.pdb', f'{base_dir}{label}/{label}.pdbqt')
# Had to manually go obabel fold_2024_08_28_hayden.pdb -o pdbqt -xr -p 7.4 > fold_2024_08_28_hayden.pdbqt
score = dock(sequence='', protein_name=label, smiles=smiles, ligand_name='ent2', residues=[100], 
            protein_dir=f'{base_dir}', ligand_dir=f'{base_dir}', output_dir=f'{base_dir}vina_output/', pH=7.4,
            method='vina', size_x=5.0, size_y=5.0, size_z=5.0)
print(score)