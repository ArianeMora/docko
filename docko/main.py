import pandas as pd
from sciutil import *
import sys
#sys.path.append('/disk1/ariane/vscode/enzymetk/')
sys.path.append('/Users/arianemora/Documents/code/docko/docko/')

from docko import *
from pathlib import Path

smiles = 'CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C=C)C)C=C)C)C(=C3CCC(=O)[O-])C)CCC(=O)[O-].[Fe]'
# smiles = 'C=CC1=C(N2[Fe]345=NC[C@]([H])(O)CO)C=C6[N+]3=C(C(CCC([O-])=O)=C6C)C=C7N4C(C(C)=C7CCC([O-])=O)=CC8=[N+]5C(C(C=C)=C8C)=CC2=C1C' #
# smiles = 'O[C@@H](CN=[N+]=[N-])CO' # nitrene
# smiles = 'C=CC(C(N1)=O)=CNC1=O' # 5 vinyluracil: C=CC(C(N1)=O)=CNC1=O

# , 270, 240
""" Test cleaning a pdb structure so that we can have """
base_dir = '/Users/arianemora/Documents/code/docko/data/casey/'
label = 'docking_heme'

#clean_one_pdb(f'{base_dir}fold_casey_04092024_model_0.cif', f'{base_dir}{label}/{label}.pdb')

# """ Test converting a pdb to a pdbqt file. """
#pdb_to_pdbqt_protein(f'{base_dir}{label}/{label}.pdb', f'{base_dir}{label}/{label}.pdbqt')

score = dock(sequence='', protein_name=label, smiles=smiles, ligand_name='heme', residues=[66, 116], 
            protein_dir=f'{base_dir}', ligand_dir=f'{base_dir}', output_dir=f'{base_dir}vina_output/', pH=7.0,
            method='vina', size_x=20.0, size_y=20.0, size_z=20.0)
print(score)