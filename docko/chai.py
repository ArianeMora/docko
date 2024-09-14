###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

"""
Wrapper around binding prediction tools to get the optimal state of a protein for predicting a reaction.
A discusting mix of hacked together elements <3 true bioinformatics style, I am going to coding hell and
apologise profusely.
"""
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
import re
import csv
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile
from tqdm import tqdm
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import PDBIO, Select, MMCIFIO
from Bio import PDB
from multiprocessing.dummy import Pool as ThreadPool
import pandas as pd
from rdkit.Chem import AllChem
import requests

import logging
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
from pathlib import Path
from Bio import SeqIO
import os
import torch
import pandas as pd
from rdkit import Chem


import argparse
from chai_lab.chai1 import run_inference

from rdkit.Chem import AllChem as Chem

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

