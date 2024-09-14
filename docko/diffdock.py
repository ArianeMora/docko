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



def dock_diffdock(output_dir: str, ligand_sdf, pdb_file: str):
    """ Dock using diffdock. """
    os.chdir(output_dir)
    print(f'conda run -n TRILL python3 trill example 1 dock DiffDock {pdb_file} {ligand_sdf}')
    os.system(f'conda run -n TRILL python3 trill example 1 dock DiffDock {pdb_file} {ligand_sdf}')
    # Also now need to move and rename it

def diffdock_df(df, output_dir, smiles, num_threads=20):
    pool = ThreadPool(num_threads)
    data = []
    num = int(len(df) / num_threads)
    smiles_sdf = smiles_to_sdf(smiles, os.path.join(output_dir, f'smiles.sdf'))
    for i in range(0, len(df), num):
        end_i = i + num if i + num < len(df) else len(df)
        sub_df = df.iloc[i: end_i]
        sub_data = [sub_df, output_dir, smiles_sdf]
        data.append(sub_data)

    # Thread it
    pool.map(run_diffdock_docking_thread, data)

def run_diffdock_docking_thread(args):
    """ Docking threading """
    # For each ID we create a new folder so that we don't have issues with the names of files
    df, output_dir, output_file, protein_dir, ligand_dir, pH = args[0], args[1], args[2], args[3], args[4], args[5], args[6]
    # we're going to read in the file and if it already has an entry for this uniprot we won't do it
    done_df = pd.read_csv(os.path.join(output_dir, output_file))
    done_ids = list(done_df['Entry'].values)
    for seq_id, seq, active_site, smiles, ligand_name in tqdm(df[['Entry', 'Sequence', 'Residue', 'Smiles', 'LigandName']].values):
        if seq_id not in done_ids:
            # 1) make a directory
            try:
                score, aux = dock(seq, seq_id.strip().replace(' ', '_'), smiles, ligand_name.strip().replace(' ', '_'), active_site.split('|'), protein_dir, ligand_dir, output_dir, pH, 'diffdock')
                with open(os.path.join(output_dir, output_file), 'a+') as fout:
                    aux = aux['affinities']
                    fout.write(f'{seq_id},{ligand_name},{aux[0]},{"|".join([str(x) for x in aux])}\n')
            except:
                print(seq_id, ligand_name)

