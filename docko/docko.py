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



def calculate_docking_affinities_across_dataset(df, output_dir, protein_dir, ligand_dir, output_file, num_threads=20, 
                                                size_x=10.0, size_y=10.0, size_z=10.0, pH=7.4, method='ad4'):
    pool = ThreadPool(num_threads)
    data = []
    with open(os.path.join(output_dir, output_file), 'w+') as fout:
        fout.write(f'uniprot_id,best_score,all_scores\n')
    num = int(len(df) / num_threads)
    for i in range(0, len(df), num):
        end_i = i + num if i + num < len(df) else len(df)
        sub_df = df.iloc[i: end_i]
        # df, output_dir, output_file, protein_dir, ligand_dir, size_x, size_y, size_z 
        sub_data = [sub_df, output_dir, output_file,  protein_dir, ligand_dir, size_x, size_y, size_z, pH, method]
        data.append(sub_data)

    # Thread it
    pool.map(run_vina_docking_thread, data)


def run_vina_docking_thread(args):
    """ Docking threading """
    # For each ID we create a new folder so that we don't have issues with the names of files
    df, output_dir, output_file, protein_dir, ligand_dir, size_x, size_y, size_z, pH, method = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9]
    # we're going to read in the file and if it already has an entry for this uniprot we won't do it
    done_df = pd.read_csv(os.path.join(output_dir, output_file))
    done_ids = list(done_df['uniprot_id'].values)
    for seq_id, seq, active_site, smiles, ligand_name in tqdm(df[['Entry', 'Sequence', 'Residue', 'Smiles', 'LigandName']].values):
        if seq_id not in done_ids:
            # 1) make a directory
            ligand_name = ligand_name.replace(' ', '_')
            try:
                affinities = dock(seq, seq_id, smiles, ligand_name, active_site.split('|'), protein_dir, ligand_dir, output_dir, pH,
                size_x=size_x, size_y=size_y, size_z=size_z, method=method)
                if method == 'vina':
                    write_vina_affinities(seq_id, ligand_name, affinities, os.path.join(output_dir, 'scores.txt'))
            except:
                print('-------------------------------------------------------------------------')
                print(seq_id, ligand_name)


def dock(sequence: str, protein_name, smiles: str, ligand_name: str, residues: list, protein_dir: str, ligand_dir: str,
            output_dir: str, pH: float, method: str, size_x=5.0, size_y=5.0, size_z=5.0, num_modes=9, exhaustivenes=32):
    """ Dock a smiles to a pdb file using vina and dockstring to automate some of the things. """

    # Step 1: Check if the protein exists, if not make it pretty and save it to the protein folder
    protein_pdbqt = format_pdb(sequence, protein_name, protein_dir, pH)

    # Step 2: Check if the ligand exists, if not make it pretty and save it to the ligand folder
    ligand_pdbqt, ligand_sdf = format_ligand(smiles, ligand_name, ligand_dir, pH)

    # step 3: dock
    if method == 'vina' or method == 'ad4':
        affinities = dock_vina(ligand_pdbqt, protein_pdbqt, ligand_name, protein_name, output_dir, residues, size_x=size_x, 
                               size_y=size_y, size_z=size_z, num_modes=num_modes, exhaustivenes=exhaustivenes, method=method)
    elif method == 'diffdock':
        dock_diffdock(output_dir, ligand_sdf, protein_pdbqt.replace('.pdbqt', '.pdb'))
    return None