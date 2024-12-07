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

from pathlib import Path
import os
from glob import glob
import torch
from ast import literal_eval
import pandas as pd
from docko.helpers import *


def run_boltz(label: str, seq: str, smiles: str, output_dir: str, cofactor_smiles:list|str = "", joinsubcofactor:bool=True) -> None:
    """
    Run CHAI on a single protein-ligand pair.

    Args:
    - label (str): label for the protein-ligand pair
    - seq (str): sequence of the protein
    - smiles (str): SMILES string of the ligand
    - output_dir (str): output directory
    - cofactor_smiles (list or str): list of SMILES strings of cofactors, default is ""
    - joinsubcofactor (bool): whether to join the substrate and cofactor in the same fasta file, default is True
    """

    # make sure output dir is dir
    output_subdir = os.path.join(output_dir, label)

    # Need to clean up the sequence
    seq = seq.strip().replace("*", "").replace(" ", "").upper()

    example_fasta = f">A|protein\n{seq}\n"

    letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    letter_position = 1
    if cofactor_smiles != "":
        # convert if cofactor_smiles to a list if it is a string
        if isinstance(cofactor_smiles, str):
            # use ast.literal_eval to convert string to list
            try:
                cofactor_smiles = literal_eval(cofactor_smiles)
            except:
                cofactor_smiles = [cofactor_smiles]

        # add cofactor SMILES to the fasta
        for cofactor_smile in cofactor_smiles:
            # Need to do new chains for each one
            example_fasta += f">{letters[letter_position]}|smiles\n{cofactor_smile}\n"
            letter_position += 1

    # now add substrate
    if smiles:
        smiles = canonicalize_smiles(smiles)
        example_fasta += f">{letters[letter_position]}|smiles\n{smiles}\n"

    if not os.path.exists(output_subdir):
        os.system(f"mkdir {output_subdir}")
        print(output_subdir)
        output_subdir = Path(output_subdir)
        # CHeck this is OK
        fasta_path = Path(f"{output_subdir}/{label}.fasta")
        fasta_path.write_text(example_fasta)
        os.system(f'boltz predict {f"{output_subdir}/{label}.fasta"} --use_msa_server --accelerator gpu --out_dir {output_subdir} ')
    else:
        print(f"Output directory exists: {output_subdir}")

    # get name of the output cif or pdb files
    output_strcut_files = glob(f"{output_subdir}/*.cif") + glob(f"{output_subdir}/*.pdb")
    
    # rename the output files cif or pdb files
    for output_strcut_file in output_strcut_files:
        os.rename(output_strcut_file, output_strcut_file.replace("pred.model_idx", label))

    # for npz files do the same
    output_scores_files = glob(f"{output_subdir}/*.npz")

    for output_scores_file in output_scores_files:
        os.rename(output_scores_file, output_scores_file.replace("scores.model_idx", label))