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
import torch
import pandas as pd
from chai_lab.chai1 import run_inference
from docko.helpers import *


def run_chai(label, seq, smiles, output_dir, cofactor_smiles=""):
    # make sure output dir is dir
    output_dir = os.path.normpath(output_dir) + "/"
    # Need to
    seq = seq.strip().replace("*", "").replace(" ", "").upper()
    if cofactor_smiles != "":
        cofactor_fasta = f">ligand|{label}-cofactor\n{cofactor_smiles}\n"
    else:
        cofactor_fasta = ""

    example_fasta = f">protein|{label}\n{seq}\n>ligand|{label}-substrate\n{smiles}\n{cofactor_fasta}"

    if not os.path.isdir(f"{output_dir}{label}"):
        os.system(f"mkdir {output_dir}{label}")
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


def run_chai_df(
    output_dir,
    filename,
    entry_column="Entry",
    seq_column="Sequence",
    ligand_column="Substrate",
    cofactor_column="Cofactor",
):
    """Run on a dataframe!"""
    df = pd.read_csv(filename)

    if cofactor_column not in df.columns:
        df[cofactor_column] = ""

    for label, seq, smiles, cofactor_smiles in df[
        [entry_column, seq_column, ligand_column, cofactor_column]
    ].values:
        run_chai(
            label=label,
            seq=seq,
            smiles=smiles,
            output_dir=output_dir,
            cofactor_smiles=cofactor_smiles,
        )