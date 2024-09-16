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
from multiprocessing.dummy import Pool as ThreadPool
import os
from docko.vina import *
from docko.diffdock import *


def calculate_docking_affinities_across_dataset(df, output_dir, protein_dir, ligand_dir, output_file, num_threads=20,
                                                size_x=10.0, size_y=10.0, size_z=10.0, pH=7.4, method='vina'):
    pool = ThreadPool(num_threads)
    data = []
    with open(os.path.join(output_dir, output_file), 'w+') as fout:
        fout.write(f'uniprot_id,best_score,all_scores\n')
    num = int(len(df) / num_threads)
    for i in range(0, len(df), num):
        end_i = i + num if i + num < len(df) else len(df)
        sub_df = df.iloc[i: end_i]
        # df, output_dir, output_file, protein_dir, ligand_dir, size_x, size_y, size_z 
        sub_data = [sub_df, output_dir, output_file, protein_dir, ligand_dir, size_x, size_y, size_z, pH, method]
        data.append(sub_data)

    # Thread it
    pool.map(run_vina_docking_thread, data)


def dock(sequence: str, protein_name, smiles: str, ligand_name: str, residues: list, protein_dir: str, ligand_dir: str,
         output_dir: str, pH: float, method: str, size_x=5.0, size_y=5.0, size_z=5.0, num_modes=9, exhaustivenes=32):
    """ Dock a smiles to a pdb file using vina and dockstring to automate some of the things. """

    # Step 1: Check if the protein exists, if not make it pretty and save it to the protein folder
    protein_name, protein_pdb_file, protein_pdbqt = format_pdb(sequence, protein_name, protein_dir, pH)

    # Step 2: Check if the ligand exists, if not make it pretty and save it to the ligand folder
    ligand_pdbqt, ligand_sdf = format_ligand(smiles, ligand_name, ligand_dir, pH)

    # step 3: dock
    if method == 'vina' or method == 'ad4':
        affinities = dock_vina(ligand_pdbqt, protein_pdbqt, ligand_name, protein_name, output_dir, residues,
                               size_x=size_x,
                               size_y=size_y, size_z=size_z, num_modes=num_modes, exhaustivenes=exhaustivenes,
                               method=method)
    elif method == 'diffdock':
        dock_diffdock(output_dir, ligand_sdf, protein_pdbqt.replace('.pdbqt', '.pdb'))
    return None
