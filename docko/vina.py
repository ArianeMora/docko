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

import requests

from rdkit.Chem import AllChem as Chem


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

def make_config_for_vina(pdb_file, ligand_file, residues, output_file: str, size_x=5.0, size_y=5.0, size_z=5.0,
                         num_modes=9, exhaustivenes=32):
    """ Make the config file for Vina based on coordinates for a specific residue (or multiple). """
    coords = []
    for position in residues:
        chain_id, coordinates = get_coordinates_without_chain_id(pdb_file, int(position))
        coords.append(coordinates)
    print(coords)
    coords = calculate_centroid(coords)
    # Now we want to save this to a config file that has everything we need
    with open(f'{output_file}', 'w+') as fout:
        #receptor = enzyme.pdbqt
        fout.write(f'receptor = {pdb_file}\n')
        fout.write(f'ligand = {ligand_file}\n')
        fout.write(f'center_x = {coords[0]}\n')
        fout.write(f'center_y = {coords[1]}\n')
        fout.write(f'center_z = {coords[2]}\n')
        fout.write(f'size_x = {size_x}\n')
        fout.write(f'size_y = {size_y}\n')
        fout.write(f'size_z = {size_z}\n')
        fout.write(f'num_modes = {num_modes}\n')
        fout.write(f'exhaustiveness = {exhaustivenes}\n')


def dock_vina(ligand_pdbqt: str, protein_pdbqt: str, ligand_name: str, protein_name: str, output_dir: str, 
         residues: str, size_x=5.0, size_y=5.0, size_z=5.0, num_modes=9, exhaustivenes=32, method='vina',
         num_cpus: Optional[int] = None, seed=42) -> Tuple[Optional[float], Dict[str, Any]]:
        """
        Given a molecule, this method will return a docking score against the current target.

        :param smiles: SMILES string of ligand
        :param pH: pH at which the docking should take place (default: 7.4, don't change unless you know what you are doing)
        :param num_cpus: number of CPUs cores available to AutoDock Vina
        :param seed: random seed for conformation generation and docking
        :param verbose: increase verbosity of log messages
        :return: docking score and dictionary containing all poses and binding free energies
        """
        conf_path = os.path.join(output_dir, f'{protein_name}-{ligand_name}_conf.txt')

        # Step 1: make config file
        make_config_for_vina(protein_pdbqt, ligand_pdbqt, residues, conf_path, size_x=size_x,
                            size_y=size_y, size_z=size_z, num_modes=num_modes, exhaustivenes=exhaustivenes)

        # Auxiliary files
        vina_logfile = os.path.join(output_dir, f'{protein_name}-{ligand_name}_log.txt')

        docked_ligand_pdb = os.path.join(output_dir, f'{protein_name}-{ligand_name}.pdb')
                # Step 2: Dock using our nice new baby's # def dock_pdbqt(ligand_pdbqt, protein_pdbqt, conf_path, log_path, out_path, seed, num_cpus: Optional[int] = None) -> None:

        if method == 'ad4':
            vina_logfile = dock_ad4_pdbqt(ligand_pdbqt, protein_pdbqt, vina_logfile, output_dir, protein_name, ligand_name)
        elif method == 'vina':
            dock_autodock_pdbqt(conf_path, vina_logfile, docked_ligand_pdb, seed=seed, num_cpus=num_cpus)

        # Step 5: convert back to a pdb structure since whotf uses pdbqt
        convert_pdbqt_to_pdb(pdbqt_file=docked_ligand_pdb, pdb_file=docked_ligand_pdb, disable_bonding=True)
        
        # Step 6: Parse scores
        return vina_logfile

def pdb_to_pdbqt_protein(input_filename, output_filename=None, ph=7.4):
    """
    Convert a pdb file to a pdbqt file.
    """
    output_filename = output_filename if output_filename else input_filename.replace('.pdb', '.pdbqt')
    os.system(f'obabel {input_filename} -xr -p {ph} --partialcharge gasteiger -O {output_filename}')
    # Now we also want to be cheeky and remove any secondary model parts from the file
    # This is a hacky way to keep a bound heme or something, seems to work fine.
    lines = []
    with open(output_filename, 'r+') as fin:
        for line in fin:
            if line.split(' ')[0] not in ['MODEL', 'TER', 'ENDMDL', 'REMARK']:
                lines.append(line)
    with open(output_filename, 'w+') as fout:
        for line in lines:
            if 'ENDMDL' not in line:
                fout.write(line)
        fout.write('TER\n')


def write_vina_affinities(protein_name, ligand_name, log_file_path, output_csv_path):
        # Read the content of the log file
    with open(log_file_path, 'r') as file:
        log_content = file.readlines()

    # Parse the relevant lines
    data_lines = []
    add_lines = False
    for line in log_content:
        if line.strip() and line.replace(' ', '')[0] == '1':
            add_lines = True
        # Pretty much you add all the lines
        if add_lines:
            data_lines.append(line.strip())

    # Write the parsed data to the CSV file
    with open(output_csv_path, 'a+') as fout:
        # Write the data rows
        for row in data_lines:
            print(row)
            fout.write(f'{protein_name}\t{ligand_name}\t{row}\n')


def run_mmff94_opt(mol: Chem.Mol, max_iters: int) -> Chem.Mol:
    """
    Optimize molecular structure with MMFF94 force field.

    :param mol: molecular structure to be optimized
    :param max_iters: maximum number of structure optimization iterations
    :return: optimized molecule
    """
    Chem.MMFFSanitizeMolecule(mol)
    opt_result = Chem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94', maxIters=max_iters)
    return opt_result


def run_uff_opt(mol: Chem.Mol, max_iters: int) -> Chem.Mol:
    """
    Optimize molecular structure with UFF.

    :param mol: molecular structure to be optimized
    :param max_iters: maximum number of structure optimization iterations
    :return: optimized molecule
    """
    opt_result = Chem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    return opt_result

def refine_mol_with_ff(mol, max_iters=1000) -> Chem.Mol:
    """
    Optimize molecular structure. Try MMFF94 first, use UFF as a backup.

    :param mol: molecular structure to be optimized
    :param max_iters: maximum number of structure optimization iterations
    :return: optimized molecule
    """
    if Chem.MMFFHasAllMoleculeParams(mol):
        try:
            opt_mol = run_mmff94_opt(mol, max_iters=max_iters)
        except Chem.rdchem.KekulizeException as exception:
            logging.info(f'Ligand optimization with MMFF94 failed: {exception}, trying UFF')
            opt_mol = run_uff_opt(mol, max_iters=max_iters)
    elif Chem.UFFHasAllMoleculeParams(mol):
        opt_mol = run_uff_opt(mol, max_iters=max_iters)
        return opt_mol
    # Otherwise it was not able to be optimised...
    return mol


def convert_pdbqt_to_pdb(pdbqt_file: str, pdb_file: str, disable_bonding=False) -> None:
    """
    Convert a PDBQT file to a PDB file with Open Babel.

    :param pdbqt_file: path to the PDBQT input file
    :param pdb_file: path to the PDB output file
    :param disable_bonding: disable automatic bonding with Open Babel
    """
    cmd_args = [
        'obabel',
        '-ipdbqt', pdbqt_file,
        '-opdb',
        '-O', pdb_file,
    ]
    if disable_bonding:
        # "a" = read option
        # "b" = disable automatic bonding
        cmd_args += ['-ab']

    cmd_return = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = cmd_return.stdout.decode('utf-8')
    logging.debug(stdout)

real_number_pattern = r'[-+]?[0-9]*\.?[0-9]+(e[-+]?[0-9]+)?'
score_re = re.compile(rf'REMARK VINA RESULT:\s*(?P<affinity>{real_number_pattern})')

def dock_ad4_pdbqt(ligand_pdbqt, protein_pdbqt, logfile, output_dir, protein_name, ligand_name) -> None:
    """ Run AD4.

    ../../software/x86_64Linux2/autogrid4 -p UYO78372.gpf -l UYO78372.glg
    pythonsh ../../enzymetk/prepare_dpf4.py -l ligand.pdbqt -r UYO78372.pdbqt -o UYO78372.dpf
    ../../software/x86_64Linux2/autodock4 -p UYO78372.dpf -l docking_results_UYO78372.dlg

    """
    gpf = os.path.join(output_dir, f'{protein_name}_{ligand_name}.gpf')
    glg = os.path.join(output_dir, f'{protein_name}_{ligand_name}.glg')
    dpf = os.path.join(output_dir, f'{protein_name}_{ligand_name}.dpf')
    dlg = os.path.join(output_dir, f'{protein_name}_{ligand_name}.dlg')

    os.chdir(output_dir)
    os.system(f'cp {protein_pdbqt} {os.path.join(output_dir, protein_name + ".pdbqt")}')
    os.system(f'cp {ligand_pdbqt} {os.path.join(output_dir, ligand_name + ".pdbqt")}')
    protein_pdbqt = protein_name + ".pdbqt"
    ligand_pdbqt = ligand_name + ".pdbqt"

    print(output_dir)
    # Step 1 prepare GPF
    cmd_list = ['/disk1/ariane/vscode/enzymetk/software/x86_64Linux2/mgltools_x86_64Linux2_1.5.7/bin/pythonsh', '/disk1/ariane/vscode/enzymetk/enzymetk/prepare_gpf.py', 
                '-l', ligand_pdbqt, 
                '-r', protein_pdbqt, 
                '-o', gpf]
    print(' '.join(cmd_list))
    os.system(' '.join(cmd_list))

    #cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # output = cmd_return.stdout.decode('utf-8')
    # logging.debug(output)

    # # Write output to the logging file
    # with open(logfile, 'w+') as fout:
    #     fout.write(output)
        
    # --------- Step 2 prepare GLG
    os.system(' '.join(['/disk1/ariane/vscode/enzymetk/software/x86_64Linux2/autogrid4', '-p', gpf, '-l', glg]))

    # cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # output = cmd_return.stdout.decode('utf-8')
    # logging.debug(output)

    # # Write output to the logging file
    # with open(logfile, 'w+') as fout:
    #     fout.write(output)
        
    # # If failure, raise DockingError
    # if cmd_return.returncode != 0:
    #     print(f'Docking with ad4 failed GLG: {output}')

     #     pythonsh ../../enzymetk/prepare_dpf4.py -l ligand.pdbqt -r UYO78372.pdbqt -o UYO78372.dpf
   # ../../software/x86_64Linux2/autodock4 -p UYO78372.dpf -l docking_results_UYO78372.dlg

    # --------- Step 3 prepare DPF
    cmd_list = ['/disk1/ariane/vscode/enzymetk/software/x86_64Linux2/mgltools_x86_64Linux2_1.5.7/bin/pythonsh', 
                '/disk1/ariane/vscode/enzymetk/enzymetk/prepare_dpf4.py',
                '-l', ligand_pdbqt, 
                '-r', protein_pdbqt,
                '-o', dpf]
    os.system(' '.join(cmd_list))

    # cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # output = cmd_return.stdout.decode('utf-8')
    # logging.debug(output)

    # # Write output to the logging file
    # with open(logfile, 'w+') as fout:
    #     fout.write(output)
        
    # --------- FINALLY RUN AD4
    cmd_list = ['/disk1/ariane/vscode/enzymetk/software/x86_64Linux2/autodock4',
                '-p', dpf,
                '-l', dlg]

    os.system(' '.join(cmd_list))

    # cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # output = cmd_return.stdout.decode('utf-8')
    # logging.debug(output)

    # # Write output to the logging file
    # with open(logfile, 'w+') as fout:
    #     fout.write(output)

    # They can get the results from here.
    return dlg
        

def dock_autodock_pdbqt(conf_path, log_path, out_path, seed, num_cpus: Optional[int] = None) -> None:
    """
    Run AutoDock Vina.

    :param pdbqt_path: path to PDBQT file
    :param log_path: path to log file
    :param out_path: path to output file
    :param seed: random seed
    :param num_cpus: number of CPU cores available to AutoDock Vina
    """
    cmd_list = [
        'vina', # Needs to be installed as vina.
        '--config', conf_path,
        '--out', out_path,
        '--seed', str(seed)
    ]
    # ToDo: add in scoring function for ad4
    print(f'vina --config {conf_path} --out {out_path}')
    if num_cpus is not None:
        cmd_list += ['--cpu', str(num_cpus)]

    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')
    logging.debug(output)

    # Write output to the logging file
    with open(log_path, 'w+') as fout:
        fout.write(output)
        
    # If failure, raise DockingError
    if cmd_return.returncode != 0:
        print(f'Docking with Vina failed: {output}')

