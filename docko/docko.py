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

# Function to convert SMILES to SDF
def smiles_to_sdf(smiles: str, output_sdf_path: str):
    # Convert SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Compute 3D coordinates
    AllChem.EmbedMolecule(mol)
    
    # Write to SDF file
    writer = Chem.SDWriter(output_sdf_path)
    writer.write(mol)
    writer.close()
    print(f'SDF file saved to {output_sdf_path}')
    return output_sdf_path

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


def dock_diffdock(output_dir: str, ligand_sdf, pdb_file: str):
    """ Dock using diffdock. """
    os.chdir(output_dir)
    print(f'trill example 1 dock DiffDock {pdb_file} {ligand_sdf}')
    os.system(f'trill example 1 dock DiffDock {pdb_file} {ligand_sdf}')
    # Also now need to move and rename it

def dock(sequence: str, protein_name, smiles: str, ligand_name: str, residues: list, protein_dir: str, ligand_dir: str,
            output_dir: str, pH: float, method: str, size_x=5.0, size_y=5.0, size_z=5.0, num_modes=9, exhaustivenes=32):
    """ Dock a smiles to a pdb file using vina and dockstring to automate some of the things. """

    # Step 1: Check if the protein exists, if not make it pretty and save it to the protein folder
    protein_pdbqt = format_pdb(sequence, protein_name, protein_dir, pH)

    # Step 2: Check if the ligand exists, if not make it pretty and save it to the ligand folder
    ligand_pdbqt, ligand_sdf = format_ligand(smiles, ligand_name, ligand_dir, pH)

    # step 3: dock
    if method == 'vina':
        affinities = dock_vina(ligand_pdbqt, protein_pdbqt, ligand_name, protein_name, output_dir, residues, size_x=size_x, 
                               size_y=size_y, size_z=size_z, num_modes=num_modes, exhaustivenes=exhaustivenes, method=method)
    elif method == 'diffdock':
        dock_diffdock(output_dir, ligand_sdf, protein_pdbqt.replace('.pdbqt', '.pdb'))
    return None

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

def get_coordinates_without_chain_id(pdb_file, seq_position):
    # Create a PDB parser
    parser = PDB.PDBParser(QUIET=True)

    # Parse the PDB file
    structure = parser.get_structure("protein", pdb_file)

    # Iterate over all chains and residues to find the matching sequence position
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == seq_position:
                    print(residue)
                    # Extract the coordinates of the alpha carbon (CA) atom
                    ca_atom = residue['CA']
                    return str(residue), ca_atom.get_coord()

    return None, None


def calculate_centroid(coords):
    x_coords = [point[0] for point in coords]
    y_coords = [point[1] for point in coords]
    z_coords = [point[2] for point in coords]

    centroid = (
        sum(x_coords) / len(coords),
        sum(y_coords) / len(coords),
        sum(z_coords) / len(coords)
    )

    return centroid

"""
------------------------------------------------------------------------------------
This is from: https://github.com/luwei0917/DynamicBind Thanks <3
------------------------------------------------------------------------------------
"""
def get_alphafold_structure(uniprot_accession, output_file):
    """ get an alphafold structure if it exists! """
    if not os.path.isfile(output_file):

        url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
        response = requests.get(url)

        if response.status_code != 200:
            print(f"Failed to retrieve data: {response.status_code}")
            return

        data = response.json()

        if not data:
            print("No structure data found for the given UniProt accession.")
            return

        pdb_url = data[0]['pdbUrl']
        pdb_response = requests.get(pdb_url)

        if pdb_response.status_code != 200:
            print(f"Failed to retrieve PDB file: {pdb_response.status_code}")
            return

        with open(output_file, 'w') as pdb_file:
            pdb_file.write(pdb_response.text)

        print(f"Structure saved to {output_file}")


def remove_hydrogen_pdb(pdbFile, toFile):
    parser = MMCIFParser(QUIET=True) if pdbFile[-4:] == ".cif" else PDBParser(QUIET=True)
    s = parser.get_structure("x", pdbFile)

    class NoHydrogen(Select):
        def accept_atom(self, atom):
            if atom.element == 'H' or atom.element == 'D':
                return False
            return True

    io = MMCIFIO() if toFile[-4:] == ".cif" else PDBIO()
    io.set_structure(s)
    io.save(toFile, select=NoHydrogen())

def save_clean_protein(s, toFile, keep_chain='A', keep_all_protein_chains=True):
    """
    First remove everything before adding everything back in.
    """
    class MySelect(Select):
        def accept_residue(self, residue, keep_chain=keep_chain):
            pdb, _, chain, (hetero, resid, insertion) = residue.full_id
            if keep_all_protein_chains or (chain == keep_chain):
                # if hetero == ' ' or hetero == 'H_MSE':
                if hetero == 'H_MSE':
                    residue.id = (' ', resid, insertion)
                    return True
                elif hetero == ' ':
                    return True
                else:
                    return False
            else:
                return False

        def accept_atom(self, atom):
            # remove altloc atoms.
            return not atom.is_disordered() or atom.get_altloc() == "A"

    if toFile[-4:] == ".cif":
        io = MMCIFIO()
    elif toFile[-4:] == ".pdb":
        io = PDBIO()
    else:
        print("toFile should end with .cif or .pdb")
    io.set_structure(s)
    io.save(toFile, MySelect())


def clean_one_pdb(proteinFile, toFile, keep_chain='keep_all'):
    """ Clean and then remove stuff from a pdb file and then read-add in missing things. """
    if proteinFile[-4:] == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    s = parser.get_structure(proteinFile, proteinFile)
    if keep_chain == 'keep_all':
        save_clean_protein(s, toFile, keep_all_protein_chains=True)
    else:
        save_clean_protein(s, toFile, keep_chain=keep_chain, keep_all_protein_chains=False)

    pdbFile = toFile
    fixed_pdbFile = toFile
    # Remove then readd hydrogens
    remove_hydrogen_pdb(pdbFile, fixed_pdbFile)
    fixer = PDBFixer(filename=fixed_pdbFile)
    fixer.removeHeterogens()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    if pdbFile[-3:] == 'pdb':
        PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdbFile, 'w'),
                          keepIds=True)
    elif pdbFile[-3:] == 'cif':
        PDBxFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdbFile, 'w'),
                           keepIds=True)
    else:
        raise 'protein is not pdb or cif'

"""
------------------------------------------------------------------------------------
This is adapted from docstring: https://github.com/dockstring/dockstring/blob/main/dockstring/utils.py
I had to edit it as I was getting issues. 
------------------------------------------------------------------------------------
"""

def read_mol_from_pdb(pdb_file) -> Chem.Mol:
    """
    Read RDKit Mol object from PDB file

    :param pdb_file: path to PDB input file
    :return: RDKit Mol object
    """
    mol = Chem.MolFromPDBFile(str(pdb_file))
    # Try running open babel and minimising the energy
    if not mol or mol.GetNumConformers() == 0:
        os.system(f'obabel {pdb_file} -O {pdb_file} --minimize --gen3d --addh')
        mol = Chem.MolFromPDBFile(str(pdb_file))
    return mol

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


def protonate_smiles(smiles: str, pH: float) -> str:
    """
    Protonate SMILES string with OpenBabel at given pH

    :param smiles: SMILES string of molecule to be protonated
    :param pH: pH at which the molecule should be protonated
    :return: SMILES string of protonated structure
    """

    # cmd list format raises errors, therefore one string
    cmd = f'obabel -:"{smiles}" -ismi -ocan -p{pH}'
    cmd_return = subprocess.run(cmd, capture_output=True, shell=True)
    output = cmd_return.stdout.decode('utf-8')
    logging.debug(output)

    if cmd_return.returncode != 0:
        print('WARNING! COULD NOT PROTONATE')
        return None

    return output.strip()

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


def format_ligand(smiles: str, name: str, ligand_dir: str, pH: float):
    """
    Check if the ligand exists already, if so, just read it in, otherwise format the smiles.
    """
    this_ligand_dir = os.path.join(ligand_dir, name)

    # CHeck if the ligand folder exists otherwise make one
    if not os.path.exists(this_ligand_dir):
        os.mkdir(this_ligand_dir)

    ligand_pdbqt_file = os.path.join(this_ligand_dir, name + '.pdbqt')
    ligand_sdf_file = os.path.join(this_ligand_dir, name + '.sdf')
    ligand_mol_file = os.path.join(this_ligand_dir, name + '.mol')
    if os.path.isfile(ligand_pdbqt_file):
        return ligand_pdbqt_file, ligand_sdf_file

    # Make sure user input is standardized
    smiles = Chem.CanonSmiles(smiles, useChiral=True)

    # Protonate ligand
    protonated_smiles = protonate_smiles(smiles, pH=pH)
    protonated_mol = Chem.MolFromSmiles(protonated_smiles, sanitize=True)

    uncharger = Uncharger()
    # Remove charges
    mol = uncharger.uncharge(protonated_mol)

    # Embed ligand
    # Add hydrogen atoms in order to get a sensible 3D structure
    mol = Chem.AddHs(mol)
    Chem.EmbedMolecule(mol)
    sterochemisty_number = refine_mol_with_ff(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True)

    # Dock
    if mol.GetNumConformers() < 1:
        print('For conversion to MDL MOL format a conformer is required')

    Chem.MolToMolFile(mol, filename=ligand_mol_file)

    # yapf: disable
    cmd_list = [
        'obabel',
        '-imol', ligand_mol_file,
        '-opdbqt',
        '-O', ligand_pdbqt_file,
        '--partialcharge', 'gasteiger'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')
    logging.debug(output)
    # Also save the mol as a sdf 
    # Write to SDF file
    writer = Chem.SDWriter(ligand_sdf_file)
    writer.write(mol)

    return ligand_pdbqt_file, ligand_sdf_file
    
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

import requests

def get_pdb_structure(pdb_id, file_path):
    """
    Retrieve a PDB structure from RCSB PDB and save it to a file.

    :param pdb_id: The PDB ID of the structure to retrieve.
    :param file_path: The path to save the retrieved PDB file.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(file_path, 'w') as file:
            file.write(response.text)
        print(f"PDB structure {pdb_id} saved to {file_path}")
    else:
        print(f"Failed to retrieve PDB structure {pdb_id}")

#     protein_pdbqt = format_pdb(sequence, protein_name, protein_dir, pH)

def format_pdb(sequence, name, protein_dir, pH):
    """ 
    Check if we have a structure, if we don't try get one from Alpha fold. If we can't get one from open Fold.
    Clean this guy for docking and other random shit.
    """ 
    this_protein_dir = os.path.join(protein_dir, name)
    protein_pdbqt_file = os.path.join(this_protein_dir, name + '.pdbqt')
    protein_pdb_file = os.path.join(this_protein_dir, name + '.pdb')

    if os.path.isfile(protein_pdbqt_file):
        return protein_pdbqt_file
    
    if not os.path.exists(this_protein_dir):
        os.system(f'mkdir {this_protein_dir}')

    # Step 1: try get PDB if PDB is a PDB id
    get_pdb_structure(name, protein_pdb_file)

    # Step 2: Try alpha fold structure (NOTE we expect name to be a uniprot ID this is the only way this will work.)
    if not os.path.isfile(protein_pdb_file):
        get_alphafold_structure(name, protein_pdb_file)

    # Step 3: Denovo make the structure using open fold
    # if not os.path.isfile(protein_pdb_file): 

    # Step 4: Clean em and return em!
    clean_one_pdb(protein_pdb_file, protein_pdb_file)

    # Step 5: Convert to pdbqt for vina
    pdb_to_pdbqt_protein(protein_pdb_file, protein_pdbqt_file, ph=pH)

    return protein_pdbqt_file

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
