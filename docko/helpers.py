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
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import PDBIO, Select, MMCIFIO
from Bio import PDB
from rdkit.Chem import AllChem
import requests
from pathlib import Path


import logging
import subprocess
import os
from rdkit import Chem
from rdkit.Chem import AllChem as Chem


def canonicalize_smiles(smiles_string):
    molecule = Chem.MolFromSmiles(smiles_string)
    if molecule:
        canonical_smiles = Chem.MolToSmiles(molecule, canonical=True)
        return canonical_smiles


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
    print(fixed_pdbFile)
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
    ligand_mol_file = os.path.join(this_ligand_dir, name + '.pdb')
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

    # Chem.MolToPDBFile(mol, filename=ligand_mol_file)
    # # mk_prepare_ligand.py
    # # yapf: disable
    # cmd_list = [
    #     'obabel',
    #     '-imol', ligand_mol_file,
    #     '-opdbqt',
    #     '-O', ligand_pdbqt_file,
    #     '--partialcharge', 'gasteiger'
    # ]
    # # yapf: enable
    # cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # output = cmd_return.stdout.decode('utf-8')

    # Also save the mol as a sdf 
    # Write to SDF file
    # Assuming 'mol' is the molecule you're working with
    if mol is None:
        print("The molecule is None. Check the input.")
    elif mol.GetNumAtoms() == 0:
        print("The molecule has no atoms. It might be invalid.")
    else:
        print(f"The molecule has {mol.GetNumAtoms()} atoms.")
        writer = Chem.SDWriter(ligand_sdf_file)
        writer.write(mol)
        writer.close()


    print('------- WARN YOU NEED VINA ENV if ya dont this next step will probs fail -----------')
    print(f'conda run -n vina mk_prepare_ligand.py -i {ligand_sdf_file} -o {ligand_pdbqt_file}')
    os.system(f'conda run -n vina mk_prepare_ligand.py -i {ligand_sdf_file} -o {ligand_pdbqt_file}')
    # Clean the ligand one more time removing any extras
    lines = []
    with open(ligand_pdbqt_file, 'r+') as fin:
        for line in fin:
            #if line.split(' ')[0] not in ['ENDBRANCH', 'BRANCH', 'ROOT', 'ENDROOT', 'MODEL', 'ENDMDL', 'TORSDOF', 'REMARK', 'CONECT'] and 'ROOT' not in line:
            lines.append(line.replace('UNL', 'LIG').replace('HETATM', 'ATOM')) # Replace any unknowns with GLC bit of a hack but yolo
                # ToDo: fix this shit.

    with open(ligand_pdbqt_file, 'w+') as fout:
        for line in lines:
            fout.write(line)

    return ligand_pdbqt_file, ligand_sdf_file


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


def pdb_to_pdbqt_protein(input_filename, output_filename=None, ph=7.4):
    """
    Convert a pdb file to a pdbqt file.
    """
    # Need to first remove stuff that is sometimes added by
    lines = []
    with open(input_filename, 'r+') as fin:
        for line in fin:
            if line.split(' ')[0] not in ['ENDBRANCH', 'BRANCH', 'ROOT', 'ENDROOT'] and 'Fe' not in line: # Add in the removal of the Iron bit
                lines.append(line)
    with open(input_filename, 'w+') as fout:
        for line in lines:
            fout.write(line)

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


def format_pdb(sequence, name, protein_dir, pH):
    """ 
    Check if we have a structure, if we don't try get one from Alpha fold. If we can't get one from open Fold.
    Clean this guy for docking and other random shit.
    """
    # First check if they passed a filename or the sequence name
    if os.path.isfile(name):
        # basically this means we have a pdb file
        # Step 4: Clean em and return em!
        protein_input_file = os.path.join(name)
        name = name.split('/')[-1].replace('.cif', '').replace('.pdb', '')

        this_protein_dir = os.path.join(protein_dir, name)

        if not os.path.exists(this_protein_dir):
            os.system(f'mkdir {this_protein_dir}')

        protein_pdb_file = os.path.join(this_protein_dir, name + '.pdb')
        protein_pdbqt_file = os.path.join(this_protein_dir, name + '.pdbqt')

        # Now run
        clean_one_pdb(protein_input_file, protein_pdb_file)

        # Step 5: Convert to pdbqt for vina
        pdb_to_pdbqt_protein(protein_pdb_file, protein_pdbqt_file, ph=pH)
        return name, protein_pdb_file, protein_pdbqt_file
    else:
        this_protein_dir = os.path.join(protein_dir, name)
        protein_pdbqt_file = os.path.join(this_protein_dir, name + '.pdbqt')
        protein_pdb_file = os.path.join(this_protein_dir, name + '.pdb')

        if os.path.isfile(protein_pdbqt_file):
            return name, protein_pdb_file, protein_pdbqt_file

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

        return name, protein_pdb_file, protein_pdbqt_file


def checkNgen_folder(folder_path: str) -> str:

    """
    Check if the folder and its subfolder exists
    create a new directory if not

    Args:
    - folder_path: str, the folder path
    """

    # if input path is file
    if bool(os.path.splitext(folder_path)[1]):
        folder_path = os.path.dirname(folder_path)

    split_list = os.path.normpath(folder_path).split("/")
    for p, _ in enumerate(split_list):
        subfolder_path = "/".join(split_list[: p + 1])
        if not os.path.exists(subfolder_path):
            print(f"Making {subfolder_path} ...")
            os.mkdir(subfolder_path)
    return folder_path
