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

import os
import re
from docko.helpers import *
import pandas as pd
from docko.vina import *
from tqdm import tqdm


def run_vina_docking_thread(args):
    """ Docking threading """
    # For each ID we create a new folder so that we don't have issues with the names of files
    df, output_dir, output_file, protein_dir, ligand_dir, size_x, size_y, size_z, pH, method = args[0], args[1], args[
        2], args[3], args[4], args[5], args[6], args[7], args[8], args[9]
    # we're going to read in the file and if it already has an entry for this uniprot we won't do it
    done_df = pd.read_csv(os.path.join(output_dir, output_file))
    done_ids = list(done_df['uniprot_id'].values)
    for seq_id, seq, active_site, smiles, ligand_name in tqdm(
            df[['Entry', 'Sequence', 'Residue', 'Smiles', 'LigandName']].values):
        if seq_id not in done_ids:
            # 1) make a directory
            ligand_name = ligand_name.replace(' ', '_')
            try:
                affinities = dock(seq, seq_id, smiles, ligand_name, active_site.split('|'), protein_dir, ligand_dir,
                                    output_dir, pH,
                                    size_x=size_x, size_y=size_y, size_z=size_z, method=method)
                if method == 'vina':
                    write_vina_affinities(seq_id, ligand_name, affinities, os.path.join(output_dir, 'scores.txt'))
            except Exception as e:
                print('-------------------------------------------------------------------------')
                print(seq_id, ligand_name)



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

def make_config_for_vina(pdb_file, ligand_file, residues, output_file: str, size_x=5.0, size_y=5.0, size_z=5.0,
                         num_modes=9, exhaustivenes=32):
    """ Make the config file for Vina based on coordinates for a specific residue (or multiple). """
    coords = []
    for position in residues:
        chain_id, coordinates = get_coordinates_without_chain_id(pdb_file, int(position))
        if coordinates is not None:
            coords.append(coordinates)
    print(coords)
    coords = calculate_centroid(coords)
    # Now we want to save this to a config file that has everything we need
    with open(f'{output_file}', 'w+') as fout:
        # receptor = enzyme.pdbqt
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
              num_cpus: int = None, seed=42):
    """
    Given a molecule, this method will return a docking score against the current target.
    """
    conf_path = os.path.join(output_dir, f'{protein_name}-{ligand_name}_conf.txt')

    # Step 1: make config file
    make_config_for_vina(protein_pdbqt, ligand_pdbqt, residues, conf_path, size_x=size_x,
                         size_y=size_y, size_z=size_z, num_modes=num_modes, exhaustivenes=exhaustivenes)

    # Auxiliary files
    vina_logfile = os.path.join(output_dir, f'{protein_name}-{ligand_name}_log.txt')

    docked_ligand_pdb = os.path.join(output_dir, f'{protein_name}-{ligand_name}.pdb')

    # Step 2: Dock using our nice new baby's
    if method == 'ad4':
        dock_ad4_pdbqt(ligand_pdbqt, protein_pdbqt, vina_logfile, output_dir, protein_name, ligand_name)
    elif method == 'vina':
        dock_autodock_pdbqt(conf_path, vina_logfile, docked_ligand_pdb, seed=seed, num_cpus=num_cpus)

    # Step 5: convert back to a pdb structure since whotf uses pdbqt
    convert_pdbqt_to_pdb(pdbqt_file=docked_ligand_pdb, pdb_file=docked_ligand_pdb, disable_bonding=True)

    # Step 6: Parse scores
    return vina_logfile


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
    package_root = Path(__file__).resolve().parent.parent

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
    cmd_list = [f'{package_root}/docko/deps/x86_64Linux2/mgltools_x86_64Linux2_1.5.7/bin/pythonsh',
                f'{package_root}/docko/deps/prepare_gpf.py',
                '-l', ligand_pdbqt,
                '-r', protein_pdbqt,
                '-o', gpf]
    print(' '.join(cmd_list))
    os.system(' '.join(cmd_list))

    # --------- Step 2 prepare GLG
    os.system(' '.join([f'{package_root}/docko/deps/x86_64Linux2/autogrid4', '-p', gpf, '-l', glg]))

    # --------- Step 3 prepare DPF
    cmd_list = [f'{package_root}/docko/deps/x86_64Linux2/mgltools_x86_64Linux2_1.5.7/bin/pythonsh',
                f'{package_root}/docko/deps/prepare_dpf4.py',
                '-l', ligand_pdbqt,
                '-r', protein_pdbqt,
                '-o', dpf]
    os.system(' '.join(cmd_list))

    # --------- FINALLY RUN AD4
    cmd_list = [f'{package_root}/docko/deps/x86_64Linux2/autodock4',
                '-p', dpf,
                '-l', dlg]

    os.system(' '.join(cmd_list))

    # They can get the results from here.
    return dlg


def dock_autodock_pdbqt(conf_path, log_path, out_path, seed, num_cpus: int = None) -> None:
    """
    Run AutoDock Vina.

    :param conf_path: config
    :param log_path: path to log file
    :param out_path: path to output file
    :param seed: random seed
    :param num_cpus: number of CPU cores available to AutoDock Vina
    """
    cmd_list = [
        'vina',  # Needs to be installed as vina.
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
