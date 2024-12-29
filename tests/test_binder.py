import unittest

import pandas as pd
from sciutil import *
import sys
from docko.docko import *
from docko.chai import *

u = SciUtil()

base_dir = '/disk1/ariane/vscode/docko/tests/'
base_dir = ''


class TestBinder(unittest.TestCase):

    def test_clean_pdb(self):
        """ Test cleaning a pdb structure so that we can have """
        clean_one_pdb(f'{base_dir}data/UYO78372.pdb', f'{base_dir}data/UYO78372_cleaned.pdb')

    def test_convert_pdb(self):
        """ Test converting a pdb to a pdbqt file. """
        pdb_to_pdbqt_protein(f'{base_dir}data/UYO78372.pdb', f'{base_dir}data/UYO78372_out.pdbqt')
        assert os.path.isfile(f'{base_dir}data/UYO78372_out.pdbqt') == True

    def test_get_coordinates(self):
        """ Tests getting the coordinates of a specific residue from a PDB file. """
        chain_id, coordinates = get_coordinates_without_chain_id(f'{base_dir}data/UYO78372.pdb', 133)
        print(chain_id, coordinates)
        assert 'GLY' in chain_id
        assert int(coordinates[0]*100) == -31  # Doesn't like the floats eww

    def test_docking(self):
        # def dock(sequence: str, protein_name, smiles: str, ligand_name: str, residues: list, protein_dir: str, ligand_dir: str,
        # output_dir: str, pH: float, method: str, size_x=5.0, size_y=5.0, size_z=5.0, num_modes=9, exhaustivenes=32):
        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'
        # Ensure the protein is in the folder
        os.system(f'mkdir {base_dir}data/UYO78372/')
        os.system(f'cp {base_dir}data/UYO78372.pdb {base_dir}data/UYO78372/UYO78372.pdb')
        dock(sequence='', protein_name='UYO78372', smiles=smiles, ligand_name='dehp', residues=[113], 
             protein_dir=f'{base_dir}data/', ligand_dir=f'{base_dir}data/', output_dir=f'{base_dir}data/',
             pH=7.4, method='vina')
        # Check output was created
        assert os.path.isfile(f'{base_dir}data/UYO78372-dehp_log.txt') == True

    def test_docking_uniprot(self):
        # def dock(sequence: str, protein_name, smiles: str, ligand_name: str, residues: list, protein_dir: str,
        # ligand_dir: str, output_dir: str, pH: float, method: str, size_x=5.0, size_y=5.0, size_z=5.0, num_modes=9,
        # exhaustivenes=32):

        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'
        dock(sequence='', protein_name='A0A0H2V871', smiles=smiles, ligand_name='DEHP', residues=[113], 
                          protein_dir=f'{base_dir}data/', ligand_dir=f'{base_dir}data/', output_dir=f'{base_dir}data/', pH=7.4, method='vina')
        assert os.path.isfile(f'{base_dir}data/A0A0H2V871-DEHP_log.txt') == True

    def test_get_structure(self):
        """ test getting the structure for alpha fold. """
        get_alphafold_structure('Q8VDG7', f'{base_dir}data/Q8VDG7.pdb')
        assert os.path.isfile(f'{base_dir}data/Q8VDG7.pdb') == True

    def test_run_para(self):
        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'
        # For a list of uniprot identifiers check their docking abilities
        df = pd.read_csv(f'{base_dir}data/active_sites.csv')
        df = df.head(2)
        df['Smiles'] = smiles
        df['LigandName'] = 'DEHP'

        # df, output_dir, protein_dir, ligand_dir, output_file,
        calculate_docking_affinities_across_dataset(df, base_dir, f'{base_dir}/data/tmp/',
                                                    f'{base_dir}/data/tmp/',  'docking_scores.csv',
                                                    num_threads=1)
        
    def test_diffdock(self):
        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'

        # For a list of uniprot identifiers check their docking abilities
        # def calculate_docking_affinities_across_dataset(df, output_dir, output_file, smiles, num_threads=20):
        df = pd.read_csv(f'{base_dir}active_sites.csv')
        df = df.head(2)
        diffdock_df(df, base_dir, smiles, num_threads=1)
        
    def test_dock_af3(self):
        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'
        # , 270, 240
        """ Test cleaning a pdb structure so that we can have """
        base_dir = f'{base_dir}data'
        label = 'fold_qhh21706_model_0.cif'
        #clean_one_pdb(f'{base_dir}{label}.cif', f'{base_dir}{label}.pdb')

        """ Test converting a pdb to a pdbqt file. """
        #pdb_to_pdbqt_protein(f'{base_dir}data/{label}{label}.pdb', f'{base_dir}data/{label}/{label}.pdbqt')

        score = dock(sequence='', protein_name=label, smiles=smiles, ligand_name='dehp', residues=[86, 89], 
                    protein_dir=f'{base_dir}data/', ligand_dir=f'{base_dir}data/', output_dir=f'{base_dir}data/',
                    pH=7.4, method='diffdock')
        print(score)

    def test_ad4(self):
        ## Testing for Autodock 4 with that forcefield, this only works on a linux computer
        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'
        #base_dir = '/disk1/ariane/vscode/docko/tests/'
        dock(sequence='', protein_name='A0A0H2V871', smiles=smiles, ligand_name='DEHP', residues=[113], 
                        protein_dir=f'{base_dir}/', ligand_dir=f'{base_dir}/', output_dir=f'{base_dir}/', pH=7.4, method='ad4')
        os.path.isfile(f'{base_dir}A0A0H2V871-DEHP_log.txt')


    def test_existing_pdb(self):
        ## Testing for Autodock 4 with that forcefield, this only works on a linux computer
        smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'
        base_dir = '/disk1/ariane/vscode/docko/tests/'
        dock(sequence='', protein_name=f'{base_dir}data/test_existing.pdb', smiles=smiles, ligand_name='n', residues=[113], 
                        protein_dir=f'{base_dir}/', ligand_dir=f'{base_dir}/', output_dir=f'{base_dir}/', pH=7.4, method='vina')
        os.path.isfile(f'{base_dir}test_existing-example_log.txt')

                    
    def test_smiles(self):
        smiles = 'CCCCCCCCCCCC(=O)Oc1ccc([N+](=O)[O-])cc'
        
        base_dir = '/disk1/ariane/vscode/docko/tests/'
        dock(sequence='', protein_name=f'{base_dir}data/test_existing.pdb', smiles=smiles, ligand_name='asljaklsd', residues=[113], 
                        protein_dir=f'{base_dir}/', ligand_dir=f'{base_dir}/', output_dir=f'{base_dir}/', pH=7.4, method='vina')
        os.path.isfile(f'{base_dir}test_existing-example_log.txt')


    def test_chia(self):
        run_chai('ihateyou', # name
         'MSIEKIPGYEVLEDQVEEILDTWYVESVPNINYRYLVAFIYPITATIKPFLARKGHTSEEVEKMHQAWFKATVLQVALWSYPYVKQGDF', # sequence
         'CCCCC(CC)COC(=O)C(CC)CCCC', # ligand as smiles
         base_dir,
        )
        

if __name__ == '__main__':
    unittest.main()
