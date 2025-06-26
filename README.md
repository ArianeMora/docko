# Docko
Docking for ligands a discusting mix of code combining the latest from Chai with old scool bro vina.

Made this for myself but others wanted to use. Love only pls. Take it as it is <3 but if you notice bugs, please submit an 
issue, be a g.

## Install
Make sure you have vina installed: https://autodock-vina.readthedocs.io/en/latest/installation.html

I have not found it to work with pip needs the executable.

Works on mac and liunx, you need big power tho for Chai so would rec linux.

Note you need `openmm==8.0` (if you get some random modeller error this should fix it, it arrises from an incompatibility with openmm and pdbfixer gah)

```
conda  create --name docko python=3.10 -y
conda activate docko
conda install -c conda-forge pdbfixer -y
pip install openmm==8.0
conda config --env --add channels conda-forge
pip install git+https://github.com/chaidiscovery/chai-lab.git
```

### install docko now
```
conda activate docko
pip install docko
```
### Maybe if you're feeling wild install boltz
```
pip install boltz
``` 

### Lucky last since vina is a b
You need to make a second environment just to prepare the ligand, I came across this issue when making all my stuff.
```
conda create --name vina python=3.9.7 -y
conda activate vina
conda install -c conda-forge numpy openbabel pdbfixer scipy rdkit -y

pip install meeko
```


## Quick start

#### Use case 1: you have a sequence and you want to bind it
Here you're best bet is using Chai, this will automatically handle everything for you:

Example:
```
from docko import *

base_dir = 'some_folder' # A folder on your computer

run_chai('A0A0E3LLD2_METBA', # name
         'MSIEKIPGYTYGKTESMSPLNLEDLKLLKDSVMFTEEDEKYLKKAGEVLEDQVEEILDTWYGFVGSHPHLLYYFTSPDGTPNEEYLAAVRKRFSKWILDTCNRNYDQAWLDYQYEIGLRHHRTKKNRTDNVESVPNINYRYLVAFIYPITATIKPFLARKGHTSEEVEKMHQAWFKATVLQVALWSYPYVKQGDF', # sequence
         'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC', # ligand as smiles
         base_dir
        )
```
The outputs will now be in `base_dir`.

Say you have a csv of these and you want to make bound structures for all of them:
```
run_chai_df(output_dir, filename, entry_column='Entry', seq_column='Sequence', ligand_column='Substrate')
```
This runs Chai on a csv that contains your sequneces, ligands (Substartes) and the entry name (Entry) 
and makes a new folder using the entry name (this would mean you ideally don't want dumb characters in there.) And puts 
all these new folders in `output_dir`.

#### Running on Boltz
About the same as chai - several options if you don't want the affinity prediction run the below:

```
import sys
from docko.boltz import *

base_dir = '/XXXXXXX/XXXXXXX/vscode/docko/'
run_boltz('boltzproteinexample', 
         'MSIEKIPGYTYGKTESMSPLNLEDLKLLKDSVMFTEEDEKYLKKAGEVLEDQVEEILDTWYGFVGSHPHLLYYFTSPDGTPNEEYLAAVRKRFSKWILDTCNRNYDQAWLDYQYEIGLRHHRTKKNRTDNVESVPNINYRYLVAFIYPITATIKPFLARKGHTSEEVEKMHQAWFKATVLQVALWSYPYVKQGDF', 
         'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC',
         base_dir
        )
# now your results are in `/XXXXXXX/XXXXXXX/vscode/docko/boltzproteinexample`
```

This will save the results including the confidence model and the `plldt.npz` file.

For affinity prediction run:
```
import sys
from docko.boltz import *

base_dir = '/XXXXXXX/XXXXXXX/vscode/docko/'
run_boltz_affinity('boltzproteinexample', 
         'MSIEKIPGYTYGKTESMSPLNLEDLKLLKDSVMFTEEDEKYLKKAGEVLEDQVEEILDTWYGFVGSHPHLLYYFTSPDGTPNEEYLAAVRKRFSKWILDTCNRNYDQAWLDYQYEIGLRHHRTKKNRTDNVESVPNINYRYLVAFIYPITATIKPFLARKGHTSEEVEKMHQAWFKATVLQVALWSYPYVKQGDF', 
         'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC',
         base_dir, 
         None # or a cofactor 
        )
# now your results are in `/XXXXXXX/XXXXXXX/vscode/docko/boltzproteinexample`
```

#### Use case 2: you have a uniprot ID and you want to get the structure and bind a ligand with vina

Here you got told "oh wow physics informed models are the best, I don't trust ML!" this will typically arise
from someone over the age of 40. Here to humour them you can also run `vina`, you'll need to have it installed.

The smiles is your ligand as smiles, `base_dir` is where you want your data to be output. Note given we are passing 
the `protein_name='A0A0H2V871'` which is a uniprot ID it will automatically get the structure for us. If we weren't 
we would need to pre-download the PDB structure, or fold it using an online server such as AF3 or Chia (you could run 
Chia and then remove the ligand as well - my fave option).

```
from docko import *

smiles = 'CCCCC(CC)COC(=O)C1=CC=CC=C1C(=O)OCC(CC)CCCC'

base_dir = 'some_folder' # A folder on your computer

dock(sequence='', 
    protein_name='A0A0H2V871', # Or the name/path of the file on your computer as pdb or cif.
    smiles=smiles, 
    ligand_name='DEHP', # name of your chemical no funny characters, the ligand will be made in a folder named this
    residues=[113, 114], # Resiudes of your active site, we position the ligand within here (I find the centroid of these guys)
    protein_dir=f'{base_dir}/', # Folder to save the proteins to
    ligand_dir=f'{base_dir}/', # Folder to save the input ligand to
    output_dir=f'{base_dir}/', # output folder with the docked ligand and config file
    pH=7.4, # pH to run docking at
    method='vina', # method can be vina, ad4, or diffdock
    size_x=5.0, # How far in x is alowed think of this as a cloud around your residues or residue centroid
    size_y=5.0, 
    size_z=5.0,
    num_modes=9, # Dunno check vina docks using the defaut
    exhaustivenes=32 ) # higher is better but slower, this is a default

# Just checks the output was logged --> this has your "energy data" about how good the docking was
os.path.isfile(f'{base_dir}A0A0H2V871-DEHP_log.txt')
```
e.g. if you wanted to run it on some file, just change it to your path to your downloaded PDB file e.g.:

```
protein_name=f'{base_dir}data/test_existing.pdb',
```
This will then make the name of your directory `test_existing` and then save the resulst in there. I guess just 
again don't have funny characters in your filename.

As above, you can also run with the option `ad4` it makes all these other random files, and again was something
that someone asked me to do, was seriously painful and I don't wish it on anyone else so have made it available. 
Basically uses some rando forcefield that makes in some cases vina dock better. Who knows. LMK if you have an opinion.


#### Use case 3: you want to use diffdock and use up all the space on your computer

```
dock(sequence='', 
    protein_name='A0A0H2V871', # Or the name/path of the file on your computer as pdb or cif.
    smiles=smiles, 
    ligand_name='DEHP', # name of your chemical no funny characters, the ligand will be made in a folder named this
    residues=[113, 114], # Resiudes of your active site, we position the ligand within here (I find the centroid of these guys)
    protein_dir=f'{base_dir}/', # Folder to save the proteins to
    ligand_dir=f'{base_dir}/', # Folder to save the input ligand to
    output_dir=f'{base_dir}/', # output folder with the docked ligand and config file
    method='diffdock', # As above just change to diffdock
    )

```
Basically exactly as above, except you need to specify the method is `diffdock`. 
Note you need to have TRILL installed for this to work:

```
micromamba create -n TRILL python=3.10 ; micromamba activate TRILL
micromamba install -c pytorch -c nvidia pytorch=2.1.2 pytorch-cuda=12.1 torchdata
micromamba install -c conda-forge openbabel pdbfixer swig openmm smina fpocket vina openff-toolkit openmmforcefields setuptools=69.5.1
micromamba install -c bioconda foldseek pyrsistent
micromamba install -c "dglteam/label/cu121" dgl
micromamba install -c pyg pyg pytorch-cluster pytorch-sparse pytorch-scatter
pip install git+https://github.com/martinez-zacharya/lightdock.git@03a8bc4888c0ff8c98b7f0df4b3c671e3dbf3b1f git+https://github.com/martinez-zacharya/ECPICK.git setuptools==69.5.1
pip install trill-proteins
```

--------------------------------------------------------------------------------------------------------------
## Other info

## PDB or structure
You need to select your structure from PDB or in liu of that, use alphafold3 server (https://alphafoldserver.com/).

Alternatively, if your IDs are PDB IDs or Uniprot IDs you can just pass those and it will get teh structures for you.

If you use the alphafoldserver you'll get `cif` files and this works with that too!

## Working with heme based files
Unfortunatley since alphafold is some new stuff and we're working with autodock vina we will need to change the files a bit. 

First, if we use the AF3 docked heme, this will be automatically "cleaned" before making the pdbqt file. So we use the pipeline on the AF3 structure.

So we need to read-add it back in after also converting it manually. 

To convert it manually, we copy and paste (i know lol) the heme from the original pdb file (if you don't have this, go into a program like chimeraX and convert the .cif file to a .pdb file).

Then you open up the alpha fold pdb in a text editor and copy out the heme atoms, ommitting the last one, the `Fe`, as vina doesn't like this one.

Then, we convert this manually by using obabel: ` obabel heme.pdb -o pdbqt > heme.pdbqt` 

Once this has been converted, we realise that vina doesn't like many of the tags. So we need to then change this so that we remove all of these.
These include (but probably not limited to)

```
ENDBRANCH
ROOT
BRANCH
ENDROOT
```

Then you can run the program as per usual :D I do this automatically within the scripts but thought I would mention it inacse
things fail for you (whoever you are.)

### References

(1) Martinez, Z. A.; Murray, R. M.; Thomson, M. W. TRILL: Orchestrating Modular Deep-Learning Workflows for Democratized, Scalable Protein Analysis and Engineering. bioRxiv October 27, 2023, p 2023.10.24.563881. https://doi.org/10.1101/2023.10.24.563881.  
(2) Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. J. Chem. Inf. Model. 2021, 61 (8), 3891–3898. https://doi.org/10.1021/acs.jcim.1c00203.  
(3) Chai Discovery. https://www.chaidiscovery.com/blog/introducing-chai-1 (accessed 2024-09-15).
(4) Boltz: 
```
@article{wohlwend2024boltz1,
  author = {Wohlwend, Jeremy and Corso, Gabriele and Passaro, Saro and Reveiz, Mateo and Leidal, Ken and Swiderski, Wojtek and Portnoi, Tally and Chinn, Itamar and Silterra, Jacob and Jaakkola, Tommi and Barzilay, Regina},
  title = {Boltz-1: Democratizing Biomolecular Interaction Modeling},
  year = {2024},
  doi = {10.1101/2024.11.19.624167},
  journal = {bioRxiv}
}
```
```
@article{mirdita2022colabfold,
  title={ColabFold: making protein folding accessible to all},
  author={Mirdita, Milot and Sch{\"u}tze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin},
  journal={Nature methods},
  year={2022},
}
```
### THANKX
Lastly if you liked this, give it a star ****