# Docko
Docking for ligands a discusting mix of code combining the latest from Chai with old scool bro vina.

Made this for myself but others wanted to use. Love only pls. Take it as it is <3 

## Install
Make sure you have vina installed: https://autodock-vina.readthedocs.io/en/latest/installation.html

```
conda  create --name docko python=3.10
conda activate docko
conda install -c conda-forge pdbfixer
conda install -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme
pip install vina
conda config --env --add channels conda-forge
pip install git+https://github.com/chaidiscovery/chai-lab.git
pip install chai_lab==0.0.1
```

g else is just pip
```
pip install -r requirements.txt
```

## PDB or structure
You need to select your structure from PDB or in liu of that, use alphafold3 server (https://alphafoldserver.com/).

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

Then you can run the program as per usual :D 

### If you wish to use DiffDock (which I would not recomend)
You'll need Trill and to run in the Trill environment:
https://trill.readthedocs.io/en/latest/home.html

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

### Everythin