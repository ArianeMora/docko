# Docko
Docking for ligands made simple.

## Install
Make sure you have vina installed. 

```
conda  create --name docko python=3.10
conda activate docko
conda install -c conda-forge pdbfixer
```

### Everything else is just pip
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