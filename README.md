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

If you use the alphafoldserver you'll get 