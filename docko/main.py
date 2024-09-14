from pathlib import Path
from Bio import SeqIO
import os
import torch
import pandas as pd
from rdkit import Chem
from docko.docko import *

import argparse
from chai_lab.chai1 import run_inference

def parse_args():
    parser = argparse.ArgumentParser(description="Run Chai on a sequence")
    parser.add_argument('-out', '--out', required=True, help='Path to the output directory')
    parser.add_argument('-df', '--df', type=str, required=True, help='Fasta of the file of interest')
    return parser.parse_args()

def main():
    args = parse_args()
    run_chai(args.out, args.df)

if __name__ == "__main__":
    main()