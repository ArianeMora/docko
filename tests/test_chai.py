import os
from docko import *

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

# set all cudo devices to cuda:1

run_chai_df(
    "/disk2/fli/docko/tests/output-fzl/",
    "/disk2/fli/docko/tests/data-fzl/trpb_sub.csv",
    entry_column="Entry",
    seq_column="Sequence",
    ligand_column="Substrate",
)

run_chai_df(
    "/disk2/fli/docko/tests/output-fzl-nocofactor/",
    "/disk2/fli/docko/tests/data-fzl/trpb_sub.csv",
    entry_column="Entry",
    seq_column="Sequence",
    ligand_column="Substrate",
    cofactor_column="",
)