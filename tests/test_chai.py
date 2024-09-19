import os
import sys
from datetime import datetime

from docko import *

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

# set all cudo devices to cuda:1

if __name__ == "__main__":

    log_folder = checkNgen_folder("logs/chai")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    run_chai_df(
        "/disk2/fli/docko/tests/fzl/output-twocofactor",
        "/disk2/fli/docko/tests/fzl/input/trpb_sub.csv",
        entry_column="Entry",
        seq_column="Sequence",
        ligand_column="Substrate",
        cofactor_column="Cofactor",
    )

    # run_chai_df(
    #     "/disk2/fli/docko/tests/fzl/output-nocofactor",
    #     "/disk2/fli/docko/tests/data-fzl/trpb_sub.csv",
    #     entry_column="Entry",
    #     seq_column="Sequence",
    #     ligand_column="Substrate",
    #     cofactor_column="",
    # )

    f.close()