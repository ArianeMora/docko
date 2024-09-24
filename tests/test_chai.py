import os
import sys
from datetime import datetime

from docko import *
from docko.chai import run_chai_df

os.environ["CUDA_VISIBLE_DEVICES"] = "1"

# set all cudo devices to cuda:1

if __name__ == "__main__":

    log_folder = checkNgen_folder("logs/chai/kms_lin")

    # log outputs
    f = open(os.path.join(log_folder, f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}.out"), 'w')
    sys.stdout = f

    run_chai_df(
        "tests/fzl/output-kms_lin_exportedsmiles",
        "tests/fzl/input/kms_3.csv",
        entry_column="Entry",
        seq_column="Sequence",
        ligand_column="Substrate",
        cofactor_column="Cofactor",
    )

    run_chai_df(
        "tests/fzl/output-kms_lin_exportedsmiles2",
        "tests/fzl/input/kms_3.csv",
        entry_column="Entry",
        seq_column="Sequence",
        ligand_column="Substrate-cofactor",
        cofactor_column="",
    )
    # run_chai_df(
    #     "tests/fzl/output-twocofactor",
    #     "tests/fzl/input/trpb_sub.csv",
    #     entry_column="Entry",
    #     seq_column="Sequence",
    #     ligand_column="Substrate",
    #     cofactor_column="Cofactor",
    # )

    # run_chai_df(
    #     "tests/fzl/output-combcofactor",
    #     "tests/fzl/input/trpb_sub.csv",
    #     entry_column="Entry",
    #     seq_column="Sequence",
    #     ligand_column="Substrate-cofactor",
    #     cofactor_column="",
    # )


    f.close()