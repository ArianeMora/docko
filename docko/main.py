from docko.chai import *
import argparse


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