#!/usr/bin/env python

import logging
import sys
import pandas as pd
from Bio import SeqIO

logfile="log.txt" # can be --> snakemake.log[0] if in rule smk file
logging.basicConfig(filename=logfile,
                    level=logging.INFO,
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )


def arguments_parser():
    '''
    Arguments parsing.

    :return: Object
    '''
    parser = argparse.ArgumentParser(
        description='Do something.',
        usage='{} --input INPUT --output OUTPUT'.format(
            os.path.basename(__file__)),
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version=__version__)

    # required group
    required_arguments = parser.add_argument_group('Required Arguments')
    required_arguments.add_argument('-f', '--fasta',
                                    type=str,
                                    required=True,
                                    help="Input fasta.")
    required_arguments.add_argument('-b', '--blast',
                                    type=str,
                                    required=True,
                                    help="Input blast.")

    # output group
    output_arguments = parser.add_argument_group('Outputs')
    output_arguments.add_argument('-o', '--output',
                                  type=str,
                                  required=True,
                                  help="output filename.")

    # logger group
    logger_arguments = parser.add_argument_group('logger arguments')
    logger_arguments.add_argument('--log_level',
                                  type=str,
                                  default='INFO',
                                  help='log level',
                                  choices=['ERROR', 'error', 'WARNING',
                                           'warning', 'INFO', 'info', 'DEBUG', 'debug'],
                                  required=False)
    logger_arguments.add_argument('--log_file',
                                  type=str,
                                  help='log file (use the stderr by default). Disable the log colors.',
                                  required=False)

    return parser


# for blast outfmt 6
def blast_import(blast_out):
    cols=[
        'qaccver', 'saccver', 'stitle', 'pident', 'length',
        'score', 'mismatch', 'gapopen', 'qstart', 'qend',
        'sstrand', 'sstart', 'send', 'qlen', 'slen'
    ]
    df = pd.read_csv(blast_out,
            sep="\t",
            header=None)
    return df.set_axis(cols, axis=1)


def fasta_to_dict(fname):
    for seq_record in SeqIO.parse(fname, "fasta"):
        header = seq_record.id
        sequence = str(seq_record.seq)
    return {'header': header, 'sequence': sequence}


def write_to_fasta(out_fasta, header, seq):
    with open(out_fasta, "w") as f:
        f.write(f">{header}\n{seq}\n")


def write_to_multifasta(out_fasta, list_entry):
    with open(out_fasta, "w") as f:
        for header, seq in list_entry:
            f.write(f">{header}\n{seq}\n")


def main(fasta, blast, output):
    # load input
    fasta_dict = fasta_to_dict(fasta)
    blastn_results = blast_import(blast)

    # Output
    logging.info(f'Done.')


def entrypoint():
    # Get CLI Args values:
    arguments_parser_obj = arguments_parser()
    args = arguments_parser_obj.parse_args()

    # do something
    main(args.fasta, args.blast, args.output)

if __name__ == '__main__':
    entrypoint()