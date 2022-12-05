#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: Mickael MELKOWSKI
"""

import os
import argparse


__version__ = "1.0.0"


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
    required_arguments.add_argument('-i', '--input',
                                    type=str,
                                    required=True,
                                    help="Input variants filename.")

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


def main(input):
    return type(input)


def entrypoint():
    # Get CLI Args values:
    arguments_parser_obj = arguments_parser()
    args = arguments_parser_obj.parse_args()

    # do something
    main(args.input)


if __name__ == "__main__":
    entrypoint()
