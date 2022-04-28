#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: Mickael MELKOWSKI
"""

import os
import argparse
import logging


__version__ = "1.0.0"


def arguments_parser():
    '''
    Arguments parsing.

    :return: Object
    '''
    parser = argparse.ArgumentParser(
        description='Prepare the working directory for targeted analysis.',
        usage='{} --study_name STUDY_NAME --archive PATH/ARCHIVE.tar.gz --md5file_name FILE_NAME.md5 --metadata_file_name FILE_NAME.yaml --schema PATH/SCHEMA.json --target_directory TARGET_PATH'.format(
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


def parametersLogger(args):
    '''
    Logger setting.
    '''
    logging_std_format = '%(asctime)s [%(name)s] [%(levelname)s] %(message)s'
    logging_debug_format = '%(asctime)s [%(levelname)s] [%(threadName)s - %(name)s] %(message)s'
    log_level = args.log_level.upper()
    if log_level == 'DEBUG':
        logging_std_format = logging_debug_format
    logging_datefmt = '%Y-%m-%d %H:%M:%S'
    if args.log_file is not None:
        logging.basicConfig(format=logging_std_format,
                            datefmt=logging_datefmt,
                            filename=args.log_file,
                            filemode='w',
                            level=log_level)
    else:
        logging.basicConfig(format=logging_std_format,
                            datefmt=logging_datefmt,
                            level=log_level)
    logger = logging.getLogger(os.path.basename(__file__))


def main(arguments_parser, args):
    '''
    Checkes integrity of study archive.
    Parses the metadata file.
    Prepares the analyses directory.
    '''
    logger = logging.getLogger(main.__name__)
    logger.info('do something with args.input: {}'.format(args.input))
    try:
        logger.info('Checking integrity...')
        print(args.input)
    except IndexError as ie:
        logger.error(ie)
        arguments_parser.exit(status=1)
    finally:
        logger.info('main function done')
    arguments_parser.exit(status=0)



def entrypoint():
    # Get CLI Args values:
    arguments_parser_obj = arguments_parser()
    args = arguments_parser_obj.parse_args()

    # Set logger args
    parametersLogger(args)

    # do something
    # main(arguments_parser, args)


if __name__ == "__main__":
    entrypoint()
