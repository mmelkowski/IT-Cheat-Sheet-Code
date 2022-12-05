#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is a default script template for any python code.
See this as an example of boilerplate code.

This provide an argparse, logging and main function examples.
"""


__version__ = '1.1.0'
__author__ = 'Mister Nobody'
__laboratory__='Evil Inc - S.A.S'
__copyright__ = 'Copyright 2022, Mister Nobody, Big Boy'
__license__ = 'GPLv3'
__version_date__='14/06/2022'
__email__ = 'M.Nobody@email.com'
__status__ = 'Squared Sun'


import argparse
import logging
import os


def arguments_parser():
    '''
    Arguments parsing.

    :return: argparse.ArgumentParser object
    :rtype: class:'argparse.ArgumentParser'
    '''
    parser = argparse.ArgumentParser(description='Script description here.')

    parser.add_argument('-v', '--version', action='version', version=__version__)

    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument('--metadata',
                                    type=str,
                                    help='Metadata file path.',
                                    required=True)
    required_arguments.add_argument('--data_type',
                                    type=str,
                                    choices=["itas", "viaas", "vitas"],
                                    help='Analysis type requested in the metadata.',
                                    required=True)

    optional_arguments = parser.add_argument_group('optional arguments')
    optional_arguments.add_argument('--option',
                                    type=str,
                                    default=os.path.join(os.path.dirname(__file__)),
                                    help='configfile file name.',
                                    required=False)

    logger_arguments = parser.add_argument_group('logger arguments')
    logger_arguments.add_argument('--log_level',
                                  type=str,
                                  default='INFO',
                                  help='log level',
                                  choices=['ERROR', 'error', 'WARNING', 'warning',
                                           'INFO', 'info', 'DEBUG', 'debug'],
                                  required=False)
    logger_arguments.add_argument('--log_file',
                                  type=str,
                                  help='log file (use the stderr by default). Disable the log colors.',
                                  required=False)

    return parser


def parametersLogger(args):
    '''
    Logger setting.

    :param args: argparse.ArgumentParser object
    :type args: class:'argparse.ArgumentParser'
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


def foo(param):
    """[Summary]

    :param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
    :type [ParamName]: [ParamType](, optional)
    ...
    :raises [ErrorType]: [ErrorDescription]
    ...
    :return: [ReturnDescription]
    :rtype: [ReturnType]
    """

    """
    The sphinx docstring style is prefered. 
    See: https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html
    """

    return type(param)


def main(args):
    """Main function. 
    """
    try:
        foo(args)
    except IOError:
        # prefer argparse exit over sys.exit
        args.exit(status=0)


def entrypoint():
    # Get CLI Args values:
    args_parser = arguments_parser()
    args = args_parser.parse_args()

    # Launch main
    main(args_parser, args)


if __name__ == "__main__":
    entrypoint()
