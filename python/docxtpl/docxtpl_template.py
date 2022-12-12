#!/usr/bin/env python3
'''
Je suis une docstring de dingue !
'''

import argparse
import os
# import jsonschema
import yaml

from docx import Document
from docxcompose.composer import Composer
from docxtpl import DocxTemplate, InlineImage


def arguments_parser():
    '''
    Arguments parsing.

    :return: Object
    '''
    parser = argparse.ArgumentParser(description='',
                                     usage='{} -i <input_path> -o <output_name>'.format(
                                         os.path.basename(__file__)),
                                     formatter_class=argparse.RawTextHelpFormatter)

    group0 = parser.add_argument_group('mandatory arguments')
    group0.add_argument('-i',
                        '--input',
                        type=str,
                        help='Path of the input YAML file.',
                        required=True)
    group0.add_argument('-o',
                        '--output',
                        type=str,
                        help='Name used for final word.',
                        required=True)

    args = parser.parse_args()
    return args


def read_yaml(file_name):
    '''
    Extract data from a YAML file to dictionary.

    :entry: File name
    :return: Dictionary type
    '''
    with open(file_name, 'r') as file:
        content = file.read()
    return yaml.safe_load(content)


def apply_context_to_a_template(context, template_file_name, output_name):
    '''
    Generates a docx file from context and the template files.

    :entry: dictionary, file name, string
    :return: output_name
    '''
    tpl = DocxTemplate(template_file_name)
    tpl.render(context)
    tpl.save(output_name)
    return output_name


def main(args):
    '''
    Main function
    '''

    print("[docxtpl_template][INFO]: Generation from", args.input)

    data = read_yaml(args.input)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    templates_path = os.path.join(dir_path, "resources/template")
    block_conf_path = os.path.join(dir_path, "resources/configs/blocks.yml")
    block_conf = read_yaml(block_conf_path)

    # for all key in the input yaml create the tmp result word file.
    tmp_files = []
    for block_name, block_context in data.items():
        template_file_path = os.path.join(templates_path,
                                          block_conf[block_name]['template_name'])
        docx_name = apply_context_to_a_template(block_context,
                                                template_file_path,
                                                '{}_tmp'.format(block_name))
        tmp_files.append(docx_name)

    # From All word generate 1 fusionned word file.
    # the first document fix the style for the following ones.
    composer = Composer(Document(tmp_files[0]))
    os.remove(tmp_files[0])
    for tmp_file in tmp_files[1:]:
        block = Document(tmp_file)
        composer.append(block)
        os.remove(tmp_file)
    composer.save(args.output)
    
    print("[docxtpl_template][INFO]: Document saved as", args.output)


if __name__ == '__main__':
    args = arguments_parser()
    main(args)
