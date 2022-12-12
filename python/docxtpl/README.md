# Docxtpl & DocxCompose

Python modules to use a docx as a jinja2 template and combine them.

This template aim to give structure to generate multiple word and combine them in one.

This template is mostly based on the work of Jordan LANGLOIS.

## Structure

The program use a yaml as an input file. Each main key will match a template name and it's value the template content.

In our test with have the following yaml file:

```yml
test_template:
  content:
    foo: oof
    bar: rab
    fizz: zzif
    buzz: zzub
  list_example:
  - name: item_1
    content: "item_1 content here"
  - name: item_2
    content: "item_2 content here"
```

The main key is `test_template` which word file can be found in `resources/template`.

This key will serve as input in the config file `resources/configs/blocks.yml` in order to retreive which word file and jsonschema should be use.

The first word used will fix the style of the whole document. It is recommended to use a main word file wiche will be the only one to maintain.

> /!\\: The json schema validation is not implemented yet.

## Test

Use conda to create an environment with the necessary module:

```bash
conda env create -f environment.yml

conda activate docxtpl_test
```

Then try it for yourself:

```bash
./docxtpl_template.py -i test/test_template.yml -o test/test_template.docx
```