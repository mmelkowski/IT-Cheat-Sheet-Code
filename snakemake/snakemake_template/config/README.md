> Describe how to configure the workflow (using config.yaml and maybe additional files).All of them need to be present with example entries inside of the config folder.

# Configuration

To use this workflow each operator has to manipulate two configuration files
the first is a general parameters file in `.yaml` format the second is a tabular file describing fastq files and their metadata, each line should have a unique `sample_id` and represent a unique sample.

### Parameters : `config/config.yaml`

```yaml
run_specific:
  operator_name: "operator_name" # <str>
  operator_mail: "operator@email.com" # <str>
  project_name: "someproject" # <str>
```

### Experimental design : `config/units.tsv`

| sample_id | fq1 | fq2 | ref |
| --- | ---| --- | ---|
| `<str>` | `<path>` | `<path>` | `<path>` |
