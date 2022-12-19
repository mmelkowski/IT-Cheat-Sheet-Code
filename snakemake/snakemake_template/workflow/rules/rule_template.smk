rule template:
    input:
        "path/to/input.file",  # all paths are relative to the working directory ,
    output:
        "{outdir}/results/<wildcards>/file.<ext>",  # It's advisable to have a {outdir} wildcard to control the absolute path of all your output files it's value will have to be specified in the config.yaml and referenced in all_input()
    log:
        "{outdir}/logs/<wildcards>/<func_tool>.log",  # Some tools write to stdout, depending on that you can redirect stdout/stderr or both to your log files 1> | 2> | &> 
    benchmark:
        "{outdir}/benchmark/<wildcards>/<func_tool>.tsv"  # Record exec statistics in a .tsv file
    container:
        "docker://condaforge/mambaforge:4.14.0-0"  # Run command inside container, can be coupled with conda to create the env in the container for maximum control over software version
    conda:
        "../envs/template.yaml"  # Can accept a file name or a named environment, for production it's always better to specify a file to create the environment from
    params:
        extra=config["function"]["tool"]["extra"],  # It's usually a grood practice to leave a back door for users to add extra params/flags
    shell:
        "command {input} > {output} {params.extra} 2> {log}"
