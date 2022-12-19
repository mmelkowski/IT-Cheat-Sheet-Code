def all_input(wildcards):
    wanted_input = []
    wanted_input.extend(
        expand(
            "{outdir}/{unit.sample_id}/{unit.sample_id}__{unit.ref}.sam",
            outdir=outdir,
            unit=get_unit(),
            )
    )

    return wanted_input
