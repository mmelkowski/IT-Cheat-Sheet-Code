import subprocess
import os
import tempfile
import shutil
import pytest
import sys
import yaml
from itertools import chain

DIFF_MASTER = os.environ.get("DIFF_MASTER", "false") == "true"
DIFF_LAST_COMMIT = os.environ.get("DIFF_LAST_COMMIT", "false") == "true"

if DIFF_MASTER or DIFF_LAST_COMMIT:

    compare = "HEAD^" if DIFF_LAST_COMMIT else "origin/master"

    # check if wrapper is modified compared to master
    DIFF_FILES = set(
        subprocess.check_output(["git", "diff", compare, "--name-only"])
        .decode()
        .split("\n")
    )

CONTAINERIZED = os.environ.get("CONTAINERIZED", "false") == "true"


class Skipped(Exception):
    pass


skip_if_not_modified = pytest.mark.xfail(raises=Skipped)


def run(wrapper, cmd, check_log=None):
    origdir = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        dst = os.path.join(d, "main")
        os.makedirs(dst, exist_ok=True)
        copy = lambda pth, src: shutil.copy(
            os.path.join(pth, src), os.path.join(dst, pth)
        )

        used_wrappers = []
        wrapper_file = "used_wrappers.yaml"
        if os.path.exists(os.path.join(wrapper, wrapper_file)):
            # is meta wrapper
            with open(os.path.join(wrapper, wrapper_file), "r") as wf:
                wf = yaml.load(wf, Loader=yaml.BaseLoader)
                used_wrappers = wf["wrappers"]
        else:
            used_wrappers.append(wrapper)

        for w in used_wrappers:
            success = False
            for ext in ("py", "R", "Rmd"):
                script = "wrapper." + ext
                if os.path.exists(os.path.join(w, script)):
                    os.makedirs(os.path.join(dst, w), exist_ok=True)
                    copy(w, script)
                    success = True
                    break
            assert success, "No wrapper script found for {}".format(w)
            copy(w, "environment.yaml")

        if (DIFF_MASTER or DIFF_LAST_COMMIT) and not any(
            any(f.startswith(w) for f in DIFF_FILES)
            for w in chain(used_wrappers, [wrapper])
        ):
            raise Skipped("wrappers not modified")

        if any(
            yaml.load(open(os.path.join(w, "meta.yaml")), Loader=yaml.BaseLoader).get(
                "blacklisted"
            )
            for w in used_wrappers
        ):
            raise Skipped("wrapper blacklisted")

        testdir = os.path.join(d, "test")
        # pkgdir = os.path.join(d, "pkgs")
        shutil.copytree(os.path.join(wrapper, "test"), testdir)
        # prepare conda package dir
        # os.makedirs(pkgdir)
        # switch to test directory
        os.chdir(testdir)
        if os.path.exists(".snakemake"):
            shutil.rmtree(".snakemake")
        cmd = cmd + [
            "--wrapper-prefix",
            "file://{}/".format(d),
            "--conda-cleanup-pkgs",
            "--printshellcmds",
            "--show-failed-logs",
        ]

        if CONTAINERIZED:
            # run snakemake in container
            cmd = [
                "sudo",
                "docker",
                "run",
                "-it",
                "-v",
                "{}:{}".format(os.getcwd(), "/workdir"),
                "snakemake/snakemake",
                " ".join(cmd),
            ]

        # env = dict(os.environ)
        # env["CONDA_PKGS_DIRS"] = pkgdir
        try:
            subprocess.check_call(cmd)
        except Exception as e:
            # go back to original directory
            os.chdir(origdir)
            logfiles = [
                os.path.join(d, f)
                for d, _, files in os.walk(os.path.join(testdir, "logs"))
                for f in files
            ]
            for path in logfiles:
                with open(path) as f:
                    msg = "###### Logfile: " + path + " ######"
                    print(msg, "\n")
                    print(f.read())
                    print("#" * len(msg))
            if check_log is not None:
                for f in logfiles:
                    check_log(open(f).read())
            else:
                raise e
        finally:
            # cleanup environments to save disk space
            subprocess.check_call(
                "for env in `conda env list | grep -P '.snakemake/conda' | "
                "cut -f1 | tr -d ' '`; do conda env remove --prefix $env; done",
                shell=True,
            )
            # go back to original directory
            os.chdir(origdir)


@skip_if_not_modified
def test_metaquast():
    run(
        "bio/metaquast",
        ["snakemake", "--cores", "1", "metaquast_out", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_fastp_pe():
    run(
        "bio/fastp_pe",
        ["snakemake", "--cores", "1", "fastp_pe.json", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_fastp_se():
    run(
        "bio/fastp_se",
        ["snakemake", "--cores", "1", "fastp_se.json", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_megahit():
    run(
        "bio/megahit", ["snakemake", "--cores", "1", "megahit_out", "--use-conda", "-F"]
    )


@skip_if_not_modified
def test_blastn():
    run(
        "bio/blastn", ["snakemake", "--cores", "1", "blastn.out", "--use-conda", "-F"]
    )
