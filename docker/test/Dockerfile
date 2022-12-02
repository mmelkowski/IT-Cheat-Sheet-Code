FROM mambaorg/micromamba:latest

COPY environment.yml /tmp/conda-tmp/environment.yml

RUN /bin/micromamba install -f /tmp/conda-tmp/environment.yml

CMD micromamba list
