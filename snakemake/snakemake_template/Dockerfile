FROM continuumio/miniconda3:4.10.3p0-alpine

LABEL Name='<name>' Version=1.0
ENV USERNAME=user
ENV GROUPID=1000
ENV USERID=1000

# Declare arguments at build time
# ARG USERNAME
# ARG GROUP
# ARG USERID

################### Copy Necessary Files For Image ###################
WORKDIR /pipeline
COPY . /pipeline

################### Config Conda And Create Snakemake Environment ###################
RUN conda config --set always_yes yes --set changeps1 no --set add_pip_as_python_dependency no \
    && conda update -n base -c defaults conda \
    && conda install -c conda-forge mamba \
    # openssl version is important for pysftp to function
    && mamba create -n snakemake -c conda-forge  -c bioconda snakemake \ 
    # beware of channel priority, for example snakemake is still available on bioconda for legacy so it might cause issues trying to install a newer version
    && conda config --set channel_priority strict \ 
    # create environments at build time, nothing more will be downloaded at runtime
    && mamba run -n snakemake snakemake --cores 2 --use-conda --conda-create-envs-only \
    && mamba clean -a \
    ### tag:latest
    # && addgroup --gid ${GROUPID} ${USERNAME} \
    # && adduser --disabled-password --gecos '' --uid ${USERID} --gid ${GROUPID} ${USERNAME} \	
    ### hard-declared user credentials ###
    # && addgroup --gid 1000 user \
    # && adduser --disabled-password --gecos '' --uid 1000 --gid 1000 user \ 
    ### tag:alpine
    # put at user level for it to be effective
    && addgroup -S usergroup \
    && adduser -S user -G usergroup -D -g '' \
    && chown -R 1000:1000 /pipeline \
    && chmod -R 755 /pipeline \
    && chown -R 1000:1000 /opt/conda/ \
    && chmod -R 755 /opt/conda/ 

USER user

################### Run pipeline ###################
ENTRYPOINT [ "mamba" , "run" , "-n" ,"snakemake","snakemake"]