#!/bin/bash

RUN=$1

export CONDA_ENVS_PATH=/var/tmp
export CONDA_PKGS_DIRS=/var/tmp
export HOME=/var/tmp/

snakemake --tes "" --use-conda --envvars CONDA_PKGS_DIRS CONDA_ENVS_PATH HOME --conda-prefix /var/tmp --jobs 40 \
	-s Archer_MCulen.smk \
	--verbose --reason --conda-frontend mamba --container-image cerit.io/ceitec/snakemake:v6.10.0
