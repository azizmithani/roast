#!/bin/bash
conda install -c conda-forge openjdk --yes
conda install -c conda-forge gcc  --yes
conda install -c "bioconda/label/cf201901" blast --yes
conda install -c bioconda samtools=1.9 --yes
conda install -c bioconda bamtools --yes
conda install -c conda-forge boost --yes
