# hiv-profile-sampling

## Installation

Install Miniconda 3 with:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh

After setting up Miniconda, install dependencies from the kantorlab channel with:

    conda create -n hiv-profile-sampling -c kantorlab biopython=1.73 mafft=7.313 raxml=8.2.12 scons=3.0.1.1

To activate the environment, use:

    source activate hiv-profile-sampling

Setup your scratch directory on Oscar:

    mkdir -p /gpfs/scratch/$USER/hiv-profile-sampling
    ln -s /gpfs/scratch/$USER/hiv-profile-sampling scratch

## Running

The run order and dependencies of the scripts are specified in the SConstruct
file.  The entire analysis can be run by executing the `scons` command from the
root directory of the repo.

