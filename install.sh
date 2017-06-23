#!/usr/bin/env bash

IS_CONDA="$(which conda)"

if [ ${#IS_CONDA} -eq "0" ]
    then
        echo "It seems that you do not have conda installed. Download it from https://www.continuum.io/downloads"
        exit 1
fi

echo "Removing preexisting asecount environment if exists"
source deactivate
conda remove --name asecount --all --yes;

echo "Creating the new asecount environment from environment.yml"
conda env create -f envs/environment.yml;

if [ $? -ne 0 ]
    then
        echo "Problem creating the conda environment. I don't know what to tell you."
        exit 1
fi

echo "Downloading the test data."
rm -rf test;
https://s3-us-west-1.amazonaws.com/mdurrant/biodb/bundles/ase_read_counting_workflow/test.tar.gz;
tar -zxvf test.tar.gz;
rm test.tar.gz;

if [ $? -ne 0 ]
    then
        echo "Could not download the test data"
        exit 1
fi



echo "------------------------------------------------------------------------------------------------------"
echo "INSTALLATION COMPLETE"
echo "------------------------------------------------------------------------------------------------------"
echo "NOTE: You must specify the required information in the config.yaml file to run the workflow properly."
echo "NOTE: Once configured, enter:"
echo ""
echo "> source activate asecount"
echo "> snakemake --use-conda"
echo ""
echo "In this directory to run the workflow."
echo "------------------------------------------------------------------------------------------------------"