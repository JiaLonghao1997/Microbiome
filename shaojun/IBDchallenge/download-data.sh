#!/usr/bin/env bash

set -ev

mkdir -p data/testset_subchallenge2_files
cd data/testset_subchallenge2_files

wget https://sbvimprover-datasets.s3-eu-west-1.amazonaws.com/testset_subchallenge2_files.zip
unzip testset_subchallenge2_files.zip
chmod 644 *.txt
rm testset_subchallenge2_files.zip

