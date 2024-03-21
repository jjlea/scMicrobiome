#!/bin/bash

GWAS_file=$1
myproject_file=$2

cd myproject_file

for k in *.tsv
do
  Rscript $myproject_file/removeMHC.R $k
done

