#!/bin/bash

GWAS_file=$1
outfile1=$2
magmadir=$3


cd $GWAS_file
for k in *_removeMHC.tsv
do
cat $GWAS_file/${k} | awk -v OFS=" " '{print $4,$2,$3,$11,$10}' > $outfile1/${k/.tsv/.magmainput}


$magmadir/soft/magma/magma \
--annotate window=10,10 \
--snp-loc $outfile1/${k/.tsv/.magmainput} \
--out $outfile1/${k/.tsv/} \
--gene-loc $mydir/magma/locations/NCBI37.3.gene.loc
done

