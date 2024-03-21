#!/bin/bash

GWAS_file=$1
outfile1=$2
magmadir=$3
myproject_file=$4
scfile=$5
## scfile is the output file.



cd $GWAS_file
for k in  *_removeMHC.tsv
do
$magmadir/soft/magma/magma \
--bfile $magmadir/magma/1000G_data/g1000_eur \
--pval $outfile1/${k/.tsv/.magmainput}  use='rsID,P.weightedSumZ' ncol='N' \
--gene-annot $outfile1/${k/.tsv/.genes.annot} \
--out $outfile1/${k/.tsv/} 
done

## extract z-scores
cd $outfile1
for k in *_removeMHC.genes.out
do 
cat $outfile1/$k | awk -v OFS="\t" '{print $1,$8}' > $scfile/${k/.out/.zscore}
done

## merge z-scores
cd $scfile
python $myproject_file/main-codes/merge_tables.py *.genes.zscore > zscore_merged.txt

## geneid
R $myproject_file/main-codes/geneid_trans.R $scfile $magmadir


# make .gs file
scdrs munge-gs \
    --out-file zscore_merged.gs \
    --zscore-file  zscore_merged2.txt  \
    --weight zscore \
    --n-max 1000


