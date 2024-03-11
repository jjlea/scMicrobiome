cd $MyDIR/LDSC-celltype
# ldsc path
ldsc=$MyDIR/soft/ldsc
# gwas path
gwasfile=$MyDIR




# munge_sumstats.py
for k in *.tsv
do
/share2/pub/lijj/lijj/soft/ldsc/munge_sumstats.py \
	--sumstats $MyDIR/$k\
	--snp id \
	--N-col num \
	--a1 ref \
	--a2 alt \
	--p pval \
	--info imputationInfo \
	--signed-sumstats beta,0 \
	--chunksize 500000 \
	--out ./LDSC_With_merge-alleles/${k//.tsv/} \
	--merge-alleles $MyDIR/w_hm3.snplist > ./LDSC_With_merge-alleles/${k//.tsv/}.munge_sumstats.log 2>&1 
done

# gwas file name
gwas=$MyDIR/gwas_removeMHC.sumstats.gz
outfile=$MyDIR/LDSC-celltype
ref=$MyDIR/LDSC-celltype
cts_name=Multi_tissue_gene_expr

${ldsc}/ldsc.py \
    --h2-cts ${gwasfile}/${gwas} \
    --ref-ld-chr ${ref}/1000G_EUR_Phase3_baseline/baseline. \
    --out ${outfile}/gwas_removeMHC.sumstats_${cts_name} \
    --ref-ld-chr-cts ${ref}/$cts_name.ldcts \
    --w-ld-chr ${ref}/weights_hm3_no_hla/weights.
