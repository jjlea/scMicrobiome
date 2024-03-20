## Linking microbial genetic signals and Single-Cell Transcriptomics to Unravel Host-microbe Interactions across Human Tissues



### Introduction

In this study, we introduce a computational framework (scBPS) to incorporate microbial GWAS summary data with human single-cell transcriptomic data for the discovery of critical cellular contexts that are highly associated with gut microbial taxon. 

![image](https://github.com/jjlea/scMicrobiome/assets/73264824/77cbab3f-4ca5-445f-87c8-0d3bbedf4d9d)


### Documentation

**1**, GWAS summary data are preprocessed to generate z-scores for microbe-gene associations. We rank the relevance of genes based on z-scores and choose top 1,000 genes for each taxon as microbe-associated genes, analogue to previous studies. 

**2**, we collapse the expression of the putative microbe-associated genes in a given cell to calculate a raw Bacteria Polygenic Score (BPS) for each cell; to enhance the statistical power, these putative microbe-associated genes are weighted by corresponding MAGMA-estimated z-scores and inversely weighted by gene-specific technical noise scores, which are estimated by the mean-variance relationship across genes. The raw BPSs are subsequently subjected to gene set-wise and cell-wise standardization for obtaining a normalized BPS, which measures the strength of association between each taxon and each host cell. 

**3**, to distinguish taxon-relevant cell types or tissue types, we estimate the relevance of each cell population using a BPS<sub>AUC</sub> metric, which is determined by calculating the area under the recovery curve (AUC) of cells in the predefined population across the ranking of all cells’ BPSs in the single-cell atlas 35. To assess statistical significance, scBPS generates 1,000 sets of control BPS<sub>AUC</sub> by randomly permuting the rank of all the cells in the single-cell atlas. Monte Carlo (MC) P value for each “taxon-cell population” pair is computed based on the empirical distribution of the BPS<sub>AUC</sub> across the corresponding 1,000 control BPS<sub>AUC</sub> (see Methods). 

**4**, the output of our framework includes: (1) the BPS values of each individual host cell for each microbial taxon; (2) BPS<sub>AUC</sub> scores representing the strength of microbe-cell type associations and the empirical P values. 


### Data analysis

- Preprocess of microbial GWAS summary data:
  
``` bash
bash run_removeMHC.sh <GWAS_file> <myproject_file>
```

- GWAS data SNP annotation:

`bash annotation.sh <GWAS_file> > <magmafile> <magmadir>`

- gene z-scores:

```bash
mkdir scfile
bash gene_anno.sh <GWAS_file> > <magmafile> <magmadir> <myproject_file> <scfile>
```


- Top 1000 putative microbial genes:



- Calculating BPS scores:



- Calculating BPS<sub>AUC</sub> scores:



- Monte-Carlo p-values:






### Resources

The microbial GWAS summary data of the Dutch microbiome project were downloaded at `https://dutchmicrobiomeproject.molgeniscloud.org`. 

The Tabula Sapiens human single-cell transcriptome data were downloaded at `https://tabula-sapiens-portal.ds.czbiohub.org/`. 

The GWAS summary data of MiBioGen project were downloaded at `https://www.mibiogen.org/`. 

The GWAS summary data of 10 liver-associated diseases were downloaded from the FinnGen database at `https://r8.risteys.finngen.fi/`. 

All codes used for the analyses and visualization are provided in the Github repository at `https://github.com/jjlea/scMicrobiome`. 




