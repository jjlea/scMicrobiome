## scBPS: Mapping genome-wide polygenic signals of gut microbiota to human single cell transcriptome profiles uncovers taxon-relevant cellular features and health implications


### Introduction


Microbial genome-wide association studies (GWAS) have uncovered numerous host genetic variants significantly associated with gut microbiota. However, the host genetic control of gut microbiome through key biological pathways within specific cellular context remains largely unknown. Here, we introduce a computational framework for integrating microbial GWAS and single-cell RNA-sequencing profiles of 24 human organs to distinguish critical host tissues and cell types relevant to gut microbes. In this study, we introduce a computational framework (scBPS) to incorporate microbial GWAS summary data with human single-cell transcriptomic data for the discovery of critical cellular contexts that are highly associated with gut microbial taxon. 

![image](https://github.com/jjlea/scMicrobiome/assets/73264824/77cbab3f-4ca5-445f-87c8-0d3bbedf4d9d)


By leveraging the largest metagenome-based GWAS data from the Lifelines Dutch Microbiome Project 12 and the high-quality multi-organ single-cell transcriptomic atlas from the Tabula Sapiens Consortium 31, we comprehensively explored the associations of 207 gut microbial taxa with 24 different human tissues and 254 cell types. We exhibit that scBPS performs well in identifying microbe-relevant cellular features at both tissue and cell-type levels. Notably, our results showcase the advantage of our framework by focusing on a hepatocyte-interacting genus Collinsella, illustrating its interaction with central-veinal hepatocyte subpopulation via modulating the activity of cholesterol metabolism. Finally, our study delves into the cellular mechanisms influenced by gut microbes that underlie 10 liver diseases.

### Documentation

**First**, GWAS summary data are preprocessed to generate z-scores for microbe-gene associations. We rank the relevance of genes based on z-scores and choose top 1,000 genes for each taxon as microbe-associated genes, analogue to previous studies. 

**Second**, we collapse the expression of the putative microbe-associated genes in a given cell to calculate a raw Bacteria Polygenic Score (BPS) for each cell; to enhance the statistical power, these putative microbe-associated genes are weighted by corresponding MAGMA-estimated z-scores and inversely weighted by gene-specific technical noise scores, which are estimated by the mean-variance relationship across genes. The raw BPSs are subsequently subjected to gene set-wise and cell-wise standardization for obtaining a normalized BPS, which measures the strength of association between each taxon and each host cell. 

**Third**, to distinguish taxon-relevant cell types or tissue types, we estimate the relevance of each cell population using a BPSAUC metric, which is determined by calculating the area under the recovery curve (AUC) of cells in the predefined population across the ranking of all cells’ BPSs in the single-cell atlas 35. To assess statistical significance, scBPS generates 1,000 sets of control BPS~AUC~ by randomly permuting the rank of all the cells in the single-cell atlas. Monte Carlo (MC) P value for each “taxon-cell population” pair is computed based on the empirical distribution of the BPSAUC across the corresponding 1,000 control BPS~AUC~ (see Methods). 

**Finally**, the output of our framework includes: (1) the BPS values of each individual host cell for each microbial taxon; (2) BPS~AUC~ scores representing the strength of microbe-cell type associations and the empirical P values. 


### Citation



### Resources

The microbial GWAS summary data of the Dutch microbiome project were downloaded at https://dutchmicrobiomeproject.molgeniscloud.org. The Tabula Sapiens human single-cell transcriptome data were downloaded at https://tabula-sapiens-portal.ds.czbiohub.org/. The GWAS summary data of MiBioGen project were downloaded at https://www.mibiogen.org/. The GWAS summary data of 10 liver-associated diseases were downloaded from the FinnGen database at https://r8.risteys.finngen.fi/. All codes used for the analyses and visualization are provided in the Github repository at https://github.com/jjlea/scMicrobiome. 



