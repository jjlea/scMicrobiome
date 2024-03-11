library(data.table)
library(TwoSampleMR)
library(data.table)
library("MendelianRandomization")
setwd("MyDIR")

# bacteria

bac <-fread('MyDIR/GWAS_k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Coriobacteriales.f__Coriobacteriaceae.g__Collinsella.tsv',header=T)
bac$phenotype <- 'Collinsella'


# hypercholesterol
chol <- fread("MyDIR/finngen_R8_E4_HYPERCHOL.gz", header = T)
chol$phenotype <- 'HYPERCHOL' 




# bacteria -> chol

bac_mr_chol <- function(bac, chol){
    bac_expo <- format_data(
                    dat=bac,
                    type = "exposure",
                    phenotype_col = "phenotype",
                    snp_col = "id",
                    beta_col = "beta",
                    se_col = "SE",
                    effect_allele_col ="alt",
                    other_allele_col = "ref",
                    eaf_col = "AF_Allele2",
                    pval_col = "pval"
                    
                     )
    bac_expo2 <- bac_expo[which(bac_expo$pval.exposure < 1e-5),]
    bac_expo2 <- clump_data2(bac_expo2, clump_r2 = 0.1)
    chol_out <- format_data(
                     dat=chol,
                     type = "outcome",
                     snps = bac_expo2$SNP,
                     header = TRUE,
                     phenotype_col = "phenotype",
                     snp_col = "rsids",
                     beta_col = "beta",
                     se_col = "sebeta",
                     effect_allele_col = "alt",
                     other_allele_col = "ref",
                     eaf = "af_alt",
                     pval_col = "pval"
                     
                    )            
    mydata_bac_chol <- harmonise_data(
           exposure_dat=bac_expo2,
           outcome_dat=chol_out,
           action= 2
           )
 #   res_bac_chol <- mr(mydata_bac_chol,method_list=c("mr_ivw_mre"))
    return(list( bac_expo=bac_expo2, chol_out=chol_out, harm=mydata_bac_chol))
}



# chol -> bacteria
chol_mr_bac <- function(chol, bac){
    chol_exp <- format_data(
                     dat=chol,
                     type = "exposure",
                     snp_col = "rsids",
                     beta_col = "beta",
                     se_col = "sebeta",
                     effect_allele_col = "alt",
                     other_allele_col = "ref",
                     eaf = "af_alt",
                     pval_col = "pval",
                     phenotype_col = "phenotype"
                        )
    chol_exp2 <- chol_exp[which(chol_exp$pval.exposure < 1e-5),]
    chol_exp2 <- clump_data2(chol_exp2, clump_r2 = 0.1)
    bac_out <- format_data(
                    dat=bac,
                    type = "outcome",
                    snps = chol_exp2$SNP,
                    header = TRUE,
                    phenotype_col = "phenotype",
                    snp_col = "id",
                    beta_col = "beta",
                    se_col = "SE",
                    effect_allele_col ="alt",
                    other_allele_col = "ref",
                    eaf_col = "AF_Allele2",
                    pval_col = "pval"
                    )
        mydata_chol_bac <- harmonise_data(
               exposure_dat=chol_exp2,
               outcome_dat=bac_out,
               action= 2
               )
  #      res_chol_bac <- mr(mydata_chol_bac)
        return(list(chol_exp=chol_exp2, bac_out=bac_out, harm=mydata_chol_bac))
}

clump_data2 <- function (dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, 
  clump_p2 = 1, pop = "EUR") 
{
  pval_column <- "pval.exposure"
  if (!is.data.frame(dat)) {
    stop("Expecting data frame returned from format_data")
  }
  if ("pval.exposure" %in% names(dat) & "pval.outcome" %in% 
    names(dat)) {
    message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
  }
  else if (!"pval.exposure" %in% names(dat) & "pval.outcome" %in% 
    names(dat)) {
    message("pval.exposure column not present, using pval.outcome column for clumping.")
    pval_column <- "pval.outcome"
  }
  else if (!"pval.exposure" %in% names(dat)) {
    message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
    dat$pval.exposure <- 0.99
  }
  else {
    pval_column <- "pval.exposure"
  }
  if (!"id.exposure" %in% names(dat)) {
    dat$id.exposure <- random_string(1)
  }
  d <- data.frame(rsid = dat$SNP, pval = dat[[pval_column]], 
    id = dat$id.exposure)
  out <- ieugwasr::ld_clump(d, clump_kb = clump_kb, clump_r2 = clump_r2, 
    clump_p = clump_p1,  bfile="MyDIR/9.R-TwoSampleMR/EUR",
                           plink_bin="plink")
  keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, 
    out$id)
  return(dat[keep, ])
}



harm_Coll_chol <- bac_mr_chol(bac, chol)
harm_chol_Coll <- chol_mr_bac( chol, bac)
lst <-  c("mr_ivw","mr_ivw_mre","mr_ivw_radial","mr_two_sample_ml","mr_egger_regression","mr_penalised_weighted_median",
                                                         "mr_weighted_median","mr_weighted_mode")
res_Coll_chol <-  mr(harm_Coll_chol$harm, method_list = lst)
res_chol_Coll <-  mr(harm_chol_Coll$harm, method_list = lst)

res <- rbind(res_Coll_chol, res_chol_Coll)
res$ID <- paste(paste( res$exposure,res$outcome, sep="_"), res$method, sep="_")
res$lower <- res$b-1.96*res$se
res$upper <- res$b+1.96*res$se
res$method <- factor(res$method, levels =unique(res$method) )
res <- res[order(res$exposure,res$outcome, res$method),]
res$b2 <- paste0(round(res$b,2)," (", round(res$lower,2),", ", round(res$upper,2),")")
res$mthd <- paste(res$exposure, res$outcome, sep="_")
res[order(res$exposure, res$outcome, res$method),]
#save(res, file="res_Collinsella.RData")