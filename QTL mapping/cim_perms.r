library(qtl)
library(rlecuyer)
library(snow)
library(stringr)

setwd("/mnt/ge-jbod/qiangchenggen/datapool/f2_data/cim")
cross <- read.cross(format ="csvr", file="f2.bin_5k_final.csv", estimate.map=F, genotypes=c("a", "h", "b"),
                  alleles=c("a", "b"),na.strings=c("-", "NA"))
cross<-calc.genoprob(cross, map.function="kosambi")

phenotype <- colnames(cross$pheno)
perms_lt <- list()
for(phe in phenotype){
  phe_col <-  which(phenotype==phe)
  perms <- cim(cross, pheno.col=phe_col, n.marcovar=3, window=10, method="hk", map.function="kosambi", n.perm=1000)
  perms_lt[[phe]] <- perms
}
 
save.image("cim_perms.RData") 
 
