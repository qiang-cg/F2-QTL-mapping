v.test.2 <- function(p1,p2,f,c=1,h=1){
  var.p1 <- var(p1, na.rm=T)
  np1 <- length(p1)-length(which(is.na(p1)))
  var.p2 <- var(p2, na.rm=T)
  np2 <- length(p2)-length(which(is.na(p2)))
  var.f <- var(f, na.rm=TRUE)
  p_var <- var(c(mean(p1, na.rm=T), mean(p2, na.rm=T)))
  v <- (p_var-var.p1/(4*np1)-var.p2/(4*np2))/(var.f*c*h)
  #v <- round(v,3)
  ff <- length(f)-length(which(is.na(f)))
  p <- pf(v, df1=1, df2=(ff-1), lower.tail=FALSE)
  #p <- round(p, 4)
  df <- data.frame(v.value=v, p.value=p, H2=h)
  return(df)
}

###########################################################
###使用sommer估算H2###
library(sommer)

pheno <- read.csv("230528_F2_pheno.csv")
geno  <- read.csv("230528_F2_genotype.csv", row.names=1)
markp <- read.csv("230528_F2_maker_position.csv")
par <- read.csv("parent_phe.csv", head=T)

names(pheno)[1] <- "id"

pheno$idd <- pheno$id; pheno$ide <- pheno$id
row.names(pheno) <- pheno$id

geno[geno=="a"] <- -1
geno[geno=="h"] <- 0
geno[geno=="b"] <- 1

rn <- row.names(geno)
geno <- as.data.frame(lapply(geno, as.numeric))
row.names(geno) <- rn
geno <- as.matrix(geno)

A <- A.mat(geno)
D <- D.mat(geno)
E <- E.mat(geno)

#FH.ADE <- mmer(FD~1, random=~vsr(id, Gu=A)+vsr(idd,Gu=D),
#               rcov=~units, nIters=3, data=pheno)
#(summary(FH.ADE)$varcomp)			   
#vpredict(FH.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) )			   
			   
################
phe <- names(pheno)[c(2:17)]

H2_al <- list()
for(p in phe){
  pheno_t <- pheno
  names(pheno_t)[which(names(pheno_t)==p)] <- "pht"
  ADE <- mmer(pht~1, random=~vsr(id, Gu=A)+vsr(idd,Gu=D),
               rcov=~units, nIters=20, data=pheno_t)
  H2 <- vpredict(ADE, h2 ~ (V1+V2) / ( V1+V2+V3) )
  H2_al[[p]] <- H2
}
names(H2_al)[1] <- "FH"
names(H2_al)[2] <- "CL"

v.list <- list()
phe <- names(f2_phe)
for(i in phe){
  p1 <- subset(N1_phe, select=i)[,1]
  p2 <- subset(R1_phe, select=i)[,1]
  f2 <- subset(f2_phe, select=i)[,1]
  H2 <- H2_al[[i]][1,1]
  v.list[[i]] <- v.test.2(p1, p2, f2, c=2, h=H2)
}
 v.test.df <- do.call(rbind, v.list)

v.test.df$v.value <- round(v.test.df$v.value, 3)
v.test.df$p.value <- signif(v.test.df$p.value, 3)

write.csv(v.test.df, file="230530_f2_vtest.csv", row.names=T, quote=F)

save.image("230530_F2_vtest.RData")


  
  
