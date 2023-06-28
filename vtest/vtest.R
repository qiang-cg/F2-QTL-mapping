v.test <- function(p1,p2,f,c=1,h=1){
  var.p1 <- var(p1, na.rm=T)
  np1 <- length(p1)-length(which(is.na(p1)))
  var.p2 <- var(p2, na.rm=T)
  np2 <- length(p2)-length(which(is.na(p2)))
  var.p <- var(c(p1,p2), na.rm=T)
  var.f <- var(f, na.rm=TRUE)
  v <- (var.p-var.p1/(4*np1)-var.p2/(4*np2))/(var.f*c*h)
  ff <- length(f)-length(which(is.na(f)))
  p <- pf(v, df1=1, df2=(ff-1), lower.tail=FALSE)
  df <- data.frame(v.value=v, p.value=p)
  return(df)
}

v.test.2 <- function(p1,p2,f,c=1,h=1){
  var.p1 <- var(p1, na.rm=T)
  np1 <- length(p1)-length(which(is.na(p1)))
  var.p2 <- var(p2, na.rm=T)
  np2 <- length(p2)-length(which(is.na(p2)))
  #var.p <- var(c(p1,p2), na.rm=T)
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

#hert.t <- function(p1,p2,f){
#  #var.p <- var(c(p1,p2), na.rm=T)
#  var.p1 <- var(p1, na.rm=T)
#  var.p2 <- var(p2, na.rm=T)
#  var.f <- var(f, na.rm=TRUE)
#  h2 <- (var.f-var.p1-var.p2)/var.f
#  return(h2)
#}
  
setwd("D:\\2018-2022\\F2群体\\v test")

par <- read.csv("parent_phe.csv", head=T)
names(par)[21] <- "GWE"
N1_phe <- subset(par, Parent=="N")
R1_phe <- subset(par, Parent=="R")
f2_phe <- read.csv("f2_phe.csv", head=T)

#fh.h2 <- hert.t(N1_phe$FH, R1_phe$FH, f2$FH)

v.list <- list()
phe <- names(f2_phe)
for(i in phe){
  p1 <- subset(N1_phe, select=i)[,1]
  p2 <- subset(R1_phe, select=i)[,1]
  f2 <- subset(f2_phe, select=i)[,1]
  #h2 <- hert.t(p1,p2,f2)
  v.list[[i]] <- v.test.2(p1, p2, f2, c=2, h=0.8)
  #v.list[[i]]$h2 <- h2
}
 v.test.df <- do.call(rbind, v.list)

v.test.df$v.value <- round(v.test.df$v.value, 3)
v.test.df$p.value <- signif(v.test.df$p.value, 3)

write.csv(v.test.df, file="230510_f2_vtest.csv", row.names=T, quote=F)

#################
(66.656-1.759/(4*25)-1.801/(4*26))/(4.678*1*0.52)

(66.656-1.759/(25)-1.801/(26))/(4.678*1*0.52) 

p1_mean <- mean(p1, na.rm=T)
p2_mean <- mean(p2, na.rm=T)
p_var <- var(c(p1_mean, p2_mean))

v <- (p_var-var.p1/(4*np1)-var.p2/(4*np2))/(var.f*c*h)

##################

shapiro.test(f2_phe$ANL)
qqnorm(f2_phe$ANL)
qqline(f2_phe$ANL)

library(car)
tra <- powerTransform(f2_phe$FH)
trans <- bcPower(f2_phe$ANL, tra$lambda)

########
hist(f2_phe$FH)
FH <- f2_phe$FH
bc <- boxcox(FH~1, lambda=seq(-5,5,1/10))
lambda <- bc$x[which.max(bc$y)]
FH_bc <- (f2_phe$FH^lambda - 1) / lambda
hist(FH_bc)
var(FH_bc)

FH_tf <- powerTransform(FH)$lambda
FH_trans <- bcPower(FH, FH_tf)

FH_t <- FH^3
hist(FH_t)

v.test.2(p1, p2, FH_bc , c=2, h=0.8)

###########################################################
###使用sommer估算H2###
library(sommer)

setwd("D:\\2018-2022\\F2群体\\v test")
pheno <- read.csv("230528_F2_pheno.csv")
geno  <- read.csv("230528_F2_genotype.csv", row.names=1)
markp <- read.csv("230528_F2_maker_position.csv")

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


#####################################################
phe_diff <- phe[c(1:3,5:17)]
t_sta <- c(); p_value <- c()
for(p in phe_diff){
  test <- t.test(R1_phe[,p], N1_phe[,p], na.rm=T)
  t_sta <- c(t_sta, test$statistic)
  p_value <- c(p_value, test$p.value)
}

t_res <- data.frame(t_statistic=round(t_sta,3), p_value=p_value)

write.csv(t_res, file="230613_par_t-test.csv", quote=F, row.names=F)



















  
  
