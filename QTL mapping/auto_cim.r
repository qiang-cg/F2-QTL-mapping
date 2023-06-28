library(qtl)
#library(rlecuyer)
#library(snow)
library(stringr)

setwd("/mnt/ge-jbod/qiangchenggen/datapool/f2_data/cim")
load("/mnt/ge-jbod/qiangchenggen/datapool/f2_data/sw/sc1_perms.RData")
sink(file="auto_cim.log", append=TRUE, type="output", split=TRUE)

###CIM模型使用IM置换检验结果得到阈值，CIM模型加入协变量后，会导致结果分布有偏，LOD值较高的位点聚集在一起，导致根据分位数取阈值时较高

#输入x为CIM扫描结果
select.cim <- function(x, threshold=3.5, mar.dist=20){
  th <-  which(x$lod >= threshold)
  l <- length(th)
  j=1
  qtl <- list()
  qtl[[j]] <- th[1]
  for (i in 2:l){
      if(th[i] - th[i-1] <= mar.dist) {     
        qtl[[j]] <- c(qtl[[j]],th[i])
      } else {
        j=j+1
        qtl[[j]] <- c(th[i])}
  }
  peak <- lapply(qtl,function(y){max(x[y,])})
  peak.n <- do.call(rbind, lapply(peak, data.frame, stringsAsFactors=FALSE))
  return(peak.n)
}

#n.marcovar从1开始，到n.marcovar大于或等于QTL数目，输出CIM结果
cim.auto.scan <- function(cross=cross, pheno.col=1, threshold=4){
    qtl.num <- 1
	n.marcovar <- 0
	i <- 1
	while(qtl.num > n.marcovar){
	    n.marcovar <- qtl.num
		out <- paste("CIM scaning", i, sep=" ")
		print(out)
		i<-i+1
	    cim.scan <- cim(cross, pheno.col=pheno.col, n.marcovar=n.marcovar, window=10, method="hk", map.function="kosambi")
		cim.qtl <- select.cim(cim.scan, threshold=threshold)
		qtl.num <- nrow(cim.qtl)
	}
	print("Intial QTL scaned by CIM")
    print(cim.qtl)	
	return(cim.scan)
}

###########根据CIM扫描结果输出QTL区间##########
cim.qtl.int <- function(cross=cross, pheno.col=1, cim=cim, threshold=threshold){
    cim_qtl <- select.cim(cim, threshold=threshold) 
    qtl <- makeqtl(cross, chr=cim_qtl$chr, pos=cim_qtl$pos, what="prob")
	qtl_Q <- paste("Q", 1:qtl$n.qtl, sep="")
	qtl_form_1 <- as.formula(paste("y~", paste(qtl_Q, collapse="+")))
    int_id <- c()  
	if(qtl$n.qtl > 1){
	    qtl_int <- addint(cross, pheno.col=pheno.col, qtl=qtl, formula=qtl_form_1, method="hk")
	    qtl_comb <- as.data.frame(t(combn(qtl_Q,2)))
		int_id <- which(qtl_int[[7]] <= 0.05)
        if(length(int_id) > 0){
		    print(paste(length(int_id), "pair of interactions detected"))		
	        for(i in int_id){
              qtl_Q <- c(qtl_Q, paste(qtl_comb[i,],collapse="*"))
            }
		qtl_form_2 <- as.formula(paste("y~", paste(qtl_Q, collapse="+")))
        rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl, formula=qtl_form_2, 
                              method="hk", model="normal", verbose=F)
		fit_f <- fitqtl(cross, pheno.col=pheno.col, qtl=rqtl, formula=qtl_form_2, method="hk", 
                              get.ests=T)
        }else{
	        rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl, formula=qtl_form_1, 
                              method="hk", model="normal", verbose=F)
            fit_f <- fitqtl(cross, pheno.col=pheno.col, qtl=qtl, formula=qtl_form_1,
                              method="hk", get.ests=T)
		}
    }else{
        rqtl <- refineqtl(cross, pheno.col=pheno.col, qtl=qtl, formula=qtl_form_1, 
                              method="hk", model="normal", verbose=F)
        fit_f <- fitqtl(cross, pheno.col=pheno.col, qtl=qtl, formula=qtl_form_1,
                              method="hk", get.ests=T)
    }
	left_ran <- c()
	right_ran <- c()
	pos <- c()
	genet_left <- c()
	Q <- c()
	start_cM <- c(); end_cM <- c(); peak_cM <- c()
	for(i in 1:rqtl$n.qtl){
	    interval<-lodint(rqtl, qtl.index=i, drop=1.5, expandtomarkers=TRUE)
        mar <- row.names(interval)
        left_ran <- c(left_ran, unlist(strsplit(mar[1], split="_"))[2])
        right_ran <- c(right_ran, unlist(strsplit(mar[3], split="_"))[3])
		pos <- c(pos, paste(unlist(strsplit(mar[2], split="_"))[2:3], collapse="_"))
		Q <- c(Q, paste("Q",i,sep=""))
		start_cM <- c(start_cM, interval[1,2]); end_cM <- c(end_cM, interval[3,2])
        peak_cM <- c(peak_cM, interval[2,2])
	}
    left_ran <- as.numeric(left_ran)
    right_ran <- as.numeric(right_ran)
    int_size <- right_ran - left_ran	
	fit <- summary(fit_f)
    add <- fit$ests[seq(2,2*rqtl$n.qtl+1,2)]
	dom <- fit$ests[seq(3,2*rqtl$n.qtl+1,2)]
	add <- round(add, 4)
	dom <- round(dom, 4)
	start_cM <- round(start_cM, 2)
    end_cM <- round(end_cM, 2)
    peak_cM <- round(peak_cM, 2)
	if(qtl$n.qtl == 1){
	    LOD <- fit$result.full[1,4]
        PVE <- fit$result.full[1,5]
    }else{
        PVE <- c(t(fit$result.drop[1:rqtl$n.qtl,4]))
	    LOD <- c(t(fit$result.drop[1:rqtl$n.qtl,3]))
    } 
    LOD <- round(LOD, 2)
    PVE <- round(PVE, 4)	
	qtl_int <- data.frame(Q=Q, chr=rqtl$chr, start_cM=start_cM, end_cM=end_cM, gent_peak=peak_cM, start=left_ran, end=right_ran, peak_pos=pos,  
	                      int_size=int_size, LOD=LOD, PVE=PVE, ADD=add, DOM=dom)
	###information of epstasis
	if(length(int_id)>0){
	    row_epi <- qtl$n.qtl+length(int_id)
		qtl_stat <- as.data.frame(fit_f$result.drop[1:row_epi,])
	    epi <- qtl_stat[(qtl$n.qtl+1):row_epi,]
		#epi %>% round(2)
		epi$LOD <- round(epi$LOD,2); epi$"%var" <- round(epi$"%var",4); epi$"F value" <- round(epi$"F value",2)
	}else{epi <- "none epistasis detected"}
	print(fit)
	qtl_info <- list()
	qtl_info[["qtl interval info"]] <- qtl_int
	qtl_info[["qtl epistasis"]] <- epi
	return(qtl_info)
}

phenotype <- colnames(cross$pheno)
cim_qtl_all_lt <- list()
cim_qtl_all_epi <- list()
for(phe in phenotype[17:20]){
  print(paste("Doing auto cim for", phe))
  phe_col <- which(phenotype==phe)
  #th <- summary(perms_lt[[phe]], alpha=0.05)[1,1]
  if(phe_col <= 16){perms <- subset(cross.sc1.perms_qt, lodcolumn=phe); model="normal"
    }else{perms <- subset(cross.sc1.perms_ql,lodcolumn=phe); model="binary"}
  th <- summary(perms, alpha=0.05)[1,1]
  cim_auto <- cim.auto.scan(cross, pheno.col=phe_col, threshold=th)
  cim_qtl_int <- cim.qtl.int(cross=cross, pheno.col=phe_col, cim=cim_auto, threshold=th)
  cim_qtl_all_lt[[phe]] <- cim_qtl_int[["qtl interval info"]]
  cim_qtl_all_epi[[phe]] <- cim_qtl_int[["qtl epistasis"]]
}

for(j in names(cim_qtl_all_lt)){
  cat(j, file="cim_qtl_all_traits.tsv", sep="\n", append=TRUE)
  if(which(phenotype==phe)==1){colnames=TRUE}else{colnames=FALSE}
  write.table(cim_qtl_all_lt[[j]], file="cim_qtl_all_traits.tsv", sep="\t", append=TRUE, row.names=TRUE, col.names=colnames, quote=FALSE)
}

for(k in names(cim_qtl_all_epi)){
  cat(k, file="qtl_epi_all_traits.tsv", sep="\n", append=TRUE)
  if(which(phenotype==phe)==1){colnames=TRUE}else{colnames=FALSE}
  write.table(cim_qtl_all_epi[[k]], file="qtl_epi_all_traits.tsv", sep="\t", append=TRUE, row.names=TRUE, col.names=colnames, quote=FALSE)
}

save.image("F2_cim_qtl.RData")
 
####2 part model
PS_2p <- scanone(cross, pheno.col=11, model="2part")
summary(PS_2p, perms=cross.sc1.perms_2p[,1:3], alpha=0.05, pvaluse=TRUE, format="allpeaks")

PE_2p <- scanone(cross, pheno.col=12, model="2part")
summary(PE_2p, perms=cross.sc1.perms_2p[,4:6], alpha=0.05, pvaluse=TRUE, format="allpeaks")
summary(PE_2p, perms=cross.sc1.perms_2p[,4:6], alpha=0.05, pvaluse=TRUE)
