setwd("D:\\2018-2022\\F2群体\\220110-遗传图谱分析\\rqtl\\230501_cim_final_results\\phenotype_distri")

f2_phe <- read.csv("230528_F2_pheno.csv", head=T)
f2_pheno$ANL <- f2_pheno$ANL*10

par_phe <- read.csv("parent_phe.csv", head=T)
par_phe$ANL <- par_phe$ANL*10
n_phe <- subset(par_phe, Parent=="N")
r_phe <- subset(par_phe, Parent=="R")

library(ggplot2)
library(cowplot)
library(patchwork)
library(svglite)

pheno_name <- c("ANL","FH","PE","AWL","CD","CH","CL","FLA","FLL","FLW","GL",
                "GWE","GWI","PL","PS","SN")
f2_pheno <- f2_phe[,pheno_name]

sz=12
lysz=0.5
#################
library(plotrix)
sd_list <- list()
for(p in pheno_name){  
  n.se <- sd(n_phe[,p], na.rm=T)
  r.se <- sd(r_phe[,p], na.rm=T)
  n.mean <- mean(n_phe[,p], na.rm=T)
  r.mean <- mean(r_phe[,p], na.rm=T)
  sd_list[[p]] <- round(data.frame(N_SD=n.se, N_AVG=n.mean, R_SD=r.se, R_AVG=r.mean),4)
}  

df <- sd_list$ANL
anl <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=ANL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=148, yend=162, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=155, yend=155, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=148, yend=162, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=155,yend=155, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$ANL), max(f2_pheno$ANL), length.out = 4),1)) + 
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FH
fh <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=FH), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=100, yend=110, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=105, yend=105, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=100, yend=110, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD,y=105, yend=105, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$FH), max(f2_pheno$FH), length.out = 4),0)) + 
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$PE
pe <-  ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=PE), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=124, yend=136, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=130, yend=130, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=124, yend=136, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=130, yend=130, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$PE), max(f2_pheno$PE), length.out = 4),0)) + 
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$AWL
awl <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=AWL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=73, yend=81, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=77, yend=77, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=73, yend=81, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=77, yend=77, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$AWL), max(f2_pheno$AWL), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$CD
cd <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=CD), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=145, yend=161, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=153, yend=153, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=167, yend=183, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=175, yend=175, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$CD), max(f2_pheno$CD), length.out = 4), 2)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$CH
ch <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=CH), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=153, yend=167, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=160, yend=160, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=153, yend=167,, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=160, yend=160, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$CH), max(f2_pheno$CH), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$CL
cl <-  ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=CL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=78, yend=86, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=82, yend=82, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=78, yend=86, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=82, yend=82, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$CL), max(f2_pheno$CL), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FLA
fla <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=FLA), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=89, yend=99, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=94, yend=94, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=89, yend=99, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=94, yend=94, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$FLA), max(f2_pheno$FLA), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FLL
fll <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=FLL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=90, yend=102, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=96, yend=96, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=104, yend=116, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=110, yend=110, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$FLL), max(f2_pheno$FLL), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$FLW
flw <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=FLW), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=176, yend=192, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=184, yend=184, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=176, yend=192, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=184, yend=184, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$FLW), max(f2_pheno$FLW), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$GL
gl <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=GL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=97, yend=107, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=102, yend=102, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=97, yend=107, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=102, yend=102, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$GL), max(f2_pheno$GL), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$GWE
gwe <-  ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=GWE), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=96, yend=108, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=102, yend=102, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=96, yend=108, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=102, yend=102, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$GWE,na.rm=T), max(f2_pheno$GWE,na.rm=T), length.out = 4), 1)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$GWI
gwi <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=GWI), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=96, yend=108, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=102, yend=102, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=96, yend=108, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=102, yend=102, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$GWI), max(f2_pheno$GWI), length.out = 4), 2)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$PL
pl <-  ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=PL), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=96, yend=108, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=102, yend=102, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=96, yend=108, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=102, yend=102, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$PL), max(f2_pheno$PL), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$PS
ps <-  ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=PS), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=280, yend=310, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=295, yend=295, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=280, yend=310, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=295, yend=295, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$PS), max(f2_pheno$PS), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))

df <- sd_list$SN
sn <- ggplot() + 
  geom_histogram(data=f2_pheno, aes(x=SN), bins=20, fill="gray70", colour="gray30") +
  xlab("") + ylab("") +
  annotate("segment", x=df$N_AVG, xend=df$N_AVG, y=78, yend=86, size=lysz) + 
  annotate("segment", x=df$N_AVG-df$N_SD, xend=df$N_AVG+df$N_SD, y=82, yend=82, size=lysz) + 
  annotate("segment", x=df$R_AVG, xend=df$R_AVG, y=92, yend=100, size=lysz) + 
  annotate("segment", x=df$R_AVG-df$R_SD, xend=df$R_AVG+df$R_SD, y=96, yend=96, size=lysz) +
  scale_x_continuous(breaks = round(seq(min(f2_pheno$SN), max(f2_pheno$SN), length.out = 4), 0)) +
  theme_classic() +
  theme(axis.text=element_text(size=sz, family="serif"))


#plot_grid(nrow=4, ncol=4, anl,fh,pe,awl,cd,ch,cl,fla,fll,flw,gl,gwe,gwi,pl,ps,sn,
         aling="v")

#ggsave("230609_F2_distri.svg", width=10, height=7)

(anl/fh/pe/awl/cd/ch/cl/fla/fll/flw/gl/gwe/gwi/pl/ps/sn) +
  plot_layout(ncol = 4, nrow=4,  heights = rep(1,16))
ggsave("230612_F2_distri.svg", width=10, height=7)

svglite("230612_F2_distri.svg", width=10, height=7)
(anl/fh/pe/awl/cd/ch/cl/fla/fll/flw/gl/gwe/gwi/pl/ps/sn) +
  plot_layout(ncol = 4, nrow=4,  heights = rep(1,16))
dev.off()




