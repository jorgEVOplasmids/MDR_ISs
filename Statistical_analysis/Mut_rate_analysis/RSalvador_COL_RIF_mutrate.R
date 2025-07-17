
# Set working directory containing S. Table 2
setwd("~/Downloads")

# Load packages
library(magicaxis)
library(rsalvador)
library(ggplot2)
library(scales)
library(tidyverse)
library(xlsx)
library(ggpubr)
library(car)
library(scales)

# Load Colistin resistant mutants data of all the strains
Mutants <- read.xlsx("supplementary_table_2.xlsx", sheetName = "mutants_COL")

# Or load Rifampicin resistant mutants data of all the strains
Mutants <- read.xlsx("supplementary_table_2.xlsx", sheetName = "mutants_COL")

# Load viables data of all the strains
Viables <- read.xlsx("supplementary_table_2.xlsx", sheetName = "mutants_COL")

# Get the strain of interest to analyze (change depending on the strain under analysis, for example, for KPN08: select(c("KPN08", "KPN08p", "KPN08pΔΔIS1")))

Mutants <- Mutants %>% select(c("KPN01", "KPN01p")) # pOXA-48- strain in first position
Viables <- Viables %>% select(c("KPN01", "KPN01p")) # pOXA-48- strain in first position

x=1 # WT column

mu=seq(1,ncol(Mutants),1)
mut_rat=seq(1,ncol(Viables),1)
int95=matrix(NA,nrow=ncol(Mutants),ncol=2)
int_Mut_rat95=matrix(NA,nrow=ncol(Mutants),ncol=2)
int84=matrix(NA,nrow=ncol(Mutants),ncol=2)
int_Mut_rat84=matrix(NA,nrow=ncol(Mutants),ncol=2)
pvalue=matrix(NA, nrow=ncol(Mutants),ncol=ncol(Mutants))
test_statistic=matrix(NA, nrow=ncol(Mutants),ncol=ncol(Mutants))


for (i in 1:ncol(Mutants)){
  
  mu[i]<-newton.LD(Mutants[,i][complete.cases(Mutants[,i])], show.iter = T)
  mut_rat[i]=mu[i]/(mean(Viables[,i][complete.cases(Viables[,i])]))
  int95[i,]<-confint.LD(Mutants[,i][complete.cases(Mutants[,i])])
  int_Mut_rat95[i,]<-int95[i,]/(mean(Viables[,i][complete.cases(Viables[,i])]))
  
  int84[i,]<-confint.LD(Mutants[,i][complete.cases(Mutants[,i])],alpha = .16)
  int_Mut_rat84[i,]<-int84[i,]/(mean(Viables[,i][complete.cases(Viables[,i])]))
  
}

pval<-matrix(ncol=2, nrow=35) # Swap ncol = 3 in the case of KPN08, in which we compare 3 genotypes
stat=seq(1,ncol(Mutants),1)
z=1
j=1

for(x in seq(1, length(Mutants), by=2)){
  for (h in 1:length(Mutants)) {
    
    print(paste0('WT=',colnames(Mutants[x])))
    print(paste0('h=',h))
    print(paste0('Z=',z))
    Mutants1<-Mutants[,x][complete.cases(Mutants[,x])]
    Viables1<-mean(Viables[,x][complete.cases(Viables[,x])])
    Mutants2<-Mutants[,h][complete.cases(Mutants[,h])]
    Viables2<-mean(Viables[,h][complete.cases(Viables[,h])])
    R=Viables2/Viables1
    
    lorr<- LRT.LD(Mutants1,Mutants2,R=R )
    
    pval[j,h]=lorr[2]
    stat[z]=lorr[1]
    z=z+1
  }
  j=j+1
}

colnames(pval)<-colnames(Mutants[1:ncol(Mutants)])
rownames(pval)<- colnames(Mutants[ seq(1, length(Mutants), by=2)])

Result<-data.frame(mu, mut_rat,int_Mut_rat95, int_Mut_rat84, colSums(!is.na(Mutants)),row.names = as.character(colnames(Mutants)))
colnames(Result)<-c("Mutations per Culture", "Mutation Rate","95% IC Lower Limit","95% IC Upper Limit","84% IC Lower Limit","84% IC Upper Limit",
                    "Number of replicates")
print(Result)

Result$Genotype<- rep(c("pOXA-48-", "pOXA-48+")) # Add pOXA-48ΔΔIS1 in the case of KPN08

### Save results in a xlsx file per strain, substitute strain with the specific strain under analysis

write.xlsx(Result, "mut_rate_KPN01.xlsx") # These results are summarized in S. Table 2 "Results" sheet, merged with all the strains
