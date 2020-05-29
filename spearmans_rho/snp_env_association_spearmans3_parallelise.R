### PSingh MARCH 2020 ####
### pooja.singh09@gmail.com ###
### CoAdapTree #####
### this script takes the unstitched scaffold and position ###
###### SNP ENV Association without correcting for population structure ####
###### : This script can be parallelised, so please set numCores!!!


#needed libraries

library(dplyr)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)
library(data.table)
library(doSNOW)
# set number of cores to use

numCores=10

cl <- makeCluster(numCores)
registerDoSNOW(cl) 

#set iterations (here it is the number of environments

iterations <- 19

# set progress bar


pb <- txtProgressBar(min = 1, max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


# read in AF table from GATK pipeline

data1 <- read.table("/data/projects/pool_seq/pangenome/JP_pangenome/JP_pooled/snpsANDindels/03_maf-p05_RD-recalculated/JP_pooled-varscan_all_bedfiles_SNP_maf_RD-recalculated.txt", header=T, sep="\t")


### select columns with FREQ and order the populations and conver FREQ % to decimal

data2 <- data1  %>% select(contains(".FREQ"))
colhead <- sub(".FREQ", "", colnames(data2))
colnames(data2) <- colhead
data3 <- data2[ , order(colnames(data2))]
data4 <- data3 %>% mutate_each(funs(as.numeric(gsub("%", "", ., fixed = TRUE))/100))

### prepare row header for final SNP table

header1 <- data1  %>% select(contains("unstitched_locus"))
header <- as.data.frame(str_replace(header1$unstitched_locus, ">", ""))
colnames(header) <- c("ID")


### paste header and FREQs and finalise


data5 <- cbind(header, data4)

snps <- data5[2:41]
rownames(snps) <- data5$ID



#### read in env data and order header the same as the SNPfile #####

env_var <- read.table("jp_std_env-19variables.txt", header=T)
env_var1 <- env_var[c(6:24)]
rownames(env_var1) <- env_var$our_id
env_var2 <- t(env_var1)
env_var3 <- env_var2[ , order(colnames(env_var2))]



###### Parallelise correlation loop through each SNP versus ENV and output file ####

scafpos <- as.matrix(rownames(snps))
envname <- as.matrix(rownames(env_var3))

args <- commandArgs(trailingOnly = TRUE) 

start1 <- 1
end1 <- nrow(snps)
end2 <- nrow(env_var3)


# Create class which holds multiple results for each loop iteration.
# Each loop iteration populates four properties: $result1 and $result2 and so on

multiResultClass <- function(result1=NULL,result2=NULL,result3=NULL,result4=NULL)
{
  me <- list(result1 = result1,result2 = result2, result3 = result3, result4 = result4)

  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}


# set counter


count <- 0

# set sys time

system.time(

# loop through envs and parallelise

output <- foreach (j = 1:end2, .options.snow=opts, .packages="foreach", .combine=rbind) %dopar% {

# loop through SNPs

foreach (i = 1:end1, .combine=rbind) %do% {
 count <- count + 1
 result <- multiResultClass()
 result$result1 <- scafpos[i,1]
 result$result2 <- envname[j,1]
 result$result3 <- cor.test(as.numeric(snps[i,]), as.numeric(env_var3[j,]), method = "spearman", exact=FALSE, use = "pairwise.complete.obs")$estimate
 result$result4 <- cor.test(as.numeric(snps[i,]), as.numeric(env_var3[j,]), method = "spearman", exact=FALSE, use = "pairwise.complete.obs")$p.value
 setTxtProgressBar(pb, i)
 return(result)
}
}
)

output1 <- data.table(output)
colnames(output1) <- c("snp", "env", "spearmansrho", "pvalue")
outname_p <- paste("snp_env_spearmans_rho_parallel",".txt",sep = "")
fwrite(output1, outname_p, sep="\t", col.names = T, row.names = F, quote = F)


## close progress bar and clean up cluster
close(pb)
stopImplicitCluster()
