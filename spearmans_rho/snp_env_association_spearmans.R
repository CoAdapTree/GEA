### PSingh OCT 2019 ####
### pooja.singh09@gmail.com ###
### CoAdapTree #####

memory.limit(size=8000)

###### SNP ENV Association without correcting for population structure ####

#needed libraries

library(dplyr)
library(tidyr)

# read in AF table from GATK pipeline

data1 <- read.table("JP_pooled-varscan_all_bedfiles_SNP_maf.txt", header=T, sep="\t")


### select columns with FREQ and order the populations and conver FREQ % to decimal

data2 <- data1  %>% select(contains(".FREQ"))
data3 <- data2[ , order(names(data2))]
data4 <-  data3 %>% mutate_each(funs(as.numeric(gsub("%", "", ., fixed = TRUE))/100))

### prepare row header for final SNP table

header <- unite(data1[1:2], ID, sep = "_")

### paste header and FREQs and finalise


data5 <- cbind(header, data4)

snps <- data5[2:41]
rownames(snps) <- data5$ID



#### read in env data and order header the same as the SNPfile #####

env_var <- read.table("JP_ENVIRONFILE_headeridx.txt", header=T)
env_var1 <- env_var[ , order(names(env_var))]



###### correlation loop through each SNP versus ENV and output file ####

scafpos <- as.matrix(rownames(snps))
envname <- as.matrix(rownames(env_var1))

args <- commandArgs(trailingOnly = TRUE) 

start1 <- 1
end1 <- nrow(snps)
results_out <- array (NA, c((nrow(snps) * nrow(env_var1)),4))

count <- 0

system.time(

# loop through SNP

for (i in 1:nrow (snps)){

# loop through env
for (j in 1:nrow (env_var1)){
 count <- count + 1
 results_out[count,4] <- cor.test(as.numeric(snps[i,]), as.numeric(env_var1[j,]), method = "spearman", exact=FALSE, use = "pairwise.complete.obs")$p.value
 results_out[count,3] <- cor.test(as.numeric(snps[i,]), as.numeric(env_var1[j,]), method = "spearman", exact=FALSE, use = "pairwise.complete.obs")$estimate
 results_out[count,2] <- envname[j,1]
 results_out[count,1] <- scafpos[i,1]
 
}
}
)

colnames(results_out) <- c("snp", "env", "spearmansrho", "pvalue")
outname_p <- paste("snp_env_spearmans_rho",".txt",sep = "")
write.table (results_out, outname_p, sep="\t", col.names = T, row.names = F, quote = F)
