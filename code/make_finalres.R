# code to make fullres.RDS file for figs

library(dplyr)
library(data.table)

iCVs <- readRDS("summary_res.RDS")
setwd("/PAINTOR_V3.0/OutDirectory")

files <- list.files(pattern="*.results",full=TRUE)

res_log <- lapply(files, function(ch) grep("Log", ch))
# which vectors contain "Log"
res_log <- sapply(res_log, function(x) length(x) > 0)

# remove those that arent proper result files
files_final <- files[!res_log]

matches <- regmatches(files_final, gregexpr("[[:digit:]]+", files_final))

locus <- as.numeric(unlist(lapply(matches, function(x) x[1])))

prop <- as.numeric(paste0(unlist(lapply(matches, function(x) x[2])),".",unlist(lapply(matches, function(x) x[3]))))

RES <- lapply(files_final, function(x) read.table(x, header = TRUE))

credset <- function(pp, thr) {
  
  cumpp = cumsum(pp[order(pp, decreasing = TRUE)])  # cum sums of ordered pps
  
  wh = which(cumpp > thr)[1]  # how many needed to exceed thr
  
  names(wh) = NULL
  
  return(wh)
}

PPs <- rep(NA, length(RES))
ranks <- rep(NA, length(RES))
ld <- rep(NA, length(RES))
OR <- rep(NA, length(RES))
NN <- rep(NA, length(RES))
max_Z <- rep(NA, length(RES))

missing <- which(!seq(1,3000,1) %in% iCVs$locus)

for(i in 1:length(RES)){
  if(locus[i] %in% missing){
    next
  } else {
    PPs[i] <- RES[[i]]$Posterior_Prob[iCVs[which(iCVs$locus==locus[i]),]$iCV]
    ranks[i] <- rank(-RES[[i]]$Posterior_Prob)[iCVs[which(iCVs$locus==locus[i]),]$iCV]
    ld[i] <- iCVs[which(iCVs$locus==locus[i]),]$ld
    OR[i] <- iCVs[which(iCVs$locus==locus[i]),]$OR
    NN[i] <- iCVs[which(iCVs$locus==locus[i]),]$NN
    max_Z[i] <- iCVs[which(iCVs$locus==locus[i]),]$max_Z
  }
} 

cs_size <- lapply(RES, function(x) credset(x$Posterior_Prob, thr = 0.95)) %>% unlist

res <- data.frame(locus, prop, PPs, ranks, ld, OR, NN, max_Z, cs_size)

# res_final <- res[-which(res$PPs==1 & res$ranks != 1),]

saveRDS(res, "fullres.RDS")
