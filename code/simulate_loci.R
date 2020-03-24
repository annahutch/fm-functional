# code to simulate loci for which to perform fine-mapping
# Locus, ld and annotation files generated for input to PAINTOR

# prerequisite: PAINTOR_V3.0 must be downloaded - see https://github.com/gkichaev/PAINTOR_V3.0#installation

# navigate to one directory above PAINTOR_V3.0 directory

library(simGWAS)
library(corrcoverage)
library(data.table)

tmp <- readRDS("medium.RDS")

h <- tmp$h

nloci <- 100

for(i in 1:nloci){
  
  start <- sample(c(1:513), 1) # choose random starting point to select 200 SNPs from
  
  h_cut <- h[,seq(start+1,start+200,1)]
  maf <- colMeans(h_cut)
  LD <- tmp$LD[seq(start+1,start+200,1), seq(start+1,start+200,1)]
  
  beta <- sample(c(log(1.05),log(1.1),log(1.2)), 1) # log OR
  
  # corrplot(LD, "color", "upper", tl.pos = "n")
  
  freq <- as.data.frame(h_cut+1)
  freq$Probability <- 1/nrow(freq)
  snps <- colnames(freq)[-ncol(freq)]
  NN <- sample(c(2000,5000,10000),1) # vary N0=N1
  
  cvtype <- sample(c('friendly','medium','lonely'),1)
  
  nfriends <- apply(LD^2>0.5,1,sum)
  iCV <- switch(cvtype,
                friendly=sample(which(nfriends>10 & pmin(maf,1-maf)>0.05), 1),
                medium=sample(which(nfriends<=10 & nfriends>2 & pmin(maf,1-maf)>0.05), 1),
                lonely=sample(which(nfriends<=2 & pmin(maf,1-maf)>0.05), 1))
  
  CV <- snps[iCV]
  
  varbeta <- Var.data.cc(maf, N=2*NN, 0.5) # variance of beta
  
  z0 = simulated_z_score(N0=NN, # number of controls
                         N1=NN, # number of cases
                         snps=snps,
                         W=CV, # causal variants, subset of snps
                         gamma.W=beta, # log odds ratios
                         freq=freq
  )
  
  locus_file <- data.frame("ZSCORE" = z0)
  
  write.table(locus_file, file = paste0("PAINTOR_V3.0/RunDirectory/Locus",i), row.names = FALSE, sep = " ", quote = FALSE)
  
  ###################
  
  colnames(LD) <- NULL
  
  write.table(LD, file = paste0("PAINTOR_V3.0/RunDirectory/Locus",i,".ld"), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
  
  #################
  
  annot <- rbinom(length(z0), 1, 0.05)
  
  annot_file <- data.frame("annot" = annot)
  
  write.table(annot_file, file = paste0("PAINTOR_V3.0/RunDirectory/Locus",i,".annotations"), row.names = FALSE, sep = " ", quote = FALSE)
  
  ####
  
  misc <- data.frame("iCV" = iCV, "CV_annot" = annot[iCV], "locus" = i, "NN" = NN, "OR" = exp(beta), "ld" = cvtype, "max_Z" = max(abs(z0)))
  
  saveRDS(misc, paste0(i,".RDS"))
}

files <- list.files(pattern="*.RDS",full=TRUE)

RES <- lapply(files,function(x) readRDS(x))
sims1 <- do.call("rbind",RES)

rownames(sims1) <- NULL

saveRDS(sims1, "summary_res.RDS")
