# code to make various combinations of loci for input into PAINTOR 
# such that the probability that the CV has the annotation varies

x <- readRDS("summary_res.RDS")
setwd("/PAINTOR_V3.0")

props <- c(0.05,seq(0.1, 0.9, 0.1),0.95)

for(i in props){
  
  for(j in 1:50){
    
    locus1 <- x[sample(which(x$CV_annot==1),100*i),'locus']
    
    locus0 <- x[sample(which(x$CV_annot==0),100-(100*i)),'locus']
    
    write.table(c(paste0("Locus",locus1), paste0("Locus",locus0)), file = paste0("RunDirectory/input.file",i,"_",j), col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
  }
}
