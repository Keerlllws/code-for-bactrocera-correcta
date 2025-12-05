#options(expressions = 500000)
library(lfmm)
library(data.table)
library(qvalue)
library(RSpectra)
#Y<-read.table("out.raw")
Y <- fread("/home/data/lws2/pop/final/lfmm/haiba/output.raw",sep = " ")
#part1 <- fread("part_aa")
#part2 <- fread("part_ab")
#part3 <- fread("part_ac")
#colnames(part2) <- colnames(part1)
#colnames(part3) <- colnames(part1)
#Y <- rbind(part1, part2, part3)
#library(readr)
#Y <- read_delim("out.raw", delim = " ")
Y<-Y[,-c(1:6)]
#X<-read.csv("out.csv",header = T,stringsAsFactors = F)
X <- fread("bio5.txt") 
#X <- read_delim("out.csv"), delim = ","
library(tictoc)
K=4
mod.lfmm <- lfmm::lfmm_ridge(Y=Y, X=X,K=K)
pv<-lfmm::lfmm_test(Y=Y,X=X,lfmm=mod.lfmm,calibrate="gif")
pvalues<-pv$calibrated.pvalue
my_qvalue <- function(x) {
  q <- qvalue::qvalue(x)
  q <- q$qvalues
  return(q)
}


qvalues<-apply(pvalues,2,my_qvalue)

#write.csv(qvalues_mat, "qvalues_results.csv", row.names = FALSE,quote=F)
write.table(qvalues, file = "bio5_qvalues.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
