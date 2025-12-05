library(vegan)
library(psych)
library(data.table)
library(qvalue)

# ----------- Custom Functions -----------
rdadapt <- function(rda, K) {
  zscores <- rda$CCA$v[, 1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
}

outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

# ----------- Set Paths and Read Data -----------
setwd("/home/lws2/workspace/pop/BM/lfmm/gf5/RDA/all/3")

gen <- fread("/home/lws2/workspace/pop/BM/lfmm/gf5/pop.RDA.geno", header = TRUE)
POP <- read.table("/home/lws2/workspace/pop/BM/lfmm/gf5/popID.txt", header = FALSE)
LOCI <- read.table("/home/lws2/workspace/pop/BM/lfmm/gf5/ID.txt", header = FALSE)
env <- read.csv("/home/lws2/workspace/pop/BM/lfmm/gf5/RDA/env.csv", header = TRUE)

rownames(gen) <- as.character(POP$V1)
colnames(gen) <- as.character(LOCI$V1)
env$individual <- as.character(env$individual)

# ----------- Select Environmental Variables -----------
pred <- env[, c("bio2", "bio3", "bio5", "bio13", "bio19")]
colnames(pred) <- c("BIO2", "BIO3", "BIO5", "BIO13", "BIO19")

# ----------- Main RDA Analysis -----------
pk.rda <- rda(gen ~ ., data = pred, scale = TRUE)
saveRDS(pk.rda, file = "pk_rda_model.rds")

load.rda <- scores(pk.rda, choices = 1:3, display = "species")
write.table(load.rda, file = "result/loaduniq_rda.txt", quote = FALSE)

# ----------- Screen Outliers for the First Three Axes (z = 3.5) -----------
cand1 <- outliers(load.rda[, 1], 3.5)
cand2 <- outliers(load.rda[, 2], 3.5)
cand3 <- outliers(load.rda[, 3], 3.5)

cand1 <- cbind.data.frame(rep(1, length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis", "snp", "loading")

write.table(cand1, file = "result/cand1.txt", quote = FALSE)
write.table(cand2, file = "result/cand2.txt", quote = FALSE)
write.table(cand3, file = "result/cand3.txt", quote = FALSE)

# ----------- Merge All Candidate Loci -----------
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# ----------- Calculate Correlation Between SNPs and Environmental Variables -----------
foo <- matrix(nrow = nrow(cand), ncol = ncol(pred))
colnames(foo) <- colnames(pred)
gen <- as.data.frame(gen)

for (i in 1:nrow(cand)) {
  nam <- cand[i, "snp"]
  snp.gen <- gen[, nam]
  foo[i, ] <- apply(pred, 2, function(x) cor(x, snp.gen))
}

cand <- cbind.data.frame(cand, foo)
cand <- cand[!duplicated(cand$snp), ]

# ----------- Add Predictor with Max Correlation and Correlation Coefficient -----------
for (i in 1:nrow(cand)) {
  bar <- cand[i, ]
  cand[i, "predictor"] <- names(which.max(abs(bar[4:8])))
  cand[i, "correlation"] <- max(abs(bar[4:8]))
}

write.table(cand, file = "result/cand.txt", quote = FALSE, row.names = FALSE)

# ----------- Plot Loading Distribution for the First Three RDA Axes -----------
pdf("RDA.pdf")
opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
hist(load.rda[, 1], main = "Loadings on RDA1")
hist(load.rda[, 2], main = "Loadings on RDA2")
hist(load.rda[, 3], main = "Loadings on RDA3")
par(opar)
dev.off()