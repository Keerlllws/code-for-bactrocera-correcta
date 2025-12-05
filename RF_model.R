library(dplyr) #package used to manipulate data see 
library(tidyr) #package used to manipulate data see 
library(stringr) #package used to change data. This package can be used 
library(randomForest) #package required to run Random Forest Analysis see
library(rpart) #package to build a single decision tree
library(party) #package used to plotting the decisions tree
library(partykit) #2nd package used to plotting the decisions tree
library(SNPRelate) #to read in genotypes from VCF
library(vcfR) #to convert vcf to genlight
library(adegenet) #for dapc
library(ggplot2) #for plotting
library(foreach)
library(doParallel)
library(readr)
library(forcats)
library(caret)


#### set up dataset ####
setwd("")
snpgdsVCF2GDS("all.recode.vcf", "Bcglobal_all.gds", method="biallelic.only")
genofile = snpgdsOpen("Bcglobal_all.gds")
# get genotypes coded as 0,1,2
genotypes <- snpgdsGetGeno(genofile, snp.id = NULL, snpfirstdim=TRUE) 
# check how many times 0,1,2 and NA occur in snps
table(c(genotypes), exclude=NULL) 

groups <- read.table("group.txt", header = TRUE)

# make dataframe with all SNPs(2)
groups.vec <- groups$group[match(read.gdsn(index.gdsn(genofile, "sample.id")), groups$sample)]
sample.vec <- read.gdsn(index.gdsn(genofile, "sample.id"))
df <- data.frame(cbind(sample.vec,groups.vec,t(genotypes)))

# df_init <- cbind(read.gdsn(index.gdsn(genofile, "sample.id")),
# str_extract((read.gdsn(index.gdsn(genofile, "sample.id"))), "[^_]+")) #the region of each sample is pulled from the sample name before the underscore, adapt regex to your needs
# df <- data.frame(df_init, t(genotypes))
#(3)
colnames(df) <- c("ID","Region", 
                  paste0("Chr0",read.gdsn(index.gdsn(genofile, "snp.chromosome")),"_",read.gdsn(index.gdsn(genofile, "snp.position"))))
dim(df)
df[1:10,1:10] #sanity check


# let's create a test set first with a reduced amount of snps
nr.snps =15000 

testdf <- cbind(df$Region, df[,sample(3:dim(df)[2], nr.snps)]) #ID column is not necessary anymore
df <- testdf #for convenience
colnames(df)[1] <- "Region"
region.vec <- df$Region
# make everything a factor (may take very long if you include many snps)
df[,1:dim(df)[2]] <- lapply(df[,1:dim(df)[2]], as.factor)


#### Prepare random forest training and validation set ####
# create training and validation set
train_index <- createDataPartition(df$Region, p = 0.7, list = FALSE, times = 1)
df.train <- df[train_index,] 
dim(df.train)
df.validate<-df[-train_index,] 
dim(df.validate)
# grow a tree as example
tree <- rpart(Region~., data=df.train)
plot(as.party(tree))
tree$cptable
# We want to prune the tree so that we have the fewest number of splits and the lowest xerror
pruned_tree <- prune(tree, cp=tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"])
# plot pruned tree
par(mar=c(5,5,5,5))
plot(as.party(pruned_tree))

#---- METHOD: Random Forest Analysis using parallel computing using the caret package (Preferred)  ----
amp_up_models <- function(){
  library(parallel)
  library(doParallel)
  no_cores <- parallel::detectCores() - 1
  #Leave one core available for Operating system
  cluster <- makePSOCKcluster(no_cores)
  registerDoParallel(cluster)
  cat("Model amped and ready to go with:", no_cores, "cores. \n")
}
amp_up_models()

#Feature selection algorithm using random forest
subsets <- c(5,10,15,20,30,40,50,80,100,130,150,180,200) #how many snps do you want to use in each subset?
set.seed(123)
rfeCtrl <- rfeControl(functions = rfFuncs,
                      method = "cv",
                      verbose = FALSE,
                      number = 10,
                      repeats = 10,
                      p = 0.7,
                      allowParallel = TRUE)
rfeCtrl
rfProfile <- rfe(x = df.train[, !names(df.train) %in% "Region"],
                 y = df.train$Region, 
                 sizes = subsets,
                 rfeControl = rfeCtrl)
plot(rfProfile, type = c("g", "o"), xlim = subsets)
rfProfile

best.features <- predictors(rfProfile)
MeanDecreaseGini <- rfProfile$fit$importance[,dim(rfProfile$fit$importance)[2]]
gini <- cbind(best.features, MeanDecreaseGini) %>% data.frame()
colnames(gini) <- c("snp","MeanDecreaseGini")
sorted_gini <- gini[order(-as.numeric(gini$MeanDecreaseGini)),]
norm_gini <- sorted_gini %>% mutate(norm = as.numeric(.[,2])/max(as.numeric(.[,2])))
norm_gini$snp <- rownames(norm_gini)
dev.off()
plot(norm_gini$norm, ylab = "Importance")

diagn_snps <- norm_gini %>% mutate(chrom = str_replace(snp, "_\\d+$", ""),
                                   pos = str_replace(snp, ".*_(\\d+)$", "\\1")) %>% dplyr::select(c(chrom, pos))

dim(diagn_snps)
write.table(diagn_snps, "rf_purged_diagn.snps_410.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

#Run random forest using only extracted features
df_purged <- cbind(df[,1], df[,best.features])
colnames(df_purged)[1] <- "Region"
df_purged[,1:dim(df_purged)[2]] <- lapply(df_purged[,1:dim(df_purged)[2]], as.factor)

train.p_index <- createDataPartition(df_purged$Region, p = 0.7, list = FALSE, times = 1)
df.train.purged <- df_purged[train.p_index,] 
dim(df.train.purged)
df.validate.purged<-df_purged[-train.p_index,] 
dim(df.validate.purged)

rf_model_purged <- randomForest(x = df.train.purged[,2:ncol(df.train.purged)],
                                y = df.train.purged[,1], ntree = 10000, importance = TRUE, ntry = 2)
purged.pred <- predict(rf_model_purged, df.validate.purged) #classifies new cases
purged.perf <- table(df.validate.purged[,1], purged.pred, dnn=c("Actual", "Predicted"))
print(purged.perf)
table.value(purged.perf, col.lab=paste("cluster", 1:12), row.lab=sort(unique(df.validate.purged[,1])))

#save and load model
save(rf_model_purged, file = "test_rf_model")
load("test_rf_model")

library(dplyr) #package used to manipulate data see 
library(tidyr) #package used to manipulate data see 
library(stringr) #package used to change data. This package can be used 
library(randomForest) #package required to run Random Forest Analysis see
library(SNPRelate) #to read in genotypes from VCF
library(foreach)
library(doParallel)
library(readr)
library(forcats)
library(caret)
library("rqdatatable")

setwd("")
        snpgdsVCF2GDS("all.recode.vcf", "query.gds", method="biallelic.only")
querygeno = snpgdsOpen("query.gds")
# get genotypes coded as 0,1,2
genotypes <- snpgdsGetGeno(querygeno, snp.id = NULL, snpfirstdim=TRUE)
# check how many times 0,1,2 and NA occur in snps (Just to check that you don't have an excessive amount of NA)
table(c(genotypes), exclude=NULL)
# make dataframe with all SNPs
df_init <- cbind(read.gdsn(index.gdsn(querygeno, "sample.id")))
df <- data.frame(df_init, t(genotypes))
colnames_vec <- c(paste0(read.gdsn(index.gdsn(querygeno, "snp.chromosome")),"_",read.gdsn(index.gdsn(querygeno, "snp.position"))))

## Load the trained model ##
load("test_rf_model2")
rf_model_purged # this is the model name that appears in your R environment to the right
predictors(rf_model_purged)

# for VCF data the colnames need to be change like: 01 to Chr01
colnames_vec2 <- paste0("Chr",colnames_vec)
colnames(df) <- c("sample", colnames_vec2)

# fill the full diagn snp dataset with NA, then add 0,1 or 2 where we have that info
NA_df <- as.data.frame(matrix("NA", ncol = length(predictors(rf_model_purged)), nrow = nrow(df)))
NA_df <- cbind(df$sample, NA_df)
colnames(NA_df) <- c(colnames(df)[1],predictors(rf_model_purged))
Q_df <- natural_join(df, NA_df,
                     by = "sample",
                     jointype = "INNER")
Q_df <- Q_df[, colnames(NA_df)]

# make the prediction
par(mfrow = c(2,4)) # make sure the plotting area is large enough by dragging the borders
for (i in 1:8) {
  predict(rf_model_purged, as.factor(Q_df[i,-1]), type = "prob") %>% barplot(main = Q_df[i,1],cex.names = 1.1,font=2, col = "lightgrey", cex.axis = 1.326)
}


