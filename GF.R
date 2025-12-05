#options(expressions = 500000)
library(gradientForest)
library(tidyverse)
library(readr)
library(data.table)
path<-"/home/lws2/workspace/pop/BM/lfmm/gf4/"

# Use fread to read the files
#Gwifl <- fread(paste0(path, 'modified_fixed_updated.csv'))
#Ewifl <- fread(paste0(path, 'out2.csv'))
Ewifl <- fread(paste0(path, 'out2.csv'), drop = 1)
Gwifl <- fread(paste0(path, 'output.csv'), drop = 1)


#Gwifl <- read_delim(paste0(path, 'modified_transposed_combined_population_freq_no_column2.csv'), delim = ",")
#Ewifl <- read_delim(paste0(path, 'out2.csv'), delim = ",")
#Gwifl<-read_csv(paste0(path,'modified_transposed_combined_population_freq_no_column2.csv'))

#Ewifl<-read_csv(paste0(path,'out2.csv'))

nSites<-26
lev<-floor(log2(nSites*0.368/2))
lev
wiflforest<-gradientForest(cbind(Ewifl,Gwifl),
predictor.vars=Ewifl %>% colnames(),
response.vars=Gwifl %>% colnames(),
ntree=500,
transform=NULL,
compact=T,
nbin=201,
maxLevel=lev,
corr.threshold=0.5,
trace=T)

# Set the output filename and path
pdf("/home/lws2/workspace/pop/BM/lfmm/gf4/Overall_Importance.pdf",
    width = 8,  # Width (inches, default is 7)
    height = 6  # Height (inches, default is 7)
)

# Plot the graph
plot(wiflforest, plot.type = "Overall.Importance")

# Close the graphics device (must be executed to save the file)
dev.off()