## Read data file into R and do basic preprcessing including,
setwd("../../../Dropbox/Diabetes/")

## 1. load metabolite measurements
data = read.csv("Data/merge_data.csv")
indx.ref= grep("Ref", data$Sample.Identification)

## 2. remove reference samples
data.samples = data[-indx.ref, ]
metabolites = colnames(data)[(10:239)*2-1]

# # 3. remove below LOD measurements
# setNAs = function(i, dat){
#   mapply(function(x, y){if(y=="<LOD"|is.na(y)) x=NA; return(x)},dat[,i*2-1], dat[,i*2])
# }
# data.samples[,metabolites] = sapply(10:239, setNAs, dat = data.samples)

## 4. load phenotype data
pheno = read.csv("Data/sc.csv")
data.samples = merge(data.samples, pheno)

## 5. load the list valid metabolites
valid_measures = read.csv("Valid measures.csv",stringsAsFactors=F, header = F)
valid_measures = valid_measures[which(valid_measures[,2]),1]
