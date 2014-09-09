setwd("../../../Dropbox/Diabetes/")

data = read.csv("Data/merge_data..csv")
indx.ref= grep("Ref", data$Sample.Identification)

data.samples = data[-indx.ref, ]
metabolites = colnames(data)[(10:239)*2-1]

setNAs = function(i, dat){
  mapply(function(x, y){if(y=="<LOD"|is.na(y)) x=NA; return(x)},dat[,i*2-1], dat[,i*2])
}
data.samples[,metabolites] = sapply(10:239, setNAs, dat = data.samples)

pheno = read.csv("Data/sc.csv")
data.samples = merge(data.samples, pheno)

##Calculate the mean and sd of the metabolites
metabo.mean = sapply(data.samples[, metabolites], 
                     function(x){
                       tapply(x, INDEX = interaction(data.samples$time.point,data.samples$Group), mean, na.rm = T)
                     })
metabo.sd = sapply(data.samples[, metabolites], 
                   function(x){
                     tapply(x, INDEX = interaction(data.samples$time.point,data.samples$Group), sd, na.rm = T)
                   })

write.csv(rbind(metabo.mean, metabo.sd), file = "metabolite mean and sd.csv")
