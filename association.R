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

## Create a ID variable
data.samples$Sample.Identification = as.character(data.samples$Sample.Identification)
data.samples$ID = sapply(data.samples$Sample.Identification, 
       function(x) {
         strsplit(x, split = "_")[[1]][1]
       })
data.samples$ID = as.numeric(data.samples$ID)

## log transformation and standardization
transform = function(x){
  return(scale(log(x)))
}
data.samples[,metabolites] = sapply(data.samples[, metabolites], transform)

## remove measurements with 5sd
outlier = function(x){
  m = mean(x, na.rm = T)
  sd = sd(x, na.rm = T)
  if (is.na(m)|is.na(sd)){return(x)}
  else {x[which(x>m+5*sd | x<m-5*sd)]=NA; return(x)}
}
data.samples[,metabolites] = sapply(data.samples[,metabolites], outlier)
                                                        
                                                        
## Apply GEE model
data.samples$Group = factor(data.samples$Group, levels = c("NGT", "IGT"))

if(!require(gee)) install.packages("gee")
require(gee)

rst = NULL
for(m in metabolites){
  data.samples$m = data.samples[,m]
  if(sum(is.na(data.samples$m))<0.5*nrow(data.samples)){
    
    model = gee((as.numeric(Group)-1) ~ m #+ as.factor(SEX) + AGE 
                #+ BMI + Syst.0 + HDLMG + HBA1C 
                #+ BZN0 + INS0_neu
                ,id = ID, data = data.samples, 
                na.action=na.omit, family = binomial, contrasts = "exchangeable")
    
    model.coef = coef(summary(model))[2,]
    pvalues = 2 * pnorm(abs(model.coef)[5], lower.tail = FALSE)
    rst = rbind(rst, c(model.coef, pvalue = pvalues))
  }
  else rst = rbind(rst, NA)
}
rownames(rst) = metabolites

write.csv(rst, file = "unadjusted model.csv")

