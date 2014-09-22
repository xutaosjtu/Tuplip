source("dataLoading.R")

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
    
    model = gee(m ~ time.point*Group #+ as.factor(SEX) + AGE 
                #+ BMI + Syst.0 + HDLMG + HBA1C 
                #+ BZN0 + INS0_neu
                ,id = ID, data = data.samples, 
                na.action=na.omit,  contrasts = "exchangeable")
    model.coef = coef(summary(model))
    pvalues = 2 * pnorm(abs(model.coef)[,5], lower.tail = FALSE)
    model.coef = cbind(model.coef, pvalues)
    
    model.coef = coef(summary(model))[2,]
    pvalues = 2 * pnorm(abs(model.coef)[5], lower.tail = FALSE)
    rst = rbind(rst, c(model.coef, pvalue = pvalues))
  }
  else rst = rbind(rst, NA)
}
rownames(rst) = metabolites

write.csv(rst, file = "unadjusted model.csv")

## Apply mixed effect model
require("lme4")
require("lmerTest")
require(nlme)
model = lmer(Asp ~  time.point*Group + (time.point|ID) , 
   data = data.samples)
# 
# model = lme(Asp ~ 1, 
#              data = data.samples,
#              random = ~time.point-1|Group)

data.samples$time.point = as.factor(data.samples$time.point)
data.samples$Group = factor(data.samples$Group, levels = c("NGT", "IGT"))

rst = NULL
for(m in metabolites){
  data.samples$m = data.samples[,m]
  if(m %in% valid_measures){
    
    model = lme(m ~ time.point*Group 
#                 + as.factor(SEX) + AGE 
#                 + BMI + Syst.0 + HDLMG + HBA1C 
#                 + BZN0 + INS0_neu 
                , na.action = na.omit,
                random = ~1|ID, 
                data = data.samples)
    
    tmp = summary(model)$tTable
    model.coef = c(tmp[2,], tmp[3,], tmp[nrow(tmp)-1,],tmp[nrow(tmp),])
    rst = rbind(rst, model.coef)
  }
  else rst = rbind(rst, NA)
}
rownames(rst) = metabolites

write.csv(rst, file = "full model_lme_categorized.csv")


## Association with insulin sensitivity
rst = NULL
for(m in valid_measures){

  metabo.3t = matrix(data.samples[,m],ncol = 3, byrow = T)
  colnames(metabo.3t) = c("T0", "T60", "T120")
  tmppheno = data.samples[which(data.samples$time.point==1),colnames(pheno)]
  tmpdata = data.frame(tmppheno, metabo.3t)
  
  model = lm(log(ISIMATS_neu) ~ T0 + T60 + T120 
              + as.factor(SEX) + AGE 
              + BMI + Syst.0 + HDLMG 
#              + HBA1C + BZN0 + INS0_neu 
             , na.action = na.omit, 
             data = tmpdata)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = valid_measures
write.csv(rst, file = "associations between T120 metabolite leve and insulin sensitivity_multivaraite model_2.csv")


model = lm(ISIMATS_neu ~ BZN0 + INS0_neu 
           + BZ60 + INS60_neu
           + BZ120 + INS120_neu
           , na.action = na.omit, 
           data = tmpdata)
plot(tmpdata$ISIMATS_neu, exp(predict(model1)))
plot(tmpdata$ISIMATS_neu, predict(model2))

pairs(tmpdata[, c("T0", "T60", "T120", "ISIMATS_neu")])

# plot of correlations
cor1 = cor(x = data.samples$ISIMATS_neu[which(data.samples$time.point==1)], 
    y = data.samples[which(data.samples$time.point==1),valid_measures[-c(156,157)]],
    use = "pair", method = "spearman")
cor2 = cor(x = data.samples$ISIMATS_neu[which(data.samples$time.point==2)], 
           y = data.samples[which(data.samples$time.point==2),valid_measures[-c(156,157)]],
           use = "pair", method = "spearman")
cor3 = cor(x = data.samples$ISIMATS_neu[which(data.samples$time.point==3)], 
           y = data.samples[which(data.samples$time.point==3),valid_measures[-c(156,157)]],
           use = "pair", method = "spearman")
cor.ISIMATS_metabo = t(rbind(cor1, cor2, cor3))

require(corrplot)
pdf("correlation between Insulin sensitivity and metabolties.pdf", width = 4, height = 10)
for( i in 1:8){
  if(i*20<nrow(cor.ISIMATS_metabo)){
    corrplot(cor.ISIMATS_metabo[(i-1)*20+1:20,], method = "number")
  }
  else corrplot(cor.ISIMATS_metabo[((i-1)*20+1):nrow(cor.ISIMATS_metabo),], method = "number")
}
dev.off()



