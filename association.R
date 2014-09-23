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

data.samples = data.samples[order(data.samples$ID, data.samples$time.point),]

data.samples$BZ = NA
data.samples$BZ[(1:146)*3-2] = data.samples$BZN0[(1:146)*3-2]
data.samples$BZ[(1:146)*3-1] = data.samples$BZ60[(1:146)*3-1]
data.samples$BZ[(1:146)*3] = data.samples$BZ120[(1:146)*3]
                                                        
data.samples$INS = NA
data.samples$INS[(1:146)*3-2] = data.samples$INS0_neu[(1:146)*3-2]
data.samples$INS[(1:146)*3-1] = data.samples$INS60_neu[(1:146)*3-1]
data.samples$INS[(1:146)*3] = data.samples$INS120_neu[(1:146)*3]

## Apply GEE model for the comparison between two groups
data.samples$Group = factor(data.samples$Group, levels = c("NGT", "IGT"))

if(!require(gee)) install.packages("gee")
require(gee)

rst = NULL
for(m in metabolites){
  data.samples$m = data.samples[,m]
  if(sum(is.na(data.samples$m))<0.5*nrow(data.samples)){
    
    model = gee((as.numeric(Group)-1) ~ scale(BZ)
                + as.factor(SEX) + AGE 
#                + BMI + Syst.0 + HDLMG 
#                + HBA1C  +BZN0 + INS0_neu
                ,id = ID, data = data.samples, 
                na.action=na.omit, family = binomial, contrasts = "exchangeable")
    
#     model = gee(m ~ time.point*Group #+ as.factor(SEX) + AGE 
#                 #+ BMI + Syst.0 + HDLMG + HBA1C 
#                 #+ BZN0 + INS0_neu
#                 ,id = ID, data = data.samples, 
#                 na.action=na.omit,  contrasts = "exchangeable")
#     model.coef = coef(summary(model))
#     pvalues = 2 * pnorm(abs(model.coef)[,5], lower.tail = FALSE)
#     model.coef = cbind(model.coef, pvalues)
    
    model.coef = coef(summary(model))[2,]
    pvalues = 2 * pnorm(abs(model.coef)[5], lower.tail = FALSE)
    rst = rbind(rst, c(model.coef, pvalue = pvalues))
  }
  else rst = rbind(rst, NA)
}
rownames(rst) = metabolites

write.csv(rst, file = "multivaraite_model 1.csv")

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
  
  model = lm(log(ISIMATS_neu) ~ T0 
              + as.factor(SEX) + AGE 
              + BMI + Syst.0 + HDLMG 
#              + HBA1C + BZN0 + INS0_neu 
             , na.action = na.omit, 
             data = tmpdata)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = valid_measures
write.csv(rst, file = "associations between T0 metabolite level and insulin sensitivity_multivaraite model_2.csv")

#Gly + lysoPC.a.C18.2 + Ile + PC.aa.C32.1 + PC.aa.C38.3+
model = lm(log(ISIMATS_neu) ~ scale(log(INS))+
           as.factor(SEX) + AGE 
           + BMI + Syst.0 + HDLMG 
#            + HBA1C + BZN0 + INS0_neu 
           , na.action = na.omit, 
           subset = time.point==3,
           data = data.samples)
summary(model)
plot(model$model$'log(ISIMATS_neu)', predict(model))


# model = lm(ISIMATS_neu ~ BZN0 + INS0_neu 
#            + BZ60 + INS60_neu
#            + BZ120 + INS120_neu
#            , na.action = na.omit, 
#            data = tmpdata)
# plot(tmpdata$ISIMATS_neu, exp(predict(model1)))
# plot(tmpdata$ISIMATS_neu, predict(model2))
# 
# pairs(tmpdata[, c("T0", "T60", "T120", "ISIMATS_neu")])



## Correlations between metabolites and insulin sensitivity
cor.test.2 = function(x, y, ...){
  corr = cor.test(x, y, ...)
  return(c(corr$estimate, corr$p.value))
}

cor.ISIMATS_metabo = NULL
for(i in 1:3){
  cor.tmp = sapply(data.samples[which(data.samples$time.point==i),valid_measures], 
                 cor.test.2, data.samples$ISIMATS_neu[which(data.samples$time.point==1)], 
                 use = "pair", method = "spearman")
  cor.ISIMATS_metabo = cbind(cor.ISIMATS_metabo, t(cor.tmp))
}

# plot the correlations
require(corrplot)
pdf("correlation between Insulin sensitivity and metabolties.pdf", width = 4, height = 10)
for( i in 1:8){
  if(i*20<nrow(cor.ISIMATS_metabo)){
    corrplot(cor.ISIMATS_metabo[(i-1)*20+1:20,], method = "number")
  }
  else corrplot(cor.ISIMATS_metabo[((i-1)*20+1):nrow(cor.ISIMATS_metabo),], method = "number")
}
dev.off()


## correlation between H1, glucose and insulin sensitivity
# correlation between H1 and blood glucose
(cor_0=cor.test(~H1+BZN0, data = data.samples, subset = time.point==1, method = "spearman"))
(cor_60=cor.test(~H1+BZ60, data = data.samples, subset = time.point==2, method = "spearman"))
(cor_120=cor.test(~H1+BZ120, data = data.samples, subset = time.point==3, method = "spearman"))
par(mfrow = c(2,2))
plot(H1~BZN0, data = data.samples, subset = time.point==1)
legend("topleft", legend = paste("Correlation:", round(cor_0$estimate,3), ", p=", format(cor_0$p.value, scientific = T, digits = 3)))
plot(H1~BZ60, data = data.samples, subset = time.point==2)
legend("topleft", legend = paste("Correlation:", round(cor_60$estimate,3), "p=", format(cor_60$p.value, scientific = T, digits = 3)))
plot(H1~BZ120, data = data.samples, subset = time.point==3)
legend("topleft", legend = paste("Correlation:", round(cor_120$estimate,3), "p=", format(cor_120$p.value, scientific = T, digits = 3)))

# correlation between blood glucose and insulin sensitivity
(cor_0 = cor.test(~ISIMATS_neu+BZN0, data = pheno, subset = time.point==1, method = "spearman"))
(cor_60 = cor.test(~ISIMATS_neu+BZ60, data = pheno, subset = time.point==2, method = "spearman"))
(cor_120 = cor.test(~ISIMATS_neu+BZ120, data = pheno, subset = time.point==3, method = "spearman"))
par(mfrow = c(2,2))
plot(ISIMATS_neu~BZN0, data = data.samples, subset = time.point==1, log = "xy")
legend("topleft", legend = paste("Correlation:", round(cor_0$estimate,3), ", p=", format(cor_0$p.value, scientific = T, digits = 3)))
plot(ISIMATS_neu~BZ60, data = data.samples, subset = time.point==2, log = "xy")
legend("topleft", legend = paste("Correlation:", round(cor_60$estimate,3), "p=", format(cor_60$p.value, scientific = T, digits = 3)))
plot(ISIMATS_neu~BZ120, data = data.samples, subset = time.point==3, log = "xy")
legend("topleft", legend = paste("Correlation:", round(cor_120$estimate,3), "p=", format(cor_120$p.value, scientific = T, digits = 3)))

# correlation between H1 and insulin sensitivity
(cor_0 = cor.test(~ISIMATS_neu+H1, data = data.samples, subset = time.point==1, method = "spearman"))
(cor_60 = cor.test(~ISIMATS_neu+H1, data = data.samples, subset = time.point==2, method = "spearman"))
(cor_120 = cor.test(~ISIMATS_neu+H1, data = data.samples, subset = time.point==3, method = "spearman"))
par(mfrow = c(2,2))
plot(ISIMATS_neu~H1, data = data.samples, subset = time.point==1, log = "xy")
legend("topleft", legend = paste("Correlation:", round(cor_0$estimate,3), ", p=", format(cor_0$p.value, scientific = T, digits = 3)))
plot(ISIMATS_neu~H1, data = data.samples, subset = time.point==2, log = "xy")
legend("topleft", legend = paste("Correlation:", round(cor_60$estimate,3), "p=", format(cor_60$p.value, scientific = T, digits = 3)))
plot(ISIMATS_neu~H1, data = data.samples, subset = time.point==3, log = "xy")
legend("topleft", legend = paste("Correlation:", round(cor_120$estimate,3), "p=", format(cor_120$p.value, scientific = T, digits = 3)))



table(data.samples$Plate.Bar.Code, data.samples$time.point)
table(data.samples$OP, data.samples$time.point)
table(data.samples$KitBarcodeNr, data.samples$time.point)
table(data.samples$Plate.Note, data.samples$time.point)
table(data.samples$Injection.Number, data.samples$time.point)
table(data.samples$Measurement.Time, data.samples$time.point)
