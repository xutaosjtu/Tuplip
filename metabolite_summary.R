source("dataLoading.R")

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


## Average value and sd of blood sugar("BZN0","BZ60","BZ120") and insulin ("INS0_neu" ,"INS60_neu","INS120_neu")
Mean = sapply(pheno[,c("INS0_neu" ,"INS60_neu","INS120_neu")], 
              function(x){
                tapply(x, INDEX = pheno$Group, mean, na.rm = T)
                })
Mean = c(Mean[1,],Mean[2,]) 
SD = sapply(pheno[,c("INS0_neu" ,"INS60_neu","INS120_neu")], 
            function(x){
              tapply(x, INDEX = pheno$Group, sd, na.rm = T)
            })
SD = c(SD[1,],SD[2,])
tmp = data.frame(mean = Mean, 
                 sd = SD, 
                 time = rep(c(1:3),2), 
                 group = rep(c("IGT", "NGT"), each = 3)
)

pdf("Blood sugar.pdf", width = 10)
grid.newpage()
grid.rect(gp=gpar(fill="grey"))
pushViewport(viewport(layout=grid.layout(2, 2)))
p = ggplot(tmp, aes(x = time, y = mean, colour = group, group = group)) + 
  geom_errorbar(aes(ymin = mean-sd/sqrt(73), ymax = mean + sd/sqrt(73)), 
                width = .1, position = pd) + 
  geom_line(position = pd) + geom_point(position = pd) +
  ggtitle("Blood sugar")
print(p)
dev.off()



## plot
require(grid)
require(ggplot2)

## 
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("test.pdf", width = 10)
## position matrix of the plots
position = rbind(c(1,1), c(1,2), c(2,1), c(2,2)); k = 1
for(i in which(colnames(metabo.mean) %in% valid_measures)){
  
  ## for every four figures great a new page
  if(k %% 4 == 1){
    grid.newpage()
    grid.rect(gp=gpar(fill="grey"))
    pushViewport(viewport(layout=grid.layout(2, 2)))
  }
  
    tmp = data.frame(mean = metabo.mean[,i], 
                     sd = metabo.sd[,i], 
                     time = rep(c(1:3),2), 
                     group = rep(c("IGT", "NGT"), each = 3)
    )
    
#     ##Set NA values to zero to make the plot run
#     tmp[which(is.na(tmp), arr.ind = T)]=0
    
    pd <- position_dodge(.1) 
    p = ggplot(tmp, aes(x = time, y = mean, colour = group, group = group)) + 
      geom_errorbar(aes(ymin = mean-sd/sqrt(73), ymax = mean + sd/sqrt(73)), 
                    width = .1, position = pd) + 
      geom_line(position = pd) + geom_point(position = pd) +
      ggtitle(colnames(metabo.mean)[i])
  
  
  ## Push the plot to the right grid according to the position matrix
  if(k %% 4 == 0){
    print(p, vp = vplayout(x = position[4,1], y = position[4,2]))
  }
  else print(p, vp = vplayout(x = position[k %% 4,1], y = position[k %% 4,2])) 

  k = k+1
}
dev.off()

# require("gplots")
# plotmeans(m~ time.point, data = data.samples, subset = data.samples$Group=="IGT", 
#           col = "red", lab = m)
# plotmeans(m~ time.point, data = data.samples, subset = data.samples$Group=="NGT", 
#           add = T, lty = 2, col = "darkgreen")

## Boxplot
pdf("boxplot.pdf", width = 10)
## position matrix of the plots
position = rbind(c(1,1), c(1,2), c(2,1), c(2,2))
for(i in 1:length(metabolites)){
  m = metabolites[i]
   
  ## for every four figures great a new page
  if(i %% 4 == 1){
    grid.newpage()
    grid.rect(gp=gpar(fill="grey"))
    pushViewport(viewport(layout=grid.layout(2, 2)))
  }
  
  pd <- position_dodge(.1) 
  p = ggplot(data.samples, aes(x = factor(time.point), y = data.samples[,m])) + geom_boxplot(aes(fill=Group)) + ggtitle(m) + xlab("Concentration (uM)") + ylab("Time Point")
  
  ## Push the plot to the right grid according to the position matrix
  if(i %% 4 == 0){
    print(p, vp = vplayout(x = position[4,1], y = position[4,2]))
  }
  else print(p, vp = vplayout(x = position[i %% 4,1], y = position[i %% 4,2])) 
}
dev.off()




pdf("Insulin sensitivity index associated with metabolites.pdf", width = 10, height = 4)
mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="time.point") { 
    value[value=="1"] <- "T0"
    value[value=="2"]   <- "T60"
    value[value=="3"]   <- "T120"
  }
  return(value)
}
p = ggplot(data.samples, aes(log(C2), log(ISIMATS_neu))) + geom_point() + facet_grid(.~time.point, labeller = mf_labeller) + stat_smooth(method = "lm") + ggtitle("C2") +xlab("Log metabolite concentration") + ylab("Log insulin sensitivity score")
print(p)
p = ggplot(data.samples, aes(log(Gly), log(ISIMATS_neu))) + geom_point() + facet_grid(.~time.point, labeller = mf_labeller) + stat_smooth(method = "lm") + ggtitle("Glycine") +xlab("Log metabolite concentration") + ylab("Log insulin sensitivity score")
print(p)
p = ggplot(data.samples, aes(log(lysoPC.a.C18.2), log(ISIMATS_neu))) + geom_point() + facet_grid(.~time.point, labeller = mf_labeller) + stat_smooth(method = "lm") + ggtitle("LPC (18:2)") +xlab("Log metabolite concentration") + ylab("Log insulin sensitivity score")
print(p)
dev.off()

