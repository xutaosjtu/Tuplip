setwd("../../../Dropbox/Diabetes/")

data = read.csv("Data/merge_data.csv")
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


## plot
require(grid)
require(ggplot2)

## 
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf("test.pdf", width = 10)
## position matrix of the plots
position = rbind(c(1,1), c(1,2), c(2,1), c(2,2))
for(i in 1:ncol(metabo.mean)){
  
  ## for every four figures great a new page
  if(i %% 4 == 1){
    grid.newpage()
    grid.rect(gp=gpar(fill="grey"))
    pushViewport(viewport(layout=grid.layout(2, 2)))
  }
  
  
  tmp = data.frame(mean = metabo.mean[,i], 
                   sd = metabo.sd[,i], 
                   time = rep(c(1:3),2), 
                   group = rep(c("IGT", "NGT"), each = 3)
                   )
  
  ##Set NA values to zero to make the plot run
  tmp[which(is.na(tmp), arr.ind = T)]=0
  
  pd <- position_dodge(.1) 
  p = ggplot(tmp, aes(x = time, y = mean, colour = group, group = group)) + 
    geom_errorbar(aes(ymin = mean-sd/sqrt(73), ymax = mean + sd/sqrt(73)), 
                  width = .1, position = pd) + 
    geom_line(position = pd) + geom_point(position = pd) +
    ggtitle(colnames(metabo.mean)[i])
  
  ## Push the plot to the right grid according to the position matrix
  if(i %% 4 == 0){
    print(p, vp = vplayout(x = position[4,1], y = position[4,2]))
  }
  else print(p, vp = vplayout(x = position[i %% 4,1], y = position[i %% 4,2])) 
}
dev.off()

# require("gplots")
# plotmeans(m~ time.point, data = data.samples, subset = data.samples$Group=="IGT", 
#           col = "red", lab = m)
# plotmeans(m~ time.point, data = data.samples, subset = data.samples$Group=="NGT", 
#           add = T, lty = 2, col = "darkgreen")

