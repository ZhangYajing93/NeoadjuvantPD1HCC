###################################

## use VDJtools to pool samples

####################
library('tzdb','backports','vroom','immunarch', lib.loc = '/home/yajingzh/R/x86_64-pc-linux-gnu-library/4.3')

poolfiles <- list.files('/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/VDJtools/PoolSamples/group12', pattern = '.pool.aa.table.txt')
poolID <- gsub(pattern = '.pool.aa.table.txt', replacement = '', poolfiles)
poolfiles <- paste0('/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/VDJtools/PoolSamples/group12/', poolfiles)

metadata.group12 <- cbind(poolfiles, poolID)
metadata.group12$Sample <- gsub(pattern = '.txt', replacement = '', list.files('/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/VDJtools/PoolSamples/group12', pattern = '.pool.aa.table.txt'))

write.table(metadata.group12, file = '/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/VDJtools/PoolSamples/metadata.group12.txt', 
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

metadata.group12 <- read.table(file = '/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/VDJtools/PoolSamples/metadata.group12.txt',
                               header = T, sep = '\t')
tcr.pool.immunarch <- repLoad(c(metadata.group12$poolfiles,'/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/VDJtools/PoolSamples/metadata.group12.txt'))

saveRDS(tcr.pool.immunarch, file = '/pathowh01/disk1/ZYJ/Nivo_ICItherapy/TCRseq/immunarch/tcr.pool12.immunarch.rds')
####################



#################### local
## data load

library(ggplot2)
library(ggpubr)

tcr.immunarch <- readRDS(file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/tcr.immunarch.rds')
tcr.immunarch.pool12 <- readRDS(file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/tcr.pool12.immunarch.rds')
load('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/pat.inf.RData')

#### cumulative percentage unique CDR3 clones

###### 12 groups  ******************
cumFreUniqueClone <- lapply(tcr.immunarch.pool12$data, function(y){
  index1<-rep(1/(dim(y)[1]),dim(y)[1])
  index3<-cumsum(y$Proportion)
  index4<-cbind(cumUniqueClone=index1,totalCloneFre=index3)
  return(index4)
})



subFre1000 <- lapply(cumFreUniqueClone, function(y){
  index1<-unlist(lapply(as.list(1:999), function(z){
    y[round(nrow(y)/1000*z),"totalCloneFre"]
  }))
  index2<-y[nrow(y),"totalCloneFre"]
  index3<-c(index1,index2)
  return(index3)
})

subFre1000 <- c(0,subFre1000$pre.NR.pool.aa.table, 0,subFre1000$pre.R.pool.aa.table,
               0,subFre1000$post.NR.pool.aa.table, 0,subFre1000$post.R.pool.aa.table,
               0,subFre1000$FU.NR.pool.aa.table, 0,subFre1000$FU.R.pool.aa.table,
               0,subFre1000$B.NR.pool.aa.table, 0,subFre1000$B.R.pool.aa.table,
               0,subFre1000$T.NR.pool.aa.table, 0,subFre1000$T.R.pool.aa.table,
               0,subFre1000$NT.NR.pool.aa.table, 0,subFre1000$NT.R.pool.aa.table)
subFre1000 <- as.data.frame(subFre1000)

subFre1000$subPer <- rep(seq(0,1000,1), 12)
subFre1000$group <- rep(c('pre.NR', 'pre.R', 'post.NR', 'post.R', 'FU.NR', 'FU.R',
                         'B.NR', 'B.R', 'T.NR', 'T.R', 'NT.NR', 'NT.R'), each = 1001)
colnames(subFre1000) <- c('cumFre', 'subPer', 'group')


library(ggplot2)
p1<- ggplot(subFre1000,aes(x=subPer,y=cumFre,fill=group))+
  geom_line(aes(colour=as.factor(group)))+
  # scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),limits = c(0,100))+
  geom_vline(xintercept = 100,colour="grey30", linetype = 'dotdash', alpha =0.5)+
  geom_hline(yintercept = 0.6, colour = 'grey30', linetype = 'dotdash', alpha =0.5) +
  geom_abline(intercept = 0, slope = 1/1000, colour = 'grey30', linetype = 'dotdash', alpha =0.5)+
  # scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))+
  # scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits = c(0,1))+
  theme_bw()+
  theme(aspect.ratio = 0.66, panel.grid = element_line(linetype = 'dashed'), panel.grid.minor = element_blank(), legend.position = 'inside', legend.position.inside = c(0.9,0.35)) +
  xlab('Cumulative percentage of TCRB clonetypes') +
  ylab('Cumulative frequency of total TCRB clonetypes')
  # + scale_color_bmj()

ggsave(filename = "meanCumulativeFrequency.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ## use

########### cumulative frequency curve: group and each sample
################################### ################################### ################################### ###################################


################################### ###################################

cumFreUniqueCloneeach <- lapply(tcr.immunarch$data, function(y){
  index1<-rep(1/(dim(y)[1]),dim(y)[1])
  index3<-cumsum(y$Proportion)
  index4<-cbind(cumUniqueClone=index1,totalCloneFre=index3)
  return(index4)
})


subFre100each <- lapply(cumFreUniqueCloneeach, function(y){
  index1<-unlist(lapply(as.list(1:99), function(z){
    y[round(nrow(y)/100*z),"totalCloneFre"]
  }))
  index2<-y[nrow(y),"totalCloneFre"]
  index3<-c(index1,index2)
  return(index3)
})

subFre100each <- lapply(subFre100each, function(x){
  index <- c(0,x)
  names(index) <- seq(0,100,1)
  return(index)
})

subFre100each <- subFre100each[intersect(c(pat.inf$pre, pat.inf$post, pat.inf$FU,
               pat.inf$B, pat.inf$T, pat.inf$NT), names(subFre100each))]
subFre100each <- lapply(names(subFre100each), function(x){
  index <- as.data.frame(subFre100each[x])
  index$subPer <- rownames(index)
  index$sampleid <- rep(x, nrow(index))
  colnames(index) <- c('cumFre', 'subPer', 'sampleid')
  return(index)
})

subFre100each <- do.call(rbind, subFre100each)

subFre100each$group <- subFre100each$sampleid
subFre100each$group[subFre100each$group %in% pat.inf$pre] <- 'pre'
subFre100each$group[subFre100each$group %in% pat.inf$post] <- 'post'
subFre100each$group[subFre100each$group %in% pat.inf$FU] <- 'FU'
subFre100each$group[subFre100each$group %in% pat.inf$B] <- 'B'
subFre100each$group[subFre100each$group %in% pat.inf$T] <- 'T'
subFre100each$group[subFre100each$group %in% pat.inf$NT] <- 'NT'

subFre100each$group1 <- subFre100each$sampleid
subFre100each$group1[subFre100each$group1 %in% pat.inf$pre[1:12]] <- 'pre.NR'
subFre100each$group1[subFre100each$group1 %in% pat.inf$pre[13:19]] <- 'pre.R'
subFre100each$group1[subFre100each$group1 %in% pat.inf$post[1:12]] <- 'post.NR'
subFre100each$group1[subFre100each$group1 %in% pat.inf$post[13:19]] <- 'post.R'

subFre100each$group1[subFre100each$group1 %in% pat.inf$FU[1:12]] <- 'FU.NR'
subFre100each$group1[subFre100each$group1 %in% pat.inf$FU[13:19]] <- 'FU.R'

subFre100each$group1[subFre100each$group1 %in% pat.inf$B[1:12]] <- 'B.NR'
subFre100each$group1[subFre100each$group1 %in% pat.inf$B[13:19]] <- 'B.R'

subFre100each$group1[subFre100each$group1 %in% pat.inf$T[1:12]] <- 'T.NR'
subFre100each$group1[subFre100each$group1 %in% pat.inf$T[13:19]] <- 'T.R'
subFre100each$group1[subFre100each$group1 %in% pat.inf$NT[1:12]] <- 'NT.NR'
subFre100each$group1[subFre100each$group1 %in% pat.inf$NT[13:19]] <- 'NT.R'

subFre100each$response <- do.call(rbind, strsplit(subFre100each$group1,'[.]'))[,2]

subFre100each$type <- subFre100each$group
subFre100each$type[subFre100each$type %in% c('pre','post','FU')] <- 'PBMC'
subFre100each$type[subFre100each$type %in% c('B','T','NT')] <- 'Tis'

subFre100each$subPer <- as.numeric(subFre100each$subPer)

###### top1%
top1per <- subFre100each[subFre100each$subPer == '1',]

top1per %>% group_by(group1) %>% summarise(mean = mean(cumFre))
# # A tibble: 12 × 2
# group1   mean
# <chr>   <dbl>
#   1 B.NR    0.300
# 2 B.R     0.298
# 3 FU.NR   0.240
# 4 FU.R    0.250
# 5 NT.NR   0.287
# 6 NT.R    0.329
# 7 T.NR    0.276
# 8 T.R     0.468
# 9 post.NR 0.282
# 10 post.R  0.295
# 11 pre.NR  0.262
# 12 pre.R   0.320


wilcox.test(top1per$cumFre[top1per$group == 'pre'][1:12], top1per$cumFre[top1per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.3402
wilcox.test(top1per$cumFre[top1per$group == 'post'][1:12], top1per$cumFre[top1per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.7108
wilcox.test(top1per$cumFre[top1per$group == 'FU'][1:12], top1per$cumFre[top1per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.8369

wilcox.test(top1per$cumFre[top1per$group == 'B'][1:6], top1per$cumFre[top1per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 1
wilcox.test(top1per$cumFre[top1per$group == 'T'][1:12], top1per$cumFre[top1per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.01707 * 
wilcox.test(top1per$cumFre[top1per$group == 'NT'][1:12], top1per$cumFre[top1per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.7108

###### top2%
top2per <- subFre100each[subFre100each$subPer == '2',]

wilcox.test(top2per$cumFre[top2per$group == 'pre'][1:12], top2per$cumFre[top2per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.2614
wilcox.test(top2per$cumFre[top2per$group == 'post'][1:12], top2per$cumFre[top2per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.5918
wilcox.test(top2per$cumFre[top2per$group == 'FU'][1:12], top2per$cumFre[top2per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.8369

wilcox.test(top2per$cumFre[top2per$group == 'B'][1:6], top2per$cumFre[top2per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.7619
wilcox.test(top2per$cumFre[top2per$group == 'T'][1:12], top2per$cumFre[top2per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.01707 * 
wilcox.test(top2per$cumFre[top2per$group == 'NT'][1:12], top2per$cumFre[top2per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.7108


###### top1%-top2%
top1to2per <- subFre100each[subFre100each$subPer == '2','cumFre']-subFre100each[subFre100each$subPer == '1','cumFre']
top1to2per <- cbind(top1to2per, top1per[,3:7])
colnames(top1to2per) <- c('cumFre', colnames(top1per)[3:7])

wilcox.test(top1to2per$cumFre[top1to2per$group == 'pre'][1:12], top1to2per$cumFre[top1to2per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9671
wilcox.test(top1to2per$cumFre[top1to2per$group == 'post'][1:12], top1to2per$cumFre[top1to2per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.06835
wilcox.test(top1to2per$cumFre[top1to2per$group == 'FU'][1:12], top1to2per$cumFre[top1to2per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.7108

wilcox.test(top1to2per$cumFre[top1to2per$group == 'B'][1:6], top1to2per$cumFre[top1to2per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.7619
wilcox.test(top1to2per$cumFre[top1to2per$group == 'T'][1:12], top1to2per$cumFre[top1to2per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.6504
wilcox.test(top1to2per$cumFre[top1to2per$group == 'NT'][1:12], top1to2per$cumFre[top1to2per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9018


###### top5%
top5per <- subFre100each[subFre100each$subPer == '5',]

wilcox.test(top5per$cumFre[top5per$group == 'pre'][1:12], top5per$cumFre[top5per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.1673
wilcox.test(top5per$cumFre[top5per$group == 'post'][1:12], top5per$cumFre[top5per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.5358
wilcox.test(top5per$cumFre[top5per$group == 'FU'][1:12], top5per$cumFre[top5per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9018

wilcox.test(top5per$cumFre[top5per$group == 'B'][1:6], top5per$cumFre[top5per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.9143
wilcox.test(top5per$cumFre[top5per$group == 'T'][1:12], top5per$cumFre[top5per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.02834 * 
wilcox.test(top5per$cumFre[top5per$group == 'NT'][1:12], top5per$cumFre[top5per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9018


###### top2%-top5%
top2to5per <- subFre100each[subFre100each$subPer == '5','cumFre']-subFre100each[subFre100each$subPer == '2','cumFre']
top2to5per <- cbind(top2to5per, top1per[,3:7])
colnames(top2to5per) <- c('cumFre', colnames(top1per)[3:7])

wilcox.test(top2to5per$cumFre[top2to5per$group == 'pre'][1:12], top2to5per$cumFre[top2to5per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 1
wilcox.test(top2to5per$cumFre[top2to5per$group == 'post'][1:12], top2to5per$cumFre[top2to5per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.6504
wilcox.test(top2to5per$cumFre[top2to5per$group == 'FU'][1:12], top2to5per$cumFre[top2to5per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.2991

wilcox.test(top2to5per$cumFre[top2to5per$group == 'B'][1:6], top2to5per$cumFre[top2to5per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.9143
wilcox.test(top2to5per$cumFre[top2to5per$group == 'T'][1:12], top2to5per$cumFre[top2to5per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.1422
wilcox.test(top2to5per$cumFre[top2to5per$group == 'NT'][1:12], top2to5per$cumFre[top2to5per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.04493  * 


###### top10%
top10per <- subFre100each[subFre100each$subPer == '10',]

top10per %>% group_by(group1) %>% summarise(mean = mean(cumFre))
# # A tibble: 12 × 2
# group1   mean
# <chr>   <dbl>
#   1 B.NR    0.580
# 2 B.R     0.580
# 3 FU.NR   0.415
# 4 FU.R    0.430
# 5 NT.NR   0.562
# 6 NT.R    0.610
# 7 T.NR    0.584
# 8 T.R     0.719
# 9 post.NR 0.454
# 10 post.R  0.479
# 11 pre.NR  0.437
# 12 pre.R   0.493

wilcox.test(top10per$cumFre[top10per$group == 'pre'][1:12], top10per$cumFre[top10per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.2268
wilcox.test(top10per$cumFre[top10per$group == 'post'][1:12], top10per$cumFre[top10per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.5358
wilcox.test(top10per$cumFre[top10per$group == 'FU'][1:12], top10per$cumFre[top10per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9018

wilcox.test(top10per$cumFre[top10per$group == 'B'][1:6], top10per$cumFre[top10per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.7619
wilcox.test(top10per$cumFre[top10per$group == 'T'][1:12], top10per$cumFre[top10per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.03584 * 
wilcox.test(top10per$cumFre[top10per$group == 'NT'][1:12], top10per$cumFre[top10per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.7732

###### top5%-top10%
top5to10per <- subFre100each[subFre100each$subPer == '10','cumFre']-subFre100each[subFre100each$subPer == '5','cumFre']
top5to10per <- cbind(top5to10per, top1per[,3:7])
colnames(top5to10per) <- c('cumFre', colnames(top1per)[3:7])

wilcox.test(top5to10per$cumFre[top5to10per$group == 'pre'][1:12], top5to10per$cumFre[top5to10per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 1
wilcox.test(top5to10per$cumFre[top5to10per$group == 'post'][1:12], top5to10per$cumFre[top5to10per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9671
wilcox.test(top5to10per$cumFre[top5to10per$group == 'FU'][1:12], top5to10per$cumFre[top5to10per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.432

wilcox.test(top5to10per$cumFre[top5to10per$group == 'B'][1:6], top5to10per$cumFre[top5to10per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.7619
wilcox.test(top5to10per$cumFre[top5to10per$group == 'T'][1:12], top5to10per$cumFre[top5to10per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.03584 * 
wilcox.test(top5to10per$cumFre[top5to10per$group == 'NT'][1:12], top5to10per$cumFre[top5to10per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.432

###### >top10%
top10to100per <- subFre100each[subFre100each$subPer == '100','cumFre']-subFre100each[subFre100each$subPer == '10','cumFre']
top10to100per <- cbind(top10to100per, top1per[,3:7])
colnames(top10to100per) <- c('cumFre', colnames(top1per)[3:7])

wilcox.test(top10to100per$cumFre[top10to100per$group == 'pre'][1:12], top10to100per$cumFre[top10to100per$group == 'pre'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.2268
wilcox.test(top10to100per$cumFre[top10to100per$group == 'post'][1:12], top10to100per$cumFre[top10to100per$group == 'post'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.5358
wilcox.test(top10to100per$cumFre[top10to100per$group == 'FU'][1:12], top10to100per$cumFre[top10to100per$group == 'FU'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.9018

wilcox.test(top10to100per$cumFre[top10to100per$group == 'B'][1:6], top10to100per$cumFre[top10to100per$group == 'B'][7:10], alternative = 'two.sided', paired = FALSE) # p = 0.7619
wilcox.test(top10to100per$cumFre[top10to100per$group == 'T'][1:12], top10to100per$cumFre[top10to100per$group == 'T'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.03584 * 
wilcox.test(top10to100per$cumFre[top10to100per$group == 'NT'][1:12], top10to100per$cumFre[top10to100per$group == 'NT'][13:19], alternative = 'two.sided', paired = FALSE) # p = 0.7732


#####  barplot
topperplot <- cbind(top10to100per$cumFre, top5to10per$cumFre,top2to5per$cumFre,top1to2per$cumFre,  top1per$cumFre) %>% as.data.frame()
colnames(topperplot) <- c('>10% clonotypes', 'top 5-10% clonotypes','top 2-5% clonotypes','top 1-2% clonotypes','top 1% clonotypes')
rownames(topperplot) <- top1per$sampleid

topperplot.group <- cbind(c(subFre1000[subFre1000$subPer == '1000', 'cumFre'] - subFre1000[subFre1000$subPer == '100', 'cumFre']),
                          c(subFre1000[subFre1000$subPer == '100', 'cumFre'] - subFre1000[subFre1000$subPer == '50', 'cumFre']),
                          c(subFre1000[subFre1000$subPer == '50', 'cumFre'] - subFre1000[subFre1000$subPer == '20', 'cumFre']),
                          c(subFre1000[subFre1000$subPer == '20', 'cumFre'] - subFre1000[subFre1000$subPer == '10', 'cumFre']),
                          subFre1000[subFre1000$subPer == '10', 'cumFre'])
colnames(topperplot.group) <- c('>10% clonotypes', 'top 5-10% clonotypes','top 2-5% clonotypes','top 1-2% clonotypes','top 1% clonotypes')
rownames(topperplot.group) <- unique(subFre1000$group)

pdf(file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig/topperCumuFreBarplot.pdf',14,7)

barplot(t(topperplot.group),col = c('#0348A6', '#7AC5FF','#C6FDEC','#FFB433','#FF4B20'),
        legend=rownames(t(topperplot.group)),las=2, ylab = 'Cumulative frequency',
        cex.names = 0.8,args.legend = list(x=3,y=1,cex=0.6))

barplot(t(topperplot),col = c('#0348A6', '#7AC5FF','#C6FDEC','#FFB433','#FF4B20'),
        legend=rownames(t(topperplot)),las=2, ylab = 'Cumulative frequency',
        cex.names = 0.8,args.legend = list(x=25,y=1,cex=0.6))

dev.off()



#############################
## boxplot statistical analysis

## top1%
p.top1 <- ggboxplot(top1per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top1% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE,
  #                    comparisons = list(c(1,2), 
  #                                       c(3,4), c(5,6),c(7,9), c(7,11),
  #                                       c(7,8),c(8,10),c(8,12),c(9,10),c(11,12)
  #                    )) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))#  T.NR vs T.R (0.017), B.R vs T.R (0.073)

## top1-2%
p.top2 <- ggboxplot(top1to2per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top1-2% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE,
  #                    comparisons = list(c(1,2), 
  #                                       c(3,4), c(5,6),c(7,9), c(7,11),
  #                                       c(7,8),c(8,10),c(8,12),c(9,10),c(11,12)
  #                    )) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # T.NR vs T.R (0.017), B.R vs T.R (0.073)

## top2-5%
p.top3 <- ggboxplot(top2to5per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top2-5% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE,
  #                    comparisons = list(c(1,2), 
  #                                       c(3,4), c(5,6),c(7,9), c(7,11),
  #                                       c(7,8),c(8,10),c(8,12),c(9,10),c(11,12)
  #                    )) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) #  T.NR vs T.R (0.017), B.R vs T.R (0.073)

## top5-10%
p.top4 <- ggboxplot(top5to10per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top5-10% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE,
  #                    comparisons = list(c(1,2), 
  #                                       c(3,4), c(5,6),c(7,9), c(7,11),
  #                                       c(7,8),c(8,10),c(8,12),c(9,10),c(11,12)
  #                    )) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) #  T.NR vs T.R (0.017), B.R vs T.R (0.073)

## >10%
p.top5 <-ggboxplot(top10to100per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of >10% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE,
  #                    comparisons = list(c(1,2), 
  #                                       c(3,4), c(5,6),c(7,9), c(7,11),
  #                                       c(7,8),c(8,10),c(8,12),c(9,10),c(11,12)
  #                    ))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) #  T.NR vs T.R (0.017), B.R vs T.R (0.073)

ggarrange(p.top1,p.top2,p.top3, p.top4, p.top5,nrow = 2, ncol = 3)

ggsave(filename = "topperCumuFreBoxplot_edit.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=12,units="in")



p.top10 <- ggboxplot(top10per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top10% clonotypes', xlab = '',
                    add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE,
                     comparisons = list(c(1,3),c(2,4),  c(3,5),c(4,6),
                                         c(1,5),c(2,6), 
                                        c(9,11),
                                        c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE,
                     comparisons = list(c(1,2),c(3,4), c(5,6), c(7,8),c(9,10),c(11,12),
                                        c(7,9), c(8,10),c(7,11),
                                       c(8,12)
                     )) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))#  T.NR vs T.R (0.017), B.R vs T.R (0.073)

ggsave(filename = "top10perCumuFreBoxplot_ref.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")

ggboxplot(top10per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top10% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3),c(2,4),  c(3,5),c(4,6),
  #                                       c(1,5),c(2,6), 
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE,
  #                    comparisons = list(c(1,2),c(3,4), c(5,6), c(7,8),c(9,10),c(11,12),
  #                                       c(7,9), c(8,10),c(7,11),
  #                                       c(8,12)
  #                    )) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))#  T.NR vs T.R (0.017), B.R vs T.R (0.073)

ggsave(filename = "top10perCumuFreBoxplot_edit.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")



#######################

topT <- rbind(top1per[top1per$group == 'T',-2], top1to2per[top1to2per$group == 'T',], top2to5per[top2to5per$group =='T',],
              top5to10per[top5to10per$group =='T',], top10to100per[top10to100per$group == 'T',])
topT$group.top <- rep(c('top1%','top1-2%','top2-5%','top5-10%','>10%'), each=19)
topT$group.top.response <- paste0(topT$response, topT$group.top)

p.top6 <-ggboxplot(topT, x = 'group.top.response', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency', xlab = '',
                   add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),c(7,8), c(9,10)))+
  # scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6 ,0.8))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) #  T.NR vs T.R (0.017), B.R vs T.R (0.073)

ggsave(filename = "topperCumuFreBoxplot.T.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")

################################### ################################### ################################### ################################### 



#######################

p2<- ggboxplot(top1per, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency of Top10% clonetypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = list(c(9,10),
                                        c(8,10))) #  T.NR vs T.R (0.017), B.R vs T.R (0.073)

p1+p2
ggsave(filename = "CumulativeFrequency.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=6,units="in") ## use

#######################
###### CumulativeFrequencyCurveEach.pdf: cumulative frequency curves of each sample
ggplot(subFre100each,aes(x=subPer,y=cumFre,fill=sampleid))+
  geom_line(aes(colour=as.factor(group1)))+
  # scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),limits = c(0,100))+
  geom_vline(xintercept = 10,colour="grey30", linetype = 'dotdash', alpha =0.5)+
  geom_hline(yintercept = 0.6, colour = 'grey30', linetype = 'dotdash', alpha =0.5) +
  geom_abline(intercept = 0, slope = 1/100, colour = 'grey30', linetype = 'dotdash', alpha =0.5)+
  # scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0),limits = c(0,1))+
  # scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),limits = c(0,1))+
  theme_bw()+
  theme(aspect.ratio = 0.66, panel.grid = element_line(linetype = 'dashed'), panel.grid.minor = element_blank(), legend.position = 'inside', legend.position.inside = c(0.9,0.35)) +
  xlab('Cumulative percentage of TCRB clonetypes') +
  ylab('Cumulative frequency of total TCRB clonetypes')
# + scale_color_bmj()

ggsave(filename = "CumulativeFrequencyCurveEach.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")  ## use

###################################












