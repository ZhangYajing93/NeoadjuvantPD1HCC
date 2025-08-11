######## top 1% clonotypes and non-top 1% clonotypes

### top 1% clonotypes
tcr.immunarch.top1per <- tcr.immunarch
tcr.immunarch.top1per$data <- lapply(tcr.immunarch.top1per$data, function(x){
  tcrdata <- x[1:round(nrow(x)/100),]
  return(tcrdata)
})

tcr.immunarch.top1per$data <- lapply(names(tcr.immunarch.top1per$data), function(x){
  tcr.immunarch.top1per$data[[x]]$sample_id <- rep(x, nrow(tcr.immunarch.top1per$data[[x]]))
  return(tcr.immunarch.top1per$data[[x]])
})

names(tcr.immunarch.top1per$data) <- names(tcr.immunarch$data)


### non-top 1% clonotypes
tcr.immunarch.nontop1per <- tcr.immunarch
tcr.immunarch.nontop1per$data <- lapply(tcr.immunarch.nontop1per$data, function(x){
  tcrdata <- x[(round(nrow(x)/100)+1):nrow(x),]
  return(tcrdata)
})

tcr.immunarch.nontop1per$data <- lapply(names(tcr.immunarch.nontop1per$data), function(x){
  tcr.immunarch.nontop1per$data[[x]]$sample_id <- rep(x, nrow(tcr.immunarch.nontop1per$data[[x]]))
  return(tcr.immunarch.nontop1per$data[[x]])
})

names(tcr.immunarch.nontop1per$data) <- names(tcr.immunarch$data)

saveRDS(tcr.immunarch.top1per, file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/tcr.immunarch.top1per.rds')
saveRDS(tcr.immunarch.nontop1per, file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/tcr.immunarch.nontop1per.rds')


######## ######## ######## ######## 


######## intersect clonotypes among top1% clonotypes

## use
unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], function(x){length(unique(x$CDR3.aa))}))
# T809 T811 T853 T841 T781 T824 T850 T791 T833 T803 T869 T788 T789 T814 T865 T795 T831 T828 T784 
# 19   31   17   23   20   71    7   21    3    7    8   21   46   20   17   33   57   58   41 
sum(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], function(x){length(unique(x$CDR3.aa))}))) # 520

tcr.immunarch.top1per.group1 <- lapply(tcr.immunarch.top1per$data, function(x){
  x %>% group_by(CDR3.aa) %>% summarise(Clones=sum(Clones),Proportion=sum(Proportion), sample_id = unique(sample_id))
})

tcr.immunarch.top1per.group1 <- list(pre = do.call(rbind, tcr.immunarch.top1per.group1[pat.inf$pre]),
                                     post = do.call(rbind, tcr.immunarch.top1per.group1[pat.inf$post]),
                                     FU = do.call(rbind, tcr.immunarch.top1per.group1[pat.inf$FU]),
                                     B = do.call(rbind, tcr.immunarch.top1per.group1[intersect(pat.inf$B,names(tcr.immunarch.top1per.group1))]),
                                     T = do.call(rbind, tcr.immunarch.top1per.group1[pat.inf$T]),
                                     NT = do.call(rbind, tcr.immunarch.top1per.group1[pat.inf$NT]))

tcr.immunarch.top1per.group1.share1 <- lapply(tcr.immunarch.top1per.group1, function(x){
  index <- x[x$CDR3.aa %in% x$CDR3.aa[duplicated(x$CDR3.aa)],]
  return(index)
})
# pre post   FU    B    T   NT 
# 372   66   50   42   54   70 

share.CDR3aa1 <- lapply(tcr.immunarch.top1per.group1.share1, function(x){
  unique(x$CDR3.aa)
})
# pre post   FU    B    T   NT 
# 177   29   22   10   13   12 
write.table(share.CDR3aa1$T, file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig/share.CDR3aa.T.txt', sep = '\t', quote = F)
write.table(share.CDR3aa1$B, file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig/share.CDR3aa.B.txt', sep = '\t', quote = F)


share.CDR3aa.each <- lapply(tcr.immunarch.top1per$data[pat.inf$T], function(x){
  x[x$CDR3.aa %in% share.CDR3aa1$T,]
})
share.CDR3aa.each <- do.call(rbind, share.CDR3aa.each)


### non-responder: nR same T cell clonotype in all patients: 0
intersect(tcr.immunarch.top1per$data$T809$CDR3.aa, 
          intersect(tcr.immunarch.top1per$data$T811$CDR3.aa, 
                    intersect(tcr.immunarch.top1per$data$T853$CDR3.aa, 
                              intersect(tcr.immunarch.top1per$data$T841$CDR3.aa,
                                        intersect(tcr.immunarch.top1per$data$T781$CDR3.aa, 
                                                  intersect(tcr.immunarch.top1per$data$T824$CDR3.aa,
                                                            intersect(tcr.immunarch.top1per$data$T850$CDR3.aa,
                                                                      intersect(tcr.immunarch.top1per$data$T791$CDR3.aa,
                                                                                intersect(tcr.immunarch.top1per$data$T833$CDR3.aa,
                                                                                          intersect(tcr.immunarch.top1per$data$T803$CDR3.aa,
                                                                                                    intersect(tcr.immunarch.top1per$data$T869$CDR3.aa,
                                                                                                              tcr.immunarch.top1per$data$T788$CDR3.aa))))))))))) # 0
### responder: R same T cell clonotype in all patients: 0
intersect(tcr.immunarch.top1per$data$T789$CDR3.aa, 
          intersect(tcr.immunarch.top1per$data$T814$CDR3.aa, 
                    intersect(tcr.immunarch.top1per$data$T865$CDR3.aa, 
                              intersect(tcr.immunarch.top1per$data$T795$CDR3.aa,
                                        intersect(tcr.immunarch.top1per$data$T831$CDR3.aa, 
                                                  intersect(tcr.immunarch.top1per$data$T828$CDR3.aa,
                                                            tcr.immunarch.top1per$data$T784$CDR3.aa)))))) # 0





##################### top1% T clonotypes overlap

tcr.immunarch.top1per.T.share <- lapply(pat.inf$patientID, function(x){
  preT <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'pre']]]$CDR3.aa,
                    tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)
  postT <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'post']]]$CDR3.aa,
                    tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)
  FUT <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'FU']]]$CDR3.aa,
                    tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)
  BT <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'B']]]$CDR3.aa,
                    tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)
  NTT <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'NT']]]$CDR3.aa,
                    tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)
  shareT <- list(preT = preT, postT = postT, FUT = FUT,
                 BT = BT, NTT = NTT)
  return(shareT)
})

names(tcr.immunarch.top1per.T.share) <- pat.inf$patientID


tcr.immunarch.top1per.T.share.prop <- lapply(pat.inf$T, function(x){
  preT.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                         tcr.immunarch$data[[pat.inf[pat.inf$T == x,'pre']]]$CDR3.aa)
  preT <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% preT.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  postT.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                          tcr.immunarch$data[[pat.inf[pat.inf$T == x,'post']]]$CDR3.aa)
  postT <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% postT.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  FUT.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                          tcr.immunarch$data[[pat.inf[pat.inf$T == x,'FU']]]$CDR3.aa)
  FUT <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% FUT.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  BT.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$T == x,'B']]]$CDR3.aa)
  BT <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% BT.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  NTT.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$T == x,'NT']]]$CDR3.aa)
  NTT <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% NTT.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  

  shareT.prop <- c(preT = preT, postT = postT, FUT = FUT,
                 BT = BT, NTT = NTT)
  return(shareT.prop)
})


tcr.immunarch.top1per.T.share.prop <- do.call(rbind, tcr.immunarch.top1per.T.share.prop) %>% as.data.frame(.)
rownames(tcr.immunarch.top1per.T.share.prop) <- pat.inf$T

unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))



tcr.immunarch.nontop1per.T.share.prop <- lapply(pat.inf$T, function(x){
  preT.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                         tcr.immunarch$data[[pat.inf[pat.inf$T == x,'pre']]]$CDR3.aa)
  preT <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% preT.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  postT.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                          tcr.immunarch$data[[pat.inf[pat.inf$T == x,'post']]]$CDR3.aa)
  postT <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% postT.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  FUT.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$T == x,'FU']]]$CDR3.aa)
  FUT <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% FUT.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  BT.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                       tcr.immunarch$data[[pat.inf[pat.inf$T == x,'B']]]$CDR3.aa)
  BT <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% BT.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  NTT.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$T == x,'NT']]]$CDR3.aa)
  NTT <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% NTT.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  shareT.prop <- c(preT = preT, postT = postT, FUT = FUT,
                   BT = BT, NTT = NTT)
  return(shareT.prop)
})

tcr.immunarch.nontop1per.T.share.prop <- do.call(rbind, tcr.immunarch.nontop1per.T.share.prop) %>% as.data.frame(.)
rownames(tcr.immunarch.nontop1per.T.share.prop) <- pat.inf$patientID



top1per.nontop1per.test <- lapply(names(tcr.immunarch.nontop1per.T.share.prop), function(x){
  all <- wilcox.test(tcr.immunarch.top1per.T.share.prop[,x], tcr.immunarch.nontop1per.T.share.prop[,x], paired = TRUE)$p.value
  nR <- wilcox.test(tcr.immunarch.top1per.T.share.prop[1:12,x], tcr.immunarch.nontop1per.T.share.prop[1:12,x], paired = TRUE)$p.value
  R <- wilcox.test(tcr.immunarch.top1per.T.share.prop[13:19,x], tcr.immunarch.nontop1per.T.share.prop[13:19,x], paired = TRUE)$p.value
  index <- c(all, nR, R)
  return(index)
})

top1per.nontop1per.test <- do.call(rbind, top1per.nontop1per.test) %>% as.data.frame(.)
rownames(top1per.nontop1per.test) <- names(tcr.immunarch.nontop1per.T.share.prop)
colnames(top1per.nontop1per.test) <- c('all', 'nR', 'R')


wilcox.test(tcr.immunarch.top1per.T.share.prop$BT[tcr.immunarch.top1per.T.share.prop$BT !=0], 
            tcr.immunarch.nontop1per.T.share.prop$BT[tcr.immunarch.nontop1per.T.share.prop$BT != 0], paired = TRUE) # 0.001953 **
wilcox.test(tcr.immunarch.top1per.T.share.prop$BT[tcr.immunarch.top1per.T.share.prop$BT !=0][1:6], 
            tcr.immunarch.nontop1per.T.share.prop$BT[tcr.immunarch.nontop1per.T.share.prop$BT != 0][1:6], paired = TRUE) # 0.03125 *
wilcox.test(tcr.immunarch.top1per.T.share.prop$BT[tcr.immunarch.top1per.T.share.prop$BT !=0][7:10], 
            tcr.immunarch.nontop1per.T.share.prop$BT[tcr.immunarch.nontop1per.T.share.prop$BT != 0][7:10], paired = TRUE) # 0.125

top1per.nontop1per.test['BT',] <- c(0.001953, 0.03125, 0.125)

#                all           nR        R
# preT  7.896423e-04 0.0268554688 0.031250
# postT 1.907349e-05 0.0024414063 0.015625
# FUT   2.098083e-04 0.0209960938 0.015625
# BT    1.953000e-03 0.0312500000 0.125000
# NTT   3.814697e-06 0.0004882813 0.015625


top1per.nontop1per.test.plot <- c(tcr.immunarch.top1per.T.share.prop$preT, tcr.immunarch.top1per.T.share.prop$postT, tcr.immunarch.top1per.T.share.prop$FUT,
                                  tcr.immunarch.top1per.T.share.prop$BT[tcr.immunarch.top1per.T.share.prop$BT !=0], tcr.immunarch.top1per.T.share.prop$NTT,
                                  tcr.immunarch.nontop1per.T.share.prop$preT, tcr.immunarch.nontop1per.T.share.prop$postT, tcr.immunarch.nontop1per.T.share.prop$FUT,
                                  tcr.immunarch.nontop1per.T.share.prop$BT[tcr.immunarch.nontop1per.T.share.prop$BT !=0], tcr.immunarch.nontop1per.T.share.prop$NTT)
top1per.nontop1per.test.plot <- cbind(prop = top1per.nontop1per.test.plot, group1 = c(rep(c('pre.top1','post.top1', 'FU.top1'), each=19), rep('B.top1',10), 
                                                                                      rep(c('NT.top1','pre.nontop1','post.nontop1','FU.nontop1'),each=19),
                                                                                      rep('B.nontop1',10),rep('NT.nontop1',19)),
                                      response = c(rep(c(rep('NR',12),rep('R',7)),3),rep('NR',6),rep('R',4),
                                                   rep(c(rep('NR',12),rep('R',7)),4),rep('NR',6),rep('R',4),rep('NR',12),rep('R',7))) %>% as.data.frame(.)

top1per.nontop1per.test.plot$prop <- as.numeric(top1per.nontop1per.test.plot$prop)
top1per.nontop1per.test.plot$group1 <- factor(top1per.nontop1per.test.plot$group1, 
                                              levels = c('pre.top1', 'pre.nontop1', 'post.top1','post.nontop1',
                                                         'FU.top1', 'FU.nontop1','B.top1','B.nontop1','NT.top1','NT.nontop1'))

p.top1nontop1 <- ggboxplot(top1per.nontop1per.test.plot, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion of shared clonotypes', xlab = '',
                    add.params = list(size=0.5)) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8),c(9,10))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # *

ggsave(filename = "top1perT.nontop1perT.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use



######### R vs nR  use
wilcox.test(tcr.immunarch.top1per.T.share.prop$preT[1:12], tcr.immunarch.top1per.T.share.prop$preT[13:19], paired = FALSE) # 0.6419
wilcox.test(tcr.immunarch.top1per.T.share.prop$postT[1:12], tcr.immunarch.top1per.T.share.prop$postT[13:19], paired = FALSE) # 0.02834
wilcox.test(tcr.immunarch.top1per.T.share.prop$FUT[1:12], tcr.immunarch.top1per.T.share.prop$FUT[13:19], paired = FALSE) # 0.03108
wilcox.test(tcr.immunarch.top1per.T.share.prop$BT[c(4:6,8,10,12)], tcr.immunarch.top1per.T.share.prop$BT[c(13,14,16,19)], paired = FALSE) # 0.9143 use 
wilcox.test(tcr.immunarch.top1per.T.share.prop$NTT[1:12], tcr.immunarch.top1per.T.share.prop$NTT[13:19], paired = FALSE) # 0.06286

######### intra-group compare: preICI.PBMC vs postICI.PBMC
wilcox.test(tcr.immunarch.top1per.T.share.prop$preT[1:12], tcr.immunarch.top1per.T.share.prop$postT[1:12], paired = TRUE) # 0.07593
wilcox.test(tcr.immunarch.top1per.T.share.prop$preT[13:19], tcr.immunarch.top1per.T.share.prop$postT[13:19], paired = TRUE) # 0.01563,  * use




######### R vs nR  use
wilcox.test(tcr.immunarch.nontop1per.T.share.prop$preT[1:12], tcr.immunarch.nontop1per.T.share.prop$preT[13:19], paired = FALSE) # 0.8369
wilcox.test(tcr.immunarch.nontop1per.T.share.prop$postT[1:12], tcr.immunarch.nontop1per.T.share.prop$postT[13:19], paired = FALSE) # 0.1422
wilcox.test(tcr.immunarch.nontop1per.T.share.prop$FUT[1:12], tcr.immunarch.nontop1per.T.share.prop$FUT[13:19], paired = FALSE) # 0.1673
wilcox.test(tcr.immunarch.nontop1per.T.share.prop$BT[c(4:6,8,10,12)], tcr.immunarch.nontop1per.T.share.prop$BT[c(13,14,16,19)], paired = FALSE) # 0.7619
wilcox.test(tcr.immunarch.nontop1per.T.share.prop$NTT[1:12], tcr.immunarch.nontop1per.T.share.prop$NTT[13:19], paired = FALSE) # 0.8369



tcr.immunarch.top1per.B.share.prop <- lapply(intersect(names(tcr.immunarch$data), pat.inf$B), function(x){

  preB.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                         tcr.immunarch$data[[pat.inf[pat.inf$B == x,'pre']]]$CDR3.aa)
  preB <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% preB.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  postB.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                          tcr.immunarch$data[[pat.inf[pat.inf$B == x,'post']]]$CDR3.aa)
  postB <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% postB.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  FUB.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$B == x,'FU']]]$CDR3.aa)
  FUB <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% FUB.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  TB.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                       tcr.immunarch$data[[pat.inf[pat.inf$B == x,'T']]]$CDR3.aa)
  TB <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% TB.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  NTB.comm <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$B == x,'NT']]]$CDR3.aa)
  NTB <- nrow(tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% NTB.comm, ])/nrow(tcr.immunarch.top1per$data[[x]])
  
  shareB.prop <- c(preB = preB, postB = postB, FUB = FUB,
                   TB = TB, NTB = NTB)
  return(shareB.prop)
})

tcr.immunarch.top1per.B.share.prop <- do.call(rbind, tcr.immunarch.top1per.B.share.prop) %>% as.data.frame(.)
rownames(tcr.immunarch.top1per.B.share.prop) <- intersect(names(tcr.immunarch$data), pat.inf$B)


######### R vs nR
wilcox.test(tcr.immunarch.top1per.B.share.prop$preB[1:6], tcr.immunarch.top1per.B.share.prop$preB[7:10], paired = FALSE) # 0.3524 ## use
wilcox.test(tcr.immunarch.top1per.B.share.prop$postB[1:6], tcr.immunarch.top1per.B.share.prop$postB[7:10], paired = FALSE) # 0.7619
wilcox.test(tcr.immunarch.top1per.B.share.prop$FUB[1:6], tcr.immunarch.top1per.B.share.prop$FUB[7:10], paired = FALSE) # 0.3359
wilcox.test(tcr.immunarch.top1per.B.share.prop$TB[1:6], tcr.immunarch.top1per.B.share.prop$TB[7:10], paired = FALSE) # 0.3524
wilcox.test(tcr.immunarch.top1per.B.share.prop$NTB[1:6], tcr.immunarch.top1per.B.share.prop$NTB[7:10], paired = FALSE) # 0.1143


####### B top 1% shared proportion with pre-treatment in NR vs R
####### T top 1% shared proportion with post-treatment and FU in NR vs R


plot <- cbind(prop = c(tcr.immunarch.top1per.B.share.prop$preB, tcr.immunarch.top1per.T.share.prop$preT, tcr.immunarch.top1per.T.share.prop$postT, 
                       tcr.immunarch.top1per.T.share.prop$FUT, tcr.immunarch.top1per.T.share.prop$BT[c(4:6,8,10,12:14,16,19)], tcr.immunarch.top1per.T.share.prop$NTT),
              group = c(rep('preB',10),rep('preT',19),rep('postT',19),rep('FUT',19),rep('BT',10),rep('NTT',19)),
              response = c(rep('NR',6),rep('R',4), rep(c(rep('NR',12),rep('R',7)), 3), rep('NR',6),rep('R',4),rep('NR',12),rep('R',7))) %>%as.data.frame(.)
plot$prop <- as.numeric(plot$prop)
plot$group1 <- paste0(plot$group,'.',plot$response)


ggboxplot(plot, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion of shared top 1% clonotypes', xlab = '',
                      add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10),c(11,12))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

ggsave(filename = "top1perTorB.RvsNR.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use



tcr.immunarch.top1per.B.share.prop.plot <- cbind(reshape2::melt(tcr.immunarch.top1per.B.share.prop), response = rep(c(rep('NR',6),rep('R',4)),5))
tcr.immunarch.top1per.B.share.prop.plot$group1 <- paste0(tcr.immunarch.top1per.B.share.prop.plot$variable,'.',tcr.immunarch.top1per.B.share.prop.plot$response)

ggboxplot(tcr.immunarch.top1per.B.share.prop.plot, x = 'group1', y = 'value', add = 'dotplot', ylab = 'Proportion of shared top 1% pre-intratumoral clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

ggsave(filename = "top1perB.RvsNR.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use


tcr.immunarch.nontop1per.B.share.prop.plot <- cbind(reshape2::melt(tcr.immunarch.nontop1per.B.share.prop), response = rep(c(rep('NR',6),rep('R',4)),5))
tcr.immunarch.nontop1per.B.share.prop.plot$group1 <- paste0(tcr.immunarch.nontop1per.B.share.prop.plot$variable,'.',tcr.immunarch.nontop1per.B.share.prop.plot$response)

ggboxplot(tcr.immunarch.nontop1per.B.share.prop.plot, x = 'group1', y = 'value', add = 'dotplot', ylab = 'Proportion of shared non-top 1% pre-intratumoral clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

ggsave(filename = "nontop1perB.RvsNR.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use






tcr.immunarch.nontop1per.B.share.prop <- lapply(intersect(names(tcr.immunarch$data), pat.inf$B), function(x){
  
  preB.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                         tcr.immunarch$data[[pat.inf[pat.inf$B == x,'pre']]]$CDR3.aa)
  preB <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% preB.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  postB.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                          tcr.immunarch$data[[pat.inf[pat.inf$B == x,'post']]]$CDR3.aa)
  postB <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% postB.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  FUB.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$B == x,'FU']]]$CDR3.aa)
  FUB <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% FUB.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  TB.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                       tcr.immunarch$data[[pat.inf[pat.inf$B == x,'T']]]$CDR3.aa)
  TB <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% TB.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  NTB.comm <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                        tcr.immunarch$data[[pat.inf[pat.inf$B == x,'NT']]]$CDR3.aa)
  NTB <- nrow(tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% NTB.comm, ])/nrow(tcr.immunarch.nontop1per$data[[x]])
  
  shareB.prop <- c(preB = preB, postB = postB, FUB = FUB,
                   TB = TB, NTB = NTB)
  return(shareB.prop)
})

tcr.immunarch.nontop1per.B.share.prop <- do.call(rbind, tcr.immunarch.nontop1per.B.share.prop) %>% as.data.frame(.)
rownames(tcr.immunarch.nontop1per.B.share.prop) <- intersect(names(tcr.immunarch$data), pat.inf$B)


top1per.nontop1per.B.test <- lapply(names(tcr.immunarch.nontop1per.B.share.prop), function(x){
  all <- wilcox.test(tcr.immunarch.top1per.B.share.prop[,x], tcr.immunarch.nontop1per.B.share.prop[,x], paired = TRUE)$p.value
  nR <- wilcox.test(tcr.immunarch.top1per.B.share.prop[1:6,x], tcr.immunarch.nontop1per.B.share.prop[1:6,x], paired = TRUE)$p.value
  R <- wilcox.test(tcr.immunarch.top1per.B.share.prop[7:10,x], tcr.immunarch.nontop1per.B.share.prop[7:10,x], paired = TRUE)$p.value
  index <- c(all, nR, R)
  return(index)
})

top1per.nontop1per.B.test <- do.call(rbind, top1per.nontop1per.B.test) %>% as.data.frame(.)
rownames(top1per.nontop1per.B.test) <- names(tcr.immunarch.nontop1per.B.share.prop)
colnames(top1per.nontop1per.B.test) <- c('all', 'nR', 'R')


#               all      nR     R
# preB  0.001953125 0.03125 0.125
# postB 0.001953125 0.03125 0.125
# FUB   0.001953125 0.03125 0.125
# TB    0.001953125 0.03125 0.125
# NTB   0.001953125 0.03125 0.125

plot <- cbind(prop = c(tcr.immunarch.nontop1per.B.share.prop$preB, tcr.immunarch.nontop1per.T.share.prop$preT, tcr.immunarch.nontop1per.T.share.prop$postT, 
                       tcr.immunarch.nontop1per.T.share.prop$FUT, tcr.immunarch.nontop1per.T.share.prop$BT[c(4:6,8,10,12:14,16,19)], tcr.immunarch.nontop1per.T.share.prop$NTT),
              group = c(rep('preB',10),rep('preT',19),rep('postT',19),rep('FUT',19),rep('BT',10),rep('NTT',19)),
              response = c(rep('NR',6),rep('R',4), rep(c(rep('NR',12),rep('R',7)), 3), rep('NR',6),rep('R',4),rep('NR',12),rep('R',7))) %>%as.data.frame(.)
plot$prop <- as.numeric(plot$prop)
plot$group1 <- paste0(plot$group,'.',plot$response)


ggboxplot(plot, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion of shared non-top 1% clonotypes', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10),c(11,12))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

ggsave(filename = "nontop1perTorB.RvsNR.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use sup

top1per.nontop1per.B.test.plot <- c(tcr.immunarch.top1per.B.share.prop$preB, tcr.immunarch.top1per.B.share.prop$postB, tcr.immunarch.top1per.B.share.prop$FUB,
                                  tcr.immunarch.top1per.B.share.prop$TB, tcr.immunarch.top1per.B.share.prop$NTB,
                                  tcr.immunarch.nontop1per.B.share.prop$preB, tcr.immunarch.nontop1per.B.share.prop$postB, tcr.immunarch.nontop1per.B.share.prop$FUB,
                                  tcr.immunarch.nontop1per.B.share.prop$TB, tcr.immunarch.nontop1per.B.share.prop$NTB)
top1per.nontop1per.B.test.plot <- cbind(prop = top1per.nontop1per.B.test.plot, group1 = c(rep(c('pre.top1','post.top1', 'FU.top1','T.top1','NT.top1',
                                                                                                'pre.nontop1','post.nontop1','FU.nontop1','T.nontop1','NT.nontop1'), each=10)),
                                      response = rep(c(rep('NR',6),rep('R',4)),10)) %>% as.data.frame(.)

top1per.nontop1per.B.test.plot$prop <- as.numeric(top1per.nontop1per.B.test.plot$prop)
top1per.nontop1per.B.test.plot$group1 <- factor(top1per.nontop1per.B.test.plot$group1, 
                                              levels = c('pre.top1', 'pre.nontop1', 'post.top1','post.nontop1',
                                                         'FU.top1', 'FU.nontop1','T.top1','T.nontop1','NT.top1','NT.nontop1'))

p.top1nontop1.B <- ggboxplot(top1per.nontop1per.B.test.plot, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion of shared top 1% clonotypes', xlab = '',
                           add.params = list(size=0.5)) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8),c(9,10))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # 

ggsave(filename = "top1perB.nontop1perB.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use





############# ITC top 1% clonotype cumulative frequency and tumor necrosis 

cor.test(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], function(x){sum(x$Proportion)})), pat.inf$necrosis, method = 'pearson') # cor 0.6591546, p 0.002143 use
cor.test(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$B[c(4:6,8,10,12:14,16,19)]], function(x){sum(x$Proportion)})), pat.inf$necrosis[c(4:6,8,10,12:14,16,19)], method = 'pearson') # cor -0.05399232 , p 0.8408 unuse
cor.test(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$NT], function(x){sum(x$Proportion)})), pat.inf$necrosis, method = 'pearson') # cor 0.3843953, p 0.1042 
cor.test(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$pre], function(x){sum(x$Proportion)})), pat.inf$necrosis, method = 'pearson') # cor 0.159737 , p 0.5136 use
cor.test(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$post], function(x){sum(x$Proportion)})), pat.inf$necrosis, method = 'pearson') # cor 0.1196951 , p 0.6255 use
cor.test(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$FU], function(x){sum(x$Proportion)})), pat.inf$necrosis, method = 'pearson') # cor 0.1007407 , p 0.6815 use

########################## use ########################## 
cor.top1perITCnecrosis <- cbind(pat.inf, top1perITCcumufre = unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], function(x){sum(x$Proportion)})))

ggscatter(cor.top1perITCnecrosis, 'top1perITCcumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency of top 1% ITC', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggsave(filename = "cor.top1perITCnecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")


cor.Btop1perITCnecrosis <- cbind(necrosis= pat.inf$necrosis[c(4:6,8,10,12:14,16,19)], 
                                 top1perITCcumufre = unlist(lapply(tcr.immunarch.top1per$data[pat.inf$B[c(4:6,8,10,12:14,16,19)]], function(x){sum(x$Proportion)}))) %>% as.data.frame(.)
cor.Btop1perITCnecrosis$label <- pat.inf$B[c(4:6,8,10,12:14,16,19)]

p.corBtop1ITCnecrosis <- ggscatter(cor.Btop1perITCnecrosis, 'top1perITCcumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.2, label.sep = "\n"),
          xlab = 'Cumulative frequency of top 1% ITC', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'label') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======


cor.NTtop1perITCnecrosis <- cbind(pat.inf, top1perITCcumufre = unlist(lapply(tcr.immunarch.top1per$data[pat.inf$NT], function(x){sum(x$Proportion)})))

p.corNTtop1ITCnecrosis <- ggscatter(cor.NTtop1perITCnecrosis, 'top1perITCcumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.2, label.sep = "\n"),
          xlab = 'Cumulative frequency of top 1% ITC', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'NT') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggarrange(p.corBtop1ITCnecrosis, p.corNTtop1ITCnecrosis)
ggsave(filename = "cor.BnNTtop1perITCnecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=6,units="in")
########################## use ########################## 

