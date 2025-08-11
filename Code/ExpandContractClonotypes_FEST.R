###########  FEST (Functional Expansion of Specific T-cells): expanded and contracted clonotypes ########### 
########### ref: PMID: 29895573 Cancer Immunology Research
###########  http://www.stat-apps.onc.jhmi.edu/FEST 

library(readxl)
prepost.FEST <- lapply(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/prepost/OR0/'), function(x){
  read_excel(paste0('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/prepost/OR0/',x), sheet = 'ref_comparison_only')
})
prepost.FEST.expand <- prepost.FEST[1:19]
prepost.FEST.contract <- prepost.FEST[20:38]

prepost.FEST <- list(expand = prepost.FEST.expand, contract = prepost.FEST.contract)
names(prepost.FEST$expand) <- paste0('QM', substr(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/prepost/OR0/')[1:19],8,10))
names(prepost.FEST$contract) <- paste0('QM', substr(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/prepost/OR0/')[1:19],8,10))

prepost.FEST$expand <- prepost.FEST$expand[pat.inf$patientID]
prepost.FEST$contract <- prepost.FEST$contract[pat.inf$patientID]

prepost.FEST <- lapply(prepost.FEST, function(x){
  lapply(x, function(y){
    colnames(y) <- c('CDR3.aa', 'FDR', 'sig_comp', 'OR', 'post_abundance', 'pre_abundance', 'post_percent', 'pre_percent')
    return(y)
  })
})

unlist(lapply(prepost.FEST$expand, nrow)) # 282
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 27    32    14    13     6    13    11    21    11     8     4    38    10    10     3    15    21     7    18 
unlist(lapply(prepost.FEST$contract, nrow)) # 245 QM841无
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 3     7     3     0     5     7    17     8     1    10    18     5     3     6    37    43    20    11    41 


expand.fre <- do.call(rbind,lapply(prepost.FEST$expand, function(x){
  apply(x[,5:8], 2, sum)
}))

contract.fre <- do.call(rbind,lapply(prepost.FEST$contract[c(1:3,5:19)], function(x){
  apply(x[,5:8], 2, sum)
})) # QM841 has no contracted colontypes，add 0.
contract.fre <- rbind(contract.fre[1:3,], QM841 = matrix(rep(0,4), ncol = 4), contract.fre[4:18,])
rownames(contract.fre)[4] <- 'QM841'





postFU.FEST <- lapply(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/postFU/OR0/'), function(x){
  read_excel(paste0('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/postFU/OR0/',x), sheet = 'ref_comparison_only')
})
postFU.FEST.expand <- postFU.FEST[1:19]
postFU.FEST.contract <- postFU.FEST[20:38]

postFU.FEST <- list(expand = postFU.FEST.expand, contract = postFU.FEST.contract)
names(postFU.FEST$expand) <- paste0('QM', substr(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/postFU/OR0/')[1:19],7,9))
names(postFU.FEST$contract) <- paste0('QM', substr(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/postFU/OR0/')[1:19],7,9))

postFU.FEST$expand <- postFU.FEST$expand[pat.inf$patientID]
postFU.FEST$contract <- postFU.FEST$contract[pat.inf$patientID]

postFU.FEST <- lapply(postFU.FEST, function(x){
  lapply(x, function(y){
    colnames(y) <- c('CDR3.aa', 'FDR', 'sig_comp', 'OR', 'FU_abundance', 'post_abundance', 'FU_percent', 'post_percent')
    return(y)
  })
})

unlist(lapply(postFU.FEST$expand, nrow)) # 172
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 2    15     2     3     0     5    16    12     5     5    29     2     5     8    18    22    11     5     7 
unlist(lapply(postFU.FEST$contract, nrow)) # 532 
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 9    32    29    32     4    13    10    12    32    11    26    51    16    36    24    12    49    81    53 





BT.FEST <- lapply(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/BT/'), function(x){
  read_excel(paste0('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/BT/',x), sheet = 'ref_comparison_only')
})
BT.FEST.expand <- BT.FEST[11:20]
BT.FEST.contract <- BT.FEST[1:10]

BT.FEST <- list(expand = BT.FEST.expand, contract = BT.FEST.contract)
names(BT.FEST$expand) <- paste0('QM', substr(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/BT/')[11:20],3,5))
names(BT.FEST$contract) <- paste0('QM', substr(dir('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/BT/')[1:10],3,5))

BT.FEST$expand <- BT.FEST$expand[intersect(pat.inf$patientID, names(BT.FEST$expand))]
BT.FEST$contract <- BT.FEST$contract[intersect(pat.inf$patientID, names(BT.FEST$contract))]

BT.FEST <- lapply(BT.FEST, function(x){
  lapply(x, function(y){
    colnames(y) <- c('CDR3.aa', 'FDR', 'sig_comp', 'OR', 'B_abundance', 'T_abundance', 'B_percent', 'T_percent')
    return(y)
  })
})

unlist(lapply(BT.FEST$expand, nrow)) # 472
# QM841 QM781 QM824 QM791 QM803 QM788 QM789 QM814 QM795 QM784 
# 23    97    38    84    48    12     9    21   103    37 

unlist(lapply(BT.FEST$contract, nrow)) # 608 
# QM841 QM781 QM824 QM791 QM803 QM788 QM789 QM814 QM795 QM784 
# 35    21   111    45     7    83    92    45    78    91 

save(prepost.FEST, postFU.FEST, BT.FEST, file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/FEST.res.RData')



##### total unique dynamic clonotypes
unique(do.call(rbind, prepost.FEST$expand)$CDR3.aa) # 278
unique(do.call(rbind, prepost.FEST$contract[-4])$CDR3.aa) # 245
unique(do.call(rbind, postFU.FEST$expand[-5])$CDR3.aa) # 171
unique(do.call(rbind, postFU.FEST$contract)$CDR3.aa) # 531
unique(do.call(rbind, BT.FEST$expand)$CDR3.aa) # 470
unique(do.call(rbind, BT.FEST$contract)$CDR3.aa) # 592

unique(c(unique(do.call(rbind, prepost.FEST$expand)$CDR3.aa),
         unique(do.call(rbind, prepost.FEST$contract[-4])$CDR3.aa),
         unique(do.call(rbind, postFU.FEST$expand[-5])$CDR3.aa),
         unique(do.call(rbind, postFU.FEST$contract)$CDR3.aa))) # 981


## unique expanded and contracted clonotypes overlap with unique top 1% post-ITCs 
length(intersect(unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], function(x){(unique(x$CDR3.aa))})),
                 unique(c(unique(do.call(rbind, prepost.FEST$expand)$CDR3.aa),
                          unique(do.call(rbind, prepost.FEST$contract[-4])$CDR3.aa),
                          unique(do.call(rbind, postFU.FEST$expand[-5])$CDR3.aa),
                          unique(do.call(rbind, postFU.FEST$contract)$CDR3.aa))))) # 141




number.plot <- cbind(number = c(unlist(lapply(prepost.FEST$expand, nrow)), -unlist(lapply(prepost.FEST$contract, nrow))),
                     patID = rep(pat.inf$patientID,2), type = c(rep('expand', 19), rep('contract',19))) %>% as.data.frame(.)
number.plot$number <- as.numeric(number.plot$number)
number.plot['QM841.1','number'] <- 0

p.dynamic1 <- ggbarplot(number.plot, 'patID', 'number',
          fill = 'type', color = 'type', palette = "Paired",
          label = TRUE, lab.col = "black", lab.pos = "in") +
  scale_y_continuous(breaks = c(-60,-40,-20,0,20,40), limits = c(-60,40))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 90)) +
  ylab('Number of dynamic clonotypes') +
  xlab('Pretreatment to posttreatment')




number.plot <- cbind(number = c(unlist(lapply(postFU.FEST$expand, nrow)), -unlist(lapply(postFU.FEST$contract, nrow))),
                     patID = rep(pat.inf$patientID,2), type = c(rep('expand', 19), rep('contract',19))) %>% as.data.frame(.)
number.plot$number <- as.numeric(number.plot$number)
number.plot['QM781','number'] <- 0

p.dynamic2 <- ggbarplot(number.plot, 'patID', 'number',
                        fill = 'type', color = 'type', palette = "Paired",
                        label = TRUE, lab.col = "black", lab.pos = "in") +
  scale_y_continuous(breaks = c(-80,-60,-40,-20,0,20,40), limits = c(-82,40))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 90)) +
  ylab('Number of dynamic clonotypes') +
  xlab('Posttreatment to follow-up')

ggarrange(p.dynamic1, p.dynamic2, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave(filename = "numberExpandContract_FEST.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use


number.plot3 <- cbind(number = c(unlist(lapply(BT.FEST$expand, nrow)), -unlist(lapply(BT.FEST$contract, nrow))),
                     patID = rep(names(BT.FEST$expand),2), type = c(rep('expand', 10), rep('contract',10))) %>% as.data.frame(.)
number.plot3$number <- as.numeric(number.plot3$number)


p.dynamic3 <- ggbarplot(number.plot3, 'patID', 'number',
                        fill = 'type', color = 'type', palette = "Paired",
                        label = TRUE, lab.col = "black", lab.pos = "in") +
  scale_y_continuous(breaks = c(-100,-50,0,50,100), limits = c(-115,110))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 90)) +
  ylab('Number of dynamic clonotypes') +
  xlab('preICI.Bx to postICI.T')
ggsave(filename = "numberExpandContract_BT.FEST.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=4,height=6,units="in") ###### use


########### ########### ########### ########### ########### ###########
########### percentage of dynamic and non-dynamic clonotypes in ITC

dynamic.per.prepost <- lapply(pat.inf$patientID, function(x){
  dyna <- unique(c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa))
  dyna.per <- length(unique(intersect(dyna,
                                      tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(dyna)
  dyna.per <- dyna.per*100
  
  allclone <- rbind(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]],
               tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]])
  nondyna <- allclone[is.na(match(allclone$CDR3.aa, dyna)),]
  nondyna.per <- length(unique(intersect(unique(nondyna$CDR3.aa),
                                         tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(nondyna$CDR3.aa))
  nondyna.per <- nondyna.per*100
  
  
  dyna.expand <- prepost.FEST$expand[[x]]$CDR3.aa
  dyna.expand.per <- length(unique(intersect(dyna.expand,
                                             tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(dyna.expand) * 100
  nondyna.expand <- allclone[is.na(match(allclone$CDR3.aa, dyna.expand)),]
  nondyna.expand.per <- length(unique(intersect(unique(nondyna.expand$CDR3.aa),
                                                tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(nondyna.expand$CDR3.aa)) *100
  
  
  dyna.contract <- prepost.FEST$contract[[x]]$CDR3.aa
  dyna.contract.per <- length(unique(intersect(dyna.contract,
                                             tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(dyna.contract) * 100
  nondyna.contract <- allclone[is.na(match(allclone$CDR3.aa, dyna.contract)),]
  nondyna.contract.per <- length(unique(intersect(unique(nondyna.contract$CDR3.aa),
                                                tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(nondyna.contract$CDR3.aa)) *100
  
  index <- c(dyna.per = dyna.per, nondyna.per = nondyna.per, 
             dyna.expand.per = dyna.expand.per, nondyna.expand.per = nondyna.expand.per, 
             dyna.contract.per = dyna.contract.per, nondyna.contract.per = nondyna.contract.per)
  return(index)
})
names(dynamic.per.prepost)<- pat.inf$patientID

dynamic.per.prepost <- do.call(rbind, dynamic.per.prepost) %>% as.data.frame(.)


dynamic.per.prepost.plot <- cbind(reshape2::melt(dynamic.per.prepost),
                                  response = rep(c(rep('NR',12), rep('R',7)), 6))

dynamic.per.prepost.plot$group <- paste0(dynamic.per.prepost.plot$variable, '.', dynamic.per.prepost.plot$response)

p.dynamic.per1 <- ggboxplot(dynamic.per.prepost.plot, x = 'group', y = 'value',add = 'dotplot', ylab = 'Percentage of clonotypes in tumor', xlab = 'Pretreatment to posttreatment', bxp.errorbar = TRUE,
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8), c(9,10),c(11,12))) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0.05,
                     comparisons = list(c(1,3),c(2,4),c(5,7),c(6,8),
                                        c(9,11),c(10,12)))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use


p.dynamic.per1.1 <- ggboxplot(dynamic.per.prepost.plot, x = 'group', y = 'value',add = 'dotplot', ylab = 'Percentage of clonotypes in tumor', xlab = 'Pretreatment to posttreatment', bxp.errorbar = TRUE,
                            add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8), c(9,10),c(11,12))) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0.05,
                     comparisons = list(c(1,3),c(2,4),c(5,7),c(6,8),
                                        c(9,11),c(10,12), c(3,5),c(4,6),c(3,9),c(4,10),
                                        c(5,9),c(6,10)))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use



dynamic.per.postFU <- lapply(pat.inf$patientID, function(x){
  dyna <- unique(c(postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
  dyna.per <- length(unique(intersect(dyna,
                                      tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(dyna) *100
  
  allclone <- rbind(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]],
                    tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'FU']]])
  nondyna <- allclone[is.na(match(allclone$CDR3.aa, dyna)),]
  nondyna.per <- length(unique(intersect(unique(nondyna$CDR3.aa),
                                         tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(nondyna$CDR3.aa)) * 100
  
  
  dyna.expand <- postFU.FEST$expand[[x]]$CDR3.aa
  dyna.expand.per <- length(unique(intersect(dyna.expand,
                                             tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(dyna.expand) * 100
  nondyna.expand <- allclone[is.na(match(allclone$CDR3.aa, dyna.expand)),]
  nondyna.expand.per <- length(unique(intersect(unique(nondyna.expand$CDR3.aa),
                                                tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(nondyna.expand$CDR3.aa)) *100
  
  
  dyna.contract <- postFU.FEST$contract[[x]]$CDR3.aa
  dyna.contract.per <- length(unique(intersect(dyna.contract,
                                               tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(dyna.contract) * 100
  nondyna.contract <- allclone[is.na(match(allclone$CDR3.aa, dyna.contract)),]
  nondyna.contract.per <- length(unique(intersect(unique(nondyna.contract$CDR3.aa),
                                                  tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(nondyna.contract$CDR3.aa)) *100
  
  index <- c(dyna.per = dyna.per, nondyna.per = nondyna.per, 
             dyna.expand.per = dyna.expand.per, nondyna.expand.per = nondyna.expand.per, 
             dyna.contract.per = dyna.contract.per, nondyna.contract.per = nondyna.contract.per)
  return(index)
})
names(dynamic.per.postFU)<- pat.inf$patientID

dynamic.per.postFU <- do.call(rbind, dynamic.per.postFU) %>% as.data.frame(.)

dynamic.per.postFU.plot <- cbind(reshape2::melt(dynamic.per.postFU),
                                  response = rep(c(rep('NR',12), rep('R',7)), 6))

dynamic.per.postFU.plot$group <- paste0(dynamic.per.postFU.plot$variable, '.', dynamic.per.postFU.plot$response)

p.dynamic.per2 <- ggboxplot(dynamic.per.postFU.plot, x = 'group', y = 'value',add = 'dotplot', ylab = 'Percentage of clonotypes in tumor', xlab = 'Posttreatment to follow-up', bxp.errorbar = TRUE,
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8), c(9,10),c(11,12))) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0.05,
                     comparisons = list(c(1,3),c(2,4),c(5,7),c(6,8),
                                        c(9,11),c(10,12)))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use


p.dynamic <- ggarrange(p.dynamic1, p.dynamic2, nrow = 1, ncol = 2, common.legend = TRUE)
p.dynamic.per <- ggarrange(p.dynamic.per1, p.dynamic.per2, nrow = 1, ncol = 2, common.legend = TRUE)

ggarrange(p.dynamic, p.dynamic.per, nrow = 2, ncol = 1)

ggsave(filename = "PBMCDynamicITC.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=10,units="in") ###### use

p.dynamic.per2.1 <- ggboxplot(dynamic.per.postFU.plot, x = 'group', y = 'value',add = 'dotplot', ylab = 'Percentage of clonotypes in tumor', xlab = 'Posttreatment to follow-up', bxp.errorbar = TRUE,
                            add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8), c(9,10),c(11,12))) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0.05,
                     comparisons = list(c(1,3),c(2,4),c(5,7),c(6,8),
                                        c(9,11),c(10,12),c(3,5),c(4,6),c(3,9),c(4,10),
                                        c(5,9),c(6,10)))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

ggarrange(p.dynamic.per1.1, p.dynamic.per2.1, nrow = 1, ncol = 2)
ggsave(filename = "PBMCDynamicITC1.1.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=10,units="in") ###### use


dynamic.per.plot <- rbind(dynamic.per.prepost.plot[dynamic.per.prepost.plot$variable %in% c('dyna.expand.per','dyna.contract.per','nondyna.per'),],
                          dynamic.per.postFU.plot[dynamic.per.postFU.plot$variable %in% c('dyna.expand.per','dyna.contract.per','nondyna.per'),])
dynamic.per.plot$stage <- c(rep('prepost',57), rep('postFU',57))
dynamic.per.plot$group1 <- paste0(dynamic.per.plot$variable, '.', dynamic.per.plot$stage)
dynamic.per.plot$variable <- factor(dynamic.per.plot$variable, levels = c('dyna.expand.per','dyna.contract.per','nondyna.per'))
dynamic.per.plot$group1 <- factor(dynamic.per.plot$group1, levels = c('dyna.expand.per.prepost','dyna.expand.per.postFU',
                                                                      'dyna.contract.per.prepost','dyna.contract.per.postFU',
                                                                      'nondyna.per.prepost','nondyna.per.postFU'))

p.dynamic.per3 <- ggboxplot(dynamic.per.plot, x = 'group1', y = 'value', ylab = 'Percentage of clonotypes in tumor', xlab = 'Posttreatment to follow-up', bxp.errorbar = TRUE,
                            add.params = list(size=0.5), color = 'stage', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
  #                    comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8), c(9,10),c(11,12))) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0.05,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(1,4),c(2,3)))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

ggsave(filename = "PBMCDynamicITC_prepost.postFU.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=10,units="in") ###### use

dynamic.per.plot$group2 <- paste0(dynamic.per.plot$group1,'.',dynamic.per.plot$response)
dynamic.per.plot$group2 <- factor(dynamic.per.plot$group2, levels = c("dyna.expand.per.prepost.NR","dyna.expand.per.postFU.NR",
                                                                      "dyna.expand.per.prepost.R","dyna.expand.per.postFU.R",
                                                                      "dyna.contract.per.prepost.NR","dyna.contract.per.postFU.NR",
                                                                      "dyna.contract.per.prepost.R","dyna.contract.per.postFU.R",
                                                                      "nondyna.per.prepost.NR","nondyna.per.postFU.NR",
                                                                      "nondyna.per.prepost.R","nondyna.per.postFU.R"))

p.dynamic.per3.1 <- ggboxplot(dynamic.per.plot, x = 'group2', y = 'value', ylab = 'Percentage of clonotypes in tumor', xlab = 'Posttreatment to follow-up', bxp.errorbar = TRUE,
                            add.params = list(size=0.5), color = 'stage', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
  #                    comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8), c(9,10),c(11,12))) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4),c(5,6),c(7,8),c(9,10),c(11,12)))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # use

###########################  no significance, unuse ###########################  

dynamic.ITC <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% prepost.FEST$expand[[x]]$CDR3.aa,]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% prepost.FEST$contract[[x]]$CDR3.aa,]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.ITC) <- pat.inf$patientID

dynamic.ITC <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.ITC[[x]]$dyna.expandITC}),
                    contract = lapply(pat.inf$patientID, function(x){dynamic.ITC[[x]]$dyna.contractITC}))
names(dynamic.ITC$expand) <- pat.inf$patientID
names(dynamic.ITC$contract) <- pat.inf$patientID

dynamic.ITC.cumufre <- lapply(dynamic.ITC, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})


cor.test(dynamic.ITC.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # r = 0.3967475, p = 0.0926
cor.test(dynamic.ITC.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # r = 0.210203, p = 0.3877
cor.test(dynamic.ITC.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # r = 0.1062843, p = 0.665
cor.test(dynamic.ITC.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # r = 0.3289369, p = 0.1691

cor.test(dynamic.ITC.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.5626965 , p = 0.05682
cor.test(dynamic.ITC.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = 0.3468027 , p = 0.446
cor.test(dynamic.ITC.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.1797135 , p = 0.5762
cor.test(dynamic.ITC.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = -0.2738537  , p = 0.5523


dynamic.top1ITC <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% prepost.FEST$expand[[x]]$CDR3.aa,]
  dyna.contractITC <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% prepost.FEST$contract[[x]]$CDR3.aa,]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.top1ITC) <- pat.inf$patientID

dynamic.top1ITC <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.top1ITC[[x]]$dyna.expandITC}),
                    contract = lapply(pat.inf$patientID, function(x){dynamic.top1ITC[[x]]$dyna.contractITC}))
names(dynamic.top1ITC$expand) <- pat.inf$patientID
names(dynamic.top1ITC$contract) <- pat.inf$patientID

dynamic.top1ITC.cumufre <- lapply(dynamic.top1ITC, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})


cor.test(dynamic.top1ITC.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # r = 0.4169423, p = 0.07574
cor.test(dynamic.top1ITC.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # r = 0.2514726 , p = 0.299
cor.test(dynamic.top1ITC.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # r = 0.1115059 , p = 0.6495
cor.test(dynamic.top1ITC.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # r = 0.3086071 , p = 0.1986

cor.test(dynamic.top1ITC.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.5740019 , p = 0.05098
cor.test(dynamic.top1ITC.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = 0.3440335 , p =  0.4499
cor.test(dynamic.top1ITC.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.1267041  , p = 0.6948
cor.test(dynamic.top1ITC.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = -0.2880126 , p = 0.5311



dynamic.top1ITC.1 <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% intersect(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa),]
  dyna.contractITC <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% intersect(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa),]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.top1ITC.1) <- pat.inf$patientID

dynamic.top1ITC.1 <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.top1ITC.1[[x]]$dyna.expandITC}),
                        contract = lapply(pat.inf$patientID, function(x){dynamic.top1ITC.1[[x]]$dyna.contractITC}))
names(dynamic.top1ITC.1$expand) <- pat.inf$patientID
names(dynamic.top1ITC.1$contract) <- pat.inf$patientID

dynamic.top1ITC.1.cumufre <- lapply(dynamic.top1ITC.1, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})


cor.test(dynamic.top1ITC.1.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.top1ITC.1.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # ns
cor.test(dynamic.top1ITC.1.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.top1ITC.1.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # ns

cor.test(dynamic.top1ITC.1.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.5762699  , p = 0.04986
cor.test(dynamic.top1ITC.1.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
cor.test(dynamic.top1ITC.1.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.top1ITC.1.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
###########################  no significance, unuse ###########################  












###########################  
dynamic.ITC.postFU <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% postFU.FEST$expand[[x]]$CDR3.aa,]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% postFU.FEST$contract[[x]]$CDR3.aa,]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.ITC.postFU) <- pat.inf$patientID

dynamic.ITC.postFU <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.ITC.postFU[[x]]$dyna.expandITC}),
                    contract = lapply(pat.inf$patientID, function(x){dynamic.ITC.postFU[[x]]$dyna.contractITC}))
names(dynamic.ITC.postFU$expand) <- pat.inf$patientID
names(dynamic.ITC.postFU$contract) <- pat.inf$patientID

dynamic.ITC.postFU.cumufre <- lapply(dynamic.ITC.postFU, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

cor.test(dynamic.ITC.postFU.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # ns
cor.test(dynamic.ITC.postFU.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # r = 0.7683098, p = 0.0001219 use
cor.test(dynamic.ITC.postFU.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # r = 0.5883925 , p = 0.008049

cor.test(dynamic.ITC.postFU.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = ns
cor.test(dynamic.ITC.postFU.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.5521211  , p = 0.0627
cor.test(dynamic.ITC.postFU.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = 0.8295468   , p = 0.02097 use
###########################  


###########################  use ###########################  
dynamic.ITC.prepostFU <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% c(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa),]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% c(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa),]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})

dynamic.ITC.prepostFU <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% unique(c(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa)),]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% unique(c(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa)),]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})

names(dynamic.ITC.prepostFU) <- pat.inf$patientID

dynamic.ITC.prepostFU <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.ITC.prepostFU[[x]]$dyna.expandITC}),
                           contract = lapply(pat.inf$patientID, function(x){dynamic.ITC.prepostFU[[x]]$dyna.contractITC}))
names(dynamic.ITC.prepostFU$expand) <- pat.inf$patientID
names(dynamic.ITC.prepostFU$contract) <- pat.inf$patientID

dynamic.ITC.prepostFU.cumufre <- lapply(dynamic.ITC.prepostFU, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

cor.test(dynamic.ITC.prepostFU.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # r = 0.753969 , p = 0.0001926 use
# cor.test(dynamic.ITC.prepostFU.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # r = 0.6306091 , p = 0.003796
cor.test(dynamic.ITC.prepostFU.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # ns
# cor.test(dynamic.ITC.prepostFU.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # ns

cor.test(dynamic.ITC.prepostFU.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.5820178  , p = 0.0471
cor.test(dynamic.ITC.prepostFU.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # r = 0.6615614  , p = 0.1056
cor.test(dynamic.ITC.prepostFU.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.prepostFU.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns

cor.expandITCnecrosis <- cbind(pat.inf, clonality = dynamic.ITC.prepostFU.cumufre$expand$V1)

ggscatter(cor.expandITCnecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency of postICI_PBMC expanded clonotypes', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggsave(filename = "cor.cumufrepostexpandnecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")

###########################  use ###########################  


###########################  unuse ###########################  
dynamic.ITC.prepostFU.1 <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% intersect(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa),]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% intersect(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa),]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.ITC.prepostFU.1) <- pat.inf$patientID

dynamic.ITC.prepostFU.1 <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.ITC.prepostFU.1[[x]]$dyna.expandITC}),
                              contract = lapply(pat.inf$patientID, function(x){dynamic.ITC.prepostFU.1[[x]]$dyna.contractITC}))
names(dynamic.ITC.prepostFU.1$expand) <- pat.inf$patientID
names(dynamic.ITC.prepostFU.1$contract) <- pat.inf$patientID

dynamic.ITC.prepostFU.1.cumufre <- lapply(dynamic.ITC.prepostFU.1, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

cor.test(dynamic.ITC.prepostFU.1.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # r = 0.380826 , p = 0.1077 ns
# cor.test(dynamic.ITC.prepostFU.1.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # r = 0.2014079 , p = 0.4083 ns
cor.test(dynamic.ITC.prepostFU.1.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # ns
# cor.test(dynamic.ITC.prepostFU.1.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # ns

cor.test(dynamic.ITC.prepostFU.1.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # r = 0.5545002  , p = 0.06134 ns
cor.test(dynamic.ITC.prepostFU.1.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
cor.test(dynamic.ITC.prepostFU.1.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.prepostFU.1.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns



dynamic.ITC.postFU.1 <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]]$CDR3.aa %in% postFU.FEST$expand[[x]]$CDR3.aa,]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]]$CDR3.aa %in% postFU.FEST$contract[[x]]$CDR3.aa,]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.ITC.postFU.1) <- pat.inf$patientID

dynamic.ITC.postFU.1 <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.ITC.postFU.1[[x]]$dyna.expandITC}),
                             contract = lapply(pat.inf$patientID, function(x){dynamic.ITC.postFU.1[[x]]$dyna.contractITC}))
names(dynamic.ITC.postFU.1$expand) <- pat.inf$patientID
names(dynamic.ITC.postFU.1$contract) <- pat.inf$patientID

dynamic.ITC.postFU.1.cumufre <- lapply(dynamic.ITC.postFU.1, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

cor.test(dynamic.ITC.postFU.1.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.1.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # ns
cor.test(dynamic.ITC.postFU.1.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.1.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # ns

cor.test(dynamic.ITC.postFU.1.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.1.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.1.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.1.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
###########################  

###########################  
dynamic.ITC.postFU.2 <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]]$CDR3.aa %in% intersect(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa),]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'post']]]$CDR3.aa %in% intersect(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa),]
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC)
  return(dyna.ITC)
})
names(dynamic.ITC.postFU.2) <- pat.inf$patientID

dynamic.ITC.postFU.2 <- list(expand = lapply(pat.inf$patientID, function(x){dynamic.ITC.postFU.2[[x]]$dyna.expandITC}),
                             contract = lapply(pat.inf$patientID, function(x){dynamic.ITC.postFU.2[[x]]$dyna.contractITC}))
names(dynamic.ITC.postFU.2$expand) <- pat.inf$patientID
names(dynamic.ITC.postFU.2$contract) <- pat.inf$patientID

dynamic.ITC.postFU.2.cumufre <- lapply(dynamic.ITC.postFU.2, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

cor.test(dynamic.ITC.postFU.2.cumufre$expand$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.2.cumufre$expand$V1, pat.inf$necrosis, method = 'spearman') # ns
cor.test(dynamic.ITC.postFU.2.cumufre$contract$V1, pat.inf$necrosis, method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.2.cumufre$contract$V1, pat.inf$necrosis, method = 'spearman') # ns

cor.test(dynamic.ITC.postFU.2.cumufre$expand$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.2.cumufre$expand$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.2.cumufre$contract$V1[1:12], pat.inf$necrosis[1:12], method = 'pearson') # ns
cor.test(dynamic.ITC.postFU.2.cumufre$contract$V1[13:19], pat.inf$necrosis[13:19], method = 'pearson') # ns
###########################  
###########################  unuse ###########################




########### ########### ########### ########### ########### ########### 
########### dynamic gene AV line plot ###########

turnclone <- lapply(pat.inf$patientID, function(x){
  index1 <- intersect(prepost.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa)
  index2 <- intersect(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa)
  return(c(length(index1),length(index2)))
}) 

continueclone <- lapply(pat.inf$patientID, function(x){
  index1 <- intersect(prepost.FEST$expand[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa)
  index2 <- intersect(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa)
  return(c(length(index1),length(index2)))
}) 




prepost.postFU.same <- lapply(pat.inf$patientID, function(x){
  intersect(c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa), 
            c(postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
})
names(prepost.postFU.same)  <- pat.inf$patientID

unlist(lapply(prepost.postFU.same,length)) # 共有 237 clonotypes
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 5    26    11    11     1     6    19    11     3     4     7    28     6     9    18    20    20     8    24 

dyn.duration.data <- lapply(pat.inf$patientID, function(x){
  pre <- tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]]$CDR3.aa %in% prepost.postFU.same[[x]],]
  pre <- split(pre, pre$CDR3.aa)
  pre <- lapply(pre, function(y){
    fresum <- sum(y$Proportion)
    countsum <- sum(y$Clones)
    index1 <- y[1,]
    index1[,'Proportion'] <- fresum
    index1[,'Clones'] <- countsum
    return(index1)
  })
  pre <- do.call(rbind, pre)
  # return(pre)
  
  post <- tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]]$CDR3.aa %in% prepost.postFU.same[[x]],]
  post <- split(post, post$CDR3.aa)
  post <- lapply(post, function(z){
    fresum <- sum(z$Proportion)
    countsum <- sum(z$Clones)
    index1 <- z[1,]
    index1[,'Proportion'] <- fresum
    index1[,'Clones'] <- countsum
    return(index1)
  })
  post <- do.call(rbind, post)
  
  FU <- tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'FU']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'FU']]]$CDR3.aa %in% prepost.postFU.same[[x]],]
  FU <- split(FU, FU$CDR3.aa)
  FU <- lapply(FU, function(n){
    fresum <- sum(n$Proportion)
    countsum <- sum(n$Clones)
    index1 <- n[1,]
    index1[,'Proportion'] <- fresum
    index1[,'Clones'] <- countsum
    return(index1)
  })
  FU <- do.call(rbind, FU)
  
  index <- list(pre = pre, post = post, FU = FU)
  return(index)
  
}) 
names(dyn.duration.data) <- pat.inf$patientID


  
dyn.duration.data <- list(pre.dur = lapply(pat.inf$patientID, function(x){dyn.duration.data[[x]]$pre}),
                          post.dur = lapply(pat.inf$patientID, function(x){dyn.duration.data[[x]]$post}),
                          FU.dur = lapply(pat.inf$patientID, function(x){dyn.duration.data[[x]]$FU}))
names(dyn.duration.data$pre.dur) <- pat.inf$patientID
names(dyn.duration.data$post.dur) <- pat.inf$patientID
names(dyn.duration.data$FU.dur) <- pat.inf$patientID


expdata <- cbind(pre = do.call(rbind,dyn.duration.data$pre.dur)$Proportion, 
                 post = do.call(rbind,dyn.duration.data$post.dur)$Proportion,
                 FU = do.call(rbind,dyn.duration.data$FU.dur)$Proportion) %>% as.data.frame(.)

expdata <- cbind(Time = c(0,1,2), t(expdata)) %>% as.data.frame(.)

expdata.zscore <- apply(t(expdata[,2:ncol(expdata)]), 1, function(dataset){
  datastd <- sd(dataset)*sqrt((length(dataset)-1)/(length(dataset)))
  datamean <- mean(dataset)
  Z_score <- (dataset - datamean) / datastd
  return(Z_score)
}) %>% as.data.frame(.)

expdata.zscore$Time <- expdata$Time
expdata.zscore <- expdata.zscore[,c(ncol(expdata.zscore),1:(ncol(expdata.zscore)-1))]

expdata.plot <- reshape2::melt(expdata.zscore,id=1)

p.dyn.duration <- ggplot(expdata.plot,aes(x=as.character(Time),y=value,group=variable))+
  # geom_line(aes(colour=color))+
  geom_line(alpha=0.3, colour="grey")+
  geom_point(aes(colour=value), alpha=0.5)+
  scale_colour_gradient2(low = "blue", mid = 'grey', high = "red")+
  scale_x_discrete(labels=c('0'='preICI_PBMC','1'='postICI_PBMC','2'='FollowUp_PBMC'))+
  # scale_x_discrete(labels=c('1'='preICI_biopsy','2'='postICI_tumor','3'='postICI_nontumor'))+
  xlab(' ')+
  ylab('Frequency of duration dynamic clonotypes')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.66) ## use

####################################

dyn.duration.ITC.data <- lapply(pat.inf$patientID, function(x){
  pre.ITC <- dyn.duration.data$pre.dur[[x]][dyn.duration.data$pre.dur[[x]]$CDR3.aa %in% tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'T']]]$CDR3.aa,]
  post.ITC <- dyn.duration.data$post.dur[[x]][dyn.duration.data$post.dur[[x]]$CDR3.aa %in% tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'T']]]$CDR3.aa,]
  FU.ITC <- dyn.duration.data$FU.dur[[x]][dyn.duration.data$FU.dur[[x]]$CDR3.aa %in% tcr.immunarch$data[[pat.inf[pat.inf$patientID == x,'T']]]$CDR3.aa,]
  ITC <- list(pre.ITC = pre.ITC, post.ITC = post.ITC, FU.ITC = FU.ITC)
  return(ITC)
})
names(dyn.duration.ITC.data) <- pat.inf$patientID
dyn.duration.ITC.data <- list(pre.durITC = lapply(pat.inf$patientID, function(x){dyn.duration.ITC.data[[x]]$pre.ITC}),
                              post.durITC = lapply(pat.inf$patientID, function(x){dyn.duration.ITC.data[[x]]$post.ITC}),
                              FU.durITC = lapply(pat.inf$patientID, function(x){dyn.duration.ITC.data[[x]]$FU.ITC}))
names(dyn.duration.ITC.data$pre.durITC) <- pat.inf$patientID
names(dyn.duration.ITC.data$post.durITC) <- pat.inf$patientID
names(dyn.duration.ITC.data$FU.durITC) <- pat.inf$patientID

unlist(lapply(dyn.duration.ITC.data$pre.durITC, nrow)) # 共有 186 clonotypes
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 5    17     9     6     1     6    11     8     2     2     4    17     6     9    16    20    20     8    19 

# unlist(lapply(prepost.postFU.same,length)) # 共有 237 clonotypes
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 5    26    11    11     1     6    19    11     3     4     7    28     6     9    18    20    20     8    24 



dyn.duration.top1ITC.data <- lapply(pat.inf$patientID, function(x){
  pre.ITC <- dyn.duration.data$pre.dur[[x]][dyn.duration.data$pre.dur[[x]]$CDR3.aa %in% tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x,'T']]]$CDR3.aa,]
  post.ITC <- dyn.duration.data$post.dur[[x]][dyn.duration.data$post.dur[[x]]$CDR3.aa %in% tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x,'T']]]$CDR3.aa,]
  FU.ITC <- dyn.duration.data$FU.dur[[x]][dyn.duration.data$FU.dur[[x]]$CDR3.aa %in% tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x,'T']]]$CDR3.aa,]
  ITC <- list(pre.ITC = pre.ITC, post.ITC = post.ITC, FU.ITC = FU.ITC)
  return(ITC)
})
names(dyn.duration.top1ITC.data) <- pat.inf$patientID
dyn.duration.top1ITC.data <- list(pre.durtop1ITC = lapply(pat.inf$patientID, function(x){dyn.duration.top1ITC.data[[x]]$pre.ITC}),
                              post.durtop1ITC = lapply(pat.inf$patientID, function(x){dyn.duration.top1ITC.data[[x]]$post.ITC}),
                              FU.durtop1ITC = lapply(pat.inf$patientID, function(x){dyn.duration.top1ITC.data[[x]]$FU.ITC}))
names(dyn.duration.top1ITC.data$pre.durtop1ITC) <- pat.inf$patientID
names(dyn.duration.top1ITC.data$post.durtop1ITC) <- pat.inf$patientID
names(dyn.duration.top1ITC.data$FU.durtop1ITC) <- pat.inf$patientID

unlist(lapply(dyn.duration.top1ITC.data$pre.durtop1ITC, nrow)) # 36 clonotypes
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 1     3     1     0     1     3     0     0     1     0     1     2     1     1     2     6     9     1     3 





expdata.ITC <- cbind(pre = do.call(rbind,dyn.duration.ITC.data$pre.durITC)$Proportion, 
                 post = do.call(rbind,dyn.duration.ITC.data$post.durITC)$Proportion,
                 FU = do.call(rbind,dyn.duration.ITC.data$FU.durITC)$Proportion) %>% as.data.frame(.)

expdata.ITC <- cbind(Time = c(0,1,2), t(expdata.ITC)) %>% as.data.frame(.)

expdata.ITC.zscore <- apply(t(expdata.ITC[,2:ncol(expdata.ITC)]), 1, function(dataset){
  datastd <- sd(dataset)*sqrt((length(dataset)-1)/(length(dataset)))
  datamean <- mean(dataset)
  Z_score <- (dataset - datamean) / datastd
  return(Z_score)
}) %>% as.data.frame(.)

expdata.ITC.zscore$Time <- expdata$Time
expdata.ITC.zscore <- expdata.ITC.zscore[,c(ncol(expdata.ITC.zscore),1:(ncol(expdata.ITC.zscore)-1))]

expdata.ITC.plot <- reshape2::melt(expdata.ITC.zscore,id=1)

p.dyn.durationITC <- ggplot(expdata.ITC.plot,aes(x=as.character(Time),y=value,group=variable))+
  # geom_line(aes(colour=color))+
  geom_line(alpha=0.3, colour="grey")+
  geom_point(aes(colour=value), alpha=0.5)+
  scale_colour_gradient2(low = "blue", mid = 'grey', high = "red")+
  scale_x_discrete(labels=c('0'='preICI_PBMC','1'='postICI_PBMC','2'='FollowUp_PBMC'))+
  # scale_x_discrete(labels=c('1'='preICI_biopsy','2'='postICI_tumor','3'='postICI_nontumor'))+
  xlab(' ')+
  ylab('Frequency of duration dynamic ITC')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.66) ## use



expdata.top1ITC <- cbind(pre = do.call(rbind,dyn.duration.top1ITC.data$pre.durtop1ITC)$Proportion, 
                     post = do.call(rbind,dyn.duration.top1ITC.data$post.durtop1ITC)$Proportion,
                     FU = do.call(rbind,dyn.duration.top1ITC.data$FU.durtop1ITC)$Proportion) %>% as.data.frame(.)

expdata.top1ITC <- cbind(Time = c(0,1,2), t(expdata.top1ITC)) %>% as.data.frame(.)

expdata.top1ITC.zscore <- apply(t(expdata.top1ITC[,2:ncol(expdata.top1ITC)]), 1, function(dataset){
  datastd <- sd(dataset)*sqrt((length(dataset)-1)/(length(dataset)))
  datamean <- mean(dataset)
  Z_score <- (dataset - datamean) / datastd
  return(Z_score)
}) %>% as.data.frame(.)

expdata.top1ITC.zscore$Time <- expdata.top1ITC$Time
expdata.top1ITC.zscore <- expdata.top1ITC.zscore[,c(ncol(expdata.top1ITC.zscore),1:(ncol(expdata.top1ITC.zscore)-1))]

expdata.top1ITC.plot <- reshape2::melt(expdata.top1ITC.zscore,id=1)

p.dyn.durationtop1ITC <- ggplot(expdata.top1ITC.plot,aes(x=as.character(Time),y=value,group=variable))+
  # geom_line(aes(colour=color))+
  geom_line(alpha=0.3, colour="grey")+
  geom_point(aes(colour=value), alpha=0.5)+
  scale_colour_gradient2(low = "blue", mid = 'grey', high = "red")+
  scale_x_discrete(labels=c('0'='preICI_PBMC','1'='postICI_PBMC','2'='FollowUp_PBMC'))+
  # scale_x_discrete(labels=c('1'='preICI_biopsy','2'='postICI_tumor','3'='postICI_nontumor'))+
  xlab(' ')+
  ylab('Frequency of duration dynamic top 1% ITC')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.66) ## use

ggarrange(p.dyn.duration, p.dyn.durationITC, p.dyn.durationtop1ITC, nrow = 2, ncol = 2)

ggsave(filename = "DurationDynamicClonotypesExp.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=10,units="in") ###### use




###################### dynamic clonotyps ##############################
dyn.data <- lapply(pat.inf$patientID, function(x){
  predyna <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]]$CDR3.aa, c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa,
                                                                                               postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
  postdyna <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]]$CDR3.aa, c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa,
                                                                                                postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
  FUdyna <- intersect(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'FU']]]$CDR3.aa, c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa,
                                                                                              postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
  prepostFUcomm <- intersect(predyna, intersect(postdyna, FUdyna))
  
  pre <- tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]]$CDR3.aa %in% prepostFUcomm,]
  pre <- split(pre, pre$CDR3.aa)
  pre <- lapply(pre, function(y){
    fresum <- sum(y$Proportion)
    countsum <- sum(y$Clones)
    index1 <- y[1,]
    index1[,'Proportion'] <- fresum
    index1[,'Clones'] <- countsum
    return(index1)
  })
  pre <- do.call(rbind, pre)
  # return(pre)
  
  post <- tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]]$CDR3.aa %in% prepostFUcomm,]
  post <- split(post, post$CDR3.aa)
  post <- lapply(post, function(z){
    fresum <- sum(z$Proportion)
    countsum <- sum(z$Clones)
    index1 <- z[1,]
    index1[,'Proportion'] <- fresum
    index1[,'Clones'] <- countsum
    return(index1)
  })
  post <- do.call(rbind, post)
  
  FU <- tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'FU']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'FU']]]$CDR3.aa %in% prepostFUcomm,]
  FU <- split(FU, FU$CDR3.aa)
  FU <- lapply(FU, function(n){
    fresum <- sum(n$Proportion)
    countsum <- sum(n$Clones)
    index1 <- n[1,]
    index1[,'Proportion'] <- fresum
    index1[,'Clones'] <- countsum
    return(index1)
  })
  FU <- do.call(rbind, FU)
  
  index <- list(pre = pre, post = post, FU = FU)
  return(index)
  
}) 
names(dyn.data) <- pat.inf$patientID



dyn.data <- list(pre = lapply(pat.inf$patientID, function(x){dyn.data[[x]]$pre}),
                          post = lapply(pat.inf$patientID, function(x){dyn.data[[x]]$post}),
                          FU = lapply(pat.inf$patientID, function(x){dyn.data[[x]]$FU}))
names(dyn.data$pre) <- pat.inf$patientID
names(dyn.data$post) <- pat.inf$patientID
names(dyn.data$FU) <- pat.inf$patientID


expdata.dyn <- cbind(pre = do.call(rbind,dyn.data$pre)$Proportion, 
                 post = do.call(rbind,dyn.data$post)$Proportion,
                 FU = do.call(rbind,dyn.data$FU)$Proportion) %>% as.data.frame(.)

expdata.dyn <- cbind(Time = c(0,1,2), t(expdata.dyn)) %>% as.data.frame(.)

expdata.dyn.zscore <- apply(t(expdata.dyn[,2:ncol(expdata.dyn)]), 1, function(dataset){
  datastd <- sd(dataset)*sqrt((length(dataset)-1)/(length(dataset)))
  datamean <- mean(dataset)
  Z_score <- (dataset - datamean) / datastd
  return(Z_score)
}) %>% as.data.frame(.)

expdata.dyn.zscore$Time <- expdata.dyn$Time
expdata.dyn.zscore <- expdata.dyn.zscore[,c(ncol(expdata.dyn.zscore),1:(ncol(expdata.dyn.zscore)-1))]

expdata.dyn.plot <- reshape2::melt(expdata.dyn.zscore,id=1)

p.dyn <- ggplot(expdata.dyn.plot,aes(x=as.character(Time),y=value,group=variable))+
  # geom_line(aes(colour=color))+
  geom_line(alpha=0.3, colour="grey")+
  geom_point(aes(colour=value), alpha=0.5)+
  scale_colour_gradient2(low = "blue", mid = 'grey', high = "red")+
  scale_x_discrete(labels=c('0'='preICI_PBMC','1'='postICI_PBMC','2'='FollowUp_PBMC'))+
  # scale_x_discrete(labels=c('1'='preICI_biopsy','2'='postICI_tumor','3'='postICI_nontumor'))+
  xlab(' ')+
  ylab('Frequency of dynamic clonotypes')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.66) ## use

ggarrange(p.dyn,p.dyn.duration, p.dyn.durationITC, p.dyn.durationtop1ITC, nrow = 2, ncol = 2, common.legend = TRUE)

ggsave(filename = "DurationDynamicClonotypesExp1.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=10,units="in") ###### use

#################################### use













