# Because peripherally dynamic clones are thought to be trafficking to
# the tumor, we tested if intratumoral clonal occupancy of clones with
# differential peripheral patterns differentiates MPRs from non-MPRs.
# ITCs were stratified as peripherally contracted, expanded, nondynamic, or tumor-resident only based on their different dynamic
# patterns during treatment. 

load('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/FEST/FEST.res.RData')


########## ITC classification: peripherally expanded, contracted, nondynamic, or tumor-resident

ITCdynprepost.Class1 <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% c(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa),]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% c(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa),]
  
  allclone <- rbind(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]],
                    tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]])
  
  dyna <- unique(c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
  nondyna <- allclone[is.na(match(allclone$CDR3.aa, dyna)),]
  
  nondynaITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa %in% unique(nondyna$CDR3.aa),]
  
  T.resident <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]][is.na(match(tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'T']]]$CDR3.aa, unique(allclone$CDR3.aa))),]
  
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC, 
                   nondynaITC = nondynaITC, T.resident = T.resident)
  return(dyna.ITC)
})
names(ITCdynprepost.Class1) <- pat.inf$patientID

ITCdynprepost.Class1 <- list(expand = lapply(pat.inf$patientID, function(x){ITCdynprepost.Class1[[x]]$dyna.expandITC}),
                            contract = lapply(pat.inf$patientID, function(x){ITCdynprepost.Class1[[x]]$dyna.contractITC}),
                            nondyna = lapply(pat.inf$patientID, function(x){ITCdynprepost.Class1[[x]]$nondynaITC}),
                            T.resident = lapply(pat.inf$patientID, function(x){ITCdynprepost.Class1[[x]]$T.resident}))
names(ITCdynprepost.Class1$expand) <- pat.inf$patientID
names(ITCdynprepost.Class1$contract) <- pat.inf$patientID
names(ITCdynprepost.Class1$nondyna) <- pat.inf$patientID
names(ITCdynprepost.Class1$T.resident) <- pat.inf$patientID


ITCdynprepost.Class1.cumufre <- lapply(ITCdynprepost.Class1, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})



# a <- tcr.immunarch$data$T809[is.na(match(tcr.immunarch$data$T809$CDR3.aa, unique(rbind(tcr.immunarch$data$pre809, tcr.immunarch$data$post809)$CDR3.aa))),]

unlist(lapply(ITCdynprepost.Class1$expand,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 21    30    28    26    13    43     8    19     6     5    14    44    35    38    18    28   135   151   104 

unlist(lapply(ITCdynprepost.Class1$contract,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 2     8     3     2     6    16    14    14     5    10    12     4     8     5    32    55    44    19    43

unlist(lapply(ITCdynprepost.Class1$nondyna,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 613   898   560   656   720  2404   122   809    99   309   171   414  1033   509   619  1266  2133  1841  1550

unlist(lapply(ITCdynprepost.Class1$T.resident,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 1297  2116  1326  1799  1226  4878   527  1241   210   414   561  1761  3546  1420  1081  1929  3557  4198  2642 


################ figure plot ################ use
##### cumulative frequency test between NR and R： boxplot

ITCdynprepost.Class1.cumufre.plot <- rbind(ITCdynprepost.Class1.cumufre$expand, ITCdynprepost.Class1.cumufre$contract,
                                          ITCdynprepost.Class1.cumufre$nondyna, ITCdynprepost.Class1.cumufre$T.resident)
ITCdynprepost.Class1.cumufre.plot <- cbind(ITCdynprepost.Class1.cumufre.plot, group = rep(c('expand','contract','nondyna','T.resident'), each = 19),
                                          response = rep(c(rep('NR',12),rep('R',7)), 4)) %>% as.data.frame(.)
ITCdynprepost.Class1.cumufre.plot$group1 <- paste0(ITCdynprepost.Class1.cumufre.plot$group,'.',ITCdynprepost.Class1.cumufre.plot$response)
ITCdynprepost.Class1.cumufre.plot$V1 <- as.numeric(ITCdynprepost.Class1.cumufre.plot$V1)

p.ITCdynprepost.Class1.cumufre <- ggboxplot(ITCdynprepost.Class1.cumufre.plot, x = 'group1', y = 'V1', add = 'dotplot', ylab = 'Cumulative frequency', xlab = '',
                                           add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use

ggsave(filename = "postexpandITC.boxplot_edit.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")


cor.test(pat.inf$necrosis, ITCdynprepost.Class1.cumufre$expand$V1) # cor: 0.753969 p: 0.0001926
cor.test(pat.inf$necrosis, ITCdynprepost.Class1.cumufre$contract$V1) # ns
cor.test(pat.inf$necrosis, ITCdynprepost.Class1.cumufre$nondyna$V1) # ns
cor.test(pat.inf$necrosis, ITCdynprepost.Class1.cumufre$T.resident$V1) # cor: -0.6583894 p: 0.002177

cor.postexpandITCnecrosis <- cbind(CYT = pat.inf$necrosis, postExpCumuFre = ITCdynprepost.Class1.cumufre$expand$V1) %>% as.data.frame(.)
cor.postexpandITCnecrosis$ID <- pat.inf$T
cor.postexpandITCnecrosis$responsive <- pat.inf$responsive

p.cor.postexpandITCnecrosis <- ggscatter(cor.postexpandITCnecrosis, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency of postICI.PBMC expanded clonotypes', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======


cor.TresidentITCnecrosis <- cbind(CYT = pat.inf$necrosis, postExpCumuFre = ITCdynprepost.Class1.cumufre$T.resident$V1) %>% as.data.frame(.)
cor.TresidentITCnecrosis$ID <- pat.inf$T
cor.TresidentITCnecrosis$responsive <- pat.inf$responsive

p.cor.TresidentITCnecrosis <- ggscatter(cor.TresidentITCnecrosis, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n", label.y = 40),
          xlab = 'Cumulative frequency of postICI_PBMC expanded clonotypes', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======



cor.test(CYT7G.Tis.G[pat.inf$surgery_tumor], ITCdynprepost.Class1.cumufre$expand$V1) # cor: 0.547794 p: 0.01518
cor.test(CYT7G.Tis.G[pat.inf$surgery_tumor], ITCdynprepost.Class1.cumufre$contract$V1)
cor.test(CYT7G.Tis.G[pat.inf$surgery_tumor], ITCdynprepost.Class1.cumufre$nondyna$V1)
cor.test(CYT7G.Tis.G[pat.inf$surgery_tumor], ITCdynprepost.Class1.cumufre$T.resident$V1) # cor: -0.7701962 p: 0.0001145


cor.postexpandITCCYT <- cbind(CYT = CYT7G.Tis.G[pat.inf$surgery_tumor], postExpCumuFre = ITCdynprepost.Class1.cumufre$expand$V1) %>% as.data.frame(.)
cor.postexpandITCCYT$ID <- pat.inf$T
cor.postexpandITCCYT$responsive <- pat.inf$responsive

p.cor.postexpandITCCYT <- ggscatter(cor.postexpandITCCYT, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n", label.y = -1),
          xlab = 'Cumulative frequency of postICI_PBMC expanded clonotypes', ylab = 'Cytolytic activity', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======


cor.TresidentITCCYT <- cbind(CYT = CYT7G.Tis.G[pat.inf$surgery_tumor], postExpCumuFre = ITCdynprepost.Class1.cumufre$T.resident$V1) %>% as.data.frame(.)
cor.TresidentITCCYT$ID <- pat.inf$T
cor.TresidentITCCYT$responsive <- pat.inf$responsive

p.cor.TresidentITCCYT <- ggscatter(cor.TresidentITCCYT, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency of postICI_PBMC expanded clonotypes', ylab = 'Cytolytic activity', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggarrange(p.cor.postexpandITCnecrosis, p.cor.TresidentITCnecrosis, 
          p.cor.postexpandITCCYT, p.cor.TresidentITCCYT, ncol = 2, nrow = 2)

ggsave(filename = "cor.cumufrepostexpandnecrosisCYT.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=12,units="in")
#######################################################################################






########## non-tumor classification: peripherally expanded, contracted, nondynamic, or tumor-resident

NTdyn.Class1 <- lapply(pat.inf$patientID, function(x){
  dyna.expandITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]]$CDR3.aa %in% c(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa),]
  dyna.contractITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]]$CDR3.aa %in% c(prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa),]
  
  allclone <- rbind(tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'pre']]],
                    tcr.immunarch$data[[pat.inf[pat.inf$patientID == x, 'post']]])
  
  dyna <- unique(c(prepost.FEST$expand[[x]]$CDR3.aa, prepost.FEST$contract[[x]]$CDR3.aa, postFU.FEST$expand[[x]]$CDR3.aa, postFU.FEST$contract[[x]]$CDR3.aa))
  nondyna <- allclone[is.na(match(allclone$CDR3.aa, dyna)),]
  
  nondynaITC <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]][tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]]$CDR3.aa %in% unique(nondyna$CDR3.aa),]
  
  T.resident <- tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]][is.na(match(tcr.immunarch$data[[pat.inf[pat.inf$patientID ==x, 'NT']]]$CDR3.aa, unique(allclone$CDR3.aa))),]
  
  dyna.ITC <- list(dyna.expandITC = dyna.expandITC, dyna.contractITC = dyna.contractITC, 
                   nondynaITC = nondynaITC, T.resident = T.resident)
  return(dyna.ITC)
})
names(NTdyn.Class1) <- pat.inf$patientID

NTdyn.Class1 <- list(expand = lapply(pat.inf$patientID, function(x){NTdyn.Class1[[x]]$dyna.expandITC}),
                             contract = lapply(pat.inf$patientID, function(x){NTdyn.Class1[[x]]$dyna.contractITC}),
                             nondyna = lapply(pat.inf$patientID, function(x){NTdyn.Class1[[x]]$nondynaITC}),
                             T.resident = lapply(pat.inf$patientID, function(x){NTdyn.Class1[[x]]$T.resident}))
names(NTdyn.Class1$expand) <- pat.inf$patientID
names(NTdyn.Class1$contract) <- pat.inf$patientID
names(NTdyn.Class1$nondyna) <- pat.inf$patientID
names(NTdyn.Class1$T.resident) <- pat.inf$patientID


NTdyn.Class1.cumufre <- lapply(NTdyn.Class1, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

# unlist(lapply(ITCdynprepost.Class1$expand,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 21    30    28    26    13    43     8    19     6     5    14    44    35    38    18    28   135   151   104 

unlist(lapply(NTdyn.Class1$expand,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 24    27    10    27    15    35    10    22    11    13    24    48    28    44    18    26    87    85    79 

unlist(lapply(NTdyn.Class1$contract,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 2     6     2     2    12    18    10    14     6    11    42     5     8     6    34    51    31     8    36

unlist(lapply(NTdyn.Class1$nondyna,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 872   321   130   348   504   840   162   817   167   703   408   292   497   428   399  1210  1664   754  1258 

unlist(lapply(NTdyn.Class1$T.resident,nrow))
# QM809 QM811 QM853 QM841 QM781 QM824 QM850 QM791 QM833 QM803 QM869 QM788 QM789 QM814 QM865 QM795 QM831 QM828 QM784 
# 1842   902   540   726  1576  1854   867  1903   600   954   762   830   855   584   635  2208  2059  1858  1445



################ figure plot ################ use
##### cumulative frequency test between NR and R： boxplot

NTdyn.Class1.cumufre.plot <- rbind(NTdyn.Class1.cumufre$expand, NTdyn.Class1.cumufre$contract,
                                   NTdyn.Class1.cumufre$nondyna, NTdyn.Class1.cumufre$T.resident)
NTdyn.Class1.cumufre.plot <- cbind(NTdyn.Class1.cumufre.plot, group = rep(c('expand','contract','nondyna','T.resident'), each = 19),
                                           response = rep(c(rep('NR',12),rep('R',7)), 4)) %>% as.data.frame(.)
NTdyn.Class1.cumufre.plot$group1 <- paste0(NTdyn.Class1.cumufre.plot$group,'.',NTdyn.Class1.cumufre.plot$response)
NTdyn.Class1.cumufre.plot$V1 <- as.numeric(NTdyn.Class1.cumufre.plot$V1)

p.NTdyn.Class1.cumufre <- ggboxplot(NTdyn.Class1.cumufre.plot, x = 'group1', y = 'V1', add = 'dotplot', ylab = 'Cumulative frequency', xlab = '',
                                            add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use

ggsave(filename = "postexpandinNT.boxplot_sup.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")
#######################################################################################





top1ITCprop <- lapply(pat.inf$patientID, function(x){
  expand <- length(intersect(unique(ITCdynprepost.Class1$expand[[x]]$CDR3.aa),unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa))
  contract <- length(intersect(unique(ITCdynprepost.Class1$contract[[x]]$CDR3.aa),unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa))
  nondyna <- length(intersect(unique(ITCdynprepost.Class1$nondyna[[x]]$CDR3.aa),unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa))
  T.resident <- length(intersect(unique(ITCdynprepost.Class1$T.resident[[x]]$CDR3.aa),unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa))
  index <- c(expand= expand, contract= contract, nondyna= nondyna, T.resident= T.resident)
  return(index)
})
names(top1ITCprop) <- pat.inf$patientID
top1ITCprop <- do.call(rbind,top1ITCprop)

top1ITCprop.plot <- reshape2::melt(top1ITCprop, value.name = 'prop')
top1ITCprop.plot$response <- rep(c(rep('NR',12),rep('R',7)), 4)
top1ITCprop.plot$group1 <- paste0(top1ITCprop.plot$Var2,'.',top1ITCprop.plot$response)


p.top1ITCprop <- ggboxplot(top1ITCprop.plot, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion among top 1% post-ITCs', xlab = '',
                                    add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use





top1ITCCF <- lapply(pat.inf$patientID, function(x){
  expand <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa %in% 
                                                                                 ITCdynprepost.Class1$expand[[x]]$CDR3.aa,]
  contract <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa %in% 
                                                                                   ITCdynprepost.Class1$contract[[x]]$CDR3.aa,]
  nondyna <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa %in% 
                                                                                  ITCdynprepost.Class1$nondyna[[x]]$CDR3.aa,]
  T.resident <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa %in% 
                                                                                     ITCdynprepost.Class1$T.resident[[x]]$CDR3.aa,]
  index <- list(expand= expand, contract= contract, nondyna= nondyna, T.resident= T.resident)
  return(index)
})
names(top1ITCCF) <- pat.inf$patientID


top1ITCCF <- list(expand = lapply(pat.inf$patientID, function(x){top1ITCCF[[x]]$expand}),
                             contract = lapply(pat.inf$patientID, function(x){top1ITCCF[[x]]$contract}),
                             nondyna = lapply(pat.inf$patientID, function(x){top1ITCCF[[x]]$nondyna}),
                             T.resident = lapply(pat.inf$patientID, function(x){top1ITCCF[[x]]$T.resident}))
names(top1ITCCF$expand) <- pat.inf$patientID
names(top1ITCCF$contract) <- pat.inf$patientID
names(top1ITCCF$nondyna) <- pat.inf$patientID
names(top1ITCCF$T.resident) <- pat.inf$patientID

top1ITCCF.cumufre <- lapply(top1ITCCF, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})


top1ITCCF.cumufre.plot <- rbind(top1ITCCF.cumufre$expand, top1ITCCF.cumufre$contract,
                                top1ITCCF.cumufre$nondyna, top1ITCCF.cumufre$T.resident)
top1ITCCF.cumufre.plot <- cbind(top1ITCCF.cumufre.plot, group = rep(c('expand','contract','nondyna','T.resident'), each = 19),
                                           response = rep(c(rep('NR',12),rep('R',7)), 4)) %>% as.data.frame(.)
top1ITCCF.cumufre.plot$group1 <- paste0(top1ITCCF.cumufre.plot$group,'.',top1ITCCF.cumufre.plot$response)
top1ITCCF.cumufre.plot$V1 <- as.numeric(top1ITCCF.cumufre.plot$V1)

p.top1ITCCF.cumufre <- ggboxplot(top1ITCCF.cumufre.plot, x = 'group1', y = 'V1', add = 'dotplot', ylab = 'Cumulative frequency', xlab = '',
                                            add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use


ggarrange(p.top1ITCprop, p.top1ITCCF.cumufre, common.legend = TRUE)

ggsave(filename = "postexpandtop1ITC.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")






cor.test(pat.inf$necrosis, top1ITCCF.cumufre$expand$V1) # cor: 0.739942 p: 0.0002927
cor.test(pat.inf$necrosis, top1ITCCF.cumufre$contract$V1) # ns
cor.test(pat.inf$necrosis, top1ITCCF.cumufre$nondyna$V1) # ns
cor.test(pat.inf$necrosis, top1ITCCF.cumufre$T.resident$V1) # cor: -0.3387335 p: 0.156 ns

cor.postexpandtop1ITCnecrosis <- cbind(CYT = pat.inf$necrosis, postExpCumuFre = top1ITCCF.cumufre$expand$V1) %>% as.data.frame(.)
cor.postexpandtop1ITCnecrosis$ID <- pat.inf$T
cor.postexpandtop1ITCnecrosis$responsive <- pat.inf$responsive

p.cor.postexpandtop1ITCnecrosis <- ggscatter(cor.postexpandtop1ITCnecrosis, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                                         cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
                                         xlab = 'Cumulative frequency of postICI.PBMC expanded clonotypes in top 1% post-ITCs', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======


cor.Tresidenttop1ITCnecrosis <- cbind(CYT = pat.inf$necrosis, postExpCumuFre = top1ITCCF.cumufre$T.resident$V1) %>% as.data.frame(.)
cor.Tresidenttop1ITCnecrosis$ID <- pat.inf$T
cor.Tresidenttop1ITCnecrosis$responsive <- pat.inf$responsive

p.cor.Tresidenttop1ITCnecrosis <- ggscatter(cor.Tresidenttop1ITCnecrosis, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                                        cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.12, label.sep = "\n"),
                                        xlab = 'Cumulative frequency of tumor-resident clonotypes', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======




cor.postcontracttop1ITCnecrosis <- cbind(CYT = pat.inf$necrosis, postExpCumuFre = top1ITCCF.cumufre$contract$V1) %>% as.data.frame(.)
cor.postcontracttop1ITCnecrosis$ID <- pat.inf$T
cor.postcontracttop1ITCnecrosis$responsive <- pat.inf$responsive

p.cor.postcontracttop1ITCnecrosis <- ggscatter(cor.postcontracttop1ITCnecrosis, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                                             cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.03, label.sep = "\n"),
                                             xlab = 'Cumulative frequency of postICI.PBMC expanded clonotypes in top 1% post-ITCs', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======


cor.nondyntop1ITCnecrosis <- cbind(CYT = pat.inf$necrosis, postExpCumuFre = top1ITCCF.cumufre$nondyna$V1) %>% as.data.frame(.)
cor.nondyntop1ITCnecrosis$ID <- pat.inf$T
cor.nondyntop1ITCnecrosis$responsive <- pat.inf$responsive

p.cor.nondyntop1ITCnecrosis <- ggscatter(cor.nondyntop1ITCnecrosis, 'postExpCumuFre', 'CYT',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                                            cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.12, label.sep = "\n"),
                                            xlab = 'Cumulative frequency of tumor-resident clonotypes', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'ID') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(0,5)) +
  # scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggarrange(p.cor.postexpandtop1ITCnecrosis, p.cor.postcontracttop1ITCnecrosis, 
          p.cor.nondyntop1ITCnecrosis, p.cor.Tresidenttop1ITCnecrosis)


ggsave(filename = "cor.cumufrepostexpandtop1ITCnecrosis1.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=12,units="in")



###### ###### ###### ###### ###### ###### 
###### plotROC
rocdata <- top1ITCCF.cumufre$expand
rocdata$response <- c(rep('nivo-NR',12), rep('nivo-R',7))
# rocdata$status <- c(rep('1',12), rep('0',7))

rocexp <- roc(rocdata$response, rocdata$V1) # Area under the curve: 0.8155

ggroc(rocexp, alpha = 0.5, colour = "black", legacy.axes = TRUE)+ 
  # theme_minimal() + 
  theme_classic()+
  theme(aspect.ratio = 0.8, plot.title = element_text(size = 7), 
        text = element_text(size = 7), axis.text = element_text(size = 7),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 5),
        axis.line = element_line(size = 0.5))+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed")

ggsave(filename = 'D:/HKU/Nivo_ICItherapy/writepaper/final/manuscript 3/Figure/postExpTop1ITC.ROCplot.pdf',
      width = 7, height = 8, units = 'cm')





