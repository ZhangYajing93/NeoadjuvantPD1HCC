# BT.FEST
# ITCdynprepost.Class1


unlist(lapply(BT.FEST$expand, nrow)) # 472
# QM841 QM781 QM824 QM791 QM803 QM788 QM789 QM814 QM795 QM784 
# 23    97    38    84    48    12     9    21   103    37 

unlist(lapply(BT.FEST$contract, nrow)) # 608 
# QM841 QM781 QM824 QM791 QM803 QM788 QM789 QM814 QM795 QM784 
# 35    21   111    45     7    83    92    45    78    91 



same.dyn <- lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){
  expand <- intersect(BT.FEST$expand[[x]]$CDR3.aa, c(prepost.FEST$expand[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa))
  contract <- intersect(BT.FEST$contract[[x]]$CDR3.aa, c(prepost.FEST$contract[[x]]$CDR3.aa,postFU.FEST$contract[[x]]$CDR3.aa))
  index <- list(expand=expand, contract=contract)
  return(index)
})
names(same.dyn) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]

same.dyn <- list(expand = lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){same.dyn[[x]]$expand}),
                 contract = lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){same.dyn[[x]]$contract}))
names(same.dyn$expand) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]
names(same.dyn$contract) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]





##### B to T dynamic clonotype: top 1% post-ITC
BT.top1ITC <- lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){
  expand <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa %in% 
                                                                                 BT.FEST$expand[[x]]$CDR3.aa,]
  contract <- tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]][tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa %in% 
                                                                                 BT.FEST$contract[[x]]$CDR3.aa,]
  index <- list(expand= expand, contract= contract)
})
names(BT.top1ITC) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]


BT.top1ITC <- list(expand = lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){BT.top1ITC[[x]]$expand}),
                   contract = lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){BT.top1ITC[[x]]$contract}))
names(BT.top1ITC$expand) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]
names(BT.top1ITC$contract) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]


unlist(lapply(BT.top1ITC$expand, nrow))
# QM841 QM781 QM824 QM791 QM803 QM788 QM789 QM814 QM795 QM784 
# 5    13    28    12     5     8     6    12    26    21 
unlist(lapply(BT.top1ITC$contract, nrow))
# QM841 QM781 QM824 QM791 QM803 QM788 QM789 QM814 QM795 QM784 
# 2     0     5     5     0     3     2     1     2     0 



BT.top1ITC.prop <- lapply(pat.inf$patientID[c(4:6,8,10,12:14,16,19)], function(x){
  expand <- length(intersect(unique(BT.FEST$expand[[x]]$CDR3.aa),unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa))
  contract <- length(intersect(unique(BT.FEST$contract[[x]]$CDR3.aa),unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa)))/length(unique(tcr.immunarch.top1per$data[[pat.inf[pat.inf$patientID == x, 'T']]]$CDR3.aa))
  index <- c(expand= expand, contract= contract)
  return(index)
})
names(BT.top1ITC.prop) <- pat.inf$patientID[c(4:6,8,10,12:14,16,19)]
BT.top1ITC.prop <- do.call(rbind,BT.top1ITC.prop)

BT.top1ITC.prop.plot <- reshape2::melt(BT.top1ITC.prop, value.name = 'prop')
BT.top1ITC.prop.plot$response <- rep(c(rep('NR',6),rep('R',4)), 2)
BT.top1ITC.prop.plot$group1 <- paste0(BT.top1ITC.prop.plot$Var2,'.',BT.top1ITC.prop.plot$response)


p.BT.top1ITC.prop <- ggboxplot(BT.top1ITC.prop.plot, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion among top 1% post-ITCs', xlab = '',
                           add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4))) +
  theme(aspect.ratio = 0.66)# use
  # theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use





BT.top1ITC.CF <- lapply(BT.top1ITC, function(x){
  do.call(rbind, lapply(x, function(y){
    sum(y$Proportion)
  })) %>% as.data.frame(.)
})

BT.top1ITC.CF.plot <- rbind(BT.top1ITC.CF$expand, BT.top1ITC.CF$contract)
BT.top1ITC.CF.plot <- cbind(BT.top1ITC.CF.plot, group = rep(c('expand','contract'), each = 10),
                                response = rep(c(rep('NR',6),rep('R',4)), 2)) %>% as.data.frame(.)
BT.top1ITC.CF.plot$group1 <- paste0(BT.top1ITC.CF.plot$group,'.',BT.top1ITC.CF.plot$response)
BT.top1ITC.CF.plot$V1 <- as.numeric(BT.top1ITC.CF.plot$V1)

p.BT.top1ITC.CF <- ggboxplot(BT.top1ITC.CF.plot, x = 'group1', y = 'V1', add = 'dotplot', ylab = 'Cumulative frequency', xlab = '',
                                 add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use


ggarrange(p.BT.top1ITC.prop, p.BT.top1ITC.CF, common.legend = TRUE)

ggsave(filename = "BTexpandtop1ITC.boxplot.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")



BT.top1ITC.prop.barplot <- cbind(prop = c(BT.top1ITC.prop[,1], -BT.top1ITC.prop[,2]),
                                 patID = rep(pat.inf$patientID[c(4:6,8,10,12:14,16,19)],2), type = c(rep('expand', 10), rep('contract',10))) %>% as.data.frame(.)
BT.top1ITC.prop.barplot$prop <- as.numeric(BT.top1ITC.prop.barplot$prop)


p.BT.top1ITC.prop.barplot <- ggbarplot(BT.top1ITC.prop.barplot, 'patID', 'prop',
                        fill = 'type', color = 'type', palette = "Paired",
                        label = TRUE, lab.col = "black", lab.pos = "in") +
  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4,0.6), limits = c(-0.3,0.8))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 90)) +
  ylab('Proportion of dynamic intratumoral clonotypes among top 1% post-ITCs') +
  xlab('')


BT.top1ITC.CF.barplot <- cbind(CF = c(BT.top1ITC.CF$expand$V1, -BT.top1ITC.CF$contract$V1),
                                 patID = rep(pat.inf$patientID[c(4:6,8,10,12:14,16,19)],2), type = c(rep('expand', 10), rep('contract',10))) %>% as.data.frame(.)
BT.top1ITC.CF.barplot$CF <- as.numeric(BT.top1ITC.CF.barplot$CF)


p.BT.top1ITC.CF.barplot <- ggbarplot(BT.top1ITC.CF.barplot, 'patID', 'CF',
                                       fill = 'type', color = 'type', palette = "Paired",
                                       label = TRUE, lab.col = "black", lab.pos = "in") +
  scale_y_continuous(breaks = c(-0.1,0,0.1,0.2,0.3, 0.4), limits = c(-0.15,0.45))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 90)) +
  ylab('CF of dynamic intratumoral clonotypes among top 1% post-ITCs') +
  xlab('')

ggarrange(p.BT.top1ITC.prop.barplot, p.BT.top1ITC.CF.barplot, common.legend = TRUE)

ggsave(filename = "BTdynaintop1ITC.NumbernCF.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=6,units="in") # use



cor.test(pat.inf$necrosis[c(4:6,8,10,12:14,16,19)], BT.top1ITC.CF$expand$V1) # cor: 0.7343812 p: 0.01558
cor.test(pat.inf$necrosis[c(4:6,8,10,12:14,16,19)], BT.top1ITC.CF$contract$V1) # ns



