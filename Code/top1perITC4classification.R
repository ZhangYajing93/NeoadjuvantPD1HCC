############## top 1% clonotypes classification ：4 mutually exclusive subsets

# Top 1% ITCs were then categorized into four mutually exclusive subsets: 
# ITCs shared in the pretreatment blood and the resected normal lung; 
# ITCs found only in the pretreatment blood; 
# ITCs found only in the resected normal lung; 
# and tumor-resident ITCs (not found in the resected normal lung or the pretreatment blood). ref: PMID: 31754049

#### load data
tcr.immunarch.top1per <- readRDS(file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/tcr.immunarch.top1per.rds')


# 1. top 1% ITCs shared in the posttreatment blood and the resected normal lung; 

top1perITC.share.postNT <- lapply(pat.inf$T, function(x){
  shareITCpostNT <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                              intersect(tcr.immunarch$data[[pat.inf[pat.inf$T == x,'post']]]$CDR3.aa,
                                        tcr.immunarch$data[[pat.inf[pat.inf$T == x,'NT']]]$CDR3.aa))
  shareITCpostNT <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% shareITCpostNT, ]
  
  return(shareITCpostNT)
}) 
names(top1perITC.share.postNT) <- pat.inf$T

top1perITC.share.postNT.prop <- unlist(lapply(top1perITC.share.postNT,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#      T809      T811      T853      T841      T781      T824      T850      T791      T833      T803      T869      T788      T789      T814      T865      T795      T831      T828      T784 
# 0.4210526 0.3870968 0.2631579 0.2400000 0.4500000 0.6712329 0.4285714 0.7142857 1.0000000 0.0000000 0.3750000 0.5000000 0.6739130 0.7000000 0.4705882 0.9090909 0.7966102 0.7258065 0.8139535 

top1perITC.share.postNT.cumufre <- lapply(top1perITC.share.postNT, function(x){
  sum(x$Proportion)
}) %>% unlist(.)

cor.test(top1perITC.share.postNT.cumufre, pat.inf$necrosis, method = 'pearson') # cor: 0.706766 p=0.0007165
cor.test(top1perITC.share.postNT.prop, pat.inf$necrosis, method = 'pearson')  # cor: 0.5715189  p=0.01058
cor.test(top1perITC.share.postNT.prop, CYT7G.Tis.G[20:38], method = 'pearson') # cor: 0.6221399   p=0.004451
cor.test(top1perITC.share.postNT.cumufre, CYT7G.Tis.G[20:38], method = 'pearson') # cor: 0.6772151   p=0.001447


cor.top1perITCinpostNTnecrosis <- cbind(pat.inf, cumufre = top1perITC.share.postNT.cumufre, prop = top1perITC.share.postNT.prop)

ggscatter(cor.top1perITCinpostNTnecrosis, 'cumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(-50,125)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggsave(filename = "cor.top1perITCinpostNTnecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")






# 2. top 1% ITCs found only in the posttreatment blood

top1perITC.post <- lapply(pat.inf$T, function(x){
  shareITCpost <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                            tcr.immunarch$data[[pat.inf[pat.inf$T == x,'post']]]$CDR3.aa)
  
  shareITCpost <- shareITCpost[is.na(match(shareITCpost, top1perITC.share.postNT[[x]]$CDR3.aa))]
  shareITCpost <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% shareITCpost, ]
  return(shareITCpost)
}) 
names(top1perITC.post) <- pat.inf$T

top1perITC.post.prop <- unlist(lapply(top1perITC.post,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#       T809       T811       T853       T841       T781       T824       T850       T791       T833       T803       T869       T788       T789       T814       T865       T795       T831       T828 
# 0.21052632 0.16129032 0.52631579 0.04000000 0.05000000 0.16438356 0.00000000 0.04761905 0.00000000 0.14285714 0.00000000 0.09090909 0.15217391 0.05000000 0.17647059 0.03030303 0.00000000 0.08064516 
#       T784 
# 0.09302326 

top1perITC.post.cumufre <- lapply(top1perITC.post, function(x){
  sum(x$Proportion)
}) %>% unlist(.)

cor.test(top1perITC.post.cumufre, pat.inf$necrosis, method = 'pearson') # cor: -0.2026126  p=0.4055


# 3. top 1% ITCs found only in the resected normal lung

top1perITC.NT <- lapply(pat.inf$T, function(x){
  shareITCnT <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                            tcr.immunarch$data[[pat.inf[pat.inf$T == x,'NT']]]$CDR3.aa)
  
  shareITCnT <- shareITCnT[is.na(match(shareITCnT, top1perITC.share.postNT[[x]]$CDR3.aa))]
  shareITCnT <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% shareITCnT, ]
  return(shareITCnT)
}) 
names(top1perITC.NT) <- pat.inf$T

top1perITC.NT.prop <- unlist(lapply(top1perITC.NT,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#       T809       T811       T853       T841       T781       T824       T850       T791       T833       T803       T869       T788       T789       T814       T865       T795       T831       T828 
# 0.05263158 0.12903226 0.10526316 0.28000000 0.25000000 0.06849315 0.28571429 0.23809524 0.00000000 0.57142857 0.25000000 0.22727273 0.00000000 0.20000000 0.23529412 0.03030303 0.20338983 0.09677419 
#       T784 
# 0.06976744 

top1perITC.NT.cumufre <- lapply(top1perITC.NT, function(x){
  sum(x$Proportion)
}) %>% unlist(.)

cor.test(top1perITC.NT.cumufre, pat.inf$necrosis, method = 'pearson') # cor: -0.03701137   p=0.8804


# 4. tumor-resident ITCs (not found in the resected normal lung or the posttreatment blood)

top1perITC.unique <- lapply(pat.inf$T, function(x){
  top1perunique <- tcr.immunarch.top1per$data[[x]]$CDR3.aa[is.na(match(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                                                                       c(top1perITC.share.postNT[[x]]$CDR3.aa, top1perITC.post[[x]]$CDR3.aa, top1perITC.NT[[x]]$CDR3.aa)))]
  top1perunique <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% top1perunique, ]
  return(top1perunique)
})

names(top1perITC.unique) <- pat.inf$T

top1perITC.unique.prop <- unlist(lapply(top1perITC.unique,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#       T809       T811       T853       T841       T781       T824       T850       T791       T833       T803       T869       T788       T789       T814       T865       T795       T831       T828 
# 0.31578947 0.32258065 0.10526316 0.44000000 0.25000000 0.09589041 0.28571429 0.00000000 0.00000000 0.28571429 0.37500000 0.18181818 0.17391304 0.05000000 0.11764706 0.03030303 0.00000000 0.09677419 
#       T784 
# 0.02325581 


top1perITC.unique.cumufre <- lapply(top1perITC.unique, function(x){
  sum(x$Proportion)
}) %>% unlist(.)


cor.test(top1perITC.unique.cumufre, pat.inf$necrosis, method = 'pearson') # cor: -0.588556  p=0.008027


cor.top1perITCuniquenecrosis <- cbind(pat.inf, cumufre = top1perITC.unique.cumufre, prop = top1perITC.unique.prop)

ggscatter(cor.top1perITCuniquenecrosis, 'cumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.05, label.sep = "\n"),
          xlab = 'Cumulative frequency', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(-45,125)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggsave(filename = "cor.top1perITCuniquenecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")


ggscatter(cor.top1perITCinpostNTnecrosis, 'cumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#FF4B20',
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency in postNT', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(-50,125)) +
  theme(aspect.ratio = 0.66) +
  ggscatter(cor.top1perITCuniquenecrosis, 'cumufre', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#0348A6',
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.05, label.sep = "\n"),
            xlab = 'Cumulative frequency unique', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(-45,125)) +
  theme(aspect.ratio = 0.66) +

ggscatter(cor.top1perITCinpostNTnecrosis, 'prop', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#FF4B20',
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Proportion in postNT', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(-50,125)) +
  theme(aspect.ratio = 0.66) +
  ggscatter(cor.top1perITCuniquenecrosis, 'prop', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#0348A6',
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.3, label.sep = "\n"),
            xlab = 'Proportion unique', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(-45,125)) +
  theme(aspect.ratio = 0.66)


ggsave(filename = "cor.top1perITC4classnecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=12,units="in")


cor.top1perITCinpostNTnecrosis$CYT.T <- CYT7G.Tis.G[20:38]
cor.top1perITCuniquenecrosis$CYT.T <- CYT7G.Tis.G[20:38]

ggscatter(cor.top1perITCinpostNTnecrosis, 'cumufre', 'CYT.T',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#FF4B20',
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'Cumulative frequency in postNT', ylab = 'Cytolytic score', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(-6,-3,0,3,6,9), limits = c(-6,9)) +
  theme(aspect.ratio = 0.66) +
  ggscatter(cor.top1perITCuniquenecrosis, 'cumufre', 'CYT.T',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#0348A6',
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.05, label.sep = "\n"),
            xlab = 'Cumulative frequency unique', ylab = 'Cytolytic score', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(-6,-3,0,3,6,9), limits = c(-6,9)) +
  theme(aspect.ratio = 0.66) +
  
  ggscatter(cor.top1perITCinpostNTnecrosis, 'prop', 'CYT.T',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#FF4B20',
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
            xlab = 'Proportion in postNT', ylab = 'Cytolytic score', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(-6,-3,0,3,6,9), limits = c(-6,9)) +
  theme(aspect.ratio = 0.66) +
  ggscatter(cor.top1perITCuniquenecrosis, 'prop', 'CYT.T',  add = 'reg.line', add.params = list(fill = 'lightgray'), color = '#0348A6',
            cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.3, label.sep = "\n"),
            xlab = 'Proportion unique', ylab = 'Cytolytic score', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(-6,-3,0,3,6,9), limits = c(-6,9)) +
  theme(aspect.ratio = 0.66)


ggsave(filename = "cor.top1perITC4classCYT.T.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=12,units="in")

#####  barplot

ITC.class.plot <- cbind(unique = top1perITC.unique.prop, NT = top1perITC.NT.prop, post = top1perITC.post.prop, postNT = top1perITC.share.postNT.prop)

pdf(file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig/top1perITCsharePropBarplot.pdf',14,7)
barplot(t(ITC.class.plot),col = c('#0348A6', '#7AC5FF','#FFB433','#FF4B20'),
        legend=rownames(t(ITC.class.plot)),las=2, ylab = 'Proportion among top 1% ITCs',
        cex.names = 0.8,args.legend = list(x=3,y=1,cex=0.6))
dev.off()


##### proportion test between NR and R： boxplot

top1perITCshareProp <- c(top1perITC.share.postNT.prop, top1perITC.post.prop, top1perITC.NT.prop, top1perITC.unique.prop)

top1perITCshareProp <- cbind(prop = top1perITCshareProp, group = rep(c('postNT','post','NT','unique'), each = 19),
                             response = rep(c(rep('NR',12),rep('R',7)), 4)) %>% as.data.frame(.)
top1perITCshareProp$group1 <- paste0(top1perITCshareProp$group,'.',top1perITCshareProp$response)
top1perITCshareProp$prop <- as.numeric(top1perITCshareProp$prop)

p.top1ITCprop <- ggboxplot(top1perITCshareProp, x = 'group1', y = 'prop', add = 'dotplot', ylab = 'Proportion of top 1% ITCs', xlab = '',
                              add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8),c(9,10))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use




##### cumulative frequency test between NR and R： boxplot

top1perITCshareCumufre <- c(top1perITC.share.postNT.cumufre, top1perITC.post.cumufre, top1perITC.NT.cumufre, top1perITC.unique.cumufre)

top1perITCshareCumufre <- cbind(cumFre = top1perITCshareCumufre, group = rep(c('postNT','post','NT','unique'), each = 19),
                                response = rep(c(rep('NR',12),rep('R',7)), 4)) %>% as.data.frame(.)
top1perITCshareCumufre$group1 <- paste0(top1perITCshareCumufre$group,'.',top1perITCshareCumufre$response)
top1perITCshareCumufre$cumFre <- as.numeric(top1perITCshareCumufre$cumFre)


p.top1ITCcumufre <- ggboxplot(top1perITCshareCumufre, x = 'group1', y = 'cumFre', add = 'dotplot', ylab = 'Cumulative frequency', xlab = '',
                    add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0,
                     comparisons = list(c(1,2),c(3,4), c(5,6),
                                        c(7,8),c(9,10))) +
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30))# use


p.top1ITCprop + p.top1ITCcumufre

ggsave(filename = "top1perITC.classtest.boxplot_edit.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=6,units="in")




