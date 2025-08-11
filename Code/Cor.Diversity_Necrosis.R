#########################################

testtar <- c('reads', 'diversity', 'observedDiversity_mean', 'efronThisted_mean', 'chao1_mean', 
             'd50Index_mean', 'shannonWienerIndex_mean', 'normalizedShannonWienerIndex_mean',
             'inverseSimpsonIndex_mean','clonality')


diversity.cor.T <- lapply(testtar, function(x){
  estimate <- cor.test(diversity[diversity$group == 'T',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(diversity[diversity$group == 'T',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})


diversity.cor.T <- do.call(rbind, diversity.cor.T)
colnames(diversity.cor.T) <- c('Cor', 'pvalue')
rownames(diversity.cor.T) <- testtar

#                                         Cor      pvalue
# reads                              0.4107314 0.080665764
# diversity                          0.2981541 0.215037699
# observedDiversity_mean             0.2981541 0.215037699
# efronThisted_mean                  0.3755510 0.113081028
# chao1_mean                         0.3614788 0.128349730
# d50Index_mean                      0.3605993 0.129349982
# shannonWienerIndex_mean           -0.3676354 0.121500462
# normalizedShannonWienerIndex_mean -0.6473198 0.002734673
# inverseSimpsonIndex_mean          -0.5655253 0.011618913
# clonality                          0.6473198 0.002734673 ###### use

########################## use ########################## 
cor.clonalitynecrosis <- cbind(pat.inf, clonality = diversity[diversity$group == 'T','clonality'])

ggscatter(cor.clonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 0.1, label.sep = "\n"),
          xlab = 'TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggsave(filename = "cor.clonalitynecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")

ggscatter(cor.clonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 0.1, label.sep = "\n"),
          xlab = 'TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'T') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggsave(filename = "cor.clonalitynecrosis_pearson.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")

########################## use ########################## 



########################## 

diversity.cor.B <- lapply(testtar, function(x){
  estimate <- cor.test(diversity[diversity$group == 'B',x], pat.inf[pat.inf$B %in% rownames(diversity),'necrosis'], alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(diversity[diversity$group == 'B',x], pat.inf[pat.inf$B %in% rownames(diversity),'necrosis'], alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})


diversity.cor.B <- do.call(rbind, diversity.cor.B)
colnames(diversity.cor.B) <- c('Cor', 'pvalue')
rownames(diversity.cor.B) <- testtar

cor.Bclonalitynecrosis <- cbind(necrosis = pat.inf[pat.inf$B %in% intersect(diversity$sample_id, pat.inf$B), 'necrosis'], 
                                clonality = diversity[diversity$group == 'B','clonality']) %>% as.data.frame(.)
cor.Bclonalitynecrosis$label <- intersect(diversity$sample_id, pat.inf$B)

p.corBnecrosis <-ggscatter(cor.Bclonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                           cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 0.1, label.sep = "\n"),
                           xlab = 'TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'label') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

################


########################## 

diversity.cor.NT <- lapply(testtar, function(x){
  estimate <- cor.test(diversity[diversity$group == 'NT',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(diversity[diversity$group == 'NT',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})

diversity.cor.NT <- do.call(rbind, diversity.cor.NT)
colnames(diversity.cor.NT) <- c('Cor', 'pvalue')
rownames(diversity.cor.NT) <- testtar


cor.NTclonalitynecrosis <- cbind(pat.inf, clonality = diversity[diversity$group == 'NT','clonality'])

p.corNTnecrosis <-ggscatter(cor.NTclonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                            cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 0.1, label.sep = "\n"),
                            xlab = 'TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'NT') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggarrange(p.corBnecrosis, p.corNTnecrosis)

ggsave(filename = "cor.BnNTclonalitynecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=16,height=6,units="in")




########################## 

diversity.cor.pre <- lapply(testtar, function(x){
  estimate <- cor.test(diversity[diversity$group == 'pre',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(diversity[diversity$group == 'pre',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})

diversity.cor.pre <- do.call(rbind, diversity.cor.pre)
colnames(diversity.cor.pre) <- c('Cor', 'pvalue')
rownames(diversity.cor.pre) <- testtar


cor.preclonalitynecrosis <- cbind(pat.inf, clonality = diversity[diversity$group == 'pre','clonality'])

p.corprenecrosis <-ggscatter(cor.preclonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                             cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"),
                             xlab = 'preICI.PBMC TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'pre') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======




########################## 


diversity.cor.post <- lapply(testtar, function(x){
  estimate <- cor.test(diversity[diversity$group == 'post',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(diversity[diversity$group == 'post',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})

diversity.cor.post <- do.call(rbind, diversity.cor.post)
colnames(diversity.cor.post) <- c('Cor', 'pvalue')
rownames(diversity.cor.post) <- testtar



cor.postclonalitynecrosis <- cbind(pat.inf, clonality = diversity[diversity$group == 'post','clonality'])

p.corpostnecrosis <-ggscatter(cor.postclonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                              cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 0.2, label.sep = "\n"),
                              xlab = 'postICI.PBMC TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'post') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======



########################## 


diversity.cor.FU <- lapply(testtar, function(x){
  estimate <- cor.test(diversity[diversity$group == 'FU',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(diversity[diversity$group == 'FU',x], pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})

diversity.cor.FU <- do.call(rbind, diversity.cor.FU)
colnames(diversity.cor.FU) <- c('Cor', 'pvalue')
rownames(diversity.cor.FU) <- testtar


cor.FUclonalitynecrosis <- cbind(pat.inf, clonality = diversity[diversity$group == 'FU','clonality'])

p.corFUnecrosis <-ggscatter(cor.FUclonalitynecrosis, 'clonality', 'necrosis',  add = 'reg.line', add.params = list(fill = 'lightgray'),
                            cor.coef = TRUE, cor.coeff.args = list(method = "spearman", label.x = 0.25, label.sep = "\n"),
                            xlab = 'FU.PBMC TCR repertoire clonality', ylab = 'Tumor necrosis', conf.int = TRUE, label = 'FU') +
  # stat_cor(aes(color = responsive), label.x = c(0.35,0.35), label.y = c(40,75)) +
  scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
  theme(aspect.ratio = 0.66) ######====== use ======

ggarrange(p.corprenecrosis, p.corpostnecrosis, p.corFUnecrosis, ncol = 3)

ggsave(filename = "cor.PBMCclonalitynecrosis.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig/use/section1",
       width=16,height=8,units="in")

########################## 


diversity.cor.postpre <- lapply(testtar, function(x){
  estimate <- cor.test(c(diversity[diversity$group == 'post',x]-diversity[diversity$group == 'pre',x]), pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$estimate
  pvalue <- cor.test(c(diversity[diversity$group == 'post',x]-diversity[diversity$group == 'pre',x]), pat.inf$necrosis, alternative = 'two.sided', method = 'spearman')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})

diversity.cor.postpre <- do.call(rbind, diversity.cor.postpre)
colnames(diversity.cor.postpre) <- c('Cor', 'pvalue')
rownames(diversity.cor.postpre) <- testtar




diversity.cor.FUpost <- lapply(testtar, function(x){
  estimate <- cor.test(c(diversity[diversity$group == 'FU',x]-diversity[diversity$group == 'post',x]), pat.inf$necrosis, alternative = 'two.sided', method = 'pearson')$estimate
  pvalue <- cor.test(c(diversity[diversity$group == 'FU',x]-diversity[diversity$group == 'post',x]), pat.inf$necrosis, alternative = 'two.sided', method = 'pearson')$p.value
  testres <- c(estimate, pvalue)
  return(testres)
})

diversity.cor.FUpost <- do.call(rbind, diversity.cor.FUpost)
colnames(diversity.cor.FUpost) <- c('Cor', 'pvalue')
rownames(diversity.cor.FUpost) <- testtar



