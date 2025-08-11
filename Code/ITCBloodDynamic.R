################## top 1% ITC: pre, post, FU overlap

###### pre
# top 1%
top1perITCinpre <- lapply(pat.inf$T, function(x){
  index <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                     tcr.immunarch$data[[pat.inf[pat.inf$T == x,'pre']]]$CDR3.aa)
  
  index <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% index, ]
  return(index)
}) 
names(top1perITCinpre) <- pat.inf$T

top1perITCinpre.prop <- unlist(lapply(top1perITCinpre,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#       T809       T811       T853       T841       T781       T824       T850       T791       T833       T803       T869       T788       T789       T814       T865       T795       T831       T828 
# 0.52631579 0.67741935 0.05263158 0.28000000 0.50000000 0.30136986 0.14285714 0.76190476 1.00000000 0.14285714 0.25000000 0.27272727 0.19565217 0.45000000 0.58823529 0.69696970 0.76271186 0.09677419 
#       T784 
# 0.39534884 

top1perITCinpre.cumufre <- lapply(top1perITCinpre, function(x){
  sum(x$Proportion)
}) %>% unlist(.)

top1perITCinpre.precumufre <- lapply(pat.inf$T, function(x){
  index <- tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'pre']]][tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'pre']]]$CDR3.aa %in% top1perITCinpre[[x]]$CDR3.aa,]
  print(index$Proportion)
  cumufre <- sum(index$Proportion)
  return(cumufre)
}) %>% unlist(.)
names(top1perITCinpre.precumufre) <- pat.inf$T


# non-top 1%
nontop1perITCinpre <- lapply(pat.inf$T, function(x){
  index <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                     tcr.immunarch$data[[pat.inf[pat.inf$T == x,'pre']]]$CDR3.aa)
  
  index <- tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% index, ]
  return(index)
}) 
names(nontop1perITCinpre) <- pat.inf$T

nontop1perITCinpre.prop <- unlist(lapply(nontop1perITCinpre,nrow))/unlist(lapply(tcr.immunarch.nontop1per$data[pat.inf$T], nrow))
#       T809       T811       T853       T841       T781       T824       T850       T791       T833       T803       T869       T788       T789       T814       T865       T795       T831       T828 
# 0.22048067 0.23145695 0.19683377 0.07892596 0.25952626 0.18923754 0.12518854 0.27206596 0.23343849 0.30000000 0.16822430 0.13266697 0.13763676 0.16273101 0.24393064 0.27281134 0.24470102 0.10935720 
#       T784 
# 0.24090485 
## 与8.top1perStattest.R 中 tcr.immunarch.nontop1per.T.share.prop$preT 相同

nontop1perITCinpre.cumufre <- lapply(nontop1perITCinpre, function(x){
  sum(x$Proportion)
}) %>% unlist(.)


nontop1perITCinpre.precumufre <- lapply(pat.inf$T, function(x){
  index <- tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'pre']]][tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'pre']]]$CDR3.aa %in% nontop1perITCinpre[[x]]$CDR3.aa,]
  print(index$Proportion)
  cumufre <- sum(index$Proportion)
  return(cumufre)
}) %>% unlist(.)
names(nontop1perITCinpre.precumufre) <- pat.inf$T



###### post
# top 1%
top1perITCinpost <- lapply(pat.inf$T, function(x){
  index <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                     tcr.immunarch$data[[pat.inf[pat.inf$T == x,'post']]]$CDR3.aa)
  
  index <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% index, ]
  return(index)
}) 
names(top1perITCinpost) <- pat.inf$T

top1perITCinpost.prop <- unlist(lapply(top1perITCinpost,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#      T809      T811      T853      T841      T781      T824      T850      T791      T833      T803      T869      T788      T789      T814      T865      T795      T831      T828      T784 
# 0.6315789 0.5483871 0.7894737 0.2800000 0.5000000 0.8356164 0.4285714 0.7619048 1.0000000 0.1428571 0.3750000 0.5909091 0.8260870 0.7500000 0.6470588 0.9393939 0.7966102 0.8064516 0.9069767 

top1perITCinpost.cumufre <- lapply(top1perITCinpost, function(x){
  sum(x$Proportion)
}) %>% unlist(.)


top1perITCinpost.postcumufre <- lapply(pat.inf$T, function(x){
  index <- tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'post']]]$CDR3.aa %in% top1perITCinpost[[x]]$CDR3.aa,]
  # print(index$Proportion)
  cumufre <- sum(index$Proportion)
  return(cumufre)
}) %>% unlist(.)
names(top1perITCinpost.postcumufre) <- pat.inf$T


# non-top 1%
nontop1perITCinpost <- lapply(pat.inf$T, function(x){
  index <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                     tcr.immunarch$data[[pat.inf[pat.inf$T == x,'post']]]$CDR3.aa)
  
  index <- tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% index, ]
  return(index)
}) 
names(nontop1perITCinpost) <- pat.inf$T

nontop1perITCinpost.prop <- unlist(lapply(nontop1perITCinpost,nrow))/unlist(lapply(tcr.immunarch.nontop1per$data[pat.inf$T], nrow))
#      T809      T811      T853      T841      T781      T824      T850      T791      T833      T803      T869      T788      T789      T814      T865      T795      T831      T828      T784 
# 0.2669801 0.1927152 0.2163588 0.2286412 0.2749743 0.2453895 0.1538462 0.3026188 0.2618297 0.3589041 0.1922563 0.1431168 0.1525164 0.2320329 0.3057803 0.2829840 0.3050146 0.2891782 0.3204291 

nontop1perITCinpost.cumufre <- lapply(nontop1perITCinpost, function(x){
  sum(x$Proportion)
}) %>% unlist(.)

nontop1perITCinpost.postcumufre <- lapply(pat.inf$T, function(x){
  index <- tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'post']]][tcr.immunarch$data[[pat.inf[pat.inf$T == x, 'post']]]$CDR3.aa %in% nontop1perITCinpost[[x]]$CDR3.aa,]
  # print(index$Proportion)
  cumufre <- sum(index$Proportion)
  return(cumufre)
}) %>% unlist(.)
names(nontop1perITCinpost.postcumufre) <- pat.inf$T

###### FU
# top 1%
top1perITCinFU <- lapply(pat.inf$T, function(x){
  index <- intersect(tcr.immunarch.top1per$data[[x]]$CDR3.aa,
                     tcr.immunarch$data[[pat.inf[pat.inf$T == x,'FU']]]$CDR3.aa)
  
  index <- tcr.immunarch.top1per$data[[x]][tcr.immunarch.top1per$data[[x]]$CDR3.aa %in% index, ]
  return(index)
}) 
names(top1perITCinFU) <- pat.inf$T

top1perITCinFU.prop <- unlist(lapply(top1perITCinFU,nrow))/unlist(lapply(tcr.immunarch.top1per$data[pat.inf$T], nrow))
#      T809      T811      T853      T841      T781      T824      T850      T791      T833      T803      T869      T788      T789      T814      T865      T795      T831      T828      T784 
# 0.5263158 0.5161290 0.2105263 0.2400000 0.6000000 0.4931507 0.0000000 0.2380952 1.0000000 0.0000000 0.3750000 0.4090909 0.3913043 0.8000000 0.4117647 0.7272727 0.7966102 0.7741935 0.8372093 

top1perITCinFU.cumufre <- lapply(top1perITCinFU, function(x){
  sum(x$Proportion)
}) %>% unlist(.)


# non-top 1%
nontop1perITCinFU <- lapply(pat.inf$T, function(x){
  index <- intersect(tcr.immunarch.nontop1per$data[[x]]$CDR3.aa,
                     tcr.immunarch$data[[pat.inf[pat.inf$T == x,'FU']]]$CDR3.aa)
  
  index <- tcr.immunarch.nontop1per$data[[x]][tcr.immunarch.nontop1per$data[[x]]$CDR3.aa %in% index, ]
  return(index)
}) 
names(nontop1perITCinFU) <- pat.inf$T

nontop1perITCinFU.prop <- unlist(lapply(nontop1perITCinFU,nrow))/unlist(lapply(tcr.immunarch.nontop1per$data[pat.inf$T], nrow))
#       T809       T811       T853       T841       T781       T824       T850       T791       T833       T803       T869       T788       T789       T814       T865       T795       T831       T828 
# 0.09143156 0.22052980 0.18522427 0.13059398 0.22760041 0.12799339 0.09502262 0.09408341 0.28391167 0.30136986 0.16021362 0.15402090 0.06411379 0.23151951 0.25086705 0.17632552 0.33258659 0.28250610 
#       T784 
# 0.24953358

nontop1perITCinFU.cumufre <- lapply(nontop1perITCinFU, function(x){
  sum(x$Proportion)
}) %>% unlist(.)




############## proportion (percentage) change during nivo-therapy in peripheral blood, compared with pre (pretreatment blood) using shared proportion with T (posttreatment tumor)

plot.errorbar <- c(rep(0,19), top1perITCinpost.prop-top1perITCinpre.prop, top1perITCinFU.prop-top1perITCinpre.prop, 
                   rep(0,19), nontop1perITCinpost.prop - nontop1perITCinpre.prop, nontop1perITCinFU.prop - nontop1perITCinpre.prop)*100
plot.errorbar <- cbind(plot.errorbar, point = c(rep(c(rep('pre',19), rep('postpre', 19), rep('FUpre', 19)),2)), 
                       response = c(rep(c(rep('NR', 12), rep('R',7)),6)), 
                       type = c(rep('top1per',57), rep('nontop1per',57))) %>% as.data.frame(.)
plot.errorbar$plot.errorbar <- as.numeric(plot.errorbar$plot.errorbar)
plot.errorbar$point <- factor(plot.errorbar$point, levels = c('pre', 'postpre', 'FUpre'))

p.per1 <- ggplot(plot.errorbar, aes(x = point, y = plot.errorbar, group = type, color = type)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.0))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",
               width = 0.25,position = position_dodge(0.0)) +
  geom_hline(yintercept = 0,colour="grey30", linetype = 'dotdash', alpha =0.5)+
  # scale_y_continuous(breaks = c(-10,0,10,20,30))+
  theme_classic()+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 0)) +
  ylab('Percentage change of ITC in blood') +
  xlab(' ')


t.test(c(top1perITCinpost.prop-top1perITCinpre.prop), alternative = 'two.sided') # p-value = 0.001379 **
t.test(c(top1perITCinFU.prop-top1perITCinpre.prop), alternative = 'two.sided') # p-value = 0.283 ns
t.test(c(top1perITCinFU.prop-top1perITCinpost.prop), alternative = 'two.sided') # p-value = 0.001938 **
# t.test(c(top1perITCinFU.prop), mu=mean(top1perITCinpost.prop), alternative = 'two.sided') # p-value = 0.01801 *
# t.test(c(top1perITCinpost.prop), mu=mean(top1perITCinFU.prop), alternative = 'two.sided') # p-value = 0.00544 **

t.test(c(nontop1perITCinpost.prop - nontop1perITCinpre.prop), alternative = 'two.sided') # p-value = 0.0005755 ***
t.test(c(nontop1perITCinFU.prop - nontop1perITCinpre.prop), alternative = 'two.sided') # p-value = 0.6515 ns
t.test(c(nontop1perITCinFU.prop - nontop1perITCinpost.prop), alternative = 'two.sided') # p-value = 0.001609 **




plot.errorbar.NR <- plot.errorbar[plot.errorbar$response == 'NR',]
plot.errorbar.NR$type1 <- paste0(plot.errorbar.NR$type,'.',plot.errorbar.NR$response)
p.per2 <- ggplot(plot.errorbar.NR, aes(x = point, y = plot.errorbar, group = type, color = type)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.0))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",
               width = 0.25,position = position_dodge(0.0)) +
  geom_hline(yintercept = 0,colour="grey30", linetype = 'dotdash', alpha =0.5)+
  # scale_y_continuous(breaks = c(-10,0,10,20,30))+
  scale_y_continuous(limits = c(-60,80))+
  theme_classic()+
  theme(aspect.ratio = 0.33, axis.text.x = element_text(angle = 0))  +
  ylab(' ') +
  xlab(' ')

t.test(c(top1perITCinpost.prop-top1perITCinpre.prop)[1:12], alternative = 'two.sided') # p-value = 0.04872 *
t.test(c(top1perITCinFU.prop-top1perITCinpre.prop)[1:12], alternative = 'two.sided') # p-value = 0.674 ns
t.test(c(top1perITCinFU.prop-top1perITCinpost.prop)[1:12], alternative = 'two.sided') # p-value = 0.01388 *

t.test(c(nontop1perITCinpost.prop - nontop1perITCinpre.prop)[1:12], alternative = 'two.sided') # p-value = 0.0165 *
t.test(c(nontop1perITCinFU.prop - nontop1perITCinpre.prop)[1:12], alternative = 'two.sided') # p-value = 0.1804 ns
t.test(c(nontop1perITCinFU.prop - nontop1perITCinpost.prop)[1:12], alternative = 'two.sided') # p-value = 0.01293 *


plot.errorbar.R <- plot.errorbar[plot.errorbar$response == 'R',]
plot.errorbar.R$type1 <- paste0(plot.errorbar.R$type,'.',plot.errorbar.R$response)
p.per3 <- ggplot(plot.errorbar.R, aes(x = point, y = plot.errorbar, group = type, color = type)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.0))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",
               width = 0.25,position = position_dodge(0.0)) +
  geom_hline(yintercept = 0,colour="grey30", linetype = 'dotdash', alpha =0.5)+
  # scale_y_continuous(breaks = c(-10,0,10,20,30))+
  theme_classic()+
  theme(aspect.ratio = 0.33, axis.text.x = element_text(angle = 0))  +
  ylab(' ') +
  xlab(' ')

t.test(c(top1perITCinpost.prop-top1perITCinpre.prop)[13:19], alternative = 'two.sided') # p-value = 0.01282 *
t.test(c(top1perITCinFU.prop-top1perITCinpre.prop)[13:19], alternative = 'two.sided') # p-value = 0.08923 ns
t.test(c(top1perITCinFU.prop-top1perITCinpost.prop)[13:19], alternative = 'two.sided') # p-value = 0.08288 ns

t.test(c(nontop1perITCinpost.prop - nontop1perITCinpre.prop)[13:19], alternative = 'two.sided') # p-value = 0.01849 *
t.test(c(nontop1perITCinFU.prop - nontop1perITCinpre.prop)[13:19], alternative = 'two.sided') # p-value = 0.5066 ns
t.test(c(nontop1perITCinFU.prop - nontop1perITCinpost.prop)[13:19], alternative = 'two.sided') # p-value = 0.06446 ns


p.per4 <- ggarrange(p.per2, p.per3, nrow = 2, ncol = 1, legend = FALSE)
ggarrange(p.per1, p.per4, nrow = 1, ncol = 2, common.legend = TRUE)

ggsave(filename = "ITCBloodDynamicCurve.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use



plot.errorbar.RNR <- rbind(plot.errorbar.R, plot.errorbar.NR)

p.per5 <- ggplot(plot.errorbar.RNR, aes(x = point, y = plot.errorbar, group = type1, color = type1)) +
  stat_summary(fun.y = mean, geom = "line")+
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.0))+
  stat_summary(fun.data = mean_se,geom = "errorbar", color = "black",
               width = 0.25,position = position_dodge(0.0)) +
  geom_hline(yintercept = 0,colour="grey30", linetype = 'dotdash', alpha =0.5)+
  # scale_y_continuous(breaks = c(-10,0,10,20,30))+
  theme_classic()+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 0))  +
  ylab(' ') +
  xlab(' ')

ggsave(filename = "ITCBloodDynamicCurve_ref.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use

ggarrange(p.per1, p.per5, nrow = 1, ncol = 2, common.legend = TRUE)
ggsave(filename = "ITCBloodDynamicCurve1.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in") ###### use









