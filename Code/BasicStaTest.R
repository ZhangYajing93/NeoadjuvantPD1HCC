## diversity change test

## load data
load(file = 'D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/nivo.TCR.basicstats.RData')

load('D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/OriData/pat.inf.RData')


###################################################
## pre and post treatment statistical analysis

###################################################
## nivo-NR

testtar <- c('count', 'diversity', 'mean_frequency', 'geomean_frequency', 'mean_cdr3nt_length', 'mean_insert_size', 
             'mean_ndn_size', 'convergence')

basicstats.test.NR <- lapply(testtar, function(x){
  prepost <- wilcox.test(basicstats[basicstats$group =='pre' & basicstats$response == 'NR',x],
                         basicstats[basicstats$group =='post'& basicstats$response == 'NR',x],
                         alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  postFU <- wilcox.test(basicstats[basicstats$group =='post' & basicstats$response == 'NR',x],
                        basicstats[basicstats$group =='FU' & basicstats$response == 'NR',x],
                        alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  preFU <- wilcox.test(basicstats[basicstats$group =='pre' & basicstats$response == 'NR',x],
                       basicstats[basicstats$group =='FU' & basicstats$response == 'NR',x],
                       alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  BT <- wilcox.test(basicstats[basicstats$group =='B' & basicstats$response == 'NR',x],
                    basicstats[basicstats$group =='T' & basicstats$response == 'NR',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  TNT <- wilcox.test(basicstats[basicstats$group =='T' & basicstats$response == 'NR',x],
                     basicstats[basicstats$group =='NT' & basicstats$response == 'NR',x],
                     alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  BNT <- wilcox.test(basicstats[basicstats$group =='B' & basicstats$response == 'NR',x],
                     basicstats[basicstats$group =='NT' & basicstats$response == 'NR',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(prepost, postFU, preFU, BT, TNT, BNT)
  return(testres)
})

basicstats.test.NR <- do.call(rbind, basicstats.test.NR)
colnames(basicstats.test.NR) <- c('prepost', 'postFU', 'preFU', 'BT', 'TNT', 'BNT')
rownames(basicstats.test.NR) <- testtar

## nR
#                      prepost    postFU     preFU         BT        TNT         BNT
# count              0.6220703 0.4697266 0.9096680 0.24962293 0.12939453 0.009696186
# diversity          0.2661133 0.3393555 0.6220703 0.05289808 0.42382813 0.024455936
# mean_frequency     0.1293945 0.6220703 0.9697266 0.05289808 0.67724609 0.024455936
# geomean_frequency  0.3012695 0.6220703 0.9096680 0.08306399 0.96972656 0.013466925
# mean_cdr3nt_length 0.6772461 1.0000000 0.7333984 0.89161819 0.09228516 0.024455936
# mean_insert_size   0.6772461 0.6772461 0.9697266 0.24962293 0.51855469 0.681965094
# mean_ndn_size      0.8500977 0.6772461 0.9697266 0.17970265 0.38037109 0.681965094
# convergence        0.5185547 0.4238281 0.5693359 0.68196509 0.46972656 0.290777850


## nivo-R
basicstats.test.R <- lapply(testtar, function(x){
  prepost <- wilcox.test(basicstats[basicstats$group =='pre' & basicstats$response == 'R',x],
                         basicstats[basicstats$group =='post'& basicstats$response == 'R',x],
                         alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  postFU <- wilcox.test(basicstats[basicstats$group =='post' & basicstats$response == 'R',x],
                        basicstats[basicstats$group =='FU' & basicstats$response == 'R',x],
                        alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  preFU <- wilcox.test(basicstats[basicstats$group =='pre' & basicstats$response == 'R',x],
                       basicstats[basicstats$group =='FU' & basicstats$response == 'R',x],
                       alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  BT <- wilcox.test(basicstats[basicstats$group =='B' & basicstats$response == 'R',x],
                    basicstats[basicstats$group =='T' & basicstats$response == 'R',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  TNT <- wilcox.test(basicstats[basicstats$group =='T' & basicstats$response == 'R',x],
                     basicstats[basicstats$group =='NT' & basicstats$response == 'R',x],
                     alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  BNT <- wilcox.test(basicstats[basicstats$group =='B' & basicstats$response == 'R',x],
                     basicstats[basicstats$group =='NT' & basicstats$response == 'R',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(prepost, postFU, preFU, BT, TNT, BNT)
  return(testres)
})

basicstats.test.R <- do.call(rbind, basicstats.test.R)
colnames(basicstats.test.R) <- c('prepost', 'postFU', 'preFU', 'BT', 'TNT', 'BNT')
rownames(basicstats.test.R) <- testtar

## R
#                    prepost  postFU   preFU        BT      TNT        BNT
# count               0.9375 0.93750 0.93750 0.1636364 0.031250 0.64848485
# diversity           0.8125 0.81250 0.46875 0.9272727 0.031250 0.16363636
# mean_frequency      0.8125 0.68750 0.68750 0.9272727 0.031250 0.16363636
# geomean_frequency   0.9375 0.46875 1.00000 0.3151515 0.015625 0.52727273
# mean_cdr3nt_length  0.9375 0.15625 0.03125 0.5272727 0.109375 0.07272727
# mean_insert_size    0.8125 0.93750 0.68750 0.4121212 0.812500 0.52727273
# mean_ndn_size       0.8125 0.93750 0.68750 0.4121212 0.578125 0.52727273
# convergence         1.0000 0.93750 0.93750 0.3151515 0.031250 1.00000000


## nivo-NR

testtar <- c('reads', 'diversity', 'observedDiversity_mean', 'chaoE_mean', 'efronThisted_mean', 'chao1_mean', 
             'd50Index_mean', 'shannonWienerIndex_mean', 'normalizedShannonWienerIndex_mean',
             'inverseSimpsonIndex_mean','clonality')

diversity.test.NR <- lapply(testtar, function(x){
  prepost <- wilcox.test(diversity[diversity$group =='pre' & diversity$response == 'NR',x],
                         diversity[diversity$group =='post' & diversity$response == 'NR',x],
                         alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  postFU <- wilcox.test(diversity[diversity$group =='post' & diversity$response == 'NR',x],
                        diversity[diversity$group =='FU' & diversity$response == 'NR',x],
                        alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  preFU <- wilcox.test(diversity[diversity$group =='pre' & diversity$response == 'NR',x],
                       diversity[diversity$group =='FU' & diversity$response == 'NR',x],
                       alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  BT <- wilcox.test(diversity[diversity$group =='B' & diversity$response == 'NR',x],
                    diversity[diversity$group =='T' & diversity$response == 'NR',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  TNT <- wilcox.test(diversity[diversity$group =='T' & diversity$response == 'NR',x],
                     diversity[diversity$group =='NT' & diversity$response == 'NR',x],
                     alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  BNT <- wilcox.test(diversity[diversity$group =='B' & diversity$response == 'NR',x],
                     diversity[diversity$group =='NT' & diversity$response == 'NR',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(prepost, postFU, preFU, BT, TNT, BNT)
  return(testres)
})

diversity.test.NR <- do.call(rbind, diversity.test.NR)
colnames(diversity.test.NR) <- c('prepost', 'postFU', 'preFU', 'BT', 'TNT', 'BNT')
rownames(diversity.test.NR) <- testtar


## nR
#                                     prepost     postFU     preFU         BT       TNT         BNT
# reads                             0.6220703 0.46972656 0.9096680 0.24962293 0.1293945 0.009696186
# diversity                         0.2661133 0.33935547 0.6220703 0.05289808 0.4238281 0.024455936
# observedDiversity_mean            0.2661133 0.33935547 0.6220703 0.05289808 0.4238281 0.024455936
# chaoE_mean                        0.3012695 0.33935547 0.6772461 0.10245637 0.4697266 0.041478130
# efronThisted_mean                 0.3012695 0.23339844 0.3012695 0.12464986 0.5185547 0.031997414
# chao1_mean                        0.2333984 0.15136719 0.3393555 0.10245637 0.4697266 0.041478130
# d50Index_mean                     0.1762695 0.30126953 0.8500977 0.15050636 0.4697266 0.041478130
# shannonWienerIndex_mean           0.2036133 0.09228516 0.5185547 0.15050636 0.7910156 0.179702650
# normalizedShannonWienerIndex_mean 0.7333984 0.23339844 0.4697266 0.96358543 0.9697266 0.891618186
# inverseSimpsonIndex_mean          0.8500977 0.23339844 0.1513672 0.61646197 0.6772461 0.493643611
# clonality                         0.7333984 0.23339844 0.4697266 0.96358543 0.9697266 0.891618186

## nivo-R

diversity.test.R <- lapply(testtar, function(x){
  prepost <- wilcox.test(diversity[diversity$group =='pre' & diversity$response == 'R',x],
                         diversity[diversity$group =='post' & diversity$response == 'R',x],
                         alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  postFU <- wilcox.test(diversity[diversity$group =='post' & diversity$response == 'R',x],
                        diversity[diversity$group =='FU' & diversity$response == 'R',x],
                        alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  preFU <- wilcox.test(diversity[diversity$group =='pre' & diversity$response == 'R',x],
                       diversity[diversity$group =='FU' & diversity$response == 'R',x],
                       alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  BT <- wilcox.test(diversity[diversity$group =='B' & diversity$response == 'R',x],
                    diversity[diversity$group =='T' & diversity$response == 'R',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  TNT <- wilcox.test(diversity[diversity$group =='T' & diversity$response == 'R',x],
                     diversity[diversity$group =='NT' & diversity$response == 'R',x],
                     alternative = "two.sided",paired = TRUE,var.equal=FALSE,conf.level = 0.95)$p.value
  BNT <- wilcox.test(diversity[diversity$group =='B' & diversity$response == 'R',x],
                     diversity[diversity$group =='NT' & diversity$response == 'R',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(prepost, postFU, preFU, BT, TNT, BNT)
  return(testres)
})

diversity.test.R <- do.call(rbind, diversity.test.R)
colnames(diversity.test.R) <- c('prepost', 'postFU', 'preFU', 'BT', 'TNT', 'BNT')
rownames(diversity.test.R) <- testtar

## R
#                                    prepost   postFU    preFU         BT      TNT       BNT
# reads                             0.937500 0.937500 0.937500 0.16363636 0.031250 0.6484848
# diversity                         0.812500 0.812500 0.468750 0.92727273 0.031250 0.1636364
# observedDiversity_mean            0.812500 0.812500 0.468750 0.92727273 0.031250 0.1636364
# chaoE_mean                        0.687500 0.937500 0.687500 0.78787879 0.015625 0.2303030
# efronThisted_mean                 0.578125 0.687500 0.375000 1.00000000 0.015625 0.2303030
# chao1_mean                        0.468750 0.687500 0.468750 0.92727273 0.015625 0.2303030
# d50Index_mean                     0.937500 0.218750 0.578125 0.78787879 0.296875 0.3151515
# shannonWienerIndex_mean           0.375000 0.578125 0.375000 0.16363636 0.687500 0.2303030
# normalizedShannonWienerIndex_mean 0.296875 0.578125 0.078125 0.04242424 0.156250 0.9272727
# inverseSimpsonIndex_mean          0.218750 0.375000 0.156250 0.07272727 0.046875 0.9272727
# clonality                         0.296875 0.578125 0.078125 0.04242424 0.156250 0.9272727


###################################################



###################################################
######### nivo-R vs nivo-NR unpaired test


testtar <- c('count', 'diversity', 'mean_frequency', 'geomean_frequency', 'mean_cdr3nt_length', 'mean_insert_size', 
             'mean_ndn_size', 'convergence')

basicstats.test.RnR <- lapply(testtar, function(x){
  pre <- wilcox.test(basicstats[basicstats$group =='pre' & basicstats$response == 'NR',x],
                     basicstats[basicstats$group =='pre'& basicstats$response == 'R',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  post <- wilcox.test(basicstats[basicstats$group =='post' & basicstats$response == 'NR',x],
                        basicstats[basicstats$group =='post' & basicstats$response == 'R',x],
                        alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  FU <- wilcox.test(basicstats[basicstats$group =='FU' & basicstats$response == 'NR',x],
                       basicstats[basicstats$group =='FU' & basicstats$response == 'R',x],
                       alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  B <- wilcox.test(basicstats[basicstats$group =='B' & basicstats$response == 'NR',x],
                    basicstats[basicstats$group =='B' & basicstats$response == 'R',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  T <- wilcox.test(basicstats[basicstats$group =='T' & basicstats$response == 'NR',x],
                     basicstats[basicstats$group =='T' & basicstats$response == 'R',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  NT <- wilcox.test(basicstats[basicstats$group =='NT' & basicstats$response == 'NR',x],
                     basicstats[basicstats$group =='NT' & basicstats$response == 'R',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(pre, post, FU, B, T, NT)
  return(testres)
})

basicstats.test.RnR <- do.call(rbind, basicstats.test.RnR)
colnames(basicstats.test.RnR) <- c('pre', 'post', 'FU', 'B', 'T', 'NT')
rownames(basicstats.test.RnR) <- testtar

#                          pre      post        FU         B          T         NT
# count              0.3844566 0.3402397 0.8369453 0.9142857 0.02834008 0.06834961
# diversity          0.9671350 0.3402397 0.8369453 0.7619048 0.08311503 0.16734143
# mean_frequency     0.9671350 0.3402397 0.8369453 0.7619048 0.08311503 0.16734143
# geomean_frequency  0.5918076 0.3844566 0.9018020 0.9142857 0.01706756 0.14221640
# mean_cdr3nt_length 0.5358419 0.2614114 0.2268397 0.9142857 0.65035326 0.59180757
# mean_insert_size   0.7732397 0.9018020 0.9671350 0.9142857 0.08311503 0.08311503
# mean_ndn_size      0.9018020 0.5358419 0.9671350 0.9142857 0.08311503 0.06834961
# convergence        0.5358419 0.3402397 0.9018020 1.0000000 0.08311503 0.73520489



testtar <- c('reads', 'diversity', 'observedDiversity_mean', 'chaoE_mean', 'efronThisted_mean', 'chao1_mean', 
             'd50Index_mean', 'shannonWienerIndex_mean', 'normalizedShannonWienerIndex_mean',
             'inverseSimpsonIndex_mean','clonality')

diversity.test.RnR <- lapply(testtar, function(x){
  pre <- wilcox.test(diversity[diversity$group =='pre' & diversity$response == 'NR',x],
                     diversity[diversity$group =='pre' & diversity$response == 'R',x],
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  post <- wilcox.test(diversity[diversity$group =='post' & diversity$response == 'NR',x],
                      diversity[diversity$group =='post' & diversity$response == 'R',x],
                      alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  FU <- wilcox.test(diversity[diversity$group =='FU' & diversity$response == 'NR',x],
                    diversity[diversity$group =='FU' & diversity$response == 'R',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  B <- wilcox.test(diversity[diversity$group =='B' & diversity$response == 'NR',x],
                   diversity[diversity$group =='B' & diversity$response == 'R',x],
                   alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  T <- wilcox.test(diversity[diversity$group =='T' & diversity$response == 'NR',x],
                   diversity[diversity$group =='T' & diversity$response == 'R',x],
                   alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  NT <- wilcox.test(diversity[diversity$group =='NT' & diversity$response == 'NR',x],
                    diversity[diversity$group =='NT' & diversity$response == 'R',x],
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(pre, post, FU, B, T, NT)
  return(testres)
})

diversity.test.RnR <- do.call(rbind, diversity.test.RnR)
colnames(diversity.test.RnR) <- c('pre', 'post', 'FU', 'B', 'T', 'NT')
rownames(diversity.test.RnR) <- testtar


#                                     pre      post        FU         B          T         NT
# reads                             0.3844566 0.3402397 0.8369453 0.9142857 0.02834008 0.06834961
# diversity                         0.9671350 0.3402397 0.8369453 0.7619048 0.08311503 0.16734143
# observedDiversity_mean            0.9671350 0.3402397 0.8369453 0.7619048 0.08311503 0.16734143
# chaoE_mean                        0.6503533 0.4824164 0.8369453 0.6095238 0.05564817 0.22683972
# efronThisted_mean                 1.0000000 0.3844566 0.9671350 0.6095238 0.03584187 0.19556244
# chao1_mean                        0.9018020 0.3402397 0.9671350 0.7619048 0.03584187 0.19556244
# d50Index_mean                     0.4824164 1.0000000 0.5918076 0.9142857 0.06834961 0.59180757
# shannonWienerIndex_mean           0.4319679 0.3844566 0.9018020 0.6095238 0.90180202 0.53584187
# normalizedShannonWienerIndex_mean 0.3402397 1.0000000 0.8369453 1.0000000 0.02210844 0.71084385
# inverseSimpsonIndex_mean          1.0000000 0.4824164 0.5358419 0.7619048 0.19556244 0.29911884
# clonality                         0.3402397 1.0000000 0.8369453 1.0000000 0.02210844 0.71084385


##################################
###### parameter change statistical test

testtar <- c('count', 'diversity', 'mean_frequency', 'geomean_frequency', 'mean_cdr3nt_length', 'mean_insert_size', 
             'mean_ndn_size', 'convergence')


basicstats.change.test.RnR <- lapply(testtar, function(x){
  postpre <- wilcox.test(c(basicstats[basicstats$group =='post'& basicstats$response == 'R',x] - basicstats[basicstats$group =='pre'& basicstats$response == 'R',x]),
                         c(basicstats[basicstats$group =='post'& basicstats$response == 'NR',x] - basicstats[basicstats$group =='pre'& basicstats$response == 'NR',x]),
                         alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  FUpost <- wilcox.test(c(basicstats[basicstats$group =='FU'& basicstats$response == 'R',x] - basicstats[basicstats$group =='post'& basicstats$response == 'R',x]),
                        c(basicstats[basicstats$group =='FU'& basicstats$response == 'NR',x] - basicstats[basicstats$group =='post'& basicstats$response == 'NR',x]),
                        alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  FUpre <- wilcox.test(c(basicstats[basicstats$group =='FU'& basicstats$response == 'R',x] - basicstats[basicstats$group =='pre'& basicstats$response == 'R',x]),
                        c(basicstats[basicstats$group =='FU'& basicstats$response == 'NR',x] - basicstats[basicstats$group =='pre'& basicstats$response == 'NR',x]),
                        alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  TB <- wilcox.test(c(basicstats[basicstats$group =='T'& basicstats$patient_id %in% basicstats[basicstats$group =='B'& basicstats$response == 'R','patient_id'],x] - basicstats[basicstats$group =='B'& basicstats$response == 'R',x]),
                    c(basicstats[basicstats$group =='T'& basicstats$patient_id %in% basicstats[basicstats$group =='B'& basicstats$response == 'NR','patient_id'],x] - basicstats[basicstats$group =='B'& basicstats$response == 'NR',x]), 
                       alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  TNT <- wilcox.test(c(basicstats[basicstats$group =='T'& basicstats$response == 'R',x] - basicstats[basicstats$group =='NT'& basicstats$response == 'R',x]),
                       c(basicstats[basicstats$group =='T'& basicstats$response == 'NR',x] - basicstats[basicstats$group =='NT'& basicstats$response == 'NR',x]),
                       alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(postpre, FUpost, FUpre, TB, TNT)
  return(testres)
})

basicstats.change.test.RnR <- do.call(rbind, basicstats.change.test.RnR)
colnames(basicstats.change.test.RnR) <- c('postpre', 'FUpost', 'FUpre', 'TB', 'TNT')
rownames(basicstats.change.test.RnR) <- testtar





testtar <- c('reads', 'diversity', 'observedDiversity_mean', 'chaoE_mean', 'efronThisted_mean', 'chao1_mean', 
             'd50Index_mean', 'shannonWienerIndex_mean', 'normalizedShannonWienerIndex_mean',
             'inverseSimpsonIndex_mean','clonality')

diversity.change.test.RnR <- lapply(testtar, function(x){
  postpre <- wilcox.test(c(diversity[diversity$group =='post'& diversity$response == 'R',x] - diversity[diversity$group =='pre'& diversity$response == 'R',x]),
                         c(diversity[diversity$group =='post'& diversity$response == 'NR',x] - diversity[diversity$group =='pre'& diversity$response == 'NR',x]),
                         alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  FUpost <- wilcox.test(c(diversity[diversity$group =='FU'& diversity$response == 'R',x] - diversity[diversity$group =='post'& diversity$response == 'R',x]),
                        c(diversity[diversity$group =='FU'& diversity$response == 'NR',x] - diversity[diversity$group =='post'& diversity$response == 'NR',x]),
                        alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  FUpre <- wilcox.test(c(diversity[diversity$group =='FU'& diversity$response == 'R',x] - diversity[diversity$group =='pre'& diversity$response == 'R',x]),
                       c(diversity[diversity$group =='FU'& diversity$response == 'NR',x] - diversity[diversity$group =='pre'& diversity$response == 'NR',x]),
                       alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  TB <- wilcox.test(c(diversity[diversity$group =='T'& diversity$patient_id %in% diversity[diversity$group =='B'& diversity$response == 'R','patient_id'],x] - diversity[diversity$group =='B'& diversity$response == 'R',x]),
                    c(diversity[diversity$group =='T'& diversity$patient_id %in% diversity[diversity$group =='B'& diversity$response == 'NR','patient_id'],x] - diversity[diversity$group =='B'& diversity$response == 'NR',x]), 
                    alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  
  TNT <- wilcox.test(c(diversity[diversity$group =='T'& diversity$response == 'R',x] - diversity[diversity$group =='NT'& diversity$response == 'R',x]),
                     c(diversity[diversity$group =='T'& diversity$response == 'NR',x] - diversity[diversity$group =='NT'& diversity$response == 'NR',x]),
                     alternative = "two.sided",paired = FALSE,var.equal=FALSE,conf.level = 0.95)$p.value
  testres <- c(postpre, FUpost, FUpre, TB, TNT)
  return(testres)
})

diversity.change.test.RnR <- do.call(rbind, diversity.change.test.RnR)
colnames(diversity.change.test.RnR) <- c('postpre', 'FUpost', 'FUpre', 'TB', 'TNT')
rownames(diversity.change.test.RnR) <- testtar

################## no difference between R 和 NR


########################## clonality: clonality boxplot


diversity$group1 <- paste0(diversity$group, '.', diversity$response)

ggboxplot(diversity, x = 'group1', y = 'clonality', add = 'dotplot', ylab = 'Clonality', xlab = '',
          add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  stat_compare_means(method = 'wilcox.test', paired = TRUE,
                     comparisons = list(c(1,3), c(3,5), c(1,5),
                                        c(2,4), c(4,6),c(2,6), 
                                        c(9,11),
                                        c(10,12))) +
  stat_compare_means(method = 'wilcox.test', paired = FALSE,
                     comparisons = list(c(1,2),c(3,4), c(5,6),c(7,8),c(9,10),c(11,12),
                                        c(7,9), c(7,11),
                                        c(8,10),c(8,12)
                     ))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) # 显著组别： T.NR vs T.R (0.017), 有差异趋势B.R vs T.R (0.073)

ggsave(filename = "clonality.boxplot_ref.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")

clonality.plot <- ggboxplot(diversity, x = 'group1', y = 'clonality', add = 'dotplot', ylab = 'Clonality', xlab = '',
                            add.params = list(size=0.5), color = 'response', palette = c("#00AFBB", "#E7B800")) +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE,
  #                    comparisons = list(c(1,3), c(1,5), c(3,5),
  #                                       c(2,4), c(2,6), c(4,6),
  #                                       c(9,11),
  #                                       c(10,12))) +
  # stat_compare_means(method = 'wilcox.test', paired = FALSE, step.increase = 0.1,
  #                    comparisons = list(c(9,10),c(8,10)
  #                    ))+
  theme(aspect.ratio = 0.66, axis.text.x = element_text(angle = 30)) 

ggsave(filename = "clonality.boxplot_edit.pdf",
       path = "D:/HKU/Nivo_ICItherapy/TCR-seq/DataAnalysis/Fig",
       width=8,height=6,units="in")




