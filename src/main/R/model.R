

#print(length(args))
args <- commandArgs(trailingOnly = TRUE)

#args[1] = "HumanR01_1123"
#args[2] = "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_psm.csv"
#args[3] = "autoRTHumanR01_1t3d.csv"
#args[4] = "humanR01_RTminfit_maxfit_diff_1t3d.txt"
#args[5] = "/data/humanR01/"

print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])

#setwd("/home/rstudio/data/")
setwd(args[5])


#
#library(keras)
#library(stringi)
#library(stringr)
#library(data.table)
#
##autortThreshold = 25
#autortThreshold = 25
##autortThreshold = 2
##模型个数
#baggingsize = 1#5
##cross validation
##CVs = c(seq(0:3)-1)
#CVs = c(0)
#CVslen = length(CVs)
##迭代次数
##Iteratorloop = c(seq(1:30))
#Iteratorloop = c(seq(1:2))
##Iteratorloop = c(seq(1:3))
#IteratorloopLen=length(Iteratorloop)
#best_ter = c(rep(0,CVslen))
#lasttrainResullt = 0
#trainSetNoDown = 0
#
##过程结果统计
#resultSave <- array(0,dim = c(IteratorloopLen,3,CVslen,baggingsize))
#
#
#dfAutoRTWithAll <- fread(args[2])
#
#dfUniqueAutoRT <- unique(dfAutoRTWithAll[, c("x", "y_pred")])
#
#
#write.table(dfUniqueAutoRT, file=args[3], sep="\t", row.names = FALSE, quote=FALSE)
#
#
#
#setwd("~/data")
#
library(keras)
library(stringi)
library(stringr)
library(data.table)

#autortThreshold = 25
autortThreshold = 25
#autortThreshold = 2
#模型个数
baggingsize = 5#5
#cross validation
#CVs = c(seq(0:3)-1)
CVs = c(0)
CVslen = length(CVs)
#迭代次数
#Iteratorloop = c(seq(1:30))
Iteratorloop = c(seq(1:10))
#Iteratorloop = c(seq(1:3))
IteratorloopLen=length(Iteratorloop)
best_ter = c(rep(0,CVslen))
lasttrainResullt = 0
trainSetNoDown = 0

#过程结果统计
resultSave <- array(0,dim = c(IteratorloopLen,3,CVslen,baggingsize))


#迭代时fdr的阈值 #迭代少，这里放大一些
#fdrThreshold = 0.3
fdrThreshold = 0.3
#中间结果存放数据信息
#modelInfo = "dp/1013/MouseR02-1t3dTop-RT25-"
modelInfo = paste0("dp/1109/TT-",args[1])
#modelInfo = "dp/1109/MouseR02-1t3dTop-RT25-7I35FdrlessIterTestPara-"
#saveResultInfo = paste0("rt",autortThreshold,"fdr",fdrThreshold,"_MouseR02Top")
saveResultInfo = paste0("rt",autortThreshold,"fdr",fdrThreshold,args[3])
#saveResultInfo = paste0("rt",autortThreshold,"fdr",fdrThreshold,"_MouseR02Top_lessIterTestPara")
saveResultInfo = gsub("\\.", "", saveResultInfo)


#seq_len(IteratorloopLen)


#resultSaveTest <- array(0,dim = c(5,3,4,5))

#resultSaveTest[1,2,3,4]=1


#第一个模型的最优结果存放
resultRound1 = c(1:baggingsize)


#dfOurAllmatchResultnewWithAutoRTbac <-
#  fread("newResult1012HumanR01.csv")#R02 trpsin without consider p

#dfOurAllmatchResultnewWithAutoRTbac <-
#  fread("newResult0919MouseR2OriTrailNewFeature.csv")#R02 trpsin without consider p

#dfOurAllmatchResultnewWithAutoRTbac <-
#  fread("newResult1017MouseR2OriTrailNewFeatureRemoveIL1t3d.csv")#R02 trpsin without consider p

#这是mouse最新的
#dfOurAllmatchResultnewWithAutoRTbac <-
#  fread("newResult1026NewFeatureTopNpredictIntensity.csv")#R02 trpsin without consider p

#这是人类最新的
dfOurAllmatchResultnewWithAutoRTbac <-
  fread(args[2])#"newResult1029HumanR01_1t3d.csv")#R02 trpsin without consider p

dfMSoneFeatureOne <- dfOurAllmatchResultnewWithAutoRTbac[msonefeature==-1]

dfOurAllmatchResultnewWithAutoRTbac[,`:=`(count=NULL)]

dfOurAllmatchResultnewWithAutoRT <-
  dfOurAllmatchResultnewWithAutoRTbac[order(-matchscore), count := seq(.N), by = msonefeature][msonefeature>-1  & count<11]

dfOurAllmatchResultnewWithAutoRTbac <- rbind(dfOurAllmatchResultnewWithAutoRT,dfMSoneFeatureOne)

#dfOurAllmatchResultnewWithAutoRTbac[,`:=`(decoyTargetUpdate=NULL)]

#这是mouse最新的
#dfAutoRTWithAll <- fread("autoRTMouseR021t3d.tsv")
#这是人类最新 的
dfAutoRTWithAll <- fread(args[3])#"autoRTHumanR01_1t3d.csv")
#dfAutoRTWithAll <- fread("mouseR02autoRT.tsv")
#dfAutoRTWithAll <- fread("autoRTdata3.tsv")
#View(dfAutoRTWithAll[,.N,.(round(y_pred,1))])
dfUniqueAutoRT <- unique(dfAutoRTWithAll[, c("x", "y_pred")])
#读取autort拟合数据
#这是mouse最新的
#dfAutoRTFit <- fread("autoRtfit_minfit_maxfit_diff.csv")
#这是人类最新 的
dfAutoRTFit <- fread(args[4])#"humanR01_RTminfit_maxfit_diff_1t3d.txt")
#dfAutoRTFit <- fread("autoRtfit_minfit_maxfit_diffHumanR01.csv")


#for autort
#View(unique(dfOurAllmatchResultnewWithAutoRTbac[msonefeature>-1,c("peptide","msonetime")]))
#dfmouseR02PepRT<-unique(dfOurAllmatchResultnewWithAutoRTbac[msonefeature>-1,c("peptide","msonetime")])
#names(dfmouseR02PepRT)=c("x","y")
#dfmouseR02PepRT[,.SD[1],.(x)]
#dfmouseR02PepRT[,`:=`(seqlen=nchar(x))]
#
#write.table(dfmouseR02PepRT[seqlen<=98,.SD[1],.(x)][,c(1,2)], file="peptide_mouseR02_forAutoRT1t3d.csv", sep="\t", row.names = FALSE, quote=FALSE)


#for prosit
#nrow(unique(dfOurAllmatchResultnewWithAutoRTbac[(msonefeature>-1) & (str_detect(proteinname,"DeBruijn_Decoy")),c("peptide","msoneZ")]))
#dfmouseR02PepZ<-unique(dfOurAllmatchResultnewWithAutoRTbac[(msonefeature>-1) & (str_detect(proteinname,"DeBruijn_Decoy")),c("peptide","msoneZ")])
#names(dfmouseR02PepZ)=c("modified_sequence","precursor_charge")
#dfmouseR02PepZ[,`:=`(collision_energy=27)]
#dfmouseR02PepZ= dfmouseR02PepZ[,c(1,3,2)]
#dfmouseR02PepZ[,`:=`(seqlen=nchar(modified_sequence))]
#
#write.table(dfmouseR02PepZ[seqlen>=7 & seqlen<=30  & precursor_charge<7   ,1:3], file="peptide_mouseR02_forProsit1t3d.csv", sep=",", row.names = FALSE, quote=FALSE)


#View(dfOurAllmatchResultnewWithAutoRTbac[(msonefeature>-1) & (str_detect(proteinname,"DeBruijn_Decoy")),][1:1000])

#View(dfAutoRTFit)
dfUniqueAutoRT[, y_predWith1decimal := round(y_pred, 1)]
dfUniqueAutoRT

dfUniqueAutoRT[dfAutoRTFit, on = .(y_predWith1decimal = autoRt), fitNew :=
                 (fit_min + fit_max) / 2]

#View(dfUniqueAutoRT[is.na(fitNew) ,.N,.(y_predWith1decimal)])

#View(dfUniqueAutoRT[is.na(fitNew) & y_predWith1decimal<min(dfAutoRTFit$autoRt) ,])
#some value not in top6w
dfUniqueAutoRT[is.na(fitNew) & y_predWith1decimal<min(dfAutoRTFit$autoRt) ,fitNew :=0]
#names(dfOurAllmatchResultnewWithAutoRTbac[,c(1:10)])
#names(dfOurAllmatchResultnewWithAutoRT[,c(1:10)])

#去除虚拟的precursor数据
#dfOurAllmatchResultnewWithAutoRT = dfOurAllmatchResultnewWithAutoRTbac[msonefeature >  -1]
#nrow(dfOurAllmatchResultnewWithAutoRTbac[msonefeature >  -1])
#nrow(dfOurAllmatchResultnewWithAutoRTbac[msonefeature ==  -1])

#cross window( generate. duplicate ) -remove duplicate
dfOurAllmatchResultnewWithAutoRT = dfOurAllmatchResultnewWithAutoRT[, head(.SD, 1),.(msonefeature,peptide)]
#nrow(dfOurAllmatchResultnewWithAutoRT)
#View(dfOurAllmatchResultnewWithAutoRTbac[msonefeature >  -1,.N,.(msonefeature,peptide)][N>1])
#dfdpround1UniqueCombineAfter[, head(.SD, 1), .(peptide)]
#dfOurAllmatchResultnewWithAutoRT = dfOurAllmatchResultnewWithAutoRTbac[msonefeature >  -1 & count<3]#for test

#dfOurAllmatchResultnewWithAutoRT[dfUniqueAutoRT, on = .(peptide = x), autort := y_pred]
#names(dfOurAllmatchResultnewWithAutoRT)
#View(dfOurAllmatchResultnewWithAutoRT[1:60000,c(1:12,75)])
#write.table(dfOurAllmatchResultnewWithAutoRT[1:60000,c(9,75,6)], file="humanR01forRTminfit_maxfit_diff_1t3d.csv", sep=",", row.names = FALSE, quote=FALSE)

dfOurAllmatchResultnewWithAutoRT[dfUniqueAutoRT, on = .(peptide = x), fitWith3 := i.fitNew]
#filter out psm with predict rt
dfOurAllmatchResultnewWithAutoRT[,.N,.(decoyTarget)]

dfOurAllmatchResultnewWithAutoRT = dfOurAllmatchResultnewWithAutoRT[abs(fitWith3 -
                                                                          msonetime) < autortThreshold]
nrow(dfOurAllmatchResultnewWithAutoRT)



minprositRt = 0 - min(dfOurAllmatchResultnewWithAutoRT$prositRT)
dfOurAllmatchResultnewWithAutoRT[, `:=`(
  diffRealWithPredictRTWith3 =
    pmax(0, 1 - (((msonetime - fitWith3) / autortThreshold
    ) ^ 2)),
  diffRealWithPrositRT =
    pmax(0, 1 - (((
      msonetime - (prositRT + minprositRt)
    ) / 10) ^ 2))
)]
#设置标志位，有autort匹配为1无则为0
dfOurAllmatchResultnewWithAutoRT[, `:=`(autoRTBz = 1)]
dfOurAllmatchResultnewWithAutoRT[is.na(fitWith3), `:=`(autoRTBz = 0)]
#没有匹配对原始值设置为0
dfOurAllmatchResultnewWithAutoRT[is.na(fitWith3), `:=`(fitWith3 = 0, diffRealWithPredictRTWith3 =
                                                         0)]
dfOurAllmatchResultnewWithAutoRT[, `:=`(prositBz = 1)]
dfOurAllmatchResultnewWithAutoRT[peptideLength < 7 |
                                   peptideLength > 30  | msoneZ >= 7, `:=`(prositBz = 0)]
#可以支持的peptide differ charge count with the same rt
peptideRTdiffZCountSTCondition = dfOurAllmatchResultnewWithAutoRTbac[(msonefeature >
                                                                        -1 & count == 1 &  (MS2TrailCoscinSimilarity > 4.5))
                                                                     |
                                                                       (msonefeature == -1 &
                                                                          matchscore > 1 & (MS2TrailCoscinSimilarity / peptideLength) > 0.40),
                                                                     .N, .(msonetime, msoneZ, peptide)][, MSOneTimepeptideCount := .N, .(msonetime, peptide)]
head(peptideRTdiffZCountSTCondition)
names(peptideRTdiffZCountSTCondition)[1:3] = c("MS1RT", "MS1Z", "Peptide")
dfOurAllmatchResultnewWithAutoRT[, `:=`(diffChargeCount = 0)]
dfOurAllmatchResultnewWithAutoRT[peptideRTdiffZCountSTCondition, on = .(msonetime =
                                                                          MS1RT ,
                                                                        msoneZ = MS1Z,
                                                                        peptide = Peptide), `:=`(diffChargeCount = MSOneTimepeptideCount)]

dfOurAllmatchResultnewWithAutoRT[, `:=`(autoRTAbs = abs(msonetime - fitWith3),autoRTSquare = (msonetime - fitWith3)*(msonetime - fitWith3) )]
dfOurAllmatchResultnewWithAutoRT[, `:=`(matchCountAll = bionMatchCount_C1 + yionMatchCount_C1 + bionMatchCount_C2 + yionMatchCount_C2)]
#View(dfOurAllmatchResultnewWithAutoRT[count == 1 & (matchCountAll>15 & (MS2TrailCoscinSimilarity / peptideLength)) > 0.40,.N,.(matchCountAll,decoyTarget)])


dfOurAllmatchResultnewWithAutoRT[,`:=`(decoyTargetUpdate="target")]

dfOurAllmatchResultnewWithAutoRT[str_detect(proteinname,"DeBruijn_Decoy"),`:=`(decoyTargetUpdate="decoy")]
dfOurAllmatchResultnewWithAutoRT[,.N,.(decoyTargetUpdate)]
dfOurAllmatchResultnewWithAutoRT[decoyTargetUpdate=='target',.N,.(decoyTarget)]

#View(dfOurAllmatchResultnewWithAutoRT[,.N,.(floor(abs(fitWith3-msonetime)),decoyTargetUpdate)])
#View(dfOurAllmatchResultnewWithAutoRT[decoyTargetUpdate=="decoy" & decoyTarget=="decoy",c("msonetime","fitWith3")])
#names(dfOurAllmatchResultnewWithAutoRT)
#tail(dfOurAllmatchResultnewWithAutoRTbac)
#
#nrow(dfOurAllmatchResultnewWithAutoRTbac[msonefeature>-1 & str_detect(proteinname,"DeBruijn_Decoy"),])
#View(dfOurAllmatchResultnewWithAutoRTbac[msonefeature>-1 ,.N,.(str_detect(proteinname,"DeBruijn_Decoy"))])
dfOurAllmatchResultnewWithAutoRTbac = NULL
gc()
#特征
#
#namesFeature <- names(dfOurAllmatchResultnewWithAutoRT)[c(84, 8:13, 15:19, 20:21, 24:35, 38:49, 57:58, 66:69, 72, 73, 75:76, 77)] #include  prositrt selected
namesFeature <- names(dfOurAllmatchResultnewWithAutoRT)[c(84, 8:13, 15:19, 20:21, 24:35, 38:49, 57:58, 66:69, 70:73, 76, 77, 79:80, 81)] #include  prositrt selected
#summary(dfOurAllmatchResultnewWithAutoRT)
namesFeatureAll = namesFeature
nFeatureCount <- length(namesFeatureAll) - 1#163
iLabelIndex <- length(namesFeatureAll)
#namesOriginalInfo = c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 20, 21, 30, 31, 22, 23)
#namesOriginalInfo = c(1,  4, 5, 6, 80, 2, 7, 8, 9, 11, 12, 20, 21, 30, 31, 22, 23)
#namesOriginalInfo = c(1,  4, 5, 6, 84, 2, 7, 8, 9, 11, 12, 20, 21, 30, 31, 22, 23)
namesOriginalInfo = c(1,  4, 5, 84, 2, 7, 8, 9, 11, 12, 20, 21, 30, 31, 22, 23)


testPrecursorpsm = dfOurAllmatchResultnewWithAutoRT#train dataset




testPrecursorpsm = testPrecursorpsm[sample(1:nrow(testPrecursorpsm)),]
sampleDataOriginal <- testPrecursorpsm[, ..namesOriginalInfo]
sampleDataTD <- testPrecursorpsm[, ..namesFeatureAll]
sampleDataTD$TDValue <- 1
names(sampleDataTD)
#sampleDataTD$TDValue[which(sampleDataTD$decoyTarget == 'decoy')] <- 0
sampleDataTD$TDValue[which(sampleDataTD$decoyTargetUpdate == 'decoy')] <- 0
sampleDataTD <- as.data.frame(sampleDataTD)
dataTopMatrix <- as.matrix(sampleDataTD[, 2:(iLabelIndex + 1)])
dimnames(dataTopMatrix) <- NULL

test <- dataTopMatrix[, 1:nFeatureCount]
testtraget <- dataTopMatrix[, iLabelIndex]
testlabel <- to_categorical(testtraget)

#sampleDataTD[,.N,.(TDValue)]
#sum(testlabel[,1])
#sum(testlabel[,2])
#
#testPrecursorpsm[,.N,.(decoyTargetUpdate)]
#sampleDataOriginal[,.N,.(decoyTarget)]
#sampleDataOriginaltest <- testPrecursorpsm[, ..namesOriginalInfo]
#
#sampleDataOriginaltest[,.N,.(decoyTarget)]
#
#sampleDataOriginaltest$decoyTarget = 'target'

dfOurAllmatchResultnewWithTop1RT <-
  dfOurAllmatchResultnewWithAutoRT[count == 1, ..namesFeature]

names(dfOurAllmatchResultnewWithTop1RT)[1]="decoyTarget"

dfOurAllmatchResultnewWithTop1RT$decoyTarget <-
  as.factor(dfOurAllmatchResultnewWithTop1RT$decoyTarget)
#dfOurAllmatchResultnewWithTop1RT[1:10]
summary(dfOurAllmatchResultnewWithTop1RT)
str(dfOurAllmatchResultnewWithTop1RT)

#tail(dfOurAllmatchResultnewWithTop1RT[decoyTarget=="target"])
#tail(dfOurAllmatchResultnewWithAutoRT[decoyTargetUpdate=="target" & count==1])

#head(dfOurAllmatchResultnewWithTop1RT)
mylogit <-
  glm(decoyTarget ~ ., data = dfOurAllmatchResultnewWithTop1RT, family = "binomial")
#Warning message:
#  glm.fit: fitted probabilities numerically 0 or 1 occurred
# names(mylogit$coefficients)
#mylogit$coefficients[[2]]
#  summary(mylogit)
#calculate lg value for each dfOurAllmatchResultnewWithAutoRT
nlgCoeff <- length(mylogit$coefficients)
#  mylogit$coefficients
#summary(dfOurAllmatchResultnewWithAutoRT)
dfOurAllmatchResultnewWithAutoRT[, `:=`(
  LGsumScore =
    mylogit$coefficients[[1]] +
    msonemass                        * mylogit$coefficients[[2]] +
    msonetime                        * mylogit$coefficients[[3]] +
    msoneQualityScore                * mylogit$coefficients[[4]] +
    msoneMZ                          * mylogit$coefficients[[5]] +
    msoneZ                           * mylogit$coefficients[[6]] +
    msonePeakAreaLocalRank           * mylogit$coefficients[[7]] +
    peptideLength                    * mylogit$coefficients[[8]] +
    MS1PeptideMassError              * mylogit$coefficients[[9]] +
    dPeptideMutationRate             * mylogit$coefficients[[10]]  +
    ibyComplemetaryCount             * mylogit$coefficients[[11]] +
    dCosWithPredictMSMS              * mylogit$coefficients[[12]] +
    dMS1WithMS2TrailCoscinSimilarity * mylogit$coefficients[[13]] +
    MS2TrailCoscinSimilarity         * mylogit$coefficients[[14]] +
    dlogBionPeakAreaSum_C1           * mylogit$coefficients[[15]] +
    dlogYionPeakAreaSum_C1           * mylogit$coefficients[[16]] +
    dBionMassErrorSum_C1             * mylogit$coefficients[[17]] +
    dYionMassErrorSum_C1             * mylogit$coefficients[[18]] +
    iBionPeakHalfSum_C1              * mylogit$coefficients[[19]] +
    iYionPeakHalfSum_C1              * mylogit$coefficients[[20]] +
    bionMatchCount_C1                * mylogit$coefficients[[21]] +
    yionMatchCount_C1                * mylogit$coefficients[[22]] +
    dBionRetentionTimeErrorSum_C1    * mylogit$coefficients[[23]] +
    dYionRetentionTimeErrorSum_C1    * mylogit$coefficients[[24]] +
    iBConsective_C1                  * mylogit$coefficients[[25]] +
    iYConsective_C1                  * mylogit$coefficients[[26]] +
    dlogBionPeakAreaSum_C2           * mylogit$coefficients[[27]] +
    dlogYionPeakAreaSum_C2           * mylogit$coefficients[[28]] +
    dBionMassErrorSum_C2             * mylogit$coefficients[[29]] +
    dYionMassErrorSum_C2             * mylogit$coefficients[[30]] +
    iBionPeakHalfSum_C2              * mylogit$coefficients[[31]] +
    iYionPeakHalfSum_C2              * mylogit$coefficients[[32]] +
    bionMatchCount_C2                * mylogit$coefficients[[33]] +
    yionMatchCount_C2                * mylogit$coefficients[[34]] +
    dBionRetentionTimeErrorSum_C2    * mylogit$coefficients[[35]] +
    dYionRetentionTimeErrorSum_C2    * mylogit$coefficients[[36]] +
    iBConsective_C2                  * mylogit$coefficients[[37]] +
    iYConsective_C2                  * mylogit$coefficients[[38]] +
    dBionCos60WithPredictMSMS_C1     * mylogit$coefficients[[39]] +
    dCosWithPredictMSMSMatchedPredict * mylogit$coefficients[[40]] +
    isRawdata                        * mylogit$coefficients[[41]] +
    peptideAnomiAcidFreRate          * mylogit$coefficients[[42]] +
    lastIonNumBC1                    * mylogit$coefficients[[43]] +
    lastIonNumYC1                    * mylogit$coefficients[[44]] +
    top6allCos                       * mylogit$coefficients[[45]] +
    top6c1Cos                        * mylogit$coefficients[[46]] +
    top6Bc1Cos                       * mylogit$coefficients[[47]] +
    top6Yc1Cos                       * mylogit$coefficients[[48]] +

    diffRealWithPredictRTWith3     * mylogit$coefficients[[49]] +
    diffRealWithPrositRT     * mylogit$coefficients[[50]] +
    #autoRTBz       *mylogit$coefficients[[44]] +
    prositBz       * mylogit$coefficients[[51]] +
    diffChargeCount       * mylogit$coefficients[[52]]
  +
    autoRTAbs       * mylogit$coefficients[[53]]
  #    diffRealWithPredictRTWith3     * mylogit$coefficients[[45]] +
  #    diffRealWithPrositRT     * mylogit$coefficients[[46]] +
  #    #autoRTBz       *mylogit$coefficients[[44]] +
  #    prositBz       * mylogit$coefficients[[47]] +
  #    diffChargeCount       * mylogit$coefficients[[48]]
  #  +
  #    autoRTAbs       * mylogit$coefficients[[49]]
  #  +
  #    autoRTSquare       * mylogit$coefficients[[47]]
  #                                         diffChargeCount                  *mylogit$coefficients[[41]] +
  #                                         diffRealWithPredictRTWith3     *mylogit$coefficients[[42]] +
  #                                         diffRealWithPrositRT     *mylogit$coefficients[[43]] +
  #                                         autoRTBz       *mylogit$coefficients[[44]] +
  #                                         prositBz       *mylogit$coefficients[[45]]
)]
#为countLG排序排名
dfOurAllmatchResultnewWithAutoRT <-
  dfOurAllmatchResultnewWithAutoRT[order(dfOurAllmatchResultnewWithAutoRT$LGsumScore, decreasing = TRUE), countLG := seq(.N), by =
                                     msonefeature]
dftest <- dfOurAllmatchResultnewWithAutoRT[countLG == 1]
dftest <- dftest[order(dftest$LGsumScore, decreasing = TRUE), ]
Un1 <- unique(sort(dftest$decoyTargetUpdate))
dftest$nth_test <- 0
dftest$resp_target <- 0
dftest$resp_decoy <- 0
dftest <-
  dftest[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
                                                              lapply(Un1, function(x)
                                                                cumsum(x == decoyTargetUpdate))) ,]
dftest$fdr <- dftest$resp_decoy / dftest$resp_target

dftest[max(which(dftest$fdr < 0.01))][, c(resp_target, LGsumScore)]
length(unique(dftest[LGsumScore >= dftest[max(which(dftest$fdr < 0.01))]$LGsumScore]$peptide))



#sampleTarDecTopCombineShuffle <-
#  dftest[trainTestBZ == FALSE ][1:200000]
#迭代少，这里放大一些
#sampleTarDecTopCombineShuffle <-
#  dftest[1:300000]
sampleTarDecTopCombineShuffle <-
  dftest[1:300000]

#dfTargetT10WTop <-
#  sampleTarDecTopCombineShuffle[decoyTarget == "target"]
dfTargetT10WTop <-  sampleTarDecTopCombineShuffle[decoyTargetUpdate == "target"]
#dfTargetT10WTop <-
#  sampleTarDecTopCombineShuffle[decoyTarget == "target" & fdr<0.1]

#dfTargetT10WTop <-dftest[decoyTarget == "target" & fdr<0.2 & trainTestBZ==FALSE]
#dfDecoyT10WTop <- dftest[decoyTarget == "decoy" & fdr<0.2 & trainTestBZ==FALSE]



#dfDecoyT10WTop <- sampleTarDecTopCombineShuffle[decoyTarget == "decoy"]
dfDecoyT10WTop <- sampleTarDecTopCombineShuffle[decoyTargetUpdate == "decoy"]
#get decoy samples
#dfDecoyT10WRest <- dfOurAllmatchResultnewWithAutoRT[decoyTarget ==
#                                                      "decoy" & msonefeature %in%
#                                                      dfTargetT10WTop$msonefeature][, .SD[][1:1], by = "msonefeature"]
#
#sampleTarDecTopCombineShuffle <-
#  rbind(dfTargetT10WTop[fdr<fdrThreshold,c(1:83)],dfDecoyT10WTop[,c(1:83)], dfDecoyT10WRest[,c(1:83)])

dfDecoyT10WRest <- dfOurAllmatchResultnewWithAutoRT[decoyTargetUpdate ==
                                                      "decoy" & msonefeature %in%
                                                      dfTargetT10WTop$msonefeature][, .SD[][1:1], by = "msonefeature"]

#head(dfDecoyT10WTop)
#head(dfTargetT10WTop[fdr<fdrThreshold])
#head(dfDecoyT10WTop[fdr<fdrThreshold])


dfTargetT10WRest <- dfOurAllmatchResultnewWithAutoRT[decoyTargetUpdate ==
                                                       "target" & msonefeature %in%
                                                       dfDecoyT10WTop[fdr<fdrThreshold]$msonefeature][, .SD[][1:1], by = "msonefeature"]

#names(dfTargetT10WTop)
#names(dfDecoyT10WRest)

#sampleTarDecTopCombineShuffle <-
#  rbind(dfTargetT10WTop[fdr<fdrThreshold,c(1:82)],dfDecoyT10WTop[,c(1:82)], dfDecoyT10WRest[,c(1:82)])
sampleTarDecTopCombineShuffle <-
  rbind(dfTargetT10WTop[fdr<fdrThreshold,c(1:86)],dfDecoyT10WTop[,c(1:86)], dfDecoyT10WRest[,c(1:86)], dfTargetT10WRest[,c(1:86)])
#sampleTarDecTopCombineShuffle <-
#  rbind(dfTargetT10WTop[fdr<fdrThreshold,c(1:82)],dfDecoyT10WTop[,c(1:82)], dfDecoyT10WRest[,c(1:82)])

sampleTarDecTopCombineShuffle[,.N,(decoyTarget)]
sampleTarDecTopCombineShuffle[,.N,(decoyTargetUpdate)]

#head(sampleTarDecTopCombineShuffle[decoyTargetUpdate=='target' & decoyTarget=='decoy'])

dsaveMax = 0
group =1
i_ter =1
i_CVs =1
for (i_ter in seq_len(IteratorloopLen)) {

  #the next iteration get the new trainning dataset
  #   sampleTarDecTopCombineShuffle = sampleTarDecTopCombineShufflenew

  sampleTarDecTopCombineShuffle = sampleTarDecTopCombineShuffle[sample(1:nrow(sampleTarDecTopCombineShuffle)),]
  namesFeatureAll = namesFeature
  sampleDataTD <- sampleTarDecTopCombineShuffle[, ..namesFeatureAll]
  sampleDataTD$TDValue <- 1
  #    names(sampleDataTD)
  sampleDataTD$TDValue[which(sampleDataTD$decoyTargetUpdate == 'decoy')] <- 0
  #  sampleDataTD$TDValue[which(sampleDataTD$decoyTarget == 'decoy')] <- 0
  nFeatureCount <- length(namesFeatureAll) - 1#163
  iLabelIndex <- length(namesFeatureAll)
  sampleDataTD <- as.data.frame(sampleDataTD)
  dataTopMatrix <- as.matrix(sampleDataTD[, 2:(iLabelIndex + 1)])
  dimnames(dataTopMatrix) <- NULL
  training <- dataTopMatrix[, 1:nFeatureCount]
  trainingtraget <- dataTopMatrix[, iLabelIndex]
  traininglabel <- to_categorical(trainingtraget)

  #sum(traininglabel[,1])
  #sum(traininglabel[,2])
  resultSave[i_ter,1,(CVs[i_CVs]+1),group]= nrow(training)

  if (exists("model1"))
  {
    rm(model1)
  }
  #blocks = c(1,  4)
  blocks = c(1, 2, 4)
  kernels = c(6, 3, 2)
  activation = "relu"
  input_tensor  <- layer_input(shape = c(nFeatureCount))
  output_tensor  <- input_tensor
  shortcut  <- input_tensor
  n_blocks  <- length(blocks)
  for (i in seq_len(n_blocks)) {
    output_tensor_x  <- layer_dense(output_tensor,
                                    units =  1024 / blocks[[i]]) #padding = "same"
    output_tensor_x  <- layer_batch_normalization(output_tensor_x)
    output_tensor_x  <-
      layer_activation(output_tensor_x,  activation = activation)
    output_tensor_x  <- layer_dropout(output_tensor_x,  rate = 0.5)
    output_tensor_y  <- layer_dense(output_tensor_x,
                                    units =  512 / blocks[[i]]) #padding = "same"
    output_tensor_y  <- layer_batch_normalization(output_tensor_y)
    output_tensor_y  <-
      layer_activation(output_tensor_y,  activation = activation)
    output_tensor_y  <- layer_dropout(output_tensor_y,  rate = 0.3)
    output_tensor_z  <- layer_dense(output_tensor_y,
                                    units =  256 / blocks[[i]]) #padding = "same"
    output_tensor_z  <- layer_batch_normalization(output_tensor_z)
    output_tensor_z  <- layer_dropout(output_tensor_z,  rate = 0.2)
    shortcut  <- layer_dense(shortcut, units =  256 / blocks[[i]])
    shortcut  <- layer_batch_normalization(shortcut)
    output_tensor  <- layer_add(list(shortcut,  output_tensor_z))
    output_tensor  <-
      layer_activation(output_tensor,  activation = activation)
    shortcut  <- output_tensor
  }
  output_tensor  <- layer_dense(output_tensor,
                                units = 2, activation = "softmax")
  FPPenaltyCorrect <- function(true, pred) {
    K <- backend()
    predGT =  K$greater_equal(pred, 0.5)
    predGT09 =  K$greater_equal(pred, 0.90)
    predGT08 =  K$greater_equal(pred, 0.85)
    #  predGT07 =  K$greater_equal(pred,0.75)
    predGT07 =  K$greater_equal(pred, 0.80)
    #  predGT09 =  K$greater_equal(pred,0.90)
    #  predGT08 =  K$greater_equal(pred,0.80)
    #  predGT07 =  K$greater_equal(pred,0.70)
    #predGT =  K$greater(pred,0.5)
    predGT = K$cast(predGT, K$floatx())
    predGT09 = K$cast(predGT09, K$floatx())
    predGT08 = K$cast(predGT08, K$floatx())
    predGT07 = K$cast(predGT07, K$floatx())
    ##  K$print_tensor('pred: ',pred)
    diffs = predGT - true
    diffs09 = predGT09 - true
    diffs08 = predGT08 - true
    diffs07 = predGT07 - true
    ##  K$print_tensor('diffs: ',diffs)
    #  y_ture_f = K$flatten(true)
    #  pred = K$cast(pred, K$floatx())
    #  y_pred_f = K$cast(  K$greater(pred,0.5), K$floatx())
    #  intersection = y_ture_f * y_pred_f
    #  K$print_tensor('intersection: ',intersection)
    #  TN = K$all(K$eval(true) == 0 , K$eval(y_pred_f) == 0)
    #  TP = K$all(K$eval(true) == 1 , K$eval(y_pred_f) == 1)
    #  K$print_tensor(K$eval(true))
    lesser = K$greater(diffs[, 2], 0)
    lesser09 = K$greater(diffs09[, 2], 0)
    lesser08 = K$greater(diffs08[, 2], 0)
    lesser07 = K$greater(diffs07[, 2], 0)
    equal09 = K$equal(diffs09[, 2], 0)
    TP09 = K$sum(K$cast(equal09 , K$floatx()))
    NP = K$sum(K$cast(K$equal(true[, 2], 1) , K$floatx()))
    N = K$sum(K$cast(K$equal(true[, 2], 0) , K$floatx()))
    FP = K$sum(K$cast(lesser, K$floatx()))
    FP09 = K$sum(K$cast(lesser09, K$floatx()))
    FP08 = K$sum(K$cast(lesser08, K$floatx())) - FP09
    FP07 = K$sum(K$cast(lesser07, K$floatx())) - FP09 - FP08
    100 * (FP09 / N) + 60 * (FP08 / N) + 30 * (FP07 / N) + K$categorical_crossentropy(true, pred)#+10*(FP0506/N)##128)
  }
  model1 <- keras_model(inputs = input_tensor,
                        outputs = output_tensor) %>% compile(loss = "categorical_crossentropy",#FPPenaltyCorrect,
                                                             #"categorical_crossentropy",# FPPenaltyCorrect,#"categorical_crossentropy",
                                                             optimizer = optimizer_rmsprop(),
                                                             # optimizer_sgd(lr=0.01),# optimizer_rmsprop(),#(lr = 0.0001, decay = 1e-6),
                                                             metrics = "accuracy")
  history <- model1 %>% fit(
    training,
    traininglabel,
    #    epochs = 10,
    epochs = 120,
    #epochs = 80,
    batch_size = 1024,
    #batch_size = 2048,
    #2048,#2048,#1024 [1] 42960
    validation_split = 0.2
    ,
    callbacks = list(
      # save best model after every epoch
      callback_model_checkpoint(
        file.path(
          paste0(modelInfo,"_g",group,"_",CVs[i_CVs],"Test-I",i_ter,"-dp.hdf5")
          #"dp/0809/mouse0928-023trainI5dp-{epoch:03d}-{val_loss:4f}-{accuracy:4f}-{val_accuracy:.4f}.hdf5"
        ),
        monitor = 'val_accuracy',
        save_best_only = TRUE,
        mode = 'max'
      ),
      #存系统认为的最好结果
      #callback_model_checkpoint(file.path("dp/dp_checkpoints.h5"),  monitor='val_accuracy', save_best_only = TRUE),#存系统认为的最好结果
      callback_csv_logger(file.path("dp/dp_csv_logs")),
      #存日志
      # callback_early_stopping(monitor='val_loss', patience=15, verbose=1),
      # only needed for visualising with TensorBoard
      callback_tensorboard(log_dir = file.path("dp/dp_logs"))#存日志
    )
  )
  if (exists("model1"))
  {
    rm(model1)
  }
  model1 <-
    load_model_hdf5(

      paste0(modelInfo,"_g",group,"_",CVs[i_CVs],"Test-I",i_ter,"-dp.hdf5"),

      custom_objects = list("FPPenaltyCorrect" = FPPenaltyCorrect)
    )

  #best dp round 1 for mouse r02
  #save_model_hdf5(model1, 'dp/my_model091-123trainFirstIterdpfdrgt20_1—all.hdf5')

  prob1 <- model1 %>%
    predict(test)
  #////////////////////////////////////////////
  #  ///分析结果，计算fdr
  #rm(testResultSort21)
  #rm(testResult1)
  testResult1 <- cbind(sampleDataOriginal, testlabel, prob1)
  names(testResult1)
  restcolnum = ncol(testResult1) - 3
  names(testResult1)[c(restcolnum:ncol(testResult1))] = c("decoyLabel", "targetlabel", "decoyscore", "targetscore")

  testResult1 <-testResult1[order(decoyscore, -targetscore, -matchscore, -MS2TrailCoscinSimilarity)]

  #  testResult1[,.N,.(decoyTarget)]
  #  testResult1[,.N,.(decoyTargetUpdate)]

  testResult1 <-
    testResult1[, countDP := seq(.N), by = msonefeature]


  testResult2 <- testResult1[countDP == 1]

  #  testResult2[,.N,.(decoyTarget)]
  #  testResult2[,.N,.(decoyTargetUpdate)]
  testResultSort21 <- as.data.frame(testResult2)

  testResultSort21$decoyTarget = 'target'
  testResultSort21$decoyTarget[which(testResultSort21$targetlabel == 0)] <- 'decoy'

  #  data.frame(table(testResultSort21$decoyTarget))
  #  data.frame(table(testResultSort21$decoyTargetUpdate))

  testResultSort21$decoyTarget <-  as.factor(testResultSort21$decoyTarget)

  testResultSort21$score <- log(testResultSort21$targetscore / testResultSort21$decoyscore)
  summary(testResultSort21$score)
  testResultSort21$score[which(is.infinite(testResultSort21$score) &
                                 testResultSort21$score > 0)] <-
    max(testResultSort21$score[which(is.finite(testResultSort21$score))]) +
    1
  testResultSort21$score[which(is.infinite(testResultSort21$score) &
                                 testResultSort21$score < 0)] <-
    min(testResultSort21$score[which(is.finite(testResultSort21$score))]) -
    1
  testResultSort21 <- as.data.table(testResultSort21)
  #最快的方法FDR
  dftest1 <-
    testResultSort21[order(testResultSort21$score, decreasing = TRUE), ]
  Un1 <- unique(sort(testResultSort21$decoyTarget))
  #Un1<-"decoy"
  dftest1$nth_test <- 0
  dftest1$resp_target <- 0
  dftest1$resp_decoy <- 0
  dftest1 <-
    dftest1[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
                                                                 lapply(Un1, function(x)
                                                                   cumsum(x == decoyTarget))) ,]
  dftest1$fdr <- dftest1$resp_decoy / dftest1$resp_target
  setorder(dftest1, -score)
  dftest1[max(which(dftest1$fdr < 0.01))][, c(resp_target, nth_test, targetscore, score)]
  length(unique(dftest1[score >= dftest1[max(which(dftest1$fdr < 0.01))]$score]$peptide))

  dftest2 = dftest1[, head(.SD, 1), .(peptide)]
  nrow(dftest2)
  dftest2$nth_test <- 0
  dftest2$resp_target <- 0
  dftest2$resp_decoy <- 0
  dftest2 <-
    dftest2[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
                                                                 lapply(Un1, function(x)
                                                                   cumsum(x == decoyTarget))) ,]
  dftest2$fdr <- dftest2$resp_decoy / dftest2$resp_target
  setorder(dftest2, -score)
  dftest2[max(which(dftest2$fdr < 0.01))][, c(resp_target, nth_test, targetscore, score)]
  resultSave[i_ter,2,(CVs[i_CVs]+1),group]=  length(unique(dftest2[score >= dftest2[max(which(dftest2$fdr < 0.01))]$score]$peptide))

  if(dsaveMax < resultSave[i_ter,2,(CVs[i_CVs]+1),group])
  {
    best_ter[i_CVs ] = i_ter #record the best iterator

    dsaveMax =resultSave[i_ter,2,(CVs[i_CVs]+1),group]
    if(CVs[i_CVs]==0){  dftrain0psm = testResult1 }
    if(CVs[i_CVs]==1){  dftrain1psm = testResult1 }
    if(CVs[i_CVs]==2){  dftrain2psm = testResult1 }
    if(CVs[i_CVs]==3){  dftrain3psm = testResult1 }

  }
  if (i_ter>6 ){
    if( lasttrainResullt > resultSave[i_ter,2,(CVs[i_CVs]+1),group])
    {
      break
    }
  }
  #迭代少，这里放大一些
  #  dfdpround1Top1 <-  dftest2[1:200000]
  dfdpround1Top1 <-  dftest2[1:200000]

  #  testResult1dp[,.N,.(decoyTarget)]
  #  testResult1dp[,.N,.(decoyTargetUpdate)]
  #  dfdpround1Top1[,.N,.(decoyTarget)]
  #  dfdpround1Top1[,.N,.(decoyTargetUpdate)]

  dfTargetT10WTopdpround1Top1 <- dfdpround1Top1[decoyTarget == "target" ]

  dfDecoyT10WTopdpround1Top1 <-  dfdpround1Top1[decoyTarget == "decoy"]
  #get decoy samples

  testResult1dp <- testResult1[order(testResult1$targetscore, decreasing = TRUE), ]


  dfDecoyT10WRestround1Top1 <- testResult1dp[decoyTargetUpdate ==
                                               "decoy" & msonefeature %in%
                                               dfTargetT10WTopdpround1Top1$msonefeature][, .SD[][1:1], by = "msonefeature"]
  dfTargetT10WRestround1Top1 <- testResult1dp[decoyTargetUpdate ==
                                                "target" & msonefeature %in%
                                                dfDecoyT10WTopdpround1Top1[fdr<fdrThreshold]$msonefeature][, .SD[][1:1], by = "msonefeature"]



  #  dfdpround1Top12 <-
  #    rbind(dfTargetT10WTopdpround1Top1[ fdr<fdrThreshold ,c(2:6,1,7:22)],dfDecoyT10WTopdpround1Top1[,c(2:6,1,7:22)] ,dfDecoyT10WRestround1Top1)
  dfdpround1Top12 <-
    rbind(dfTargetT10WTopdpround1Top1[ fdr<fdrThreshold ,c(2:4,5,1,6:21)],dfDecoyT10WTopdpround1Top1[,c(2:4,5,1,6:21)]
          ,dfDecoyT10WRestround1Top1,dfTargetT10WRestround1Top1)


  #  names(dfTargetT10WTopdpround1Top1)
  #  names(dfDecoyT10WRestround1Top1)
  #fast for filter data.table by other data.table with multiple fields
  fileterDP= dfOurAllmatchResultnewWithAutoRT[sort(dfOurAllmatchResultnewWithAutoRT[dfdpround1Top12, on=.(msonefeature,peptide), which=TRUE, nomatch=0])]
  #  fileterDP[,.N,.(decoyTarget)]
  #  fileterDP[,.N,.(decoyTargetUpdate)]
  #  sampleTarDecTopCombineShuffle[,.N,(decoyTarget)]
  #  sampleTarDecTopCombineShuffle[,.N,(decoyTargetUpdate)]

  sampleTarDecTopCombineShuffle = funion(sampleTarDecTopCombineShuffle,fileterDP)

  #record each result
  trainSetNoDown = sampleTarDecTopCombineShuffle
  lasttrainResullt = resultSave[i_ter,2,(CVs[i_CVs]+1),group]

  #20221020
  #sampleTarDecTopCombineShuffle = funion(sampleTarDecTopCombineShuffle[decoyTargetUpdate == "decoy"],fileterDP)

}


#获取了非下降的训练集结果，准备用目前积累的训练集重复跑机器学习多次，多次训练出来的模型去合并

trainSetNoDown[,.N,.(decoyTargetUpdate)]


sampleTarDecTopCombineShuffle = trainSetNoDown[sample(1:nrow(trainSetNoDown)),]
namesFeatureAll = namesFeature
sampleDataTD <- sampleTarDecTopCombineShuffle[, ..namesFeatureAll]
sampleDataTD$TDValue <- 1
#    names(sampleDataTD)
sampleDataTD$TDValue[which(sampleDataTD$decoyTargetUpdate == 'decoy')] <- 0
#  sampleDataTD$TDValue[which(sampleDataTD$decoyTarget == 'decoy')] <- 0
nFeatureCount <- length(namesFeatureAll) - 1#163
iLabelIndex <- length(namesFeatureAll)
sampleDataTD <- as.data.frame(sampleDataTD)
dataTopMatrix <- as.matrix(sampleDataTD[, 2:(iLabelIndex + 1)])
dimnames(dataTopMatrix) <- NULL
training <- dataTopMatrix[, 1:nFeatureCount]
trainingtraget <- dataTopMatrix[, iLabelIndex]
traininglabel <- to_categorical(trainingtraget)

#记录新的多模型结果
i_ter = 1
#resultSave[i_ter,1,(CVs[i_CVs]+2),group]= nrow(training)

#dsaveMax

#best_ter[i_CVs ]

resultSave

resultSaveRound1 <- c(1:IteratorloopLen)

#生成新的测试集数据

#testPrecursorpsm = dfOurAllmatchResultnewWithAutoRT#train dataset
#namesFeature <- names(dfOurAllmatchResultnewWithAutoRT)[c(6, 8:13, 15:19, 20:21, 24:35, 38:49, 57:58, 66:69, 72, 73, 75:76, 77)] #include  prositrt selected
namesFeature <- names(dfOurAllmatchResultnewWithAutoRT)[c(6, 8:13, 15:19, 20:21, 24:35, 38:49, 57:58, 66:69, 70:73, 76, 77, 79:80, 81)] #include  prositrt selected
testPrecursorpsm = dfOurAllmatchResultnewWithAutoRT[decoyTargetUpdate=="target"]#train dataset
namesOriginalInfo = c(1,  4, 5, 6, 2, 7, 8, 9, 11, 12, 20, 21, 30, 31, 22, 23)


#dfOurAllmatchResultnewWithAutoRT[,.N,.(decoyTarget)]
#dfOurAllmatchResultnewWithAutoRT[,.N,.(decoyTargetUpdate)]

#testPrecursorpsm[,.N,.(decoyTarget)]
namesFeatureAll = namesFeature
testPrecursorpsm = testPrecursorpsm[sample(1:nrow(testPrecursorpsm)),]
sampleDataOriginal <- testPrecursorpsm[, ..namesOriginalInfo]
sampleDataTD <- testPrecursorpsm[, ..namesFeatureAll]
sampleDataTD$TDValue <- 1
names(sampleDataTD)
#sampleDataTD$TDValue[which(sampleDataTD$decoyTarget == 'decoy')] <- 0
sampleDataTD$TDValue[which(sampleDataTD$decoyTarget == 'decoy')] <- 0


sampleDataTD <- as.data.frame(sampleDataTD)
dataTopMatrix <- as.matrix(sampleDataTD[, 2:(iLabelIndex + 1)])
dimnames(dataTopMatrix) <- NULL

test <- dataTopMatrix[, 1:nFeatureCount]
testtraget <- dataTopMatrix[, iLabelIndex]
testlabel <- to_categorical(testtraget)


iTopN = 10
#iTopN = 2


valueThreshold = 1e-5
resultCombiineSaveRound1 <- c(1:iTopN)
for (i_ter in seq_len(IteratorloopLen)) {
  #for (i_ter in c(1:9)) {

  #i_ter = 30

  if (exists("model1"))
  {
    rm(model1)
  }
  blocks = c(1, 2, 4)
  kernels = c(6, 3, 2)
  activation = "relu"
  input_tensor  <- layer_input(shape = c(nFeatureCount))
  output_tensor  <- input_tensor
  shortcut  <- input_tensor
  n_blocks  <- length(blocks)
  for (i in seq_len(n_blocks)) {
    output_tensor_x  <- layer_dense(output_tensor,
                                    units =  1024 / blocks[[i]]) #padding = "same"
    output_tensor_x  <- layer_batch_normalization(output_tensor_x)
    output_tensor_x  <-
      layer_activation(output_tensor_x,  activation = activation)
    output_tensor_x  <- layer_dropout(output_tensor_x,  rate = 0.5)
    output_tensor_y  <- layer_dense(output_tensor_x,
                                    units =  512 / blocks[[i]]) #padding = "same"
    output_tensor_y  <- layer_batch_normalization(output_tensor_y)
    output_tensor_y  <-
      layer_activation(output_tensor_y,  activation = activation)
    output_tensor_y  <- layer_dropout(output_tensor_y,  rate = 0.3)
    output_tensor_z  <- layer_dense(output_tensor_y,
                                    units =  256 / blocks[[i]]) #padding = "same"
    output_tensor_z  <- layer_batch_normalization(output_tensor_z)
    output_tensor_z  <- layer_dropout(output_tensor_z,  rate = 0.2)
    shortcut  <- layer_dense(shortcut, units =  256 / blocks[[i]])
    shortcut  <- layer_batch_normalization(shortcut)
    output_tensor  <- layer_add(list(shortcut,  output_tensor_z))
    output_tensor  <-
      layer_activation(output_tensor,  activation = activation)
    shortcut  <- output_tensor
  }
  output_tensor  <- layer_dense(output_tensor,
                                units = 2, activation = "softmax")

  model1 <- keras_model(inputs = input_tensor,
                        outputs = output_tensor) %>% compile(loss = "categorical_crossentropy",#FPPenaltyCorrect,
                                                             #"categorical_crossentropy",# FPPenaltyCorrect,#"categorical_crossentropy",
                                                             optimizer = optimizer_rmsprop(),
                                                             # optimizer_sgd(lr=0.01),# optimizer_rmsprop(),#(lr = 0.0001, decay = 1e-6),
                                                             metrics = "accuracy")
  history <- model1 %>% fit(
    training,
    traininglabel,
    #    epochs = 10,
    epochs = 120,
    #epochs = 80,
    batch_size = 1024,
    #batch_size = 2048,
    #2048,#2048,#1024 [1] 42960
    validation_split = 0.2
    ,
    callbacks = list(
      # save best model after every epoch
      callback_model_checkpoint(
        file.path(
          paste0(modelInfo,"_ReTrain_Test-I",i_ter,"-dp.hdf5")
          #"dp/0809/mouse0928-023trainI5dp-{epoch:03d}-{val_loss:4f}-{accuracy:4f}-{val_accuracy:.4f}.hdf5"
        ),
        monitor = 'val_accuracy',
        save_best_only = TRUE,
        mode = 'max'
      ),
      #存系统认为的最好结果
      #callback_model_checkpoint(file.path("dp/dp_checkpoints.h5"),  monitor='val_accuracy', save_best_only = TRUE),#存系统认为的最好结果
      callback_csv_logger(file.path("dp/dp_csv_logs")),
      #存日志
      # callback_early_stopping(monitor='val_loss', patience=15, verbose=1),
      # only needed for visualising with TensorBoard
      callback_tensorboard(log_dir = file.path("dp/dp_logs"))#存日志
    )
  )
  if (exists("model1"))
  {
    rm(model1)
  }
  model1 <-
    load_model_hdf5(
      paste0(modelInfo,"_ReTrain_Test-I",i_ter,"-dp.hdf5")
      #,
      #custom_objects = list("FPPenaltyCorrect" = FPPenaltyCorrect)
    )

  #best dp round 1 for mouse r02
  #save_model_hdf5(model1, 'dp/my_model091-123trainFirstIterdpfdrgt20_1—all.hdf5')

  prob1 <- model1 %>%
    predict(test)


  testResult1 <- cbind(sampleDataOriginal, testlabel, prob1)
  names(testResult1)
  restcolnum = ncol(testResult1) - 3
  names(testResult1)[c(restcolnum:ncol(testResult1))] = c("decoyLabel", "targetlabel", "decoyscore", "targetscore")

  testResult1 <-testResult1[order(decoyscore, -targetscore, -matchscore, -MS2TrailCoscinSimilarity)]

  testResult1 <-
    testResult1[, countDP := seq(.N), by = msonefeature]


  testResult2 <- testResult1[countDP == 1]
  testResultSort21 <- as.data.frame(testResult2)
  #  testResult2[,.N,.(decoyTarget)]
  #  testResult2[,.N,.(decoyTargetUpdate)]

  testResultSort21$decoyTarget = 'target'
  testResultSort21$decoyTarget[which(testResultSort21$targetlabel == 0)] <- 'decoy'


  # data.frame(table(testResultSort21$decoyTarget))
  testResultSort21$decoyTarget <-  as.factor(testResultSort21$decoyTarget)

  testResultSort21$score <- log(testResultSort21$targetscore / testResultSort21$decoyscore)
  summary(testResultSort21$score)
  testResultSort21$score[which(is.infinite(testResultSort21$score) &
                                 testResultSort21$score > 0)] <-
    max(testResultSort21$score[which(is.finite(testResultSort21$score))]) +
    1
  testResultSort21$score[which(is.infinite(testResultSort21$score) &
                                 testResultSort21$score < 0)] <-
    min(testResultSort21$score[which(is.finite(testResultSort21$score))]) -
    1
  testResultSort21 <- as.data.table(testResultSort21)
  #最快的方法FDR
  dftest1 <-
    testResultSort21[order(testResultSort21$score, decreasing = TRUE), ]
  Un1 <- unique(sort(testResultSort21$decoyTarget))
  #Un1<-"decoy"
  dftest1$nth_test <- 0
  dftest1$resp_target <- 0
  dftest1$resp_decoy <- 0
  dftest1 <-
    dftest1[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
                                                                 lapply(Un1, function(x)
                                                                   cumsum(x == decoyTarget))) ,]
  dftest1$fdr <- dftest1$resp_decoy / dftest1$resp_target
  setorder(dftest1, -score)
  dftest1[max(which(dftest1$fdr < 0.01))][, c(resp_target, nth_test, targetscore, score)]
  length(unique(dftest1[score >= dftest1[max(which(dftest1$fdr < 0.01))]$score]$peptide))
  dftest2 = dftest1[, head(.SD, 1), .(peptide)]
  nrow(dftest2)
  dftest2$nth_test <- 0
  dftest2$resp_target <- 0
  dftest2$resp_decoy <- 0
  dftest2 <-
    dftest2[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
                                                                 lapply(Un1, function(x)
                                                                   cumsum(x == decoyTarget))) ,]
  dftest2$fdr <- dftest2$resp_decoy / dftest2$resp_target
  setorder(dftest2, -score)
  dftest2[max(which(dftest2$fdr < 0.01))][, c(resp_target, nth_test, targetscore, score)]
  length(unique(dftest2[score >= dftest2[max(which(dftest2$fdr < 0.01))]$score]$peptide))
  #resultSave[best_ter,3,(CVs[i_CVs]+1)]=  length(unique(dftest2[score >= dftest2[max(which(dftest2$fdr < 0.01))]$score]$peptide))
  resultSaveRound1[i_ter] =  length(unique(dftest2[score >= dftest2[max(which(dftest2$fdr < 0.01))]$score]$peptide))



  #combine model

  testResult1[,`:=`(targetscoreadjust=targetscore,decoyscoreadjust=decoyscore)]
  testResult1[targetscore<valueThreshold,`:=`(targetscoreadjust=valueThreshold)]
  testResult1[decoyscore<valueThreshold,`:=`(decoyscoreadjust=valueThreshold)]
  testResult1[,`:=`(scoreBeforeSigmod=log(targetscoreadjust/(decoyscoreadjust)))]

  if(i_ter==1)
  {
    dfdpround1UniqueCombine = testResult1
  }
  if(i_ter>1)
  {

    if(i_ter>2)
    {
      dfdpround1UniqueCombine[,`:=`(newCombineScore=NULL)]
    }

    dfdpround1UniqueCombine = rbind(  dfdpround1UniqueCombine, testResult1    )

    dfdpround1UniqueCombineAfter = dfdpround1UniqueCombine[,`:=`(newCombineScore=sum(scoreBeforeSigmod)/i_ter),.(msonefeature,peptide)][order(-newCombineScore)]
    dftest2CombineAfter = dfdpround1UniqueCombineAfter[, head(.SD, 1), .(peptide)]
    nrow(dftest2CombineAfter)
    dftest2CombineAfter$nth_test <- 0
    dftest2CombineAfter$resp_target <- 0
    dftest2CombineAfter$resp_decoy <- 0
    Un1 <- unique(sort(dftest2CombineAfter$decoyTarget))

    dftest2CombineAfter <-
      dftest2CombineAfter[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
                                                                               lapply(Un1, function(x)
                                                                                 cumsum(x == decoyTarget))) ,]
    dftest2CombineAfter$fdr <- dftest2CombineAfter$resp_decoy / dftest2CombineAfter$resp_target
    setorder(dftest2CombineAfter, -newCombineScore)
    dftest2CombineAfter[max(which(dftest2CombineAfter$fdr < 0.01))][, c(resp_target, nth_test, scoreBeforeSigmod, newCombineScore)]

    resultCombiineSaveRound1[i_ter]=length(unique(dftest2CombineAfter[newCombineScore >= dftest2CombineAfter[max(which(dftest2CombineAfter$fdr < 0.01))]$newCombineScore]$peptide))
  }
}

#print results
resultSaveRound1
resultCombiineSaveRound1


write.table(dftest2CombineAfter, file=paste0(args[1],"_CombineResult.csv", sep="\t", row.names = FALSE, quote=FALSE)

#length(unique(dftest2CombineAfter[newCombineScore >= dftest2CombineAfter[max(which(dftest2CombineAfter$fdr < 0.005))]$newCombineScore]$peptide))
#diann<-fread("~/data1/diannMourseR02/conDecoyRTarget/diannResult.tsv")
#diann<-fread("diannResult.tsv")
#
#diannForDecoyTarget = diann[,1-min(Q.Value),.(Stripped.Sequence)]
#write.table(diannForDecoyTarget, file="dp/diannForDecoyTarget.csv", sep="\t", row.names = FALSE, quote=FALSE)
#
#
#diannPeptideTD<-fread("~/data1/peptideFDR.txt")
#diannPeptideTD<-fread("peptideFDRHumanR01.txt")
#diannPeptideTD[,.N,.(decoyTarget)]
#head(diannPeptideTD)
#View(diann[,min(Q.Value),.(Stripped.Sequence)])
#
#
#diannPeptideQvalueTarget = diann[,1-min(Q.Value),.(Stripped.Sequence)][diannPeptideTD,on=.(Stripped.Sequence==peptide),decoyTarget:=i.means]
#names(diannPeptideQvalueTarget)[1,2]=c("peptide","ReverseQValueScore")
#View(diannPeptideQvalueTarget)
#diannPeptideQvalueTarget = diannPeptideQvalueTarget[decoyTarget!="BothTargetDeco"][order(-ReverseQValueScore)]
#diannPeptideQvalueTarget = diannPeptideTD[decoyTarget!="BothTargetDeco"][order(-score)]
#
#diannPeptideQvalueTarget[decoyTarget=="NoFound",decoyTarget:="Decoy"]
#diannPeptideQvalueTarget[decoyTarget=="Decoy",decoyTarget:="decoy"]
#diannPeptideQvalueTarget[decoyTarget=="Target",decoyTarget:="target"]
#
#diannPeptideQvalueTarget$nth_test <- 0
#diannPeptideQvalueTarget$resp_target <- 0
#diannPeptideQvalueTarget$resp_decoy <- 0
#Un1 <- unique(sort(diannPeptideQvalueTarget$decoyTarget))
#
#diannPeptideQvalueTarget <-
#  diannPeptideQvalueTarget[,  c("nth_test", paste("resp", Un1, sep = "_")) := c(list(1:.N),
#                                                                           lapply(Un1, function(x)
#                                                                             cumsum(x == decoyTarget))) ,]
#diannPeptideQvalueTarget$fdr <- diannPeptideQvalueTarget$resp_decoy / diannPeptideQvalueTarget$resp_target
##setorder(diannPeptideQvalueTarget,-ReverseQValueScore)
#setorder(diannPeptideQvalueTarget,-score)
##diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < 0.01))][, c(resp_target, nth_test, ReverseQValueScore, ReverseQValueScore)]
#diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < 0.01))][, c(resp_target, nth_test, score, score)]
#
#ourFdrLg005 = dftest2CombineAfter[newCombineScore >= dftest2CombineAfter[max(which(dftest2CombineAfter$fdr < 0.005))]$newCombineScore]
#diannFdrLg005=diannPeptideQvalueTarget[ReverseQValueScore >= diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < 0.005))]$ReverseQValueScore]
#
#
#fdrThreadForOurAndDiann = 0.0065
#ourFdrLg005 = dftest2CombineAfter[newCombineScore >= dftest2CombineAfter[max(which(dftest2CombineAfter$fdr < fdrThreadForOurAndDiann))]$newCombineScore]
##diannFdrLg005=diannPeptideQvalueTarget[ReverseQValueScore >= diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < fdrThreadForOurAndDiann))]$ReverseQValueScore]
#diannFdrLg005=diannPeptideQvalueTarget[score >= diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < fdrThreadForOurAndDiann))]$score]
#
#nrow(ourFdrLg005)
#nrow(diannFdrLg005)
#
#nrow(ourFdrLg005[decoyTarget=="decoy"])
#nrow(diannFdrLg005[decoyTarget=="decoy"])
#
#length(unique(diannPeptideQvalueTarget[ReverseQValueScore >= diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < 0.005))]$ReverseQValueScore]$Stripped.Sequence))
#length(unique(diannPeptideQvalueTarget[ReverseQValueScore >= diannPeptideQvalueTarget[max(which(diannPeptideQvalueTarget$fdr < 0.005))]$ReverseQValueScore]$peptide))
#
#
#nrow(diannPeptideQvalueTarget)
#
#
#library("ggVennDiagram")
#x <- list(
#  Our = ourFdrLg005$peptide,
#  Diann = diannFdrLg005$peptide
##  Diann = diannFdrLg005$Stripped.Sequence
#)
#ggVennDiagram(x)
##01   34063+67066+33603
##0065 28847+64532+32473
##007  29627+65056+32861
##005  26734+63122+31425
#
##humanr01
##005 # 15823+26823+16475
##0065 # 17435+27523+15775
#
#x <- list(
#  Ourdecoy = ourFdrLg005[decoyTarget=="decoy"]$peptide,
#
# # Dianndecoy = diannFdrLg005[decoyTarget=="decoy"]$Stripped.Sequence
#  Dianndecoy = diannFdrLg005[decoyTarget=="decoy"]$peptide
#)
#ggVennDiagram(x)
##01 967+34+895
##0065 578+25+602
##007 633+25+657
##005 432+15+456
#
##humanr01
##005 # 212+0+3
##0065 290+0+3
#View(ourFdrLg005[peptide %in% setdiff(ourFdrLg005$peptide,diannFdrLg005$peptide)])
#View(ourFdrLg005[peptide %in% setdiff(ourFdrLg005$peptide,diannFdrLg005$peptide)][,.N,.(msoneZ)])
#
#
#(212+0+3)/(15823+26823+16475-(212+0+3))
#(15823+26823+16475)/43298
#(17435+27523+15775)/43298
#(290+0+3)/(17435+27523+15775-(290+0+3))
#903/(121281-903)
#(903+15)/(121281-903)
#31268+65858+33177
#
#742+28+760
#(967+34+895)/(34063+67066+33603-(967+34+895))
#(578+25+602)/(28847+64532+32473-(578+25+602))
#(633+25+657)/(29627+65056+32861-(633+25+657))
#121281/100886
#(28847+64532+32473)/100886
#26734+63122+31425 -
#447+471
#
#918/(26734+63122+31425)
#
#125852/100886
#
#ggVennDiagram(x)
##end
############################################################
#
#version
