

# Libraries ---------------------------------------------------------------

library(ggplot2)
library(pROC)
library(randomForest)

# Data Import -------------------------------------------------------------

data <- read.csv('/Users/murra668/Documents/projects/ems_serum_metabolomics/data/cleaned/all/pos/mccu0173_murra668_ems_serum_all_rplc_pos_cleaned_data.csv',
                 stringsAsFactors = F, row.names = 1)
cohort <- read.csv('/Users/murra668/Documents/projects/ems_serum_metabolomics/data/cleaned/all/pos/mccu0173_murra668_ems_serum_all_rplc_pos_cleaned_cohort.csv',
                   stringsAsFactors = F, row.names = 1)

# Remove horses without oral sugar test measurements
data <- data[-which(is.na(cohort$INS_OGT)),]
cohort <- cohort[-which(is.na(cohort$INS_OGT)),]

data <- data[-which(is.na(cohort$INS)),]
cohort <- cohort[-which(is.na(cohort$INS)),]


# Remove horses without laminitis history measrument
data <- data[-which(is.na(cohort$Lam_history)),]
cohort <- cohort[-which(is.na(cohort$Lam_history)),]

# Remove horses without BCS measurment
data <- data[-which(is.na(cohort$BCS)),]
cohort <- cohort[-which(is.na(cohort$BCS)),]

# Remove horses from the Farnley farm
data <- data[-which(cohort$Owner == 'Farnley'),]
cohort <- cohort[-which(cohort$Owner == 'Farnley'),]

# Code insulin resistance 
cohort$IR <- ifelse(cohort$INS_OGT < 45, 'SI', 'IR')
cohort$IR <- as.factor(cohort$IR)

# Clean insulin measures outside test limits
cohort[which(cohort$INS < 1),]$INS <- 1
cohort[which(cohort$INS > 200),]$INS <- 200

cohort[which(cohort$INS_OGT < 1),]$INS_OGT <- 1
cohort[which(cohort$INS_OGT > 200),]$INS_OGT <- 200

# Transform INS measurements to RISQI
cohort$logINS <- log(cohort$INS)
cohort$logINS_OGT <- log(cohort$INS_OGT)

# Remove horses with fewer than 2 horses per farm
x <- table(cohort$Owner)
x <- names(x[which(x < 2)])
data.pre <- data
cohort.pre <- cohort
data <- data.pre[-which(cohort.pre$Owner %in% x),]
cohort <- cohort.pre[-which(cohort.pre$Owner %in% x),]
x <- table(cohort$Owner)
x <- names(x[which(x < 6)])
data1 <- data.pre[-which(cohort.pre$Owner %in% x),]
cohort1 <- cohort.pre[-which(cohort.pre$Owner %in% x),]

data.a <- data[which(cohort$Sampler == 'Nichol'),]
cohort.a <- cohort[which(cohort$Sampler == 'Nichol'),]

data.b <- data[which(cohort$Sampler == 'Elaine'),]
cohort.b <- cohort[which(cohort$Sampler == 'Elaine'),]


met <- colnames(data)
icc_adj <- sapply(met, function(x) do.icc(x, data, cohort$Owner))
icc_full <- sapply(met, function(x) do.icc(x, data1, cohort1$Owner))
icc_a <- sapply(met, function(x) do.icc(x, data.a, cohort.a$Owner))
icc_b <- sapply(met, function(x) do.icc(x, data.b, cohort.b$Owner))

temp <- data.frame(icc = c(icc_full))
ggplot(data = temp, mapping = aes(x = icc)) + geom_histogram(fill = 'green', col = 'darkgreen') + theme_bw()
temp <- data.frame(icc = c( icc_full, icc_adj),
                   group = c(rep('Complete', 1552), rep('N in Farm > 6', 1552)))
ggplot(data = temp, mapping = aes(x = icc, y = ..count.., fill = group)) + geom_histogram(alpha = 0.6, position = 'identity') + theme_bw()
cor.test(icc_full, icc_adj)
wilcox.test(icc_adj, icc_full, conf.int = T)


temp <- data.frame(icc = c( icc_a, icc_b),
                   group = c(rep('Nichol', 1552), rep('Elaine', 1552)))
ggplot(data = temp, mapping = aes(x = icc, y = ..count.., fill = group)) + geom_histogram(alpha = 0.6, position = 'identity', na.rm = T) + theme_bw()
wilcox.test(icc_b, icc_a, conf.int = T)

x <- data.frame(nichol = icc_a, elaine = icc_b)
cor.test(icc_a, icc_b, use = 'na.or.complete')
ggplot(x, aes(x = nichol, y = elaine)) + geom_point() + geom_abline(slope = 1, col = 'red', lty=2, lwd = 1) + theme_bw()

t.m <- sapply(met, function(x) length(which(is.na(data[[x]]))) / length(data[[x]]))
t.test(t.m)
t.a <- sapply(met, function(x) length(which(is.na(data.a[[x]]))) / length(data.a[[x]]))
t.b <- sapply(met, function(x) length(which(is.na(data.b[[x]]))) / length(data.b[[x]]))
wilcox.test(t.a, t.b, conf.int = T)
ggplot(mapping = aes(x = t.m)) + geom_histogram(fill = 'green', col = 'darkgreen') + theme_bw() + xlab("Fraction Missing")
temp <- data.frame(x = c(t.a, t.b), group = c(rep('Nichol', 1552), rep('Elaine', 1552)))
ggplot(data = temp, mapping = aes(x = x, fill = group)) + geom_histogram(alpha = 0.6, position = 'identity', na.rm = T) + theme_bw() + xlab("Fraction Missing")
wilcox.test(t.a, t.b, conf.int = T)
ggplot(mapping = aes(x = t.a, y = t.b)) + geom_point() + geom_abline(slope = 1, col = 'red', lty = 2, lwd =1) + theme_bw() + xlab('Nichol - Fraction Missing') + ylab('Elaine - Fraction Missing')
cor.test(t.a, t.b)

s.m <- 0
for(i in 1:nrow(data)){
  s.m[i] <- length(which(is.na(data[i,]))) / length(data[i,])
}
s.a <- 0
for(i in 1:nrow(data.a)){
  s.a[i] <- length(which(is.na(data.a[i,]))) / length(data.a[i,])
}
s.b <- 0
for(i in 1:nrow(data.b)){
  s.b[i] <- length(which(is.na(data.b[i,]))) / length(data.b[i,])
}
temp <- data.frame(x = c(s.a, s.b), group = c(rep('Nichol', length(s.a)), 
                                              rep('Elaine', length(s.b))))
ggplot(temp, aes(x = x, y = ..density.., fill = group)) + geom_histogram(position = 'identity', alpha = 0.6) + theme_bw() + xlab('Fraction Missing in Sample')
t.test(s.a, s.b)

ggplot(mapping = aes(x = t.m, y = icc_adj)) + geom_point() + geom_abline(slope = 1, col = 'red', lty = 2, lwd = 1) + theme_bw() +
  xlab('Fraction Missing in Feature') + ylab('ICC - Farm')
cor.test(t.m, icc_adj)



gplots::heatmap.2(as.matrix(data),
                  col = my_palette,
                  breaks = col_breaks,
                  trace = 'none', 
                  labRow = FALSE, labCol = FALSE,
                  dendrogram = 'column',
                  ColSideColors = brewer.pal(length(unique((cohort$Owner))), name = 'Dark2'))

my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 30)
col_breaks = c(seq(min(temp2), 0.87,length=10), # for red
               seq(0.88 ,0.93,length=11),  # for yellow
               seq(0.94,1,length=10)) # for green

lam <- names(which(volcano$lam$p <= calcBH(volcano$lam$p)))
ins <- names(which(volcano$ins$p <= calcBH(volcano$ins$p)))
bcs <- names(which(volcano$bcs$p <= calcBH(volcano$bcs$p)))

lam.ins <- intersect(lam, ins)
lam.bcs <- intersect(lam, bcs)
ins.bcs <- intersect(ins, bcs)
lam.ins.bcs <- intersect(lam.ins, bcs)

tempa <- data[which(cohort$Lam_history =='y'),1]
tempb <- data[which(cohort$Lam_history == 'n'),1]

x <- data[,1:776]
x[is.na(x)] <-  1
x <- cor(x)
gplots::heatmap.2(x, col = my_palette, trace = 'none',
                  dendrogram = 'col', 
                  labRow = F, labCol = F)
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 30)
col_breaks = c(seq(min(temp2), 0.87,length=10), # for red
               seq(0.88 ,0.93,length=11),  # for yellow
               seq(0.94,1,length=10)) # for green

# Univariate Analysis -----------------------------------------------------

# Get feature names
met <- colnames(data)

volcano <- list()

# Compute difference of mean for the 3 primary phenotypes of EMS
volcano[['lam']]$p <- sapply(met, function(x) do.mnu(x, cbind(y = cohort$Lam_history, data), 'y'))
volcano[['ins']]$p <- sapply(met, function(x) do.mnu(x, cbind(y = cohort$IR, data), 'IR'))
volcano[['bcs']]$p <- sapply(met, function(x) do.mnu(x, cbind(y = cohort$Obese, data), 'obese'))

# Compute fold change between groupings
volcano[['lam']]$x <- sapply(met, function(x) do.mnu.f(x, cbind(y = cohort$Lam_history, data), 'y'))
volcano[['ins']]$x <- sapply(met, function(x) do.mnu.f(x, cbind(y = cohort$IR, data), 'IR'))
volcano[['bcs']]$x <- sapply(met, function(x) do.mnu.f(x, cbind(y = cohort$Obese, data), 'obese'))

# Compute ICC for each metabolic feature
volcano[['icc']] <- sapply(met, function(x) do.icc(x, data, cohort$Owner))

# Visualize changes
temp <- data.frame(p = c(volcano$lam$p, volcano$ins$p, volcano$bcs$p),
                   pheno = c(rep('lam', 1552), rep('ins', 1552), rep('bcs', 1552)))
ggplot(data = temp, aes(x = p, fill = pheno)) + geom_density(alpha = 0.1)

temp <- data.frame(x = c(volcano$lam$x, volcano$ins$x, volcano$bcs$x),
                   pheno = c(rep('lam', 1552), rep('ins', 1552), rep('bcs', 1552)))
ggplot(data = temp, aes(x = x, fill = pheno)) + geom_density(alpha = 0.1)

ggplot(mapping = aes(x = volcano[['icc']])) + geom_density()


# Plot volcano plots for each phenotype of EMS
plot_volcano(volcano$lam, 'Laminitis History')
plot_volcano(volcano$ins, 'Insulin Resistance')
plot_volcano(volcano$bcs, 'Obesity')

plot_icc(volcano$lam, volcano$icc, 'Laminitis History')
plot_icc(volcano$ins, volcano$icc, 'Insulin Resistance')
plot_icc(volcano$bcs, volcano$icc, 'Obesity')

ggplot(mapping = aes(x = data[['X46.POST']], fill = cohort$Lam_history)) +
  geom_density(alpha = 0.1)
ggplot(mapping = aes(x = cohort$logINS_OGT, fill = cohort$Lam_history)) +
  geom_density(alpha = 0.1)
ggplot(mapping = aes(x = data[['X201.POST']], fill = cohort$IR)) +
  geom_density(alpha = 0.1)
ggplot(mapping = aes(x = data[['X776.PRE']], fill = cohort$Obese)) +
  geom_density(alpha = 0.1)


lam <- which(volcano$lam$p <= calcBH(volcano$lam$p))
ins <- which(volcano$ins$p <= calcBH(volcano$ins$p))
bcs <- which(volcano$bcs$p <= calcBH(volcano$bcs$p))

setdiff(lam, c(ins, bcs))
setdiff(ins, c(lam, bcs))
setdiff(bcs, c(ins, lam))
intersect(lam, ins)
intersect(lam, bcs)
intersect(ins, bcs)
temp <- data[,intersect(lam, c(ins, bcs))]
ggplot(mapping = aes(x = temp[,4], fill = cohort$Lam_history)) + geom_density(alpha = 0.5) + theme_bw()
ggplot(mapping = aes(x = temp[,4], fill = cohort$IR)) + geom_density(alpha = 0.5) + theme_bw()
ggplot(mapping = aes(x = temp[,4], fill = cohort$Obese)) + geom_density(alpha = 0.5) + theme_bw()

temp[is.na(temp)] <- 0
temp$y <- as.factor(cohort$Lam_history)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
seed <- 7
metric <- "Kappa"
set.seed(seed)
mtry <- sqrt(ncol(temp))
tunegrid <- expand.grid(.mtry = sqrt(ncol(data)))
rf_default <- train(y~., data=temp, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_default)
rf_default$finalModel

temp <- cbind(y = cohort$logINS_OGT, data)
temp[is.na(temp)] <- 0
rf <- randomForest(y ~., data = temp, ntree = 1000)
rf

transparentTheme(trans = .2)
library(caret)
temp2 <- cbind(cohort[,c('logINS', 'BCS', 'Leptin')], temp)
temp2 <- temp2[which(cohort$Sampler == 'Elaine'),]
featurePlot(x = temp2[, c(1:3, sample(4:95, 7))], 
            y = as.factor(cohort[which(cohort$Sampler == 'Elaine'),]$Lam_history), 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 2))

gplots::heatmap.2(t(as.matrix(temp)), trace = 'none',
                  labRow = F, labCol = F, RowSideColors = cohort[$Owner)

#
temp <- data
temp[is.na(temp)] <- 0

# Obesity classification
x <- cbind(y = cohort$Obese, temp)
cut <- length(which(cohort$Obese == 'obese')) / nrow(data)
set.seed(10)
rf.o <- randomForest(y ~., data = x, importance = T, cutoff=c(1-cut,cut))
rf.roc.o<- pROC::roc(d$y,rf.o$votes[,2])
# plot(rf.roc.o)
# auc(rf.roc.o)

# History of Laminitis classification
x <- cbind(y = cohort$Lam_history, temp)
cut <- length(which(cohort$Lam_history == 'y')) / nrow(data)
set.seed(10)
rf.l <- randomForest(y ~., data = x, importance = T, cutoff=c(1-cut,cut))
rf.roc.l <- pROC::roc(x$y, rf.l$votes[,2])
# plot(rf.roc.l)
# auc(rf.roc.l)

# Insulin dysregulation classification
x <- cbind(y = cohort$IR, temp)
cut <- length(which(cohort$IR == 'IR')) / nrow(data)
set.seed(10)
rf.i <- randomForest(y ~., data = x, importance = T, ntree = 2000, classwt = c(1e5, 0.01))
rf.rf.roc.i <- pROC::roc(x$y,rf.i$votes[,2])
# plot(rf.roc.i)
# auc(rf.roc.i)

# Plot all 3 ROC curves simulatenous. ROC curves shows the the sensitivity and
# specificity of our classification as function of classification probability
# thresholds
ggplot() +
  geom_path(mapping = aes(x=rf.roc.o$specificities, y=rf.roc.o$sensitivities, colour = paste0('Obesity:         AUC ', format(round(auc(rf.roc.o), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=rf.roc.l$specificities, y=rf.roc.l$sensitivities, colour = paste0('Laminitis:       AUC ', format(round(auc(rf.roc.l), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=rf.roc.i$specificities, y=rf.roc.i$sensitivities, colour = paste0('Insulin:           AUC ', format(round(auc(rf.roc.i), 3), nsmall = 2)))) +  
  geom_abline(intercept = 1, linetype = 2, colour = 'gray') +
  guides(color = guide_legend("Classification Performance")) +
  scale_x_reverse(name = "\nSpecificity", lim=c(1,0), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), minor_breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01))) +
  scale_y_continuous(name = "Sensitivity\n", lim=c(0,1), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), minor_breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01)))+
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.783, 0.115),
        legend.background = element_rect(color = "black", size = 0.25, linetype = "solid"))

cohort$IR_20 <- ifelse(cohort$INS_OGT < 20, 'SI', 'IR')
cohort$IR_25 <- ifelse(cohort$INS_OGT < 25, 'SI', 'IR')
cohort$IR_30 <- ifelse(cohort$INS_OGT < 30, 'SI', 'IR')
cohort$IR_35 <- ifelse(cohort$INS_OGT < 35, 'SI', 'IR')
cohort$IR_40 <- ifelse(cohort$INS_OGT < 40, 'SI', 'IR')
cohort$IR_45 <- ifelse(cohort$INS_OGT < 45, 'SI', 'IR')
cohort$IR_50 <- ifelse(cohort$INS_OGT < 50, 'SI', 'IR')
cohort$IR_55 <- ifelse(cohort$INS_OGT < 55, 'SI', 'IR')
cohort$IR_60 <- ifelse(cohort$INS_OGT < 60, 'SI', 'IR')
cohort$IR_65 <- ifelse(cohort$INS_OGT < 65, 'SI', 'IR')
cohort$IR_70 <- ifelse(cohort$INS_OGT < 70, 'SI', 'IR')
cohort$IR_75 <- ifelse(cohort$INS_OGT < 75, 'SI', 'IR')
cohort$IR_80 <- ifelse(cohort$INS_OGT < 80, 'SI', 'IR')

do.rf <- function(data, cohort, factor){
  
  x <- cbind(y = cohort[[factor]], data)
  cut <- length(which(cohort[[factor]] == 'IR')) / nrow(data)
  set.seed(10)
  rf.i <- randomForest(y ~., data = x, importance = T, cutoff=c(1-cut,cut))
  print(paste(factor, '-', ))
  rf.roc.i <- pROC::roc(x$y,rf.i$votes[,2])
  return(rf.roc.i)
}


ir <- c('IR_20', 'IR_25', 'IR_30', 'IR_35', 'IR_40', 'IR_45', 'IR_50', 'IR_55',
        'IR_60', 'IR_65', 'IR_70', 'IR_75', 'IR_80')

ir.model <- lapply(ir, function(x) do.rf(temp, cohort, x))

ggplot() +
  geom_path(mapping = aes(x=ir.model[[1]]$specificities, y=ir.model[[1]]$sensitivities, colour = paste0('20:       AUC ', format(round(auc(ir.model[[1]]), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=ir.model[[2]]$specificities, y=ir.model[[2]]$sensitivities, colour = paste0('25:       AUC ', format(round(auc(ir.model[[2]]), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=ir.model[[3]]$specificities, y=ir.model[[3]]$sensitivities, colour = paste0('30:       AUC ', format(round(auc(ir.model[[3]]), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=ir.model[[4]]$specificities, y=ir.model[[4]]$sensitivities, colour = paste0('35:       AUC ', format(round(auc(ir.model[[4]]), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=ir.model[[5]]$specificities, y=ir.model[[5]]$sensitivities, colour = paste0('40:       AUC ', format(round(auc(ir.model[[5]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[6]]$specificities, y=ir.model[[6]]$sensitivities, colour = paste0('45:       AUC ', format(round(auc(ir.model[[6]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[7]]$specificities, y=ir.model[[7]]$sensitivities, colour = paste0('50:       AUC ', format(round(auc(ir.model[[7]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[8]]$specificities, y=ir.model[[8]]$sensitivities, colour = paste0('55:       AUC ', format(round(auc(ir.model[[8]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[9]]$specificities, y=ir.model[[9]]$sensitivities, colour = paste0('60:       AUC ', format(round(auc(ir.model[[9]]), 3), nsmall = 2)))) +
  geom_path(mapping = aes(x=ir.model[[10]]$specificities, y=ir.model[[10]]$sensitivities, colour = paste0('65:       AUC ', format(round(auc(ir.model[[10]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[11]]$specificities, y=ir.model[[11]]$sensitivities, colour = paste0('70:       AUC ', format(round(auc(ir.model[[11]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[12]]$specificities, y=ir.model[[12]]$sensitivities, colour = paste0('75:       AUC ', format(round(auc(ir.model[[12]]), 3), nsmall = 2)))) + 
  geom_path(mapping = aes(x=ir.model[[13]]$specificities, y=ir.model[[13]]$sensitivities, colour = paste0('80:       AUC ', format(round(auc(ir.model[[13]]), 3), nsmall = 2)))) + 
  geom_abline(intercept = 1, linetype = 2, colour = 'gray') +
  guides(color = guide_legend("Classification Performance")) +
  scale_x_reverse(name = "\nSpecificity", lim=c(1,0), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), minor_breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01))) +
  scale_y_continuous(name = "Sensitivity\n", lim=c(0,1), breaks = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), minor_breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01)))+
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.783, 0.115),
        legend.background = element_rect(color = "black", size = 0.25, linetype = "solid"))

library(lme4)


temp <- cbind(logINS = cohort$logINS, logINS_OGT = cohort$logINS_OGT, data)
f <- c('logINS', 'logINS_OGT', 'X201.POST')
icc <- sapply(f, function(x) do.icc(x, temp, cohort$Owner))


