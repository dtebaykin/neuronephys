# Compare electrophysiology properties to metadata variations, include neuron types.
library(ggplot2)
library(devtools)
library(ggbiplot)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
library(cluster)
library(HSAUR)
library(httr)
library(RCurl)
library(splines)
library(dplyr)
library(glmnetUtils)
library(magrittr)
library(ggbiplot)
library(ggthemes)

setwd("~/Documents/Neuroelectro documents")
set_config( config( ssl_verifypeer = 0L ))

# Load the data
ions <- c("Na", "Ca", "Mg", "Cl", "K")
rev_ions <- c("Veq_Na", "Veq_Ca", "Veq_Mg", "Veq_Cl", "Veq_K")
solns_orig <- c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Ca", "external_0_Mg", "internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Ca", "internal_0_Mg")
solns <- c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Ca", "external_0_Mg", "internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Ca", "internal_0_Mg",
           "external_0_Cs", "internal_0_Cs", "external_0_glucose", "internal_0_glucose", "external_0_HEPES", "internal_0_HEPES", "external_0_EDTA", "internal_0_EDTA",
           "external_0_EGTA", "internal_0_EGTA", "external_0_BAPTA", "internal_0_BAPTA", "external_0_ATP", "internal_0_ATP", "external_0_GTP", "internal_0_GTP")

new_solns <- c("[Na]external", "[K]external", "[Cl]external", "[Ca]external", "[Mg]external", "[Na]internal", "[K]internal", "[Cl]internal", "[Ca]internal", "[Mg]internal")
ephys_props <- c("rin","rmp","apthr","apamp","aphw","tau", "ahpamp","apwidth","cap","ahpdur","rheo","apfreq",
                 "adratio","sagratio","fahpamp","spike.peak","maxfreq","other", "spontfreq","fislope","apdelay","sahpamp","apriseslope","adpamp",
                 "sagamp","apdecayslope","apriset","fahpdur","apdecayt","accres", "sahpdur","celldiam","mahpamp","mahpdur","adpdur","surfarea")
metadata <- c("Species", "Strain", "ElectrodeType", "PrepType", "JxnPotential", "JxnOffset", "RecTemp", "AnimalAge", "AnimalWeight", "ExternalSolution", "InternalSolution")
log10_ephys <- c("cap", "rin", "tau", "aphw", "rheo", "apwidth", "maxfreq")

thresholds <- data.frame(prop.name = character(), min.range = numeric(), max.range = numeric(), stringsAsFactors = F)

### Get data from the dev neuroelectro website
URL <- "http://dev.neuroelectro.org/static/src/article_ephys_metadata_curated.csv"
ne_all_curated <- read.delim(textConnection(getURL(URL)))

### Or get data from the existing spreadsheet
ne_all_curated <- read.delim("~/Documents/Neuroelectro documents/article_ephys_metadata_curated.csv", stringsAsFactors=FALSE)

### Find out which error values are missing n's, write the corresponding table ids to a csv file
# q <- ne_all_curated[,grepl('_err', colnames(ne_all_curated)) | grepl('_n', colnames(ne_all_curated))]
# q$ArticleID <- ne_all_curated$ArticleID
# q$TableID <- ne_all_curated$TableID
# 
# seq <- sub("_err|_n", "\\1", colnames(q))
# seq <- unique(seq)
# seq <- setdiff(seq, c("ArticleID", "TableID"))
# 
# qa <- data.frame(ArticleID = numeric(0), TableID = numeric(0), err = numeric(0), n = numeric(0))
# 
# for (ep in seq) {
#   q_temp <- q[,c("ArticleID", "TableID", paste0(ep, "_err"), paste0(ep, "_n"))]
#   colnames(q_temp) <- c("ArticleID", "TableID", "err", "n")
#   qa <- rbind(qa, q_temp)
# }
# nrow(qa[which(is.na(qa$n) & !is.na(qa$err)),])
# missing_n <- qa[which(is.na(qa$n) & !is.na(qa$err)),]
# k <- table(missing_n$TableID)
# 
# k <- TableID.with.missing.n[with(TableID.with.missing.n, order(-Count)),]
# 
# write.table(k, "TableID with missing n")
# 
# rm(seq, q, qa, q_temp, ep, k, TableID.with.missing.n)

### Solutions variance overview
ne_all_filter_solns <- ne_all_curated[apply(ne_all_curated[, solns], 1, function(y) !all(is.na(y))),]

### Assign Na's: a small number to avoid 'divide by zero' issue when computing Veq
#rev_data <- subset(ne_all_filter_solns, !is.na(RecTemp))
rev_data <- ne_all_filter_solns
rev_data[,solns][is.na(rev_data[, solns])] <- 0.000001

### Constants (at rest), P's are 'typical' values from neuroscience course / textbook / the internet
Far <- 9.6485 * 10000 / 1000
K <- 273.15
R <- 8.314
P_K <- 1
P_Na <- 0.05
P_Cl <- 0.45

### Calculate reversal potentials
for (ion in ions) {
  new_col <- data.frame(a = R * (K + rev_data[,"RecTemp"]) / Far / 2 * log(rev_data[, paste("external_0_", ion, sep = "")] / rev_data[, paste("internal_0_", ion, sep = "")]))
  colnames(new_col) <- paste("Veq_", ion, sep = "")
  rev_data <- cbind(rev_data, new_col)
}
rm(new_col, ion)

### For now: we are interested in the common relations, so filter out strange solutions
rev_data_filtered <- subset(rev_data, Veq_Na > 0)
rev_data_filtered <- subset(rev_data_filtered, Veq_K < 0)
rev_data_filtered <- subset(rev_data_filtered, AnimalAge > 0)
rev_data_filtered[,log10_ephys] <- log10(rev_data_filtered[,log10_ephys])
rev_data_filtered <- subset(rev_data_filtered,  PrepType == 'in vitro' & (Species == 'Rats' | Species == 'Mice' | Species == 'Guinea Pigs'))
#rev_data_filtered <- subset(rev_data_filtered, ExternalSolution_conf == 5 & InternalSolution_conf == 5)
rev_data_filtered <- droplevels(rev_data_filtered)

# Filter out neuron types that have < 3 mentions


# Create data for models to be used in feature_selection.R
data_models <- rev_data_filtered
data_models$Species <- factor(data_models$Species)

### Remove outliers
#rev_data_filtered <- subset(rev_data_filtered, (ahpdur < 2000 | is.na(ahpdur)) & external_0_Cl < 400 & Veq_K > -100)

### Run this if adjusting for Age and temp
# x <- c("Species", "bs(log10(AnimalAge),df=3)", "bs(RecTemp, df=5)")
# 
# for (ep in ephys_props[1:16]) {
#   curr_formula = formula(paste(ep, " ~", paste(x, collapse="+"))) 
#   
#   glm_mod_fit = cv.glmnet(curr_formula, data = rev_data_filtered, nlambda = 100, alpha=.99, na.action = na.omit)
#   predicted_y = as.vector(predict(glm_mod_fit, newdata = rev_data_filtered, s = "lambda.min"))
#   
#   xdf_new = rev_data_filtered
#   xdf_new[2:nrow(xdf_new),] = xdf_new[1,]
#   xdf_new$NeuronName = rev_data_filtered$NeuronName
#   temp_corrected_y = as.vector(predict(glm_mod_fit, newdata = xdf_new, s = "lambda.min"))
#   corrected_y = rev_data_filtered[ep,] - predicted_y + temp_corrected_y
#   test[ep,] <- corrected_y
# }
# rm(x, xdf_new, curr_formula, temp_corrected_y, corrected_y, predicted_y)

### Facet view of all ephys props by all rev potentials
rev_data_filtered[,log10_ephys] <- 10 ** rev_data_filtered[,log10_ephys]
rev_data_filtered <- subset(rev_data_filtered, (is.na(aphw) | aphw < 4) & (is.na(rheo) | rheo < 500) & (is.na(apamp) | apamp > 30) &
                              (is.na(apdecayt) | apdecayt > 0))
plot_data <- melt(rev_data_filtered, measure = rev_ions, variable.name = "Veq_ion", value.name = "Veq_value")
plot_data <- melt(plot_data, measure = ephys_props, variable.name = "ephys.prop", value.name = "ephys.value")
plot_data <- subset(plot_data, complete.cases(plot_data[,c("Species", "AnimalAge", "RecTemp", "ephys.value")]))
plot_data <- droplevels(plot_data)

ephys_interest <- c("rin","rmp","apthr","apamp","aphw","tau",
                    "rheo","apfreq",
                    "spontfreq","apriseslope","adpamp",
                    "sagamp","apdecayslope","apriset","fahpdur","apdecayt")

#f_linear = formula(y ~ x + Species + log10(AnimalAge) + JxnPotential + RecTemp)
#f_splines = formula(y ~ x + Species + bs(log10(AnimalAge),df=3) + JxnPotential + bs(RecTemp, df=5))
# & NeuronName == "Hippocampus CA1 pyramidal cell"
g <- ggplot(subset(plot_data, Veq_ion == "Veq_K"  & Veq_value > -60 & Veq_value < -40 &  ephys.prop == "rmp"), aes(x = Veq_value, y = ephys.value))
g <- ggplot(subset(plot_data, Veq_ion == "Veq_Cl"  & Veq_value > -10 & Veq_value < 60 &  ephys.prop %in% ephys_props[1:16]), aes(x = Veq_value, y = ephys.value))
g <- ggplot(subset(plot_data, Veq_ion == "Veq_Na"  & Veq_value > 0 & Veq_value < 100 &  ephys.prop %in% ephys_interest), aes(x = Veq_value, y = ephys.value))
g <- ggplot(subset(plot_data, Veq_ion == "Veq_Ca" & Veq_value > 120 & Veq_value < 140  &  ephys.prop == "apamp"), aes(x = Veq_value, y = ephys.value))
g <- ggplot(subset(plot_data, Veq_ion == "Veq_Mg" & Veq_value > -30 & Veq_value < 30  &  ephys.prop %in% ephys_props[1:16]), aes(x = Veq_value, y = ephys.value))

g + geom_point(alpha = 0.3) +
   stat_smooth(method = lm, formula = y ~ x, alpha = 0.2, size = 1, colour = "blue") +
   labs(x = "Veq_K, mV", y = "Resting membrane potential, mV") +
   facet_wrap(~ephys.prop, scales = "free") +
   theme_bw(20)

outliers = function(x, zs) {
  temp <- abs(apply(x, 1, scale))
  return(x[temp > zs])
}

temp <- rev_data_filtered[,solns]
colnames(temp) <- new_solns
temp_data <- melt(temp)
ggplot(subset(temp_data, variable %in% new_solns[1:5])) + geom_density(aes(x = value, fill = variable), alpha = 0.5) + xlab("Concentration, mM") + ylab("Relative density") + theme_bw(20)
ggplot(temp_data) + geom_density(aes(x = value, fill = variable), alpha = 0.5) + xlim(-80, 200)  + xlab("Reversal potentials, mV") + ylab("Relative density") + theme_bw(20)

ggplot(plot_data) + geom_density(aes(x = Veq_value, fill = Veq_ion), alpha = 0.5) + xlim(-80, 200)  + xlab("Reversal potentials, mV") + ylab("Relative density") + theme_bw(20)

pheatmap(t(rev_data_filtered[,rev_ions]), show_colnames = F, fontsize = 14, cluster_rows = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100))

test <- data.frame(scale(rev_data_filtered[,solns[1:5]]))
test1 <- sapply(test, function(x) {abs(x) < 2})

test %<>% filter(apply(test1, 1, function(x) { all(x) }))

remove = apply(test1, 1, function(x) { all(x) })

colnames(test) <- c("External [Na]", "External [K]", "External [Cl]", "External [Ca]", "External [Mg]")
pheatmap(t(test), show_colnames = F, fontsize = 14, cluster_rows = F, clustering_method = "average",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100))


pheatmap(t(rev_data_filtered[,solns[6:10]]), show_colnames = F, fontsize = 14, cluster_rows = F, scale = "column",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100))

temp.hist <- melt(rev_data_filtered[,solns[1:5]])
ggplot(temp.hist) + geom_bar(aes(x = value, fill = variable), binwidth = 0.1) + xlim(0, 10)  + xlab("Concs, mM") + ylab("Relative density") + theme_bw(20)

### Plot Na conc wrt time (years)
ggplot(rev_data_filtered, aes(x = as.factor(PubYear), y = Veq_Na)) +
  geom_point() +
  theme_bw(20)

########### Models
### PCA of Na, K, Cl on Vm
data_pca <- subset(rev_data_filtered, complete.cases(rev_data_filtered[, c("Veq_Na", "Veq_Cl", "Veq_K", "rmp")]))
rownames(data_pca) <- 1:nrow(data_pca)
Vm_pca <- prcomp(data_pca[, c("Veq_Na", "Veq_Cl", "Veq_K", "rmp")], scale = T)

### Plot the data points on PC1-PC2 space
g <- ggbiplot(Vm_pca, obs.scale = 1, var.scale = 1)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + xlim(-2,4) + ylim(-3, 2.5)
print(g)

### Run this if adjusting for Age, Temp and Species
x <- c("Species", "bs(log10(AnimalAge),df=3)", "bs(RecTemp, df=5)")

data_models <- subset(rev_data_filtered, complete.cases(rev_data_filtered[, c("rmp", "Species", "AnimalAge", "RecTemp")]))
data_models$Species <- factor(data_models$Species)
data_models <- data_models[rowSums(is.na(data_models)) < ncol(data_models),]

### Plot mean+/-SEM RMP for Neocortex Martinotti cell and Neostriatum medium spiny neuron
plot_data <- subset(data_models, NeuronName %in% c("Neocortex Martinotti cell", "Neostriatum medium spiny neuron"))
plot_data$NeuronName = as.factor(plot_data$NeuronName)
plot_data <- plot_data[,c("NeuronName", "rmp", "rmp_err", "rmp_n", "rmp_sd", "TableID", solns)]
plot_data <- plot_data[!is.na(plot_data$rmp),]
plot_data <- plot_data[with(plot_data, order(NeuronName, -rmp)),]
plot_data$index <- 1:nrow(plot_data)
plot_data$NeuronName = droplevels(plot_data$NeuronName)
levels(plot_data$NeuronName) = c("Martinotti cells", "Medium spiny neurons")
p = ggplot(plot_data, aes(x = index, y = rmp)) + 
  geom_point(aes(colour = NeuronName)) +
  geom_pointrange(aes(ymin = rmp - rmp_err, ymax = rmp + rmp_err, colour = NeuronName)) +
  scale_color_manual(name = "Neuron type:", values=c("red", "blue"),  guide='legend') +
  labs(x = "", y = "Resting membrane potential, mV") +
  theme_few(15) +
  theme(legend.position="none")
  #theme(axis.text.x = element_text(angle = 15, hjust = 1))
print(p)
rm(plot_data)

ggsave(plot=p,height=6,width=6,dpi=200, filename= paste0(getwd(), "/Plots/example.pdf"), useDingbats=FALSE)


#data_models <- data_models[apply(data_models[, c("external_0_Na", "external_0_K", "external_0_Cl", "internal_0_Na", "internal_0_K", "internal_0_Cl")], 1, function(x) all(x > 0.000001)), ]

### Plot histograms for each patch-clamp, rats/mice with RecTemp listed article.
glist = list()
plot_data <- data_models[,c("external_0_K", "internal_0_K")]
#colnames(plot_data) <- c("External [K]", "Internal [K]")
g = ggplot(melt(plot_data), aes(x = value)) +
  geom_histogram(stat = "bin", aes(fill = variable)) +
  #geom_density(aes(color = variable, fill = variable)) +
  scale_x_continuous(breaks = pretty_breaks(n = 10)) + 
  #ggtitle("Patch-clamp, rats/mice/guinea pigs, in vitro studies") +
  xlim(c(0, 180)) +
  labs(x = "Concentration, mM", y = "Table count") +
  theme_bw(15)
print(g)
rm(plot_data)

glist = c(glist, c(g))

outlier_removal <- subset(data_models, internal_0_K < 100)

### 'Correct' ephys_interest properties using glmnet predictions
for (ep in ephys_interest) {
  ep = "rmp"
  curr_formula = formula(paste(ep, " ~", paste(x, collapse="+"))) 
  
  glm_mod_fit = cv.glmnet(curr_formula, data = data_models, nlambda = 100, alpha=.99, na.action = na.omit)
  predicted_y = as.vector(predict(glm_mod_fit, newdata = data_models, s = "lambda.min"))
  
  xdf_new = data_models
  xdf_new[2:nrow(xdf_new),] = xdf_new[1,]
  xdf_new$NeuronName = data_models$NeuronName
  temp_corrected_y = as.vector(predict(glm_mod_fit, newdata = xdf_new, s = "lambda.min"))
  data_models[,ep] = data_models[,ep] - predicted_y + temp_corrected_y
}
rm(x, xdf_new, curr_formula, temp_corrected_y, predicted_y)

### Let's try predicting RMPs (without Veq's)
rmp_predicted <- with(data_models, (RecTemp + K) * R / Far * log((P_Na * external_0_Na + P_K * external_0_K + P_Cl * internal_0_Cl) / 
                                                           (P_Na * internal_0_Na + P_K * internal_0_K + P_Cl * external_0_Cl)))

library(MASS)
par(oma = c(1,2,0,0))
plot(sort(data_models$rmp / rmp_predicted), ylab = "", yaxt = "n", las = 1, ylim = c(-1, 5))
abline(h=1, lty=2)

my.at = c( -1 : 5)
my.labels <- c(as.character(fractions(-1/(my.at[ my.at < 1 ] -2))), my.at[my.at >= 1])
axis(side = 2, at = my.at, labels = my.labels, las = 1)
mtext(side = 2, line = 3, text = 'Correction factor: actual / predicted')

data_models$rmp_predicted <- rmp_predicted

### Expect Neostriatum medium spiny neuron to have low RMP - check
### Expect Neocortex Martinotti cell to have high RMP - check => Ion Permeability matters

### Plot rmp measures vs rmp predicted (just for CA1)
ggplot(subset(data_models, NeuronName == "Hippocampus CA1 pyramidal cell"  & rmp_predicted < -40), aes(x = rmp_predicted, y = rmp, color = JxnPotential)) + geom_point(aes(size = 2)) + theme_bw(15)
summary(lm("rmp ~ rmp_predicted", subset(data_models, NeuronName == "Hippocampus CA1 pyramidal cell" & rmp_predicted < -40)))

ggplot(subset(data_models, NeuronName == "Hippocampus CA1 pyramidal cell"), aes(x = internal_0_Cl, y = rmp, color = JxnPotential)) + geom_point(aes(size = 2)) + theme_bw(15)


### Correct uncorrected jxn using formula: Vm_corr = Vm_old - |jxn| (if jxn not measured - assume 10 mV)
data_models$rmp_corr = data_models$rmp 
data_models$rmp_corr = as.numeric(apply(data_models, 1 , function(x) {
  if (x['JxnPotential'] == "Uncorrected") {
    if (is.na(x['JxnOffset'])) {
      as.numeric(x['rmp_unc']) - 10
    } else {
      as.numeric(x['rmp_unc']) - abs(as.numeric(x['JxnOffset']))
    }  
  } else {
    x['rmp_corr']
  }
}))

temp <- table(data_models$NeuronName)
data_models_filtered <- subset(data_models, NeuronName %in% c(names(temp[temp > 9]), "Neostriatum medium spiny neuron") & NeuronName != "Other" & NeuronName != "Neocortex uncharacterized cell")
rm(temp)
ggplot(data_models_filtered, aes(x = as.numeric(rmp_unc), y = rmp_predicted, colour = NeuronName)) +
  geom_point(cex = 5) +
  xlim(-90, -50) +
  ylim(-90, -50)

cor(data_models$rmp_predicted, data_models$rmp_corr)

### Try to optimizing permeability for each cell type, get the upper bound on the model explanatory strength
library("RCEIM") # "Let's shoot in the dark and hope to hit something good" library
library("lbfgs") # stats library
### Try just the default R stats package optimization method (optim)

result <- data.frame(NeuronName = character(0), P_Na = numeric(0), P_K = numeric(0), P_Cl = numeric(0), 
                                P_Na_bfgs = numeric(0), P_K_bfgs = numeric(0), P_Cl_bfgs = numeric(0),
                                P_Na_sann = numeric(0), P_K_sann = numeric(0), P_Cl_sann = numeric(0))
for (nt in unique(data_models_filtered$NeuronName)) {
  model_data <- subset(data_models, NeuronName == nt)
  optFunc <- function(x) {
    dist_matrix <- matrix(nrow = 2, ncol = 0)
    for (row in nrow(model_data)) {
      predicted_rmp <- with(model_data[row,], (RecTemp + K) * R / Far * log((x[1] * external_0_Na + x[2] * external_0_K + x[3] * internal_0_Cl) / 
                                                          (x[1] * internal_0_Na + x[2] * internal_0_K + x[3] * external_0_Cl)))
      extracted_rmp <- model_data[row, "rmp_corr"]
      dist_matrix <- cbind(dist_matrix, c(predicted_rmp, extracted_rmp))
    }
    
    dist(dist_matrix, method = "euclidean")[1]
  }
  
  sann_output <- optim(c(0.05, 1, 0.45), optFunc, method = "SANN", control = list(trace = 1))
  bfgs_output <- optim(c(0.05, 1, 0.45), optFunc, method = "L-BFGS-B", lower = 0, upper = 1, control = list(trace = 1))
  
  output <- ceimOpt(OptimFunction = "optFunc", nParams = 3, minimize = T, maxIter = 1000, epsilon = 0.05, waitGen = 20, boundaries = matrix(c(0,0,0,1,1,1), nrow = 3, ncol = 2), plotResultDistribution = T, verbose = T, parallelVersion = T)
  result <- rbind(result, data.frame(NeuronName = nt, P_Na = output$BestMember[1], P_K = output$BestMember[2], P_Cl = output$BestMember[3],
                                             P_Na_bfgs = bfgs_output$par[1], P_K_bfgs = bfgs_output$par[2], P_Cl_bfgs = bfgs_output$par[3],
                                             P_Na_sann = sann_output$par[1], P_K_sann = sann_output$par[2], P_Cl_sann = sann_output$par[3]))
}
rm(optFunc, output, model_data, sann_output, bfgs_output)

data_models_filtered$rmp_predicted_best <- apply(data_models_filtered, 1, function(x) {
  res <- subset(result, NeuronName == x["NeuronName"])
  (as.numeric(x["RecTemp"]) + K) * R / Far * log( (res$P_Na * as.numeric(x["external_0_Na"]) + res$P_K * as.numeric(x["external_0_K"]) + res$P_Cl * as.numeric(x["internal_0_Cl"])) / 
                                                  (res$P_Na * as.numeric(x["internal_0_Na"]) + res$P_K * as.numeric(x["internal_0_K"]) + res$P_Cl * as.numeric(x["external_0_Cl"])) )
})

ggplot(subset(data_models_filtered, NeuronName == "Hippocampus CA1 pyramidal cell"), aes(x = as.numeric(rmp_corr), y = rmp_predicted_best)) +
  geom_point(cex = 5) 

### Try multivariate regressions with RMP vs jxn, jxn offset, solutions (just for CA1 pyramidal cells and in general)
test <- subset(rev_data, NeuronName == "Hippocampus CA1 pyramidal cell" & JxnPotential %in% c("Corrected", "Not corrected", "Unreported"))
test_lm <- lm(rmp ~ JxnPotential, test)
summary(test_lm)

temp <- subset(test, !is.na(JxnOffset) & JxnPotential == "Not corrected")
temp$rmp_new <- temp$rmp - abs(temp$JxnOffset)

test$rmp_new <- test$rmp
test[rownames(temp), "rmp_new"] <- temp[,"rmp_new"]
test[rownames(temp), "JxnPotential"] <- "Recorrected"

ggplot(test, aes(x = JxnPotential, y = as.numeric(rmp_new), colour = JxnPotential)) + geom_point(cex = 5) + ggtitle("Hippocampus CA1 pyramidal cells") + theme_bw(15)

summary(lm(rmp ~ JxnOffset, jxn_data))

### Try making a model for Jxn potential offset prediction based on solutions
model_jxn = as.formula(paste("JxnOffset ~ internal_0_K + internal_0_Cl"))
jxn_data = subset(rev_data, !is.na(JxnOffset))
jxn_data$JxnOffset <- abs(jxn_data$JxnOffset)
summary(lm(model_jxn, data = jxn_data))

ggplot(jxn_data, aes(x = internal_0_K, y = JxnOffset)) + geom_point(size = 2) + stat_smooth(method = "lm")

summary(lm("JxnOffset ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8", ca1_corr_dat))

ggplot(test, aes(x = internal_0_Cl, y = rmp_new, color = JxnPotential)) + geom_point(size = 2)

model_rmp_ca1 = as.formula(paste("rmp ~ RecTemp + ", paste(solns, collapse = "+")))
summary(lm(model_jxn, data = test))

jxn_pca <- prcomp(jxn_data[, solns <- c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Mg", "internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Mg")], scale = T)

### Plot the data points on PC1-PC2 space
g <- ggbiplot(jxn_pca, obs.scale = 1, var.scale = 1)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + xlim(-2,3) + ylim(-2, 2.5)
print(g)

### Try just corrected and recorrected jxn potential points for the GHK model comparison for CA1 pyramidal cells
ca1_corr_dat <- test[, solns <- c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Mg", "internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Mg")]
ca1_pca <- prcomp(ca1_corr_dat, scale = T)
### Plot the data points on PC1-PC2 space
g <- ggbiplot(ca1_pca, obs.scale = 1, var.scale = 1)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

(loadings <- ca1_pca$rotation)
summary(ca1_pca)
ca1_corr_dat <- cbind(test, predict(ca1_pca, newdata = ca1_corr_dat))

ggplot(ca1_corr_dat, aes(x = PC6, y = rmp)) + geom_jitter(aes(size = 2), position = position_jitter(width = 0.3)) + stat_smooth(method = "lm")
summary(lm("rmp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + JxnPotential", ca1_corr_dat))
hist(ca1_corr_dat$internal_0_Mg)

cor(ca1_corr_dat$rmp, ca1_corr_dat$internal_0_Mg, use="complete.obs")

### Identify 'schools of thought' for creating external and internal solutions
library(cluster)

solns_data <- subset(rev_data, !duplicated(ArticleID))
solns_data <- subset(solns_data, ElectrodeType == "Patch-clamp")

# External solutions outliers
solns_data <- subset(solns_data, external_0_Cl < 200)
solns_data <- subset(solns_data, external_0_Na < 200)
solns_data <- subset(solns_data, external_0_K < 20)
solns_data <- subset(solns_data, external_0_Mg < 8.1)
solns_data <- subset(solns_data, external_0_Ca < 3)

# Internal solutions outliers
solns_data <- subset(solns_data, internal_0_Cl < 200)
solns_data <- subset(solns_data, internal_0_Na < 100)
solns_data <- subset(solns_data, internal_0_K < 200)
solns_data <- subset(solns_data, internal_0_Mg < 8.1)
solns_data <- subset(solns_data, internal_0_Ca < 2)

# for (x in c("internal_0_Na", "internal_0_K", "internal_0_Cl")) {
#   solns_data <- solns_data[which(
#     solns_data[,x] > quantile(solns_data[,x], 0.01) &
#     solns_data[,x] < quantile(solns_data[,x], 0.99)
#   ),]
# }

plot(agnes(data_models[, solns_orig]))
plot(pam(data_models[, solns_orig], k = 4))

solns_data_pca <- prcomp(solns_data[,c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Mg", "external_0_Ca")], center = T, scale. = T, tol = sqrt(.Machine$double.eps))
solns_data_pca <- prcomp(solns_data[,c("internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Mg", "internal_0_Ca")], center = T, scale. = T, tol = sqrt(.Machine$double.eps))

### Plot the data points on PC1-PC2 space
g <- ggbiplot(solns_data_pca, obs.scale = 1, var.scale = 1)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + theme_bw(20)
g

bk = c(seq(0, 5, length = 6), 10, seq(15, 30, length = 10), 50, seq(60, 150, length = 10), 300)
colors = c(colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(6), colorRampPalette(brewer.pal(n = 7, name = "Purples"))(10), colorRampPalette(brewer.pal(n = 7, name = "BrBG"))(12))
pheatmap(solns_data[, c("internal_0_Na", "internal_0_K", "internal_0_Cl", "internal_0_Mg", "internal_0_Ca")],
         color = colors, breaks = bk, legend_breaks = c(0, 10, 20, 30, 40, 50, 100, 150, 200), cluster_cols = F, show_rownames = F)

dupes = solns_data[duplicated(solns_data[,c("external_0_Na", "external_0_K", "external_0_Cl", "external_0_Mg", "external_0_Ca")]),]

# Uses plot_data from Martinotti cells vs Medium Spiny neurons, near line 223
library(pdist)

x.pdist = as.data.frame(as.matrix(pdist(plot_data[,c("internal_0_Na", "internal_0_K", "internal_0_Cl", "external_0_Na", "external_0_K", "external_0_Cl")], indices.A = 12:19, indices.B = c(5, 7, 21, 25, 31, 40, 54))))
y.pdist = as.data.frame(as.matrix(pdist(plot_data[,c("internal_0_Na", "internal_0_K", "internal_0_Cl", "external_0_Na", "external_0_K", "external_0_Cl")], indices.A = 12:19, indices.B = 12:19)))
# Not really different solutions...

ggplot(data_models, aes(x = as.character(PubYear), y = internal_0_Na)) + geom_boxplot() + labs(x = "Year published", y = "Internal [Na], mM") + theme_fivethirtyeight(20)

# To attach text to ggplot points: geom_text(aes(wt, mpg, label = rownames(mtcars))) +