######################
## RDA Fagus APPATS ##
######################

library(vegan)
library(robust)
library(qvalue)
library(ade4)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(adegenet)
library(raster)
library(pcadapt)
library(LEA)
library(RColorBrewer)
library(reshape2)
library(spatstat)
library(geosphere)
library(Imap)
library("rgdal")
library("maptools")
library("plyr")

##############
#### Data ####

# Genet
GenData<-as.data.frame(t(read.table("batch_1.plink.pcadapt")))
row.names(GenData)<- unlist(lapply(strsplit(as.character(read.table("batch_1.plink.ped", header=F)[,2]), ".", fixed = T), function(x) x[1]))
colnames(GenData) <- as.character(read.table("batch_1.plink.map", header=F)[,2])

# Removing unused population of Ste Baumes and replicates
GenData<-GenData[-c(grep(pattern = "Bis", row.names(GenData)), grep(pattern = "Ter", row.names(GenData))),]
GenData<-GenData[-c(grep(pattern = "SB", row.names(GenData))),]

# Samples information
info_inds<-read.table("Infos_Individus.csv", header=T, sep=",")
info_inds<-info_inds[match(row.names(GenData), info_inds[,1]),]

# Climatic data
fich_biocl<-list.files(path='./', pattern='.img$',full.names=T)
var_env<-projectRaster(stack(fich_biocl), crs=proj4string(raster("./bio3_16.tif")))
var_env_df<-as.data.frame(extract(var_env, info_inds[,c(6,5)]))

# Ancestry
write.geno(GenData, "genotypes.geno")
object <- snmf("genotypes.geno", K = 2, entropy = T, ploidy = 2, alpha = 100, project = "new", repetitions = 10)

Q<-Q(object, K=2, run=which.min(cross.entropy(object, K = 2)))

## Variables 
EnvData <- data.frame(info_inds[,c(6,5)], var_env_df, Q=Q[,2])
EnvData <- as.data.frame(aggregate(EnvData, by=list(info_inds$pop), mean))

# Correlation among variables
pdf("./Pairs_Fagus.pdf")
pairs(EnvData[,-1], lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = panel.hist)
dev.off()

## Imputing missing data
for (i in 1:ncol(GenData))
{
  GenData[which(GenData[,i]==9),i]<-median(GenData[-which(GenData[,i]==9),i], na.rm=TRUE)
}

## Allelic frequencies
FreqData<-as.data.frame(apply(GenData, 2, function(x) by(x, as.character(info_inds$pop), mean))/2)
freq_mean<-apply(FreqData, 2, mean)
FreqData<-FreqData[,-which(freq_mean>=0.99 | freq_mean<=0.01)]

#####################
#### partial RDA ####

## Model climate
pRDA1 <- rda(FreqData ~ etp_mean + mind_mean + prec_sum + tmax_mean + tmin_mean + Condition(X + Y + Q),  EnvData)
anova(pRDA1)
pRDA1
RsquareAdj(pRDA1)

## Model geography
pRDA2 <- rda(FreqData ~ X + Y + Condition(etp_mean + mind_mean + prec_sum + tmax_mean + tmin_mean + Q),  EnvData)
anova(pRDA2)
pRDA2
RsquareAdj(pRDA2)

## Modele ancestry
pRDA3 <- rda(FreqData ~ Q + Condition(etp_mean + mind_mean + prec_sum + tmax_mean + tmin_mean + Y + X),  EnvData)
anova(pRDA3)
pRDA3
RsquareAdj(pRDA3)

## Full model
RDAfull <- rda(FreqData ~ etp_mean + mind_mean + prec_sum + tmax_mean + tmin_mean + X + Y + Q,  EnvData)
anova(RDAfull)
RsquareAdj(RDAfull)

#################
#### RDAdapt ####

RDA_pop <- rda(FreqData ~ etp_mean + mind_mean + prec_sum + tmax_mean + tmin_mean + Condition(Q),  EnvData)

rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

pdf("./ScreePlot_RDA_Fagus.pdf", width = 4, height = 3)
ggplot() +
  geom_line(aes(x=1:length(RDA_pop$CCA$eig), y=as.vector(RDA_pop$CCA$eig)), linetype="dotted", size = 1.5) +
  geom_point(aes(x=1:length(RDA_pop$CCA$eig), y=as.vector(RDA_pop$CCA$eig)), size = 3) +
  scale_x_discrete(name = "Ordination axes", limits=c(1:10)) +
  ylab("Inertia") +
  theme_bw()
dev.off()

rdadapt_Fagus_pop<-rdadapt(RDA_pop, 3)

## Visualization
pdf("./P-values_RDA_Fagus.pdf", width = 5, height = 4)
ggplot() +
  geom_point(aes(x=as.vector(which(rdadapt_Fagus_pop$q.values<1)), y=-log10(rdadapt_Fagus_pop$p.values)), col = "gray86") +
  geom_point(aes(x=as.vector(which(rdadapt_Fagus_pop$q.values<0.0001)), y=-log10(rdadapt_Fagus_pop$p.values)[which(rdadapt_Fagus_pop$q.values<0.0001)]), colour = "orange", size = 2.5) +
  xlab("SNPs") + ylab("RDA -log(p-values)") +
  theme_bw() +
  theme(legend.position="none")
dev.off()

# Outliers
locilist.names <- colnames(FreqData)
outliers_rdadapt_pop <- locilist.names[which(rdadapt_Fagus_pop$q.values < 0.0001)]

# Projection SNPs and variables
p1 <- ggplot() +
  geom_point(aes(x=RDA_pop$CCA$v[,1], y=RDA_pop$CCA$v[,2]), col = "gray83") +
  geom_point(aes(x=RDA_pop$CCA$v[which(rdadapt_Fagus_pop$q.values < 0.0001),1], y=RDA_pop$CCA$v[which(rdadapt_Fagus_pop$q.values < 0.0001),2]), col = "orange") +
  geom_segment(aes(xend=RDA_pop$CCA$biplot[,1]/10, yend=RDA_pop$CCA$biplot[,2]/10, x=0, y=0), colour="black", size=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*RDA_pop$CCA$biplot[,1]/10, y=1.2*RDA_pop$CCA$biplot[,2]/10, label = colnames(EnvData[,c(6:10)]))) +
  xlab("RDA 1") + ylab("RDA 2") +
  theme_bw() +
  theme(legend.position="none")

p2 <- ggplot() +
  geom_point(aes(x=RDA_pop$CCA$v[,1], y=RDA_pop$CCA$v[,3]), col = "gray83") +
  geom_point(aes(x=RDA_pop$CCA$v[which(rdadapt_Fagus_pop$q.values < 0.0001),1], y=RDA_pop$CCA$v[which(rdadapt_Fagus_pop$q.values < 0.0001),3]), col = "orange") +
  geom_segment(aes(xend=RDA_pop$CCA$biplot[,1]/10, yend=RDA_pop$CCA$biplot[,3]/10, x=0, y=0), colour="black", size=0.5, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*RDA_pop$CCA$biplot[,1]/10, y=1.2*RDA_pop$CCA$biplot[,3]/10, label = colnames(EnvData[,c(6:10)]))) +
  xlab("RDA 1") + ylab("RDA 3") +
  theme_bw() +
  theme(legend.position="none")

pdf("./Projection_LociRDA_Fagus.pdf", width = 10, height = 5)
plot_grid(p1, p2, labels = "AUTO", align = "h")
dev.off()

## Outliers allelic frequency in populations
toto1<-melt(data.frame(FreqData[,outliers_rdadapt_pop], RDA = RDA_pop$CCA$u[,1]), id="RDA")
toto1$axe<-"RDA1"
toto2<-melt(data.frame(FreqData[,outliers_rdadapt_pop], RDA = RDA_pop$CCA$u[,2]), id="RDA")
toto2$axe<-"RDA2"

p1 <- ggplot(data = toto1, aes(x = variable, y = value)) +
  geom_point(aes(colour = RDA, size = RDA)) +
  scale_color_continuous(low="white", high="black") +
  labs(colour = "RDA1", size = "RDA1") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 8))

p2 <- ggplot(data = toto2, aes(x = variable, y = value)) +
  geom_point(aes(colour = RDA, size = RDA)) +
  scale_color_continuous(low="white", high="black") +
  labs(colour = "RDA2", size = "RDA2") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 8))

pdf("./Freq_outliers_Fagus.pdf", width = 8, height = 8)
ggarrange(p1,p2, ncol = 1, nrow = 2)
dev.off()

cor1<-NULL
cor2<-NULL
pv1<-NULL
pv2<-NULL
for (i in outliers_rdadapt_pop) {
  cor1<-c(cor1, summary(lm(FreqData[,i]~RDA_pop$CCA$u[,1]))$r.squared)
  pv1<-c(pv1, summary(lm(FreqData[,i]~RDA_pop$CCA$u[,1]))$coefficients[2,4])
  cor2<-c(cor2, summary(lm(FreqData[,i]~RDA_pop$CCA$u[,2]))$r.squared)
  pv2<-c(pv2, summary(lm(FreqData[,i]~RDA_pop$CCA$u[,2]))$coefficients[2,4])
}

## Delta allelic frequencies
hist(apply(FreqData, 2, function(x) max(x)-min(x)))
abline(v=apply(FreqData[, outliers_rdadapt_pop], 2, function(x) max(x)-min(x)))

delta <- apply(FreqData[,-which(colnames(FreqData)%in%outliers_rdadapt_pop)], 2, function(x) max(x)-min(x))
delta_out <- apply(FreqData[, outliers_rdadapt_pop], 2, function(x) max(x)-min(x))
delta_tot <- data.frame(delta = c(delta, delta_out), loci = c(rep("neutral", length(delta)), rep("outliers", length(delta_out))))
pdf("./Delta_Freq_Fagus.pdf")
ggplot(data=delta_tot) + 
  geom_histogram(aes(x=delta, fill=as.character(loci)), binwidth = 0.1, col='black') +
  stat_bin(aes(x=delta, y=..count.., label=..count.., fill=as.character(loci)), geom="text", binwidth = 0.1, vjust=-.5) +
  scale_fill_manual(values = c("gray83", "orange")) +
  guides(fill=guide_legend(title="Loci")) +
  xlab("Delta allelic frequencies")
dev.off()

#########################
#### RDA on outliers ####

RDA_enrichspace <- rda(FreqData[, outliers_rdadapt_pop] ~ etp_mean + mind_mean + prec_sum + tmax_mean + tmin_mean + Condition(Q),  EnvData)

# Projection SNPs and variables
pdf("./Projection_EnrichedLoci_Fagus.pdf", width = 4, height = 4)
ggplot() +
  geom_point(aes(x=RDA_enrichspace$CCA$v[,1], y=RDA_enrichspace$CCA$v[,2]), col = "orange", size = 3) +
  geom_segment(aes(xend=RDA_enrichspace$CCA$biplot[,1], yend=RDA_enrichspace$CCA$biplot[,2], x=0, y=0), colour="black", size=0.75, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(aes(x=1.2*RDA_enrichspace$CCA$biplot[,1], y=1.2*RDA_enrichspace$CCA$biplot[,2], label = row.names(RDA_enrichspace$CCA$biplot)), size = 3.5) +
  xlab("RDA 1") + ylab("RDA 2") +
  xlim(-1,1) + ylim(-1,1) +
  theme_bw() +
  theme(legend.position="none")
dev.off()


############################
#### Adaptive gradients ####

### scenario hc_hadrm3q0_hadcm3_ar3/a1b/clim5180/

# Variables in current conditions
fich_biocl<-list.files(path='./', pattern='.img$',full.names=T)
var_env<-projectRaster(stack(fich_biocl), crs=proj4string(raster("./bio3_16.tif")))
var_env<-crop(var_env, extent(projectRaster(raster("./alt_FA.tif"), crs=proj4string(var_env[[1]]))))
var_env_proj<-as.data.frame(extract(var_env[[-c(1,2,8)]], xyFromCell(var_env, 1:ncell(var_env)), cellnumbers=T))

# Variables in futur conditions
fich_biocl<-list.files(path='./', pattern='.img$',full.names=T)
var_env_fu<-projectRaster(stack(fich_biocl), var_env[[1]])
var_env_fu<-crop(var_env_fu, extent(projectRaster(raster("./alt_FA.tif"), crs=proj4string(var_env[[1]]))))
var_env_proj_fu<-as.data.frame(extract(var_env_fu[[-1]], xyFromCell(var_env_fu, 1:ncell(var_env_fu)), cellnumbers=T))

# Standardisation 
var_env_proj_tot<-as.data.frame(apply(rbind(var_env_proj, var_env_proj_fu), 2, scale))
var_env_proj <- var_env_proj_tot[1:nrow(var_env_proj),]
var_env_proj_fu <- var_env_proj_tot[(nrow(var_env_proj)+1):nrow(var_env_proj_tot),]

colpal=colorRampPalette(c("white", "black"))

cities<-crop(shapefile("./cities.shp"), extent(c(4.7, 6.4, 42.6565, 46.85202)))

# Projection present futur RDA1
RDA1_proj <-apply(var_env_proj[,-c(1)],1, function(x) sum( x * RDA_enrichspace$CCA$biplot[,1]))
rda1<-var_env[[1]]
rda1[]<-as.vector(RDA1_proj)
raster_rda1<-rasterToPoints(crop(rda1, extent(4.7, 8.199158, 42.6565, 46.85202)))
p1 <- ggplot() + 
  geom_raster(aes(x=raster_rda1[,1],y=raster_rda1[,2],fill=cut(raster_rda1[,3], breaks=c(-6,-3,0,3,6,9,12,15)))) + 
  scale_fill_manual(values=colpal(8)[-1]) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA scores")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

RDA1_proj_fu <-apply(var_env_proj_fu[,-c(1)],1, function(x) sum( x * RDA_enrichspace$CCA$biplot[,1]))
rda1_fu<-var_env[[1]]
rda1_fu[]<-as.vector(RDA1_proj_fu)
raster_rda1_fu<-rasterToPoints(crop(rda1_fu, extent(4.7, 8.199158, 42.6565, 46.85202)))
p2 <- ggplot() + 
  geom_raster(aes(x=raster_rda1_fu[,1],y=raster_rda1_fu[,2],fill=cut(raster_rda1_fu[,3], breaks=c(-6,-3,0,3,6,9,12,15,18)))) + 
  scale_fill_manual(values=colpal(9)[-1]) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA scores")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

# Projection present futur RDA2
RDA2_proj <-apply(var_env_proj[,-c(1)],1, function(x) sum( x * RDA_enrichspace$CCA$biplot[,2]))
rda2<-var_env[[1]]
rda2[]<-as.vector(RDA2_proj)
raster_rda2<-rasterToPoints(crop(rda2, extent(4.7, 8.199158, 42.6565, 46.85202)))
p3 <- ggplot() + 
  geom_raster(aes(x=raster_rda2[,1],y=raster_rda2[,2],fill=cut(raster_rda2[,3], breaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5)))) + 
  scale_fill_manual(values=colpal(8)[-1]) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA scores")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

RDA2_proj_fu <-apply(var_env_proj_fu[,-c(1)],1, function(x) sum( x * RDA_enrichspace$CCA$biplot[,2]))
rda2_fu<-var_env[[1]]
rda2_fu[]<-as.vector(RDA2_proj_fu)
raster_rda2_fu<-rasterToPoints(crop(rda2_fu, extent(4.7, 8.199158, 42.6565, 46.85202)))
p4 <- ggplot() + 
  geom_raster(aes(x=raster_rda2_fu[,1],y=raster_rda2_fu[,2],fill=cut(raster_rda2_fu[,3], breaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5)))) + 
  scale_fill_manual(values=colpal(8)[-1]) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="RDA scores")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())


pdf("./RDAindex_FuturePresent_Fagus.pdf", width = 8, height = 8)
ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = FALSE)
dev.off()

##############################################################################################################################
#### Distribution of RDA1 and RDA2 scores in current and futur conditions across the current locations of Fagus sylvatica ####

occurences <- read.csv("./OccurenceFagus.csv",header = T, sep=",", dec=".")

scores_rda1 <- extract(rda1, occurences)
scores_rda1_fu <- extract(rda1_fu, occurences)
table <- data.frame(scores=c(scores_rda1, scores_rda1_fu), cond = as.factor(c(rep("present", length(scores_rda1)), rep("futur", length(scores_rda1)))))
p1 <- ggplot() + 
  geom_density(data=table, aes(x=scores, fill=cond), alpha=.7) +
  scale_fill_manual(values = c("black", "white")) +
  xlim(c(-5,10)) +
  geom_vline(xintercept=quantile(extract(rda1, occurences), probs= c(0.01, 0.99)),  colour=c("black","black"), linetype="dashed", size=0.5) +
  xlab("Scores RDA1") +
  theme_bw()

scores_rda2 <- extract(rda2, occurences)
scores_rda2_fu <- extract(rda2_fu, occurences)
table <- data.frame(scores=c(scores_rda2, scores_rda2_fu), cond = as.factor(c(rep("present", length(scores_rda2)), rep("futur", length(scores_rda2)))))
p2 <- ggplot() + 
  geom_density(data=table, aes(x=scores, fill=cond), alpha=.7) +
  scale_fill_manual(values = c("black", "white")) +
  xlim(c(-3,2)) +
  geom_vline(xintercept=quantile(extract(rda2, occurences), probs= c(0.01, 0.99)),  colour=c("black","black"), linetype="dashed", size=0.5) +
  xlab("Scores RDA2") +
  theme_bw()


pdf("./RDAindex_FuturePresent_Fagus.pdf", width = 4, height = 4)
ggarrange(p1,p2, ncol=1, nrow=2, common.legend = TRUE, legend = "right")
dev.off()

#########################################
#### Favourable pixels present/futur ####

rda1<-var_env[[1]]
rda1[]<-as.vector(RDA1_proj)
rda2<-var_env[[1]]
rda2[]<-as.vector(RDA2_proj)

fav_pres_RDA<-var_env[[1]]
fav_pres_RDA[]<-NA
fav_pres_RDA[which(rda1[]>quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[1] & rda1[]<quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[2] & rda2[]>quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[1] & rda2[]<quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[2])]<-1

rda1_fu<-var_env[[1]]
rda1_fu[]<-as.vector(RDA1_proj_fu)
rda2_fu<-var_env[[1]]
rda2_fu[]<-as.vector(RDA2_proj_fu)

fav_fut_RDA<-var_env[[1]]
fav_fut_RDA[]<-NA
fav_fut_RDA[which(rda1_fu[]>quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[1] & rda1_fu[]<quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[2] & rda2_fu[]>quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[1] & rda2_fu[]<quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[2])]<-1

plot(fav_pres_RDA)
plot(fav_fut_RDA)

coord_fu<-as.data.frame(rasterToPoints(fav_fut_RDA))
coord_fu[,3]<-"futur"
coord_pres<-as.data.frame(rasterToPoints(fav_pres_RDA))
coord_pres[,3]<-"present"
table<-rbind(coord_pres, coord_fu)

nochange<-var_env[[1]]
nochange[]<-NA
nochange[][which(fav_pres_RDA[]==1 & fav_fut_RDA[]==1)]<-1

lost<-var_env[[1]]
lost[]<-NA
lost[][which(fav_pres_RDA[]==1 & is.na(fav_fut_RDA[]))]<-1

new<-var_env[[1]]
new[]<-NA
new[][which(is.na(fav_pres_RDA[]) & fav_fut_RDA[]==1)]<-1

country<-crop(shapefile("./country.shp"), extent(c(4.6,8,42.8,46.7)))
cities<-crop(shapefile("./cities.shp"), extent(c(4.6,8,42.8,46.7)))

pdf("./LossWin_RDAprojection.pdf", width = 5, height = 6.5)
plot(nochange, col="grey", xlim=c(4.6,8), ylim=c(42.6,46.7), legend=F)
plot(lost, add = TRUE, col="#D55E00", xlim=c(4.6,8), ylim=c(42.6,46.7), legend=F)
lines(country, col="black",lty=1, xlim=c(4.6,8), ylim=c(42.6,46.7))
points(cities, col="black", pch=16, cex=0.7)
text(cities, cities$CITY_NAME, pos=3, cex=0.8)
plot(new, add = TRUE, col="forestgreen", xlim=c(4.6,8), ylim=c(42.6,46.7), legend=F)
dev.off()

#####################################################
#### Species response to climate change for RDA1 ####

## Difference present/futur RDA1
diff_rda1<-var_env[[1]]
diff_rda1[]<-as.vector(RDA1_proj)
diff_rda1_fu<-var_env[[1]]
diff_rda1_fu[]<-as.vector(RDA1_proj_fu)
diff_rda1[which(rda1_fu[]<quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[1] | rda1_fu[]>quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[2] | rda2_fu[]<quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[1] | rda2_fu[]>quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[2])] <- NA

diff_rda1<-rasterToPoints(crop(diff_rda1_fu-diff_rda1, extent(4.7, 8.199158, 42.6565, 46.85202)))

colors <- colorRampPalette(c("#FFEDA0", "red"))

p1 <- ggplot() + 
  geom_raster(aes(x=diff_rda1[,1],y=diff_rda1[,2],fill=cut(abs(diff_rda1[,3]), breaks=c(0,1,2,3,4)))) + 
  scale_fill_manual(values=colors(5)[c(1,2,5)]) +
  geom_polygon(aes(x = country.df$long, y = country.df$lat, group=country.df$group), color="black", alpha=0) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Change of genetic composition") +
  guides(fill=guide_legend(title="S")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())


## Distance between favourable pixels in the futur and present
coord_fu<-as.data.frame(rasterToPoints(rda1_fu))
coord_fu[,3]<-"futur"
coord_pres<-as.data.frame(rasterToPoints(rda1))
coord_pres[,3]<-"present"
table<-rbind(coord_pres, coord_fu)

nb_classes <- 4
dist_classes <- list()
coord_classes <- list()
for(i in 1:nb_classes)
{
  RDA1_pres_classes<-var_env[[1]]
  RDA1_pres_classes[]<-NA
  RDA1_pres_classes[which(rda1[]>quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i] & rda1[]<quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i+1])]<-1
  
  RDA1_fut_classes<-var_env[[1]]
  RDA1_fut_classes[]<-NA
  RDA1_fut_classes[which(rda1_fu[]>quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i] & rda1_fu[]<quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i+1])]<-1
  
  coord_fu_classes<-as.data.frame(rasterToPoints(RDA1_fut_classes))
  coord_fu_classes[,3]<-"futur"
  coord_pres_classes<-as.data.frame(rasterToPoints(RDA1_pres_classes))
  coord_pres_classes[,3]<-"present"
  table_classes<-rbind(coord_pres_classes, coord_fu_classes)
  
  dist_classes[[i]] <- 6371*nndist(table_classes[,1:2], by=as.factor(table_classes[,3]), k=1)[which(table_classes[,3]=="futur"),2]*pi/180
  coord_classes[[i]] <- coord_fu_classes[,c(1,2)]
}

table <- data.frame(Distance = as.vector(unlist(dist_classes)), Coordinates = do.call(rbind, coord_classes), Cell = cellFromXY(rda1, do.call(rbind, coord_classes)))

dDIST <- rda1
dDIST[] <- NA
dDIST[table$Cell]<-table$Distance
dDIST[which(rda1_fu[]<quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[1] | rda1_fu[]>quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[2] | rda2_fu[]<quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[1] | rda2_fu[]>quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[2])] <- NA
dDIST <- crop(dDIST, extent(c(4.7, 8.199158, 42.6565, 46.85202)))

mat_dDIST<-rasterToPoints(dDIST)
p2 <- ggplot() + 
  geom_raster(aes(x=mat_dDIST[,1],y=mat_dDIST[,2],fill= mat_dDIST[,3])) + 
  scale_fill_gradient(low = "#FFEDA0", high = "red") +
  geom_polygon(aes(x = country.df$long, y = country.df$lat, group=country.df$group), color="black", alpha=0) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Distance from present equivalent genetic composition") +
  guides(fill=guide_legend(title="D")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())


## Demo-genetic index for RDA1
toto <- as.matrix(rda1)
toto_fu <- as.matrix(rda1_fu)

invDist <- matrix(0, ncol=401, nrow=401)
for(i in 1:401) {
  for(j in 1:401) {
    invDist[i,j] <- dist(rbind(c(201,201), c(i,j)))+1
  }
}
dist <- invDist

library(foreach)
library(doParallel)

value <- matrix(NA, nrow(toto), ncol(toto))
combi = expand.grid(x = 1:ncol(value), y = 1:nrow(value))

registerDoParallel(cores = 12)
RES = foreach(x = combi$x, y = combi$y, .combine = "rbind") %dopar% {
  pixel = c(x,y)
  ind_x_min = max(x - 200, 1)
  ind_x_max = min(x + 200, ncol(value))
  ind_y_min = max(y - 200, 1)
  ind_y_max = min(y + 200, nrow(value))
  diff <- toto[ind_y_min:ind_y_max, ind_x_min:ind_x_max]-toto_fu[y,x]
  
  ind_x_min = ifelse(ind_x_min==1, abs(x - 201) + 1, 1)
  ind_x_max = ifelse(ind_x_max==ncol(value), 201 + (ncol(value)- x), 401)
  ind_y_min = ifelse(ind_y_min==1, abs(y - 201) + 1, 1)
  ind_y_max = ifelse(ind_y_max==nrow(value), 201  + (nrow(value) - y), 401)
  res <- sum(diff/dist[ind_y_min:ind_y_max, ind_x_min:ind_x_max],na.rm=T)/sum(!is.na(diff))
  return(data.frame(x=x,y=y,res=res))
}

MAT = matrix(abs(RES$res), nrow = nrow(toto), ncol = ncol(toto), byrow=T)

PROB <- rda1
PROB[] <- MAT
PROB[which(rda1_fu[]<quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[1] | rda1_fu[]>quantile(extract(rda1, occurences), probs= c(0.01, 0.99))[2] | rda2_fu[]<quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[1] | rda2_fu[]>quantile(extract(rda2, occurences), probs= c(0.01, 0.99))[2])] <- NA
PROB <- crop(PROB, extent(c(4.7, 8.199158, 42.6565, 46.85202)))


mat_PROB<-rasterToPoints(PROB)
p3 <- ggplot() + 
  geom_raster(aes(x=mat_PROB[,1],y=mat_PROB[,2],fill=mat_PROB[,3])) + 
  scale_fill_gradient(low = "#FFEDA0", high = "red", limits = c(0,0.06)) +
  geom_polygon(aes(x = country.df$long, y = country.df$lat, group=country.df$group), color="black", alpha=0) +
  geom_point(aes(x= cities@coords[,1], y=cities@coords[,2]), size=1) +
  geom_text(aes(x= cities@coords[,1]+0.05, y=cities@coords[,2]+0.05, label=cities$CITY_NAME), size=3.5) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="P")) +
  ggtitle("Difficulty of adaptation") +
  theme_bw() +
  theme(panel.grid.major = element_blank())

pdf("./ClimateResponse_Fagus.pdf", width = 12, height = 5)
ggarrange(p1, p2, p3, ncol=3, nrow=1, widths = c(1,1,1), labels=c("A", "B", "C"))
dev.off()
