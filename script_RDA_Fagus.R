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
GenData<-GenData[-c(grep(pattern = "Bis", row.names(GenData)), grep(pattern = "Ter", row.names(GenData)), grep(pattern = "SB", row.names(GenData))),]

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
                              
##############
                              
                              
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

#####################
                              
                              
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
                   
#################                   
                   

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
                   
#########################
                   
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
                     
############################
            
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

##############################################################################################################################                     
                     
#############################################################################
#### Genetic Offset based on RDA scores in current and future conditions ####

colors <- colorRampPalette(c("#FFEDA0", "red"))

## Offset present/future for RDA1
proj_rda1<-var_env[[1]]
proj_rda1[]<-as.vector(RDA1_proj)
proj_rda1 <- mask(proj_rda1, range)
proj_rda1_fu<-var_env[[1]]
proj_rda1_fu[]<-as.vector(RDA1_proj_fu)
proj_rda1_fu <- mask(proj_rda1_fu, range)
diff_rda1<-rasterToPoints(crop(proj_rda1_fu-proj_rda1, extent(4.7, 8.199158, 42.6565, 46.85202)))

p1 <- ggplot() + 
  geom_polygon(data=country, aes(long, lat, group = group), colour = NA, fill = "gray95", size = 0.3) +
  geom_raster(aes(x=diff_rda1[,1],y=diff_rda1[,2],fill=cut(abs(diff_rda1[,3]), breaks=c(0,1,2,3,4)))) + 
  scale_fill_manual(values=colors(4)) +
  geom_polygon(data=country, aes(long, lat, group = group), colour = "gray25", fill = NA, size = 0.3) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Genetic offset")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

## Offset present/future for RDA2
proj_rda2<-var_env[[1]]
proj_rda2[]<-as.vector(RDA2_proj)
proj_rda2 <- mask(proj_rda2, range)
proj_rda2_fu<-var_env[[1]]
proj_rda2_fu[]<-as.vector(RDA2_proj_fu)
proj_rda2_fu <- mask(proj_rda2_fu, range)
diff_rda2<-rasterToPoints(crop(proj_rda2_fu-proj_rda2, extent(4.7, 8.199158, 42.6565, 46.85202)))

p2 <- ggplot() + 
  geom_polygon(data=country, aes(long, lat, group = group), colour = NA, fill = "gray95", size = 0.3) +
  geom_raster(aes(x=diff_rda2[,1],y=diff_rda2[,2],fill=cut(abs(diff_rda2[,3]), breaks=c(0,1,2,3,4)))) + 
  scale_fill_manual(values=colors(4)) +
  geom_polygon(data=country, aes(long, lat, group = group), colour = "gray25", fill = NA, size = 0.3) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Genetic offset")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

## Offset present/future for RDA1 + RDA2
diff_tot<-rasterToPoints(crop(proj_rda1_fu-proj_rda1, extent(4.7, 8.199158, 42.6565, 46.85202))-crop(proj_rda2_fu-proj_rda2, extent(4.7, 8.199158, 42.6565, 46.85202)))

p_offset <- ggplot() + 
  geom_polygon(data=country, aes(long, lat, group = group), colour = NA, fill = "gray95", size = 0.3) +
  geom_raster(aes(x=diff_tot[,1],y=diff_tot[,2],fill=cut(abs(diff_tot[,3]), breaks=c(1,2,3,4,5)))) + 
  scale_fill_manual(values=c(colors(4)[-4],"brown")) +
  geom_polygon(data=country, aes(long, lat, group = group), colour = "gray25", fill = NA, size = 0.3) +
  geom_point(data=EnvData, aes(x= Y, y= X), cex=0.7) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Genetic offset")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

#############################################################################

################################################################################
#### Geographic offset based on RDA scores in current and future conditions ####

coord_fu1<-as.data.frame(rasterToPoints(rda2_fu))
coord_fu1[,3]<-"futur"
coord_pres1<-as.data.frame(rasterToPoints(rda2))
coord_pres1[,3]<-"present"
table1<-rbind(coord_pres1, coord_fu1)
coord_fu2<-as.data.frame(rasterToPoints(rda2_fu))
coord_fu2[,3]<-"futur"
coord_pres2<-as.data.frame(rasterToPoints(rda2))
coord_pres2[,3]<-"present"
table2<-rbind(coord_pres2, coord_fu2)

nb_classes <- 4
dist_classes_rda1 <- list()
dist_classes_rda2 <- list()
coord_classes_rda1 <- list()
coord_classes_rda2 <- list()
for(i in 1:nb_classes)
{
  RDA1_pres_classes<-var_env[[1]]
  RDA1_pres_classes[]<-NA
  RDA1_pres_classes[which(rda1[]>quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i] & rda1[]<quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i+1])]<-1
  RDA2_pres_classes<-var_env[[1]]
  RDA2_pres_classes[]<-NA
  RDA2_pres_classes[which(rda2[]>quantile(rda2[], probs= seq(0,1,1/nb_classes), na.rm=T)[i] & rda2[]<quantile(rda2[], probs= seq(0,1,1/nb_classes), na.rm=T)[i+1])]<-1
  
  RDA1_fut_classes<-var_env[[1]]
  RDA1_fut_classes[]<-NA
  RDA1_fut_classes[which(rda1_fu[]>quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i] & rda1_fu[]<quantile(rda1[], probs= seq(0,1,1/nb_classes), na.rm=T)[i+1])]<-1
  RDA2_fut_classes<-var_env[[1]]
  RDA2_fut_classes[]<-NA
  RDA2_fut_classes[which(rda2_fu[]>quantile(rda2[], probs= seq(0,1,1/nb_classes), na.rm=T)[i] & rda2_fu[]<quantile(rda2[], probs= seq(0,1,1/nb_classes), na.rm=T)[i+1])]<-1
  
  coord_fu_classes_rda1<-as.data.frame(rasterToPoints(RDA1_fut_classes))
  coord_fu_classes_rda1[,3]<-"futur"
  coord_pres_classes_rda1<-as.data.frame(rasterToPoints(RDA1_pres_classes))
  coord_pres_classes_rda1[,3]<-"present"
  table_classes_rda1<-rbind(coord_pres_classes_rda1, coord_fu_classes_rda1)  
  coord_fu_classes_rda2<-as.data.frame(rasterToPoints(RDA2_fut_classes))
  coord_fu_classes_rda2[,3]<-"futur"
  coord_pres_classes_rda2<-as.data.frame(rasterToPoints(RDA2_pres_classes))
  coord_pres_classes_rda2[,3]<-"present"
  table_classes_rda2<-rbind(coord_pres_classes_rda2, coord_fu_classes_rda2)
  
  dist_classes_rda1[[i]] <- 6371*nndist(table_classes_rda1[,1:2], by=as.factor(table_classes_rda1[,3]), k=1)[which(table_classes_rda1[,3]=="futur"),2]*pi/180
  coord_classes_rda1[[i]] <- coord_fu_classes_rda1[,c(1,2)]
  dist_classes_rda2[[i]] <- 6371*nndist(table_classes_rda2[,1:2], by=as.factor(table_classes_rda2[,3]), k=1)[which(table_classes_rda2[,3]=="futur"),2]*pi/180
  coord_classes_rda2[[i]] <- coord_fu_classes_rda2[,c(1,2)]
}

table_rda1 <- data.frame(Distance = as.vector(unlist(dist_classes_rda1)), Coordinates = do.call(rbind, coord_classes_rda1), Cell = cellFromXY(rda1, do.call(rbind, coord_classes_rda1)))
table_rda2 <- data.frame(Distance = as.vector(unlist(dist_classes_rda2)), Coordinates = do.call(rbind, coord_classes_rda2), Cell = cellFromXY(rda2, do.call(rbind, coord_classes_rda2)))

dDIST_rda1 <- rda1
dDIST_rda1[] <- NA
dDIST_rda1[table_rda1$Cell]<-table_rda1$Distance
dDIST_rda2 <- rda2
dDIST_rda2[] <- NA
dDIST_rda2[table_rda2$Cell]<-table_rda2$Distance

dDIST <- mean(dDIST_rda1, dDIST_rda2, na.rm = T)

mat_dDIST<-rasterToPoints(dDIST)
p_dist <- ggplot() + 
  geom_polygon(data=country, aes(long, lat, group = group), colour = NA, fill = "gray95", size = 0.3) +
  geom_raster(aes(x=mat_dDIST[,1],y=mat_dDIST[,2],fill= cut(mat_dDIST[,3], breaks=c(-1,5,10,15,20,25,30,35,40,45,50)))) + 
  scale_fill_manual(values=colors(10)) +
  geom_polygon(data=country, aes(long, lat, group = group), colour = "gray25", fill = NA, size = 0.3) +
  geom_point(data=EnvData, aes(x= Y, y= X), cex=0.7) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Distance \n from similar \n current conditions")) +
  theme_bw() +
  theme(panel.grid.major = element_blank())

###########################################################################
                    
################################################
#### PAI and SGV in the sampled populations ####

mean_q <- apply(1-FreqData[, outliers_rdadapt_pop], 2, mean)

# Population adaptive index (Bonin et al. 2007) 
PAI <- lapply(1:nrow(FreqData), function(x) sum(abs((1-FreqData[x,outliers_rdadapt_pop])-mean_q), na.rm=T))
names(PAI) <- row.names(FreqData)

SVG <- lapply(lapply(1:nrow(FreqData), function(x) FreqData[x,outliers_rdadapt_pop]*(1-FreqData[x,outliers_rdadapt_pop])), function(y) mean(unlist(y)))
names(SVG) <- row.names(FreqData)

TAB <- data.frame(Name = EnvData$Group.1, Latitude = EnvData$X, Longitude = EnvData$Y, PAI = unlist(PAI), SVG = unlist(SVG))

p_SGV <- ggplot() +
  geom_polygon(data=country, aes(long, lat, group = group), colour = NA, fill = "gray95", size = 0.3) +
  geom_raster(aes(x=mat_dDIST[,1],y=mat_dDIST[,2],fill="gray75"), fill="gray75") +
  #geom_polygon(data=range, aes(long, lat, group = group), colour = NA, fill = "gray75", size = 0.3) +
  geom_polygon(data=country, aes(long, lat, group = group), colour = "gray35", fill = NA, size = 0.3) +
  geom_point(data=TAB, aes(x= Longitude, y= Latitude, size = SVG, colour = PAI)) + 
  scale_color_gradient(low="#FFEDA0", high = "red") +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_blank())

g1 <- ggplotGrob(p_offset)
g2 <- ggplotGrob(p_dist)
g3 <- ggplotGrob(p_SGV)
g = cbind(g1, g2, g3, size = "first")

pdf("./GeneticOffset_Distance_SGV_Fagus.pdf", width = 12, height = 4)
grid.newpage()
grid.draw(g)
dev.off()

################################################
                 
