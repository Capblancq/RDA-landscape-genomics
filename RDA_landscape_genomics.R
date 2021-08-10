#####################################################
##                                                 ##
##  Redundancy Analysis (RDA): a Swiss-army knife  ##
##             for landscape genomics              ##
##                                                 ##
##         Capblancq & Forester - raw code         ##
##                                                 ##
#####################################################

# Required libraries
library(pegas)
library(ggplot2)
library(raster)
library(rgdal)
library(ade4)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(qvalue)
library(robust)
library(WMDB)
library(ggVennDiagram)
library(cowplot)
library(corrplot)

# Set working directory
setwd("./")

##########################################################################################
###################################### DATA ##############################################

#### GENETIC MARKERS

# Loading genetic dataset
Genotypes <- read.table("./Data/Pine_AllNatural_GCandTotemIndivs_GWAS_SNPs_June8th2019.txt", header = T)
names_ind <- row.names(Genotypes)

# Function to transform the alleles into counts
genotypes2geno <- function(vec) {
  vec <- as.character(vec)
  lev <- levels(as.factor(vec))[-which(levels(as.factor(vec))=="00")]
  vec[vec==lev[1]] <- 0
  vec[vec==lev[2]] <- 1
  vec[vec==lev[3]] <- 2
  vec[vec=="00"] <- NA
  return(as.integer(vec))
}

# Running the function on the lodgepole pine dataset
Genotypes <- as.data.frame(apply(Genotypes, 2, genotypes2geno))
row.names(Genotypes) <- names_ind

# Loading samples metadata
InfoInd <- read.table("./Data/Pine_TotemField_AllNaturalIndivResiduals_Jan20th2017.csv", sep = ",", header = T)

# Sorting genetic data
Genotypes <- Genotypes[match(InfoInd$Internal_ID, row.names(Genotypes), nomatch = 0),]

# Estimating population allele frequencies
AllFreq <- aggregate(Genotypes, by = list(InfoInd$ProvSeedlotCode), function(x) mean(x, na.rm = T)/2)
row.names(AllFreq) <- as.character(AllFreq$Group.1)

# Filtering out the loci with a lot of missing data
na_pop <- apply(AllFreq[,-1], 2, function(x) sum(is.na(x)))
AllFreq <- AllFreq[,(which(na_pop<12)+1)]

# Imputation of missing population frequencies by the median across pop
for (i in 1:ncol(AllFreq))
{
  AllFreq[which(is.na(AllFreq[,i])),i] <- median(AllFreq[-which(is.na(AllFreq[,i])),i], na.rm=TRUE)
}

# Removing sites with very low allele frequency
freq_mean <- colMeans(AllFreq)
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]

# Ordering the loci based on their scaffold
AllFreq <- AllFreq[,order(colnames(AllFreq))]


#### CLIMATIC VARIABLES

# Coordinates of the populations
Coordinates <- read.table("./Data/PlSeedlots.csv", sep = ",", header = T, row.names = 1)
Coordinates <- Coordinates[match(row.names(AllFreq), Coordinates$id2, nomatch = 0),]
colnames(Coordinates) <- c("Population", "Latitude", "Longitude", "Elevation")

# Load ClimateNA bioclimatic variables
ras <- stack(list.files("./Data/ClimateNA/", pattern = ".img$", full.names = T))
names(ras) <- unlist(strsplit(unlist(lapply(strsplit(list.files("./Data/ClimateNA/", pattern = ".img$", full.names = T), split = "./Data/ClimateNA//"), function(x) x[2])), split = ".img"))

ras_6190 <- ras[[grep("6190", names(ras))]]
names(ras_6190) <- unlist(strsplit(names(ras_6190), split = "_6190"))
ras_2050 <- ras[[grep("2050_85", names(ras))]]
names(ras_2050) <- unlist(strsplit(names(ras_2050), split = "_2050_85"))
ras_2080 <- ras[[grep("2080_85", names(ras))]]
names(ras_2080) <- unlist(strsplit(names(ras_2080), split = "_2080_85"))

# Aligning NAs
remove.NAs.stack<-function(rast.stack){
  nom<-names(rast.stack)
  test1<-calc(rast.stack, fun=sum)
  test1[!is.na(test1)]<-1
  test2<-rast.stack*test1
  test2<-stack(test2)
  names(test2)<-nom
  return(test2)
}
ras_6190 <- remove.NAs.stack(ras_6190)
ras_2050 <- remove.NAs.stack(ras_2050)
ras_2080 <- remove.NAs.stack(ras_2080)

# Extract climatic values for each population in the present
Env <- data.frame(extract(ras_6190, Coordinates[,3:2]))

# Standardization of the variables
Env <- scale(Env)

# Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')

# Climatic table
Env <- as.data.frame(Env)
row.names(Env) <- c(Coordinates$Population)


#### PHENOTYPIC TRAITS

# Estimating mean trait value per population
traits <- aggregate(InfoInd[,19:24], by = list(InfoInd$ProvSeedlotCode), function(x) mean(x, na.rm = T))
traits <- as.data.frame(traits[match(Coordinates$Population, traits$Group.1, nomatch = 0),-1])
colnames(traits) <- c("Height","GthRate","ShootDryMass","Gth5pct","Gth95pct","ColdInjury")
  

#### POPULATION STRUCTURE

# Loading the intergenic SNPs dataset
Neutral <- read.table("./Data/Pine_AllNatural_GCandTotemIndivs_ControlSNPs_June8th2019.txt", header = T)
names_ind_neutral <- row.names(Neutral)

# Formatting genotypes
Neutral <- as.data.frame(apply(Neutral, 2, genotypes2geno))
row.names(Neutral) <- names_ind_neutral

# Sorting genetic data
Neutral <- Neutral[match(InfoInd$Internal_ID, row.names(Neutral), nomatch = 0),]

# Estimating allele frequencies for each source population
AllFreq_neutral <- aggregate(Neutral, by = list(InfoInd$ProvSeedlotCode), function(x) mean(x, na.rm = T)/2)

# Imputation of missing population frequencies by the median across the complete sampling
for (i in 2:ncol(AllFreq_neutral))
{
  AllFreq_neutral[which(is.na(AllFreq_neutral[,i])),i] <- median(AllFreq_neutral[-which(is.na(AllFreq_neutral[,i])),i], na.rm=TRUE)
}
AllFreq_neutral[,-1] <- AllFreq_neutral[,-1][,-which(is.na(colMeans(AllFreq_neutral[,-1])))]

# Running a PCA on neutral genetic markers
pca <- dudi.pca(AllFreq_neutral[,-1], scannf = FALSE, nf = 6)

# Neutral population structure table
PopStruct <- data.frame(Population = AllFreq_neutral[,1], scale(pca$li[,1:4]))
colnames(PopStruct) <- c("Population", "PC1", "PC2", "PC3", "PC4")


#### MASTER VARIABLES TABLE

# Table gathering all variables
Variables <- data.frame(Coordinates, PopStruct[,-1], Env, traits)


#### MAPPING FEATURES

# Administrative boundaries
admin <- ne_countries(scale = "medium", returnclass = "sf")

# Species range
range <- readOGR("./Data/pinucon/pinucont.shp")
crs(range) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

##########################################################################################


##########################################################################################
################################ VARIABLE SELECTION ######################################

# Null model
RDA0 <- rda(AllFreq ~ 1,  Variables) 

# Full model
RDAfull <- rda(AllFreq ~ AHM + bFFP + CMD + DD_0 + DD_18 + DD18 + DD5 + eFFP + EMT + Eref + EXT + FFP + MAP + MAR + MAT + MCMT + MSP + MWMT + NFFD + PAS + PPT_sm + PPT_wt + RH + SHM + Tave_sm + Tave_wt + TD, Variables)

# Stepwise procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T)
mod$anova

##########################################################################################


##########################################################################################
###################################### PARTIAL RDA #######################################

# Full model
pRDAfull <- rda(AllFreq ~ PC1 + PC2 + PC3 + Longitude + Latitude + MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS,  Variables)
RsquareAdj(pRDAfull)
anova(pRDAfull)

# Pure climate model
pRDAclim <- rda(AllFreq ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + Condition(Longitude + Latitude + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAclim)
anova(pRDAclim)

# Pure neutral population structure model 
pRDAstruct <- rda(AllFreq ~ PC1 + PC2 + PC3 + Condition(Longitude + Latitude + MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS),  Variables)
RsquareAdj(pRDAstruct)
anova(pRDAstruct)

# Pure geography model
pRDAgeog <- rda(AllFreq ~ Longitude + Latitude + Condition(MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAgeog)
anova(pRDAgeog)

# Correlogram 
M <- cor(Variables[, c("PC1","PC2","PC3","Longitude","Latitude","MAR","EMT","MWMT","CMD","Tave_wt","DD_18","MAP","Eref","PAS")])

# Figure correlogram
pdf("Figures-Tables/Correlogram.pdf")
corrplot(M, type="upper")
dev.off()

##########################################################################################


##########################################################################################
################################## RDA GENOME SCANS ######################################

#### RDA-based GEA

# Run RDA on population allele frequencies while accounting for population structure using three genetic PCs
RDA_env <- rda(AllFreq ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + Condition(PC1 + PC2 + PC3),  Variables)

# Screeplot
screeplot(RDA_env, main="Eigenvalues of constrained axes")

# Function rdadapt
source("./src/rdadapt.R")

# Running the function with K = 2
rdadapt_env<-rdadapt(RDA_env, 2)

# P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rdadapt_env$p.values)

# Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))

# Top hit outlier per contig
outliers <- outliers[order(outliers$contig, outliers$p.value),]

# List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

#### Importance of the axes

pdf("Figures-Tables/Inertia_genome_scans.pdf", width = 1, height = 1)
ggplot() +
  geom_col(aes(x=1:length(RDA_env$CCA$eig), y=as.vector(RDA_env$CCA$eig), fill = c("1","1", rep("0", length(RDA_env$CCA$eig)-2)))) +
  scale_fill_manual(values = c("grey", "black")) +
  theme_bw() +
  theme(legend.position = "none", panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
dev.off()

### RDA genetic space

# Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

# Biplot of RDA loci and variables scores
p1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (29.3%)") + ylab("RDA 2 (11.5%)") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

### Manhattan plot

Outliers <- rep("Neutral", length(colnames(AllFreq)))
Outliers[colnames(AllFreq)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(AllFreq)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(AllFreq)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
p2 <- ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  facet_wrap(~"Manhattan plot", nrow = 3) +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

### Figure

pdf("Figures-Tables/Genome_scans.pdf", width = 12, height = 3.8)
ggarrange(p1, p2, widths = c(1.1,2), common.legend = T, legend = "right", labels = c("A", "B"))
dev.off()

### Comparing outlier list to Mahony et al. 2020

# Loading previously identified outliers
outliers_mahony <- read.table("./Data/gea_outliers_by_env_Mahonyetal2020.txt", header = T, row.names = 1)
outliers_mahony_env <- row.names(outliers_mahony)[apply(outliers_mahony[,3:21], 1, function(x) any(as.logical(x)))]
outliers_mahony_gpa <- row.names(outliers_mahony)[apply(outliers_mahony[,22:25], 1, function(x) any(as.logical(x)))]
list_outliers <- list(RDA_best = outliers_rdadapt_env, RDA_all = as.character(outliers$Loci), Bayenv2 = outliers_mahony_env, GPA = outliers_mahony_gpa)

# Venn diagram
pdf("Figures-Tables/Venn_diagram_RDA_Bayenv2.pdf", width = 6, height = 5)
ggVennDiagram(list_outliers, category.names = c("RDA best hit", "RDA all", "Bayenv2", "Bayenv2 & GPA"), lty="solid", color="black", size=0.2) + 
  scale_fill_gradient2(low = "gray90", high = 'gray20') + guides(fill = FALSE) + theme(text = element_text(size=16, family = "Times"))
dev.off()

#### Not accounting for population structure

# Running a simple RDA model
RDA_env_unconstrained <- rda(AllFreq ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS,  Variables)

# Running the rdadapt function
rdadapt_env_unconstrained <- rdadapt(RDA_env_unconstrained, 2)

# Setting the p-value threshold
thres_env <- 0.01/length(rdadapt_env_unconstrained$p.values)

# Identifying the outliers for the simple RDA
outliers_unconstrained <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env_unconstrained$p.values<thres_env)], p.value = rdadapt_env_unconstrained$p.values[which(rdadapt_env_unconstrained$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env_unconstrained$p.values<thres_env)], split = "_"), function(x) x[1])))
outliers_unconstrained <- outliers_unconstrained[order(outliers_unconstrained$contig, outliers_unconstrained$p.value),]
outliers_rdadapt_env_unconstrained <- as.character(outliers_unconstrained$Loci[!duplicated(outliers_unconstrained$contig)])

# For all the outliers
list_outliers_RDA_all <- list(RDA_constrained = as.character(outliers$Loci), RDA_unconstrained = as.character(outliers_unconstrained$Loci))
p1 <- ggVennDiagram(list_outliers_RDA_all, category.names = c("RDA constrained all", "RDA unconstrained all"), lty="solid", color="black", size=0.2) + 
  scale_fill_gradient2(low = "gray90", high = 'gray20') + guides(fill = FALSE) + theme(text = element_text(size=16, family = "Times"))

# Only for the top hit locus per contig
list_outliers_RDA_top <- list(RDA_constrained = outliers_rdadapt_env, RDA_unconstrained = outliers_rdadapt_env_unconstrained)
p2 <- ggVennDiagram(list_outliers_RDA_top, category.names = c("RDA constrained top hits", "RDA unconstrained top hits"), lty="solid", color="black", size=0.2) + 
  scale_fill_gradient2(low = "gray90", high = 'gray20') + guides(fill = FALSE) + theme(text = element_text(size=16, family = "Times"))

pdf("Figures-Tables/Venn_diagram_RDA_constrained_Vs_unconstrained.pdf", width = 12, height = 5)
ggarrange(p1, p2)
dev.off()

# Common outliers
common_outliers_RDA_top <- Reduce(intersect, list_outliers_RDA_top)

##########################################################################################


##########################################################################################
############################## PROJECTION ADAPTIVE INDEXES ###############################

# Adaptively enriched RDA
RDA_outliers <- rda(AllFreq[,common_outliers_RDA_top] ~ MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS,  Variables)

# RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
p1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1 (67%)") + ylab("RDA 2 (21%)") +
  facet_wrap(~"Adaptively enriched genetic space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

# Function to predict the adaptive index across the landscape
source("./src/adaptive_index.R")

# Running the function for all the climatic pixels of lodgepole pine distribution range
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras_6190, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

# Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}

# Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))

p2 <- ggplot(data = TAB_RDA) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ variable) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

# Figure
pdf("Figures-Tables/Genetic_turnover_predict.pdf", width = 12, height = 3.5)
cowplot::plot_grid(p1, p2, nrow = 1, align = "h", axis = "bt", rel_widths = c(0.39,1), labels = c("A", "B"))
dev.off()

##########################################################################################


##########################################################################################
################################ GENOMIC OFFSET ##########################################

# Function to predict genomic offset from a RDA model
source("./src/genomic_offset.R")

# Running the function for 2020 and 2080
res_RDA_proj2080 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_6190, env_fut = ras_2080, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
res_RDA_proj2050 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_6190, env_fut = ras_2050, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

# Table global genetic offset predicted for 2050 and 2080
RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj2050$Proj_offset_global), rasterToPoints(res_RDA_proj2080$Proj_offset_global)), Date = c(rep("2050", nrow(rasterToPoints(res_RDA_proj2050$Proj_offset_global))), rep("2080", nrow(rasterToPoints(res_RDA_proj2080$Proj_offset_global)))))

# Projecting genomic offset on a map
colors <- c(colorRampPalette(brewer.pal(11, "Spectral")[6:5])(2), colorRampPalette(brewer.pal(11, "Spectral")[4:3])(2), colorRampPalette(brewer.pal(11, "Spectral")[2:1])(3))
p_offset <- ggplot(data = RDA_proj_offset) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = cut(Global_offset, breaks=seq(1, 8, by = 1), include.lowest = T)), alpha = 1) + 
  scale_fill_manual(values = colors, labels = c("1-2","2-3","3-4","4-5","5-6","6-7","7-8"), guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
  geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  facet_grid(~ Date) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

# Figure
pdf("./Figures-Tables/Genomic_offset_loadings.pdf", width = 9, height = 3.5)
p_offset
dev.off()

##########################################################################################


##########################################################################################
################################ RDA OFFSET VALIDATION ###################################

# Environmental variables for the period the seedlings were growing in the common garden 2012-2015
env_garden <- read.table("./Data/ClimateData_VancouverGarden.csv", sep = ",", header = T)
env_garden <- env_garden[, row.names(RDA_outliers$CCA$biplot)]
env_garden$MAR <- extract(ras_6190[[row.names(RDA_outliers$CCA$biplot)]]$MAR, c(-123.250, 49.256))[2]
env_garden <- as.data.frame(scale(env_garden, center = center_env[row.names(RDA_outliers$CCA$biplot)], scale = scale_env[row.names(RDA_outliers$CCA$biplot)]))
env_garden <- colMeans(env_garden)

# Environmental variables for the source populations between 61 and 90
env_provenance <- data.frame(extract(ras_6190[[row.names(RDA_outliers$CCA$biplot)]], Variables[,c("Longitude", "Latitude")]))
env_provenance <- as.data.frame(scale(env_provenance, center_env[row.names(RDA_outliers$CCA$biplot)], scale_env[row.names(RDA_outliers$CCA$biplot)]))
row.names(env_provenance) <- Variables[,"Population"]

# Function to estimate genetic offset between provenance and garden  
source("./src/provgar_offset.R")

# Running the provenance to garden genetic offset function
provgar_offset <- provgar_offset(RDA = RDA_outliers, K = 2, env_garden = env_garden, env_provenance = env_provenance, weights = TRUE)

# Climate distance between provenances and garden
tabmaha <- as.data.frame(rbind(env_garden, env_provenance))
mahaclim <- mahalanobis(as.matrix(tabmaha), center = as.numeric(tabmaha[1,]), cov = cov(as.matrix(tabmaha)))

# Comparison climate transfer vs. genomic offset impact on fitness
TAB_comp <- data.frame(x = rep(c(mahaclim[-1], provgar_offset[-1]), 2), y = c(rep(traits$Height, 2), rep(traits$GthRate, 2)), method = rep(rep(c("Climate", "Genetics"), each = nrow(traits)), 2), variable = rep(c("Final height", "Growth rate"), each = 2*nrow(traits)))
p_garden <- ggplot(TAB_comp, aes(x=x, y=y)) +
  geom_point(shape=19, size = 2, alpha = .7) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_grid(variable ~ method, scales = "free") +
  ylab("Scaled fitness trait value") + 
  xlab("Climate transfer distance                              Genomic offset") +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(panel.grid = element_blank(), panel.background = element_blank())

# Figure
pdf("Figures-Tables/Genetic_offset_gardens.pdf", width = 7, height = 4.5)
p_garden
dev.off()

# Quadratic regression statistical testing
summary(lm(traits$Height ~ poly(mahaclim[-1])))
summary(lm(traits$Height ~ poly(provgar_offset[-1])))
summary(lm(traits$GthRate ~ poly(mahaclim[-1])))
summary(lm(traits$GthRate ~ poly(provgar_offset[-1])))

##########################################################################################


#########################################################################################
################################ POPULATION ADAPTIVE UNIT ###############################

library(factoextra)
library(plyr)
library(geometry)

# RDA1 and RDA2 scores for all the pixels of the range
TAB <- data.frame(x = rasterToPoints(res_RDA_proj_current$RDA1)[,1], y = rasterToPoints(res_RDA_proj_current$RDA1)[,2], RDA1 = rasterToPoints(res_RDA_proj_current$RDA1)[,3], RDA2 = rasterToPoints(res_RDA_proj_current$RDA2)[,3])

# Hierarchical K-means clustering for 10,000 random pixels
TAB_sample <- TAB[sample(1:nrow(TAB), 10000),3:4]
clust <- list()
for(i in 1:15){
  clust[[i]] <- hkmeans(x = TAB_sample, k = i, hc.metric = "euclidean", hc.method = "ward.D2", iter.max = 10, km.algorithm = "Hartigan-Wong")
}

# Select the best number of clusters
within_cust_var <- lapply(clust, function(x) x$withinss/x$totss)
plot(unlist(lapply(within_cust_var, mean)))

# Extrapolate the clustering to all the pixels with 4 clusters
df <- data.frame(x = TAB_sample$RDA1, y = TAB_sample$RDA2, cluster = clust[[4]]$cluster) 
find_hull <- function(df) df[chull(df$x, df$y),]
hulls <- ddply(df, "cluster", find_hull) # Convex hulls.
cluster <- rep(NA, nrow(TAB))
for(i in 1:4){
  inhull <- inhulln(convhulln(hulls[hulls$cluster==i,1:2]), as.matrix(TAB[,3:4])) # Intersection convex hulls and points
  cluster[inhull] <- i
}
clusterNA <- cluster[which(is.na(cluster))]
for(i in 1:length(clusterNA)){
  clusterNA[i] <- which.min(pointDistance(TAB[which(is.na(cluster))[i], 3:4], clust[[4]]$centers, lonlat = F)) # Fill the NA using the distance to the hull centers
}
cluster[which(is.na(cluster))] <- clusterNA

# Biplot 
TAB$color <- as.factor(cluster)
TAB_var <- as.data.frame(RDA_outliers$CCA$biplot[,1:2])
p_unit_biplot <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB, aes(x=RDA1/20, y=RDA2/20, colour = color), size = 1) +
  scale_colour_viridis_d(alpha = 0.8, direction = -1, option = "A", begin = 0.15) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=FALSE) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

# Map 
p_unit_map <- ggplot(data = TAB) + 
  geom_sf(data = admin, fill=gray(.9), size=0) +
  geom_raster(aes(x = x, y = y, fill = color), ) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", begin = 0.15) +
  geom_sf(data = admin, fill=NA, size=0.1) +
  coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive Units")) +
  facet_grid(~ "Pinus contorta range") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

# Figure
png("Figures-Tables/Adaptive_population_units.png", width = 8, height = 3.3, units = "in", res = 300)
cowplot::plot_grid(p_unit_biplot, p_unit_map, nrow = 1, align = "h", rel_widths = c(1,1.8), labels = c("A", "B"), axis = "bt")
dev.off()

#########################################################################################


