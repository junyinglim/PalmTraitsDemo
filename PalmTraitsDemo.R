## R script to make Fig. 3 of Kissling et al.: 
## PalmTraits 1.0, a species-level functional trait database for palms worldwide
## 
## Author of script: Jun Ying Lim
## Modifications and additional annotations: W. Daniel Kissling

# CLEAN WORKSPACE ====================
rm(list=ls())

## DIRECTORIES ====================
main.dir <- "~/Dropbox/Projects/2019/palms/projects/palmTraits/repo/PalmTraitsDemo/"
fig.dir <- file.path(main.dir, "figs")
data.dir <- file.path(main.dir, "data")

## LOAD PACKAGES ====================
# You will need to install ggtree which is crucial for the plotting in this script, run the following code to ensure that all dependent packages are also installed and updated
# install.packages(c("devtools", "rgdal", "rgeos", "ape", "plyr", "wesanderson", "cowplot", "scatterpie", "ggrepel", "ggplot2", "BiocManager", "viridis", "reshape2"), dependencies = TRUE)
# devtools::install_github('GuangchuangYu/ggtree', force = TRUE)

# spatial packages for maps
library(rgdal); library(rgeos)

# tidy packages for data handling and summary statistics
library(plyr); library(reshape2)

# phylogenetic package
library(ape)

# general plotting packages
library(ggplot2);library(wesanderson); library(cowplot); library(ggrepel); library(scatterpie); library(ggtree); library(viridis)

## IMPORT DATA ====================
# Load spatial data
#     these polygons represent TDWG level 3 units ('botanical countries') as available from the
#     World Geographical Scheme for Recording Plant Distributions (WGSRPD) 
#     http://www.tdwg.org/standards/109
#     https://github.com/tdwg/wgsrpd 
shape <- readOGR(dsn = file.path(main.dir, "tdwg_level3_shp", "."), layer = "level3")
shape@data$id <- rownames(shape@data)

# "Simplify" polygons for easier plotting
shape2 <- gSimplify(spgeom = shape, tol = 0.01, topologyPreserve = TRUE) 
shape <- SpatialPolygonsDataFrame(shape2, shape@data)
shape_centroid <- gCentroid(shape, byid = TRUE)
shape@data$centroid_long <- shape_centroid$x
shape@data$centroid_lat <- shape_centroid$y
shape_fort <- fortify(shape, id = id)
shape_gg <- merge(shape_fort, shape@data, by = "id")

# Load palm trait data
traitData <- read.table(file.path(data.dir, "PalmTraits_1.0.txt"), stringsAsFactors = FALSE, sep= "\t", header = TRUE)

# Load palm occurence data
#     These are palm species occurrences within TDWG level 3 units ('botanical countries')
#     as provided by the world checklist of palms, available from the Royal Botanic Gardens Kew
#     http://wcsp.science.kew.org/
#     Govaerts, R. & Dransfield, J. (2005). World checklist of palms. Royal Botanic Gardens Kew
#     Richmond.
#     Downloaded on July 2015
occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"), stringsAsFactors = FALSE)
occ$SpecName <- gsub(occ$SpecName, pattern = "_", replacement = " ")

# Load phylogeny
#     This is a Maximum Clade Credibility (MCC) phylogenetic tree as used by Onstein et al. (2017)
#     Nature Ecology & Evolution 1: 1903-1911.
#     and by Onstein et al. (2018), Proceedings of the Royal Society B: Biological Sciences, 285,
#     20180882.
#     Available from DRYAD: https://datadryad.org/resource/doi:10.5061/dryad.cm4nm
palmPhylo <- read.nexus(file.path(data.dir, "TREE.nex"))
palmPhylo$tip.label <- gsub(palmPhylo$tip.label, pattern = "_", replacement = " ")

# Identify species that differ between the phylogeny and trait dataset
setdiff(palmPhylo$tip.label, traitData$SpecName)
setdiff(traitData$SpecName, palmPhylo$tip.label)

# Exclude those taxa from the phylogeny and the trait dataset
intersectTaxa <- intersect(traitData$SpecName, palmPhylo$tip.label)
traitDataSubset <- subset(traitData, SpecName %in% intersectTaxa)
palmPhyloSubset <- drop.tip(palmPhylo,
                            tip = palmPhylo$tip.label[!palmPhylo$tip.label %in% intersectTaxa])

################################################################
## FIGURE 2a,b from Kissling et al.: PLOTTING COMPLETENESS ONTO WORLD MAP
# Merge target traits to occurrence dataset
nm <- c("Climbing", "StemArmed", "MaxStemHeight_m", "Max_Blade_Length_m", "AverageFruitLength_cm")
occ_traitcompl <- merge(occ, traitData[c("SpecName",nm)], by = "SpecName" )

traitList <- traitListToPlot <- c("Climbing", "StemArmed", "MaxStemHeight_m", "Max_Blade_Length_m", "AverageFruitLength_cm")

# Calculate completeness and trait means for each botanical country
meanTrait <- function(x, traitList){
  ntaxa <- nrow(x)
  res <- list()
  res["ntaxa"] <- ntaxa
  
  for(i in 1:length(traitList)){
    # calculate proportion of empty values for each area for each trait
    res[paste0(traitList[i], "_completeness")] <- round(1 - ( sum(is.na(x[[traitList[i]]])) / ntaxa), 3)
    if(traitList[i] == "StemArmed"){
      res[paste0(traitList[i], "_mean")] <- round(sum(ifelse(x[[traitList[i]]]>=1,1,0), na.rm = T),3)
    } else {
      res[paste0(traitList[i], "_mean")] <- round(mean(x[[traitList[i]]], na.rm = TRUE),3)  
    }
  }
  data.frame(res)
}

traitByArea <- ddply(.data = occ_traitcompl, .variables = .(Area_code_L3), .fun = meanTrait, traitList = traitList)

# Merge with shapefile
shape@data <- merge(traitByArea, shape@data, by.x ="Area_code_L3", by.y = "LEVEL3_COD", all.y = TRUE)
shapeDF <- fortify(shape) 
finalDataSet <- merge(shapeDF, shape@data, all.x = TRUE)

finalDataSet_long <- melt(finalDataSet, id.vars = c("id", "long", "lat", "group", "order", "hole", "piece", "Area_code_L3", "LEVEL3_NAM"))
finalDataSet_long$value <- as.numeric(finalDataSet_long$value)

traitListLabels <- c("Growth form \n(climbing/non-climbing)","Stem armature", "Maximum stem height", "Maximum blade length", "Average fruit length")

# Create list of plot of trait completeness for each trait
completenessPlotList <- list()
for(i in 1:length(traitList)){
  completenessPlotList[[i]] <- 
    ggplot(data = subset(finalDataSet_long, variable == paste0(traitList[i], "_completeness") & !Area_code_L3 == "ANT")) + 
    geom_polygon(aes(y = lat, x = long, group = group, fill = value)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Longitude", y = "Latitude", title = traitListLabels[i]) + 
    theme(legend.position = "right",
          plot.title = element_text(size = 20, color = "grey20", hjust = 0),
          axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(),
          panel.background = element_blank()) +
    scale_fill_viridis(limits = c(0,1), name = "Proportion\ncomplete", direction = 1)
}

# Create list of plot of trait variation for each trait
traitMeanPlotList <- list()
legendKey <- c("Prop. of taxa", "Species\nrichness", "Mean (m)", "Mean (m)", "Mean (cm)")
q_probs = seq(0, 1.0, 0.1)

for(i in 1:length(traitList)){
  temp <- subset(finalDataSet_long, variable == paste0(traitList[i], "_mean") & !Area_code_L3 == "ANT")
  if(i == 1 | i == 2){
    temp$value2 <- temp$value
    discrete = FALSE
  } else {
    temp$value2 <- cut(temp$value, quantile(temp$value, probs = q_probs, na.rm = T), include.lowest = T)
    levels(temp$value2) <- gsub(levels(temp$value2), pattern = "\\(|\\]|\\[", replacement = "")
    levels(temp$value2) <- gsub(levels(temp$value2), pattern = ",", replacement = " - ")  
    discrete = TRUE
    breaks = levels(temp$value2)
  }
  
  traitMeanPlotList[[i]] <- 
    ggplot(data = temp) +
    geom_polygon(aes(y = lat, x = long, group = group, fill = value2)) +
    labs(x = "Longitude", y = "Latitude", title = traitListLabels[i], color = "white") + 
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = "right",
          plot.title = element_text(size = 20, color = "white", hjust = 0),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(),
          panel.background = element_blank())
  
  if(i %in% c(3,4,5)){
    traitMeanPlotList[[i]] <- traitMeanPlotList[[i]] + scale_fill_viridis(name = legendKey[i], discrete = discrete, na.value="grey", breaks = breaks, direction = 1)
  } else {
    traitMeanPlotList[[i]] <- traitMeanPlotList[[i]] + scale_fill_viridis(name = legendKey[i], discrete = discrete, na.value="grey")
  }
}


completenessPlotGrid <- plot_grid(completenessPlotList[[1]],
                                  completenessPlotList[[2]],
                                  completenessPlotList[[3]],
                                  completenessPlotList[[4]],
                                  completenessPlotList[[5]],
                                  align = "v", nrow = 5, ncol = 1)

traitMeanPlotGrid  <- plot_grid(traitMeanPlotList[[1]],
                                traitMeanPlotList[[2]],
                                traitMeanPlotList[[3]],
                                traitMeanPlotList[[4]],
                                traitMeanPlotList[[5]],
                                align = "v", nrow = 5, ncol = 1)

completenessPlotGridTitle <- ggdraw() + draw_label("Geographic completeness", size = 30)
traitMeanPlotGridTitle <- ggdraw() + draw_label("Trait variation", size = 30)

trait_compl_plot <- plot_grid(completenessPlotGridTitle,
                              completenessPlotGrid, ncol = 1, rel_heights = c(0.05, 0.95))
trait_var_plot <- plot_grid(traitMeanPlotGridTitle,
                            traitMeanPlotGrid, ncol = 1, rel_heights = c(0.05, 0.95))

fig2_combined <- plot_grid(trait_compl_plot, trait_var_plot, ncol = 2,
                           labels = c("a", "b"), label_size = 40)

ggsave(fig2_combined,
       filename = file.path(fig.dir, "fig2_traitcoverage.pdf"), width = 20, height = 20)


################################################################
## FIGURE 3 from Kissling et al.: PLOTTING COMPLETENESS ONTO PHYLOGENY
# Convert trait values into binary values (data present or absent)
rownames(traitDataSubset) <- traitDataSubset$SpecName
traitDataSubset_PA <- as.data.frame(ifelse(is.na(traitDataSubset), NA,  1))
traitDataSubset_PA[,1:3] <- traitDataSubset[,1:3] # restore the taxonomic classifications

# Plot phylogeny and label clades
cladeCol <- c(wes_palette("IsleofDogs2", 5, type = "discrete"), "grey20")
bs = 4; fs = 10; ofs = 22; cl = c("grey80", "grey20")
phyloCladePlot <- ggtree(palmPhyloSubset, layout = "circular", size = 0.2) +
  geom_cladelabel(node = 4421, label = "Calamoideae", hjust = 0,
                  offset = ofs, offset.text = 1.9, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 1, label = "Nypoideae", hjust = 0.5, extend = 0.5,
                  offset = ofs, offset.text = 6.0, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 3022, label = "Arecoideae", hjust = 1,
                  offset = ofs, offset.text = 2, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 4376, label = "Ceroxyloideae", hjust = 1,
                  offset = ofs, offset.text = 1.3, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 2529, label = "Coryphoideae", hjust = 1,
                  offset = ofs, offset.text = 2.5, barsize = bs, fontsize = fs, color = cl) +
  scale_x_continuous(limits = c(0, 200)) +
  theme(panel.background = element_blank())

ggsave(plot = phyloCladePlot, file.path(fig.dir, "phyloTaxo.pdf"), height = 10, width = 10)

# Convert trait values into binary values (i.e., data present or absent)
traitListToPlot <- c("Climbing", "Max_Blade_Length_m", "MaxStemHeight_m", "AverageFruitLength_cm", "StemArmed")
for(i in 1:length(traitListToPlot)){
  traitDataSubset_PA[traitListToPlot[i]] <- factor(ifelse(is.na(traitDataSubset_PA[traitListToPlot[i]]), NA, i))
}

# Plot traits onto phylogeny
traitListCols <- c(wes_palette("Cavalcanti1", 5, type = "discrete"))
traitCoveragePlot <- gheatmap(phyloCladePlot, traitDataSubset_PA[traitListToPlot], colnames = FALSE, width = 0.15, color = "transparent") + 
  scale_fill_manual(values = traitListCols) + #ignore error message please
  theme(legend.position = "none",
        plot.background = element_blank(),
        panel.background = element_blank())
ggsave(plot = traitCoveragePlot, file.path(fig.dir, "phyloTraitCoverage.pdf"), height = 10, width = 10)

# Generate legend separately
traitListLabels <- c("Growth form", "Maximum blade length", "Maximum stem height", "Average fruit length","Stem armature")
legendDF <- data.frame(y = factor(traitListLabels, levels = traitListLabels))
legendPlot <- ggplot(data = legendDF) + geom_bar(aes(y, fill = y)) + 
  #guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text=element_text(size = 20,color = "grey20"),
        legend.key.size = unit(3, 'lines'),
        legend.justification = "center",
        legend.spacing.x = unit(0.5,'lines'), legend.position = "bottom") + 
  scale_fill_manual(values = traitListCols)
legendGrob <- get_legend(legendPlot)

# Combine plots
phylotraitCoveragePlot <- plot_grid(traitCoveragePlot,
                                    legendGrob, 
                                    nrow = 2,
                                    rel_heights= c(0.8, 0.2),
                                    scale = c(1.5, 2.5))
phyloCoveragePlotGridTitle <- ggdraw() + draw_label("Phylogenetic coverage", size = 30)
phylo_compl_plot <- plot_grid(phyloCoveragePlotGridTitle,
                              phylotraitCoveragePlot, ncol = 1, rel_heights = c(0.05, 0.95))

ggsave(phylo_compl_plot,
       filename = file.path(fig.dir, "fig3_traitcoverage.pdf"), width = 20, height = 15)

################################################################
## FIGURE 4a from Kissling et al.: MAPPING GROWTH FORM PROPORTION ONTO WORLD MAP
# Merge palm occurrences at the TDWG unit scale with trait data
occ_trait <- merge(occ, traitData[c("SpecName", "Climbing", "Acaulescent", "Erect")], by = "SpecName", all.x = TRUE)

# Only include species with are mono-morphic (i.e., unambiguously a specific growth form )
occ_trait_subset <- subset(occ_trait, !(is.na(Climbing) & is.na(Acaulescent) & is.na(Erect)))

# Species with more than 1 growth form are coded as 0.5s so they will contribute to both growth forms categories when calculating proportions
occ_trait_subset$Climbing[occ_trait_subset$Climbing >1] <- 0.5
occ_trait_subset$Acaulescent[occ_trait_subset$Acaulescent > 1] <- 0.5
occ_trait_subset$Erect[occ_trait_subset$Erect > 1] <- 0.5

# Some minor changes to original geom_scatterpie_legend function (original source code from ggtree) so labels and size of pies can be more easily modified
geom_scatterpie_legend2 <- function (radius, x, y, n = 5, labeller) {
  if (length(radius) > n) {
    radius <- unique(sapply(seq(min(radius, na.rm = T), max(radius, na.rm = T),
                                length.out = n), round))
  }
  label <- FALSE
  if (!missing(labeller)) {
    if (!inherits(labeller, "function")) {
      stop("labeller should be a function for converting radius")
    }
    label <- TRUE
  }
  dd <- data.frame(r = radius, start = 0, end = 2 * pi, x = x,
                   y = y + radius - max(radius), maxr = max(radius))
  if (label) {
    dd$label <- labeller(dd$r)
  }
  else {
    dd$label <- dd$r
  }
  list(ggforce::geom_arc_bar(aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r,
                                  start = ~start, end = ~end), data = dd, inherit.aes = FALSE),
       geom_segment(aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y +
                           r, yend = ~y + r), data = dd, inherit.aes = FALSE),
       geom_text(aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label),
                 data = dd, hjust = "left", inherit.aes = FALSE))
}

# Calculate mean proportion of each growth form in each TDWG unit
tdwg_growthform <- ddply(.data = occ_trait_subset, .variable = .(Area_code_L3),
                         .fun = summarise,
                         Climber = mean(Climbing, na.rm = TRUE),
                         Acaulescent = mean(Acaulescent, na.rm = TRUE),
                         Erect = mean(Erect, na.rm = TRUE),
                         Nsp = length(unique(SpecName)) )

# Merge growth form proportion with botanical country polygons
tdwg_growthform2 <- merge(x = tdwg_growthform, y = shape@data, all = TRUE,
                          by.x = "Area_code_L3")

# Plot map
tdwg_growthform2$radius <- log(tdwg_growthform2$Nsp+1) # plotting parameter
traitCols <- c("navyblue",wes_palette("Zissou1", n = 5)[c(5,3)])
growthform_plot <- ggplot() +
  # Plot base world map polygons (Antarctica excluded for clarity)
  geom_polygon(aes(y = lat, x = long, group = group),
              data = subset(shape_gg, !LEVEL3_COD == "ANT"), fill = "grey40") +
  geom_scatterpie(aes(y = centroid_lat, x= centroid_long, group = Area_code_L3, r = radius),
                  data = tdwg_growthform2,
                  cols = c("Climber","Acaulescent","Erect"),
                  colour = NA, alpha = 0.9) +
  coord_equal() +
  scale_fill_manual(name = "Growth Form", values = traitCols) +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank()) 
growthform_plot_wleg <- growthform_plot + geom_scatterpie_legend2(radius = tdwg_growthform2$radius, x=-130, y=-50, n = 3, labeller = function(x) round( exp( (x / 1) - 1), digits = 0)  ) + geom_text(aes(x = -130, y = -50, label = "No. of species    "), hjust  = 1)

ggsave(growthform_plot_wleg, filename = file.path(fig.dir, "growthform_piechart.pdf"), width = 9, height = 4.8)

########################################################################
## FIGURE 4b from Kissling et al.: MAPPING GROWTH FORM ONTO PHYLOGENY
# Convert trait values into binary values for plotting
rownames(traitDataSubset) <- traitDataSubset$SpecName
traitListToPlot <- c("Acaulescent", "Climbing", "Erect")

for(i in 1:length(traitListToPlot)){
  traitDataSubset[traitListToPlot[i]] <- factor(ifelse(traitDataSubset[traitListToPlot[i]] >= 1, i, NA ))
}

# Plot phylogeny with growth forms highlighted (uncomment lines if you would like tip and node labels)
cladeCol <- c(wes_palette("IsleofDogs2", 5, type = "discrete"), "grey20")
traitListCols <- c("navyblue",wes_palette("Zissou1", n = 5)[c(5,3)])
bs = 2; fs = 3; ofs = 30; cl = c("grey90", "grey20")
phyloCladePlot <- ggtree(palmPhyloSubset, layout = "circular", size = 0.2) +
  # geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 1) +
  # geom_tiplab2(size = 1)
  xlim(c(-10, 180)) +
  geom_cladelabel(node = 4421, label = "Calamoideae", hjust = 0,
                  offset = ofs, offset.text = 1.2, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 1, label = "Nypoideae", hjust = 0.5, extend = 0.5,
                  offset = 35, offset.text = 4.5, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 3022, label = "Arecoideae", hjust = 1,
                  offset = ofs, offset.text = 1.2, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 4376, label = "Ceroxyloideae", hjust = 1,
                  offset = ofs, offset.text = 0.8, barsize = bs, fontsize = fs, color = cl) +
  geom_cladelabel(node = 2529, label = "Coryphoideae", hjust = 1,
                  offset = ofs, offset.text = 1.2, barsize = bs, fontsize = fs, color = cl) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        panel.background = element_rect(fill = "transparent"))
#ggsave("~/Desktop/palmtree.pdf", phyloCladePlot, width = 20, height = 20)

traitCoveragePlot <- gheatmap(phyloCladePlot, traitDataSubset[traitListToPlot], colnames = FALSE, width = 0.2, color = "transparent") +
  scale_fill_manual(values = traitListCols,
                    breaks = c("1", "2", "3"),
                    labels = c("Acaulescent", "Climbing", "Erect"),
                    expand = c(0,0)) + 
  theme(legend.title = element_blank())
# ignore error message, it is due to our code coercing the plotting of different variables in different colours as opposed to a plotting of factor levels across variables
ggsave(plot = traitCoveragePlot, file.path(fig.dir, "phyloTrait.pdf"), height = 10, width = 10)

########################################################################
## FIGURE 4c from Kissling et al.: Perform a principal component analysis to explore information of
## continuous traits in the context of growth forms
# Standardize variables
traitData$logBladeLength <- log(traitData$Max_Blade_Length_m)
traitData$logFruitLength <- log(traitData$AverageFruitLength_cm)
traitData$logFruitWidth <- log(traitData$AverageFruitWidth_cm)
traitData$logRachisLength <- log(traitData$Max_Rachis_Length_m)
traitData$logStemHeight <- log(traitData$MaxStemHeight_m+ 1) # acaulescent palms often have underground stems, and so heights are equals to zero

# Only include species with complete trait values
targetCol <- c("SpecName","logBladeLength", "logFruitLength", "logFruitWidth", "logRachisLength", "logStemHeight")
targetTraits <- traitData[targetCol]
targetTraitsSubset <- na.omit(targetTraits)
rownames(targetTraitsSubset) <- targetTraitsSubset$SpecName

# Perform principal component analysis
traitPca <- prcomp(targetTraitsSubset[,2:6], center = TRUE, scale = TRUE)
traitPcaCoord <- as.data.frame(traitPca$x)
traitPcaCoord$SpecName <- rownames(traitPcaCoord)
traitPcaCoordRes <- merge(traitPcaCoord, traitData, by = "SpecName", all.x = TRUE)
traitPcaAxes <- as.data.frame(traitPca$rotation)
traitPcaAxes$label <- rownames(traitPcaAxes)

# Group points by life form
traitPcaCoordRes$LifeForm <- NA
traitPcaCoordRes$LifeForm[traitPcaCoordRes$Climbing == 1] <- "Climbing"
traitPcaCoordRes$LifeForm[traitPcaCoordRes$Acaulescent == 1] <- "Acaulescent"
traitPcaCoordRes$LifeForm[traitPcaCoordRes$Erect == 1] <- "Erect"
traitPcaCoordRes$LifeForm[rowSums(traitPcaCoordRes[c("Climbing", "Acaulescent","Erect")]) > 1] <- NA

growthformPCA <- ggplot(data = subset(traitPcaCoordRes, !is.na(LifeForm))) + 
  geom_point(shape = 21, alpha = 0.7,aes(y = PC1, x = PC2, color = LifeForm), size = 2) +
  stat_ellipse(type="norm", aes(y = PC1, x = PC2, color = LifeForm), show.legend = FALSE) +
  geom_segment(aes(y = 0, x = 0, yend = PC1*5, xend = PC2*5), alpha = 0.8, data = traitPcaAxes) +
  geom_text_repel(aes(y = PC1*5.5, x = PC2*5.5, label = label), alpha = 0.8, data = traitPcaAxes) + 
  theme(panel.background = element_blank(), axis.line = element_line()) +
  scale_color_manual(name = "Growth Form",
                     values = c("navyblue",wes_palette("Zissou1", n = 5)[c(5,3)]))


ggsave(growthformPCA, filename = file.path(fig.dir, "growthformPCA.pdf"), height= 5, width = 6)

########################################################################
# Combining and annotating sub-figures into a single panelled figure
growthform_leg <- get_legend(growthform_plot_wleg + theme(legend.justification="center"))
fig_panel <- plot_grid(growthform_plot_wleg + theme(legend.position = "none"),
                       plot_grid(traitCoveragePlot + theme(legend.position = "none"),
                                 growthformPCA + theme(legend.position = "none"),
                                 labels= c("b", "c"), label_size = 20,
                                 nrow =1, rel_widths = c(0.5, 0.5), scale =c(1.3, 0.95)),
          nrow= 2, labels = "a", 
          rel_heights = c(1, 0.8), label_size = 20)

# Add legend
fig_panel_wleg <- plot_grid(fig_panel, growthform_leg, nrow = 2, rel_heights = c(0.95, 0.05))
ggsave(fig_panel_wleg, filename = file.path(fig.dir, "fig3_combined.pdf"), width = 10, height = 9)
