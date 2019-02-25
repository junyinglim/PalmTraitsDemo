## Script used to make the Scientific Data paper Figures v.1 by Renske Onstein: maps and the distribution of traits
# Author: Jun Ying Lim

## DIRECTORIES ====================
main.dir <- "~/Dropbox/Projects/2019/palms/projects/palmTraits/demo/PalmTraitsDemo/"
fig.dir <- file.path(main.dir, "figs")
data.dir <- file.path(main.dir, "data")

## LOAD PACKAGES ====================
# spatial packages for maps
library(rgdal); library(rgeos)

# tidy packages for data handling and summary statistics
library(reshape2); library(plyr)

# phylogenetic package
library(ape)

# general plotting packages
library(viridis); library(ggplot2);library(wesanderson); library(ggpubr); library(cowplot); library(viridis); library(ggrepel); library(scatterpie)

## IMPORT DATA ====================
# load spatial data 
shape <- readOGR(dsn = file.path(main.dir, "tdwg_level3_shp", "."), layer = "level3")
shape@data$id <- rownames(shape@data)

# simplify polygons for easier plotting
shape2 <- gSimplify(spgeom = shape, tol = 0.01, topologyPreserve = TRUE) 
shape <- SpatialPolygonsDataFrame(shape2, shape@data)
shape_centroid <- gCentroid(shape, byid = TRUE)
shape@data$centroid_long <- shape_centroid$x
shape@data$centroid_lat <- shape_centroid$y
shape_fort <- fortify(shape, id = id)
shape_gg <- merge(shape_fort, shape@data, by = "id")

# load palm trait data
traitData <- read.csv(file.path(data.dir, "PalmTraits_10.csv"), stringsAsFactors = FALSE)

# load palm occurence data
occ <- read.csv(file.path(data.dir, "palms_in_tdwg3.csv"))
occ$SpecName <- gsub(occ$SpecName, pattern = "_", replacement = " ")

# load phylogeny
palmPhylo <- read.nexus(file.path(data.dir, "TREE.nex"))
palmPhylo$tip.label <- gsub(palmPhylo$tip.label, pattern = "_", replacement = " ")

setdiff(palmPhylo$tip.label, traitData$SpecName)
setdiff(traitData$SpecName, palmPhylo$tip.label)

# exclude those taxa for now
intersectTaxa <- intersect(traitData$SpecName, palmPhylo$tip.label)
traitDataSubset <- subset(traitData, SpecName %in% intersectTaxa)
palmPhyloSubset <- drop.tip(palmPhylo,
                            tip = palmPhylo$tip.label[!palmPhylo$tip.label %in% intersectTaxa])

## VIGNETTE 1: MAPPING GROWTH FORM PROPORTION ONTO WORLD MAP
# Merge palm occurrences at the TDWG unit scale with trait data
occ_trait <- merge(occ, traitData[c("SpecName", "Climbing", "Acaulescence", "Errect")], by = "SpecName", all.x = TRUE)

# Only include species with are mono-morphic (i.e., unambiguously )
occ_trait_subset <- subset(occ_trait, !(Climbing > 1 | Acaulescence > 1 | Errect > 1))

# Some minor changes to original geom_scatterpie_legend function so labels and size of pies can be more easily modified
geom_scatterpie_legend2 <- function (radius, x, y, n = 5, labeller) {
  #radius = rnorm(n = 10, 1,1)
  #n = 4; round_digit = 0
  if (length(radius) > n) {
    radius <- unique(sapply(seq(min(radius), max(radius),
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

# Calculate mean proportion of each growth form, here, we can simply
tdwg_growthform <- ddply(.data = occ_trait_subset, .variable = .(Area_code_L3),
                         .fun = summarise,
                         Climber = mean(Climbing, na.rm = TRUE),
                         Acaulescent = mean(Acaulescence, na.rm = TRUE),
                         Erect = mean(Errect, na.rm = TRUE),
                         Nsp = length(unique(SpecName)) )
tdwg_growthform2 <- merge(x = tdwg_growthform, y = shape@data, all.x = TRUE,
                          by.x = "Area_code_L3",
                          by.y = "LEVEL3_COD")
tdwg_growthform2$radius <- log10(tdwg_growthform2$Nsp) * 2.5 # plotting parameter
sum(tdwg_growthform2$Nsp == 0)
traitCols <- c("navyblue",wes_palette("Zissou1", n = 5)[c(5,3)])
growthform_plot <- ggplot() +
  # Plot base world map polygons (Antarctica excluded for clarity)
  geom_polygon(aes(y = lat, x = long, group = group),
               data = subset(shape_gg, !LEVEL3_COD == "ANT"), fill = "grey40") +
  geom_scatterpie(aes(y = centroid_lat, x= centroid_long, group = Area_code_L3, r = radius),
                  data = tdwg_growthform2,
                  cols = c("Climber","Acaulescent","Erect"), colour = NA, alpha = 0.9) +
  coord_equal() +
  scale_fill_manual(name = "Growth Form", values = traitCols) +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()) 
growthform_plot_wleg <- growthform_plot + geom_scatterpie_legend2(radius = tdwg_growthform2$radius, x=-130, y=-50, n = 3, labeller = function(x) round(10^(x / 2.5), digits = 0)  ) + geom_text(aes(x = -130, y = -50, label = "No. of species    "), hjust  = 1)

ggsave(growthform_plot_wleg, filename = file.path(fig.dir, "growthform_piechart.pdf"), width = 9, height = 4.8)

## VIGNETTE 2: MAPPING GROWTH FORM ONTO PHYLOGENY ========

# Convert trait values into binary values for plotting
rownames(traitDataSubset) <- traitDataSubset$SpecName
traitDataSubset_PA <- as.data.frame(ifelse( is.na(traitDataSubset), NA,  ifelse(traitDataSubset >= 1, 1, NA )))
traitDataSubset_PA[,1:3] <- traitDataSubset[,1:3] # restore the taxonomic classifications

# Plot circular phylogeny with subfamilies highlighted
cladeCol <- c(wes_palette("IsleofDogs2", 5, type = "discrete"), "grey20")
barsize = 2; ofs = 20

traitListToPlot <- c("Acaulescence", "Climbing", "Errect")
for(i in 1:length(traitListToPlot)){
  traitDataSubset_PA[traitListToPlot[i]] <- factor(ifelse(is.na(traitDataSubset_PA[traitListToPlot[i]]), NA, i))
}
traitListCols <- c("navyblue",wes_palette("Zissou1", n = 5)[c(5,3)])
#c(wes_palette("Cavalcanti1", 3, type = "discrete"))

phyloCladePlot <- ggtree(palmPhyloSubset, layout = "circular", size = 0.2) +
  geom_cladelabel(node = 3025, barsize = barsize, color = cladeCol[c(1,6)], label = "Arecoideae", hjust = 1, offset = ofs, offset.text = 1.2) + #"Arecoideae"
  geom_cladelabel(node = 4380, barsize = barsize, color = cladeCol[c(2,6)], label = "Ceroxyloideae", hjust = 1, offset = ofs, offset.text = 0.8) +#"Ceroxyloideae"
  geom_cladelabel(node = 4425, barsize = barsize, color = cladeCol[c(3,6)], label = "Calamoideae", hjust = 0, offset = ofs, offset.text = 1.2) +#"Calamoideae"
  geom_cladelabel(node = 2531, barsize = barsize, color = cladeCol[c(4,6)], label = "Coryphoideae", hjust = 1, offset = ofs, offset.text = 1.2) + #"Coryphoideae"
  geom_cladelabel(node = 1, barsize = barsize, extend = 0.5, color = cladeCol[c(5,6)], label = "Nypoideae", hjust = 0.5, offset = ofs, offset.text = 3) #+ #"Nypoideae"

traitCoveragePlot <- gheatmap(phyloCladePlot, traitDataSubset_PA[traitListToPlot], colnames = FALSE, width = 0.15) +
  scale_fill_manual(values = traitListCols,
                    breaks = c("1", "2", "3"),
                    labels = c("Acaulescent", "Climber","Erect"))
ggsave(plot = traitCoveragePlot, file.path(fig.dir, "phyloTrait.pdf"), height = 10, width = 10)
# ignore error message, it is due to our code coercing the plotting of different variables in different colours as opposed to a plotting of factor levels across variables

## VIGNETTE 3: Perform a principal component analysis
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
traitPcaCoordRes$LifeForm[traitPcaCoordRes$Acaulescence == 1] <- "Acaulescent"
traitPcaCoordRes$LifeForm[traitPcaCoordRes$Errect == 1] <- "Erect"
traitPcaCoordRes$LifeForm[rowSums(traitPcaCoordRes[c("Climbing", "Acaulescence","Errect")]) > 1] <- NA

growthformPCA <- ggplot(data = subset(traitPcaCoordRes, !is.na(LifeForm))) + 
  geom_point(shape = 21, alpha = 0.7,aes(y = PC1, x = PC2, color = LifeForm), size = 2) +
  stat_ellipse(type="norm", aes(y = PC1, x = PC2, color = LifeForm), show.legend = FALSE) +
  geom_segment(aes(y = 0, x = 0, yend = PC1*5, xend = PC2*5), alpha = 0.8, data = traitPcaAxes) +
  geom_text_repel(aes(y = PC1*5.5, x = PC2*5.5, label = label), alpha = 0.8, data = traitPcaAxes) + 
  theme(panel.background = element_blank(), axis.line = element_line()) +
  scale_color_manual(name = "Growth Form",
                     values = c("navyblue",wes_palette("Zissou1", n = 5)[c(5,3)]))


ggsave(growthformPCA, filename = file.path(fig.dir, "growthformPCA.pdf"), height= 5, width = 6)
