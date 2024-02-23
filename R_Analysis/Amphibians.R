######################################## setwd,load libraries, source functions ####
setwd("/Users/ljpinnell/Documents/VERO/Amphibians16S/R_analysis/")

library(phyloseq); library(ggplot2); library(btools); library(dplyr);library(multcompView);library(rcompanion)
library(tidyr); library(stringr); library(randomcoloR); library(metagMisc);
library(metagenomeSeq); library(GUniFrac); library(randomForest); library(pairwiseAdonis)
library(knitr); library(readr); library(kableExtra); library(scales); library(vegan);
library(ggdendro); library(forcats); library(ggbreak); library(BiodiversityR); library(ggrepel); library(ggforce)

source("/Users/ljpinnell/Documents/R_Functions_Scripts/change16STaxaNames.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/w_unifrac.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/g_unifrac.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/uw_unifrac.R")
source("/Users/ljpinnell/Documents/R_Functions_Scripts/MergeLowAbund.R")

######################################## IMPORTING DATA ####
qiimedata <- import_biom("table-with-taxonomy.biom","tree.nwk","dna-sequences.fasta")
qiimedata # 6584 taxa and 83 samples
#write.csv(sample_names(qiimedata),"sample_names.csv")
map_file <- import_qiime_sample_data("metadata.txt")

# combining data with metadata
data <- merge_phyloseq(qiimedata, map_file)
data # 6584 taxa and 82 samples (lost the duplicate sample that failed (PQ51A))

data <- prune_taxa(taxa_sums(data) > 0, data)
data # 5583 taxa (guess the single ASV in the dropped sample was unique!)

######################################## PRE-PROCESSING ####
## renaming ranks
rank_names(data) # "Rank1" - "Rank7" not ideal, lets change em
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rank_names(data) # beauty, now they are named properly

# changing the SILVA style naming (k__Bacteria, etc.)
tax.data <- data.frame(tax_table(data)) # extract the taxonomy table as a data frame
tax.data.names <- change16Staxa(tax.data) # this gets rid of the GG format

# now to change the NAs to a better naming scheme
for (i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])} # converting all columns to characters
tax.data.names[is.na(tax.data.names)] <- "" # replacing the NAs with an empty string

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  }
}

head(tax.data.names) # great, no more NAs and no more k__
tax_table(data) <- as.matrix(tax.data.names) # re-insert the taxonomy table into the phyloseq object
head(tax_table(data), 20) # sweet, lookin good!

## removing non Bacteria/Archaea but keeing Unassigned so I know classification %
data # 6583 
data <- subset_taxa(data, Kingdom=="Bacteria" | Kingdom=="Archaea" | Kingdom=="Unassigned") # 
data # 6529 (lost 54 taxa)

### I think it's best to not include the on-exhibit "controls"
### they have different treatments and swabing locations
### so I don't think they are comparable or add anything

data <- subset_samples(data, group!="exhibit control" )
data <- prune_taxa(taxa_sums(data) > 0, data)
data # 6456 taxa, 74 samples

# some QC checks
min(sample_sums(data)) #1
max(sample_sums(data)) # 126,957
mean(sample_sums(data)) # 57,991
median(sample_sums(data)) # 57,575
sort(sample_sums(data)) # 1,1,1,455,3102,14545... cutoff at 10K

data <- prune_samples(sample_sums(data) > 10000, data)
data # left with 71 samples
data <- prune_taxa(taxa_sums(data) > 0, data)
data # 6450 taxa 71 samples remain

min(sample_sums(data)) #14,545
max(sample_sums(data)) # 126,957
mean(sample_sums(data)) # 60,391.41
median(sample_sums(data)) # 59,199
sort(sample_sums(data))


############### VARIABLE PALETTES ####
week_palette <- c("forestgreen","cornflowerblue","#d87cec","#af49c5","#843b94","#5b146a","#33003e")
######################################## ALPHA DIVERSITY ####
alpha_div1 <- estimate_richness(data, measures = c("Observed", "Shannon", "Simpson","InvSimpson")) # richness, div
alpha_div2 <- estimate_pd(data) # faith's pd

# combine metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div
alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div_meta <- cbind(data@sam_data, alpha_div)
alpha_div_meta # metadata and div metrics

count(alpha_div_meta, timepoint_cat)

### RICHNESS
ggplot(alpha_div_meta, aes(x= week_cat, y= Observed, fill = week_cat, colour = week_cat)) +
  theme_bw() +
  labs(y= "Observed ASVs", title = "RICHNESS") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(shape = 18, size =5) +
  scale_fill_manual(values = week_palette) +
  #geom_text(label = alpha_div_meta$SampleID) +
  scale_y_continuous(limits = c(0,1800), expand = c(0.001,0,0.04,0)) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"))

pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$week_cat, p.adjust.method = "BH")
# A from all others, thats it

### SHANNON
ggplot(alpha_div_meta, aes(x= week_cat, y= Shannon, fill = week_cat, colour = week_cat)) +
  theme_bw() +
  labs(y= "Shannon", title = "DIVERSITY") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(shape = 18, size =5) +
  scale_fill_manual(values = week_palette) +
  scale_y_continuous(limits = c(0,7), expand = c(0.001,0,0.04,0)) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"))

pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$week_cat, p.adjust.method = "BH")
# A from all others, B from C, F and G

### FPD
ggplot(alpha_div_meta, aes(x= week_cat, y= PD, fill = week_cat, colour = week_cat)) +
  theme_bw() +
  labs(y= "Faith's PD", title = "PHYLO. DIVERSITY") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(shape = 18, size =5) +
  scale_fill_manual(values = week_palette) +
  scale_y_continuous(limits = c(0,120), expand = c(0.001,0,0.04,0), breaks = c(0,30,60,90)) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"))

pairwise.wilcox.test(alpha_div_meta$PD, alpha_div_meta$week_cat, p.adjust.method = "BH")
# A from all others, B from F and G

  

######################################## BETADIVERSITY ####
#################### CSS NORMALIZATION & RA TRANSFORM ####
data.css <- phyloseq_transform_css(data, log = F)
data_ra <- transform_sample_counts(data.css, function(x) {x/sum(x)} * 100)


#################### TAX GLOMMING ####
data_phylum <- tax_glom(data_ra, taxrank = "Phylum", NArm = F)
write.csv(otu_table(data_phylum),"phylum_otus.csv")
write.csv(tax_table(data_phylum),"phylum_taxa.csv")

data_class <- tax_glom(data_ra, taxrank = "Class", NArm = F)
write.csv(otu_table(data_class),"class_otus.csv")
write.csv(tax_table(data_class),"class_taxa.csv")

data_order <- tax_glom(data_ra, taxrank = "Order", NArm = F)
write.csv(otu_table(data_order),"order_otus.csv")
write.csv(tax_table(data_order),"order_taxa.csv")

data_family <- tax_glom(data_ra, taxrank = "Family", NArm = F) # 517 families
write.csv(otu_table(data_family),"family_otus.csv")
write.csv(tax_table(data_family),"family_taxa.csv")

data_genus <- tax_glom(data_ra, taxrank = "Genus", NArm = F)
write.csv(otu_table(data_genus),"genus_otus.csv")
write.csv(tax_table(data_genus),"genus_taxa.csv")

#################### DISTANCE AND ORDINATIONS ####
data.df <- as(data.css@sam_data,"data.frame")
#### WEIGHTED
data.wunifrac <- wunifrac(data.css)
data.wunifrac.ord <- ordinate(data.css, method = "NMDS", distance = data.wunifrac)

#### GENERALIZED
data.gunifrac <- gunifrac(data.css)
data.gunifrac.ord <- ordinate(data.css, method = "NMDS", distance = data.gunifrac)

#### UNWEIGHTED
data.uwunifrac <- uwunifrac(data.css)
data.uwunifrac.ord <- ordinate(data.css, method = "NMDS", distance = data.uwunifrac)

######## CENTROIDS AND PLOTS ####
#### WEIGHTED
# find centroids
w_ordination_plot <- ordiplot(data.wunifrac.ord$points)
w_ordination_sitesLong <- sites.long(w_ordination_plot, data.df)
w_ordination_centroids <- envfit(data.wunifrac.ord ~ data.df$week_cat)

# make df for centroids
w_ordination_variable_col <- c("A","B","C","D","E","F","G")
w_ordiantion_NMDS1_col <- c(0.3043,-0.0826, -0.0575,-0.0581,-0.0423,-0.0292,-0.0331)
w_ordiantion_NMDS2_col <- c(0.0069,-0.0348,0.0101,0.0081,0.0065,-0.0225,0.0107)

w_centroids_df <- data.frame(w_ordination_variable_col, w_ordiantion_NMDS1_col, w_ordiantion_NMDS2_col)
w_centroids_df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = w_ordination_sitesLong, aes(x=axis1,y=axis2, colour= week_cat), alpha = 0.5, shape = 18, size = 3) +
  stat_ellipse(data = w_ordination_sitesLong, geom = "polygon", aes(x=axis1,y=axis2, colour =week_cat, fill = week_cat), alpha = c(0.1), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = w_centroids_df, aes(x=w_ordiantion_NMDS1_col, y=w_ordiantion_NMDS2_col), colour = week_palette, fill = week_palette, size = 12, shape = 18) +
  geom_text(data = centroids_df, aes(x=w_ordiantion_NMDS1_col, y=w_ordiantion_NMDS2_col), label = c("W","T","Q1","Q2","Q4","Q7","Q9"), colour = "white", size = 5, fontface = "bold") +
  scale_colour_manual(values = week_palette) +
  scale_fill_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

#### GENERALIZED
# find centroids
g_ordination_plot <- ordiplot(data.gunifrac.ord$points)
g_ordination_sitesLong <- sites.long(g_ordination_plot, data.df)
g_ordination_centroids <- envfit(data.gunifrac.ord ~ data.df$week_cat)

# make df for centroids
g_ordination_variable_col <- c("A","B","C","D","E","F","G")
g_ordiantion_NMDS1_col <- c(0.3043,-0.0826, -0.0575,-0.0581,-0.0423,-0.0292,-0.0331)
g_ordiantion_NMDS2_col <- c(0.0069,-0.0348,0.0101,0.0081,0.0065,-0.0225,0.0107)

g_centroids_df <- data.frame(g_ordination_variable_col, g_ordiantion_NMDS1_col, g_ordiantion_NMDS2_col)
g_centroids_df

ggplot() + theme_bw() +
  labs(x= "NMDS1", y= "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_point(data = g_ordination_sitesLong, aes(x=axis1,y=axis2, colour= week_cat), alpha = 0.5, shape = 18, size = 4.5) +
  stat_ellipse(data = g_ordination_sitesLong, geom = "polygon", aes(x=axis1,y=axis2, colour =week_cat, fill = week_cat), alpha = c(0.1), lty=2, level = 0.90, linewidth = 1) +
  geom_point(data = g_centroids_df, aes(x=g_ordiantion_NMDS1_col, y=g_ordiantion_NMDS2_col), colour = week_palette, fill = week_palette, size = 15, shape = 18) +
  geom_text(data = centroids_df, aes(x=g_ordiantion_NMDS1_col, y=g_ordiantion_NMDS2_col), label = c("W","T","Q1","Q2","Q4","Q7","Q9"), colour = "white", size = 6, fontface = "bold") +
  scale_colour_manual(values = week_palette) +
  scale_fill_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12, colour = "black"))

gunifrac.adonis <- pairwise.adonis2(data.gunifrac ~ timepoint_cat, data.df, nperm = 9999, p.adjust.methods = "BH")
gunifrac.adonis

## PERMDISP
gunifrac.disper <- betadisper(data.gunifrac, data.df$timepoint_cat)
plot(gunifrac.disper)
gunifrac.permdisp <- permutest(gunifrac.disper, permutations = 9999, pairwise = T)
write.csv(gunifrac.permdisp$pairwise,"permdisp_results.csv")

# wild comparisons
gunifrac.adonis$wild_vs_transport
gunifrac.adonis$`wild_vs_quarantine week 1`
gunifrac.adonis$`wild_vs_quarantine week 2`
gunifrac.adonis$`wild_vs_quarantine week 4`
gunifrac.adonis$`wild_vs_quarantine week 7`
gunifrac.adonis$`wild_vs_quarantine week 9`
#transport comparisons
gunifrac.adonis$`transport_vs_quarantine week 1`
gunifrac.adonis$`transport_vs_quarantine week 2`
gunifrac.adonis$`transport_vs_quarantine week 4`
gunifrac.adonis$`transport_vs_quarantine week 7`
gunifrac.adonis$`transport_vs_quarantine week 9`
# quarantine comparisons
gunifrac.adonis$`quarantine week 1_vs_quarantine week 2`
gunifrac.adonis$`quarantine week 1_vs_quarantine week 4`
gunifrac.adonis$`quarantine week 1_vs_quarantine week 7`
gunifrac.adonis$`quarantine week 9_vs_quarantine week 1`
gunifrac.adonis$`quarantine week 2_vs_quarantine week 4`
gunifrac.adonis$`quarantine week 2_vs_quarantine week 7`
gunifrac.adonis$`quarantine week 9_vs_quarantine week 2`
gunifrac.adonis$`quarantine week 4_vs_quarantine week 7`
gunifrac.adonis$`quarantine week 9_vs_quarantine week 4`
gunifrac.adonis$`quarantine week 9_vs_quarantine week 7`


#################### DENDOGRAM + RA PLOTS ####
#### DENDRO ####
data.hclust <- hclust(data.gunifrac, method = "ward.D2")
data.dendro <- as.dendrogram(data.hclust)
data.dendro.data <- dendro_data(data.dendro, type = "rectangle")
data.dendro_metadata <- as_tibble(data.css@sam_data)
data.dendro.data$labels <- data.dendro.data$labels %>%
  left_join(data.dendro_metadata, by = c("label" = "SampleID"))
data.dendro.data$labels
#write.csv(data.dendro.data$labels,"dendro_meta.csv")
# want to add W,T, Q1-Q9 for labels 
short_name_label <- c("W",	"W",	"W",	"W",	"W",	"W",	"W",	"W",	"Q2",	"T",	"T",	"T",	"T",	"T",	"T",	"Q2",	"Q4",	"T",	"T",	"T",	"T",	"Q4",	"Q7",	"Q9",	"Q9",	"Q9",	"Q9",	"Q9",	"Q7",	"Q9",	"Q9",	"Q9",	"Q9",	"Q7",	"Q9",	"Q7",	"W",	"Q9",	"Q2",	"Q2",	"Q2",	"Q1",	"Q1",	"Q1",	"Q1",	"Q1",	"Q1",	"Q1",	"Q1",	"Q1",	"Q1",	"Q9",	"Q4",	"Q7",	"Q9",	"Q9",	"Q9",	"Q9",	"W",	"Q9",	"Q4",	"Q7",	"Q4",	"Q2",	"Q4",	"Q4",	"Q4",	"Q4",	"Q2",	"Q2",	"Q2")
data.dendro.data$labels <- cbind(data.dendro.data$labels, short_name_label)
data.dendro.data$labels

ggplot(data.dendro.data$segments) +
  theme_minimal() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = data.dendro.data$labels, 
             aes(x=x,y=y, colour = week_cat, fill= week_cat),
             size = 8, shape=15, stroke =1.5, position = position_nudge(y=-0.15)) +
  geom_text(data = data.dendro.data$labels, aes(x=x, y=y, label = short_name_label), colour = "white", size = 4.2, position = position_nudge(y=-0.14), fontface = "bold") +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = week_palette) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

#### RA PLOT FOR UNDER DENDRO ####
dendro_sample_order <- data.dendro.data$labels$label
# GOING TO DO CLASS, THEN BREAK IT DOWN BY FAMILY
class_melt <- psmelt(data_class) # 123 classes

# make palette for class (no. colours = 123)
#ra_class_palette <- distinctColorPalette(123)
#write.csv(ra_class_palette, "class_palette.csv")
ra_class_palette <- c("#E4BC2D",	"#D2F3D0",	"#AAD4EE",	"#C9C3C1",	"#AECCB7",	"#EF5A48",	"#F52EE8",	"#593AA8",	"#E4679B",	"#9ED9B8",	"#6DCD2E",	"#93B7ED",	"#EBB7A3",	"#CAF2B8",	"#BD8370",	"#5FC27E",	"#5E1FF4",	"#929FED",	"#AC889A",	"#E6667C",	"#8BEEC7",	"#E73D77",	"#EC8BD2",	"#558BE8",	"#C09647",	"#A545B9",	"#C9B070",	"#A89E73",	"#DFF1E0",	"#8AEAAE",	"#F1A4AC",	"#748F65",	"#919EAC",	"#DDA7BE",	"#9C4958",	"#47BF9B",	"#36C4CF",	"#65ADCF",	"#51D3C7",	"#CAEF4D",	"#54EDBF",	"#7D6640",	"#CF40E3",	"#EDEC5C",	"#E68FB1",	"#E6C9A2",	"#BEEDE7",	"#A3EB95",	"#8FF582",	"#DEBF55",	"#AF76DC",	"#495888",	"#E8EB31",	"#F4329D",	"#C0AFCA",	"#E6AEE0",	"#E58AED",	"#C8EC8C",	"#419B65",	"#7AB996",	"#4F98DC",	"#A9B73C",	"#A7F42A",	"#E466E7",	"#F0CCEE",	"#317568",	"#ABB980",	"#CDD4A9",	"#B2E75E",	"#B7CBF1",	"#9BC586",	"#DA5BAF",	"#B28BC8",	"#6EAFB5",	"#A14233",	"#F7E27C",	"#5950EA",	"#89D3CB",	"#66F03A",	"#EDEDB9",	"#6E6EE7",	"#DDE870",	"#55C9EE",	"#777FE3",	"#E9EFED",	"#846AB2",	"#D0E2ED",	"#F3E7D0",	"#A95097",	"#868DB7",	"#955B7C",	"#F1AC73",	"#52F1EA",	"#E8C6C4",	"#49E3CA",	"#F3E0EB",	"#586B76",	"#66F26C",	"#CAD5C1",	"#48CF63",	"#B7F7E0",	"#D9DAF1",	"#65E3F4",	"#C0B4EA",	"#4ABBF5",	"#9AAD55",	"#42EB9E",	"#F38E8A",	"#EAEEA2",	"#BE98E7",	"#77C44E",	"#A0BAC5",	"#E5A1EC",	"#423480",	"#A6A097",	"#9E4FE8",	"#89F7F0",	"#EBD888",	"#E141BE",	"#BF2FF3",	"#6823C7",	"#EC8C2F",	"#A4DBE4")

ggplot(class_melt, aes(x= SampleID, y= Abundance, fill= Class)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(limits = dendro_sample_order, expand = c(0.03,0,0.015,0)) +
  scale_fill_manual(values = ra_class_palette) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", size = 0.75),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

#################### LET'S LOOK AT THE 5 CLASSES OVER 1% AT THE FAMILY LEVEL ####
######## GAMMAPROTEO
gammaproteo_class <- subset_taxa(data_class, Class=="Gammaproteobacteria")
#write.csv(otu_table(gammaproteo_class),"gammaproteo_class_ra.csv")
#write.csv(data.df, "metadata_phylo.csv")
gammaproteo_class_melt <- psmelt(gammaproteo_class)
gammaproteo_family <- subset_taxa(data_family, Class=="Gammaproteobacteria") # 72 families
write.csv(tax_table(gammaproteo_family),"gammaproteo_family.csv")
#write.csv(otu_table(gammaproteo_family),"gamma_otus.csv")
gammaproteo_family_melt <- psmelt(gammaproteo_family)
#gammaproteo_palette <- distinctColorPalette(72)
#write.csv(gammaproteo_palette,"gammaproteo_palette.csv")
gammaproteo_palette <- c("#51E9BF",	"#4C9DE0",	"#E1EAED",	"#A3EC97",	"#AD8FD7",	"#C4D0F1",	"#9FC787",	"#A9D5E7",	"#E7A2BA",	"#D5F0B5",	"#6CF5EC",	"#ED574B",	"#BB42A4",	"#6C6DE3",	"#E9C1A4",	"#E4BD2F",	"#B25585",	"#6D983B",	"#6FEE39",	"#E0A9DD",	"#5ADCEA",	"#582DE4",	"#9348DE",	"#E356DB",	"#EDCAED",	"#94B7EE",	"#93825D",	"#E0D7A6",	"#6A8BE8",	"#A2AAAD",	"#EADEEA",	"#DDBD56",	"#CDED8E",	"#E2C4BF",	"#E77FC2",	"#DF9080",	"#A04337",	"#E9EBD1",	"#4F7E72",	"#52D38A",	"#ADBC45",	"#6E94AC",	"#E9EA44",	"#8587B7",	"#DEAE6F",	"#B9B0E5",	"#E6EF71",	"#5BBFE5",	"#58D9C5",	"#503690",	"#ED7992",	"#D832E9",	"#81CEC8",	"#5DDD5E",	"#A7AF78",	"#8FEEBD",	"#9C7380",	"#B1E53E",	"#B6CCB8",	"#E68D39",	"#C0EDE7",	"#A3EC72",	"#C2AAC2",	"#E49DF0",	"#E9448C",	"#4FF29D",	"#AD75DB",	"#B9ECCD",	"#E57DE7",	"#56B082",	"#EFE18E",	"#4B5881")

ggplot() +
  coord_flip() + theme_bw() +
  labs(y= "Relative Abundance (%)", title = "GAMMAPROTEOBACTERIA") +
  geom_bar(data = gammaproteo_family_melt, aes(x= week_cat, y= Abundance, fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = gammaproteo_class_melt, aes(x= week_cat, y= Abundance), stat = "summary", width = 0.45) +
  scale_fill_manual(values = gammaproteo_palette) +
  scale_x_discrete(limits = c("G","F","E","D","C","B","A"), labels = c("WEEK 9","WEEK 7","WEEK 4","WEEK 2","WEEK 1","TRANSPORT","WILD")) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,0.001,0)) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(size = 32),
        plot.margin = unit(c(0.5,1,0.5,0.5),"lines"),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.ticks = element_line(colour = "black", size = 0.8))

gammaproteo_anova <- pairwise.wilcox.test(gammaproteo_class_melt$Abundance, gammaproteo_class_melt$timepoint_cat, p.adjust.method = "BH") # wild from all, transport from all, no others

gammaproteo_pvalues <- gammaproteo_anova$p.value %>%
  fullPTable()
multcompLetters(gammaproteo_pvalues, reversed = T)

#### TOP 5 MOST ABUNDANT FAMILIES
## Moraxellaceae, Pseudomonadaceae, Comamonadaceae, Alcaligenaceae, Yersiniaceae, Enterobacteriaceae

### Moraxellaceae
moraxellaceae <- subset_taxa(data_family, Family=="Moraxellaceae") %>%
  psmelt()

ggplot(moraxellaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "MORAXELLACEAE") +
  geom_boxplot(linewidth = 1, fill = "#E9EBD1", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,81.7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

moraxellaceae_anova <- pairwise.wilcox.test(moraxellaceae$Abundance, moraxellaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
morax_pvalues <- moraxellaceae_anova$p.value %>%
  fullPTable()
multcompLetters(morax_pvalues, reversed = T)

### Pseudomonadaceae
pseudomonadaceae <- subset_taxa(data_family, Family=="Pseudomonadaceae") %>%
  psmelt()

ggplot(pseudomonadaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "PSEUDOMONADACEAE") +
  geom_boxplot(linewidth = 1, fill = "#503690", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,81.7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

pseudomonadaceae_anova <- pairwise.wilcox.test(pseudomonadaceae$Abundance, pseudomonadaceae$timepoint_cat, p.adjust.method = "BH") # anova
pseudomonadaceae_pvalues <- pseudomonadaceae_anova$p.value %>%
  fullPTable() # grab the p-values and use that fullPTable to convert from the half to full version of a matrix
multcompLetters(pseudomonadaceae_pvalues, reversed = T) # get the letters, and put on figure in photoshop

### Comamonadaceae
comamonadaceae <- subset_taxa(data_family, Family=="Comamonadaceae") %>%
  psmelt()

ggplot(comamonadaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "COMAMONADACEAE") +
  geom_boxplot(linewidth = 1, fill = "#6FEE39", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,81.7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

comamonadaceae_anova <- pairwise.wilcox.test(comamonadaceae$Abundance, comamonadaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
comamonadaceae_pvalues <- comamonadaceae_anova$p.value %>%
  fullPTable()
multcompLetters(comamonadaceae_pvalues, reversed = T)

### Alcaligenaceae
Alcaligenaceae <- subset_taxa(data_family, Family=="Alcaligenaceae") %>%
  psmelt()

ggplot(Alcaligenaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "ALCALIGENACEAE") +
  geom_boxplot(linewidth = 1, fill = "#A3EC97", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,81.7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Alcaligenaceae_anova <- pairwise.wilcox.test(Alcaligenaceae$Abundance, Alcaligenaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Alcaligenaceae_pvalues <- Alcaligenaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Alcaligenaceae_pvalues, reversed = T)

### Yersiniaceae
Yersiniaceae <- subset_taxa(data_family, Family=="Yersiniaceae") %>%
  psmelt()

ggplot(Yersiniaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "YERSINIACEAE") +
  geom_boxplot(linewidth = 1, fill = "#EFE18E", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,81.7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Yersiniaceae_anova <- pairwise.wilcox.test(Yersiniaceae$Abundance, Yersiniaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Yersiniaceae_pvalues <- Yersiniaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Yersiniaceae_pvalues, reversed = T)

### Enterobacteriaceae
Enterobacteriaceae <- subset_taxa(data_family, Family=="Enterobacteriaceae") %>%
  psmelt()

ggplot(Enterobacteriaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "ENTEROBACTERIACEAE") +
  geom_boxplot(linewidth = 1, fill = "#582DE4", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,81.7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Enterobacteriaceae_anova <- pairwise.wilcox.test(Enterobacteriaceae$Abundance, Enterobacteriaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Enterobacteriaceae_pvalues <- Enterobacteriaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Enterobacteriaceae_pvalues, reversed = T)

######## BACTEROIDIA
Bacteroidia_class <- subset_taxa(data_class, Class=="Bacteroidia")
#write.csv(otu_table(Bacteroidia_class),"bacteroidia_class_ra.csv")
Bacteroidia_class_melt <- psmelt(Bacteroidia_class)
Bacteroidia_family <- subset_taxa(data_family, Class=="Bacteroidia") # 38 families
#write.csv(tax_table(Bacteroidia_family),"Bacteroidia_family.csv")
#write.csv(otu_table(Bacteroidia_family),"Bacteroidia_otus.csv")
Bacteroidia_family_melt <- psmelt(Bacteroidia_family)
Bacteroidia_palette <- distinctColorPalette(38)
#write.csv(Bacteroidia_palette,"Bacteroidia_palette.csv")
Bacteroidia_palette <- c("#8E8960",	"#69ECAE",	"#D88F7E",	"#DBBDB6",	"#DEE83E",	"#E5C553",	"#706C95",	"#ABBEDF",	"#E1DDE6",	"#E198E2",	"#E66EE9",	"#58A3E0",	"#E08840",	"#CE68AF",	"#8AEC42",	"#D531EC",	"#9F75D6",	"#EB3EB6",	"#83B143",	"#73CDE6",	"#B8EAB1",	"#AB50D4",	"#5F81DF",	"#E1C88D",	"#B2E5DC",	"#C1879E",	"#6AE8DC",	"#6BB980",	"#EBBFDF",	"#DEE8C9",	"#76A39F",	"#E05057",	"#B8A9E8",	"#C2E582",	"#7731E0",	"#E66997",	"#5A5AE5",	"#67E272")

ggplot() +
  coord_flip() + theme_bw() +
  labs(y= "Relative Abundance (%)", title = "BACTEROIDIA") +
  geom_bar(data = Bacteroidia_family_melt, aes(x= week_cat, y= Abundance, fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = Bacteroidia_class_melt, aes(x= week_cat, y= Abundance), stat = "summary", width = 0.45) +
  scale_fill_manual(values = Bacteroidia_palette) +
  scale_x_discrete(limits = c("G","F","E","D","C","B","A"), labels = c("WEEK 9","WEEK 7","WEEK 4","WEEK 2","WEEK 1","TRANSPORT","WILD")) +
  scale_y_continuous(limits = c(0,100), expand = c(0.001,0,0.001,0)) +
  theme(#legend.position = "bottom",
        legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(size = 32),
        plot.margin = unit(c(0.5,1,0.5,0.5),"lines"),
        panel.grid.major.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(colour = "black", size = 14),
        axis.ticks = element_line(colour = "black", size = 0.8))

Bacteroidia_class_melt_anova <- pairwise.wilcox.test(Bacteroidia_class_melt$Abundance, Bacteroidia_class_melt$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Bacteroidia_class_melt_pvalues <- Bacteroidia_class_melt_anova$p.value %>%
  fullPTable()
multcompLetters(Bacteroidia_class_melt_pvalues, reversed = T)

##### 6 most abundant families
### Sphingobacteriaceae, Weeksellaceae, Flavobacteriaceae, Chitinophagaceae, Spirosomaceae, Bacteroidaceae

## Sphingobacteriaceae
Sphingobacteriaceae <- subset_taxa(data_family, Family=="Sphingobacteriaceae") %>%
  psmelt()

ggplot(Sphingobacteriaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "SPHINGOBACTERIACEAE") +
  geom_boxplot(linewidth = 1, fill = "#EBBFDF", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,41.09)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Sphingobacteriaceae_anova <- pairwise.wilcox.test(Sphingobacteriaceae$Abundance, Sphingobacteriaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Sphingobacteriaceae_pvalues <- Sphingobacteriaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Sphingobacteriaceae_pvalues, reversed = T)

## Weeksellaceae
Weeksellaceae <- subset_taxa(data_family, Family=="Weeksellaceae") %>%
  psmelt()

ggplot(Weeksellaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "WEEKSELLACEAE") +
  geom_boxplot(linewidth = 1, fill = "#5A5AE5", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,41.09)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Weeksellaceae_anova <- pairwise.wilcox.test(Weeksellaceae$Abundance, Weeksellaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Weeksellaceae_pvalues <- Weeksellaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Weeksellaceae_pvalues, reversed = T)

## Flavobacteriaceae
Flavobacteriaceae <- subset_taxa(data_family, Family=="Flavobacteriaceae") %>%
  psmelt()

ggplot(Flavobacteriaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "FLAVOBACTERIACEAE") +
  geom_boxplot(linewidth = 1, fill = "#8AEC42", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,41.09)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Flavobacteriaceae_anova <- pairwise.wilcox.test(Flavobacteriaceae$Abundance, Flavobacteriaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Flavobacteriaceae_pvalues <- Flavobacteriaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Flavobacteriaceae_pvalues, reversed = T)

## Chitinophagaceae
Chitinophagaceae <- subset_taxa(data_family, Family=="Chitinophagaceae") %>%
  psmelt()

ggplot(Chitinophagaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "CHITINOPHAGACEAE") +
  geom_boxplot(linewidth = 1, fill = "#ABBEDF", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,41.09)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Chitinophagaceae_anova <- pairwise.wilcox.test(Chitinophagaceae$Abundance, Chitinophagaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Chitinophagaceae_pvalues <- Chitinophagaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Chitinophagaceae_pvalues, reversed = T)

## Spirosomaceae
Spirosomaceae <- subset_taxa(data_family, Family=="Spirosomaceae") %>%
  psmelt()

ggplot(Spirosomaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "SPIROSOMACEAE") +
  geom_boxplot(linewidth = 1, fill = "#DEE8C9", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,41.09)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Spirosomaceae_anova <- pairwise.wilcox.test(Spirosomaceae$Abundance, Spirosomaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Spirosomaceae_pvalues <- Spirosomaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Spirosomaceae_pvalues, reversed = T)

## Bacteroidaceae
Bacteroidaceae <- subset_taxa(data_family, Family=="Bacteroidaceae") %>%
  psmelt()

ggplot(Bacteroidaceae, aes(x= week_cat, y= Abundance)) + 
  theme_bw() +
  labs(y= "Relative Abundance (%)", title = "BACTEROIDACEAE") +
  geom_boxplot(linewidth = 1, fill = "#DBBDB6", colour = "black") + geom_point(size = 3) +
  scale_x_discrete(labels = c("W","T","Q1","Q2","Q4","Q7","Q9")) +
  scale_y_continuous(limits = c(0,41.09)) +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"lines"),
        plot.title = element_text(size = 24),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 24, angle = 45, hjust = .95, vjust = 0.95),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 14, colour = "black")) 

Bacteroidaceae_anova <- pairwise.wilcox.test(Bacteroidaceae$Abundance, Bacteroidaceae$timepoint_cat, p.adjust.method = "BH") # wild from all others, that's it
Bacteroidaceae_pvalues <- Bacteroidaceae_anova$p.value %>%
  fullPTable()
multcompLetters(Bacteroidaceae_pvalues, reversed = T)

##### ANTIFUNGAL DB ALIGNMENT ####
antifungal_ASVs <- read.table("antifungalDB/all_antifungal_ASVs.txt", header = F) # bring in the ASVs that aligned to the antifungal DB
head(antifungal_ASVs)
str(antifungal_ASVs) #328 deduped ASV IDs
antifungal_taxa <- antifungal_ASVs$V1 # create character array for pulling these ASVs from phyloseq object
str(antifungal_taxa)

antifungal_ra <- prune_taxa(antifungal_taxa, data_ra)
antifungal_ra_kingdom <- tax_glom(antifungal_ra, taxrank = "Kingdom", NArm = F) %>%
  psmelt()

antifungal_ra_family <- tax_glom(antifungal_ra, taxrank = "Family", NArm = F) %>%
  psmelt()
length(unique(antifungal_ra_family$Family)) # 50 families

#antifungal_family_palette <- distinctColorPalette(50)
antifungal_family_palette <- c("#D1C8B8",	"#6EB9E4",	"#B762E3",	"#4E8741",	"#53AE89",	"#6FD8E8",	"#B4EE8F",	"#DAE198",	"#AE8EE2",	"#AF964D",	"#9bafd7",	"#de7b32",	"#EAAAE1",	"#E9C3D6",	"#4500dd",	"#DC9F2D",	"#568EA1",	"#E17B82",	"#7bed33",	"#7DD044",	"#DB9E89",	"#E7DE74",	"#DF7FD8",	"#89C0B3",	"#c39be5",	"#C5F339",	"#B4CDE1",	"#BCEBB7",	"#e4e936",	"#E13E71",	"#3e237d",	"#C8F06C",	"#66F03C",	"#686FDF",	"#e5afd7",	"#A58DB1",	"#d7e4bd",	"#B0F0DF",	"#E54BB6",	"#81516C",	"#D870A7",	"#DF6B47",	"#D93CE1",	"#944EA4",	"#E3EAEC",	"#E8DF3A",	"#6A91D9",	"#473ede",	"#6FEDC3",	"#ebdc7c")

ggplot(antifungal_ra_family, aes(x= week_cat, y=Abundance)) +
  coord_flip() + theme_bw() +
  labs(y= "Relative Abundance (%)", title = "ANTIFUNGAL TAXA") +
  geom_bar(data = antifungal_ra_family, aes(x= week_cat, y= Abundance, fill = Family), stat = "summary", colour = "black") +
  geom_errorbar(data = antifungal_ra_kingdom, aes(x= week_cat, y= Abundance), stat = "summary", width = 0.45) +
  scale_fill_manual(values = antifungal_family_palette) +
  scale_x_discrete(limits = c("G","F","E","D","C","B","A"), labels = c("WEEK 9","WEEK 7","WEEK 4","WEEK 2","WEEK 1","TRANSPORT","WILD")) +
  scale_y_continuous(limits = c(0,100), expand = c(0.0005,0,0.0005,0)) +
  theme(#legend.position = "bottom",
    legend.position = "none",
    panel.border = element_rect(colour = "black", linewidth = 1),
    plot.title = element_text(size = 40),
    plot.margin = unit(c(0.5,1,0.5,0.5),"lines"),
    panel.grid.major.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 24),
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.ticks = element_line(colour = "black", size = 0.7))

antifungal_RA_test <- pairwise.wilcox.test(antifungal_ra_kingdom$Abundance, antifungal_ra_kingdom$week_cat, p.adjust.method = "BH")
antifungal_RA_pvalues <- antifungal_RA_test$p.value %>%
  fullPTable()
multcompLetters(antifungal_RA_pvalues, reversed = F)

## antifungal alpha-div
antifungal_data <- prune_taxa(antifungal_taxa, data)

antifungal_alpha <- estimate_richness(antifungal_data)
antifungal_alpha2 <- estimate_pd(antifungal_data)
antifungal_alpha_meta <- cbind(data@sam_data, antifungal_alpha, antifungal_alpha2)

ggplot(antifungal_alpha_meta, aes(x= week_cat, y= Observed, fill = week_cat, colour = week_cat)) +
  theme_bw() +
  labs(y= "Richness") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(shape = 18, size =5) +
  scale_fill_manual(values = week_palette) +
  #geom_text(label = alpha_div_meta$SampleID) +
  scale_y_continuous(limits= c(0,84), expand = c(0.0015,0,0.08,0)) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"))

antifungal_richness_test <- pairwise.wilcox.test(antifungal_alpha_meta$Observed, antifungal_alpha_meta$week_cat, p.adjust.method = "BH")
antifungal_richness_pvalues <- antifungal_richness_test$p.value %>%
  fullPTable()
multcompLetters(antifungal_richness_pvalues, reversed = F)

ggplot(antifungal_alpha_meta, aes(x= week_cat, y= Shannon, fill = week_cat, colour = week_cat)) +
  theme_bw() +
  labs(y= "Shannon") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(shape = 18, size =5) +
  scale_fill_manual(values = week_palette) +
  #geom_text(label = alpha_div_meta$SampleID) +
  scale_y_continuous(limits= c(0,4), expand = c(0.0015,0,0.08,0)) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"))

antifungal_shannon_test <- pairwise.wilcox.test(antifungal_alpha_meta$Shannon, antifungal_alpha_meta$week_cat, p.adjust.method = "BH")
antifungal_shannon_pvalues <- antifungal_shannon_test$p.value %>%
  fullPTable()
multcompLetters(antifungal_shannon_pvalues, reversed = F)

ggplot(antifungal_alpha_meta, aes(x= week_cat, y= PD, fill = week_cat, colour = week_cat)) +
  theme_bw() +
  labs(y= "Faith's PD") +
  geom_boxplot(alpha = 0.5, linewidth = 1) + geom_point(shape = 18, size =5) +
  scale_fill_manual(values = week_palette) +
  #geom_text(label = alpha_div_meta$SampleID) +
  scale_y_continuous(limits= c(0,5), expand = c(0.0015,0,0.08,0)) +
  scale_colour_manual(values = week_palette) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 32),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 32),
        axis.text.y = element_text(size = 14, colour = "black"))

antifungal_FPD_test <- pairwise.wilcox.test(antifungal_alpha_meta$PD, antifungal_alpha_meta$week_cat, p.adjust.method = "BH")
antifungal_FPD_pvalues <- antifungal_FPD_test$p.value %>%
  fullPTable()
multcompLetters(antifungal_FPD_pvalues, reversed = F)
