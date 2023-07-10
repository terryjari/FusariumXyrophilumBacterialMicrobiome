ps_16s <- readRDS("~/Downloads/ps_16s.rds")

## Remove unexpected ASVs (i.e., non bacterial ASVs)
library(phyloseq)
ps_16s_bacteria <- subset_taxa(ps_16s, domain == "Bacteria")
ps_16s_bacteria <- subset_taxa(ps_16s_bacteria, phylum != "Cyanobacteria")
ps_16s_bacteria_wo_mitochondria <- subset_taxa(ps_16s_bacteria, (family!="Mitochondria") | is.na(family))
library(dplyr)
## Extract fasta file 
ps_16s_bacteria %>%
  refseq() %>%
  Biostrings::writeXStringSet("asv.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

## Calculatetotal count of ASVs in each samples 
ps_16s_bacteria_asv <- as.data.frame(rowSums(otu_table(ps_16s_bacteria_wo_mitochondria)))
write.csv(ps_16s_bacteria_asv, "ps_16s_asv_bacteria.csv")

## Check which ASVs assigend to the mitochondria
library(devtools)
source_url("https://github.com/umerijaz/microbiomeSeq/blob/master/R/taxa_level.R")

mitochondria <- subset_taxa(ps_16s_bacteria, (family=="Mitochondria"))
mitochondria_matrix <- taxa_level(mitochondria, which_level = "family")
write.csv(x = otu_table(mitochondria_matrix), "mitochondria_matrix.csv")

## Export ASVs table and taxonomy table as .csv
asv <- otu_table(ps_16s_bacteria_wo_mitochondria)
tax <- tax_table(ps_16s_bacteria_wo_mitochondria)

write.csv(as.data.frame(asv), "ASVs_16s_terry.csv")
write.csv(as.data.frame(tax), "Taxonomy_16s_terry.csv")

## Convert phyloseq object as family level
phyloseq_family <- tax_glom(ps_16s_bacteria_wo_mitochondria, taxrank = "family", NArm= FALSE)
asv_family <- otu_table(phyloseq_family)
tax_family <- tax_table(phyloseq_family)

write.csv(as.data.frame(asv_family), "Family_asvs_16s_terry.csv")
write.csv(as.data.frame(tax_family), "Family_taxonomy_16s_terry.csv")
phyloseq_family %>%
  refseq() %>%
  Biostrings::writeXStringSet("family.fna", append=FALSE,
                              compress=FALSE, compression_level=NA, format="fasta")

# Compositional anaylsis of overall microbiome
### 1. stacked bar plot
### convert phyloseq object into designated level of taxonomy (family)
library(phyloseq)
rank_names(ps_16s_bacteria_wo_mitochondria) # check the taxonomy structure
stacked_bar_plot_function <- function(phyloseq_object, taxrank = "family", NArm = TRUE){
  phyloseq <- tax_glom(phyloseq_object, taxrank = "family", NArm= TRUE)
  phyloseq_rel <- transform_sample_counts(phyloseq, function(x) x/sum(x)) # convert data into relative abundance
  phyloseq_melt <- psmelt(phyloseq_rel) # melt data
  phyloseq_melt[,taxrank] <- as.character(phyloseq_melt[,taxrank])
  phyloseq_melt$original_family <- phyloseq_melt$family
  phyloseq_melt[,taxrank][phyloseq_melt$Abundance < 0.1] <- "Other"
  write.csv(phyloseq_melt, "psmelt.csv")
  # Check number of families 
  Count = length(unique(phyloseq_melt[,taxrank]))
  if (Count > 25) {
    print("Warning: You have more than 25 taxa to plot, consider using higher cut off")
    print(Count)
  } else {
    print(Count)
  } 
  
  # Getting color code for plotting - randomly assigned color
  library(RColorBrewer)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category =='qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col=sample(col_vector, Count)
  
  #plot stacked barplot
  library(ggplot2)
  plot <- ggplot(data=phyloseq_melt, aes(x=Sample, y=Abundance, fill = family)) + 
    geom_bar(stat="identity", position="stack") + 
    facet_grid(.~Host + Collection_site + Type_of_sample, scales= "free", space = "free") +
    guides(fill=guide_legend(nrow=5)) + 
    scale_fill_manual(values=col)+
    theme(legend.position = "bottom", axis.text.x = element_text(angle=90), 
          axis.ticks.x = element_blank(), 
          panel.background = element_blank(),
          legend.title = element_blank()) +
    xlab("Samples") + 
    ylab("Relative abundance") 
  plot
  
  return(plot)
}


### Stacked bar plot WITH NA in family level
plot_without_numbers <- stacked_bar_plot_function(ps_16s_bacteria_wo_mitochondria, 
                          taxrank = "family", 
                          NArm=TRUE)

plot_without_numbers

ggsave("stacked_bar_plot_without_NA.jpg", plot=plot_without_numbers, width = 15, height= 5, units = "in", dpi = 600)
ggsave("stacked_bar_plot_without_NA.svg", plot=plot_without_numbers, width = 12.5, height= 5, units = "in", dpi = 600, device="svg")
ggsave("stacked_bar_plot_without_NA.pdf", plot=plot_without_numbers, width = 12.5, height= 5, units = "in", dpi = 600, device="pdf")

## Diversity analysis
## 1. Alpha Diversity
library(ggpubr)
library(FSA)
library(reshape2)
library(stringr)
library(rcompanion)
install_github("cran/multcompView")
library(multcompView)
alpha_diversity <- estimate_richness(ps_16s_bacteria_wo_mitochondria, measures= c("Observed", "Shannon", "InvSimpson")) # calculate alpha diversity of selected measures
alpha_comparison <- cbind(alpha_diversity, sample_data(ps_16s_bacteria_wo_mitochondria)) # add metadata to the alpha diversity output
alpha_comparison$plot_x <- NULL
alpha_comparison$plot_x <- paste0(alpha_comparison[,"Host"], "\n",
                                  alpha_comparison[,"Collection_site"], "\n",
                                  alpha_comparison[,"Type_of_sample"])

melt_plot <- melt(alpha_comparison) # melt data for plotting

# Seperate data into different indices
melt_plot_observed <- subset(melt_plot, variable == "Observed")
melt_plot_shannon <- subset(melt_plot, variable == "Shannon")
melt_plot_InvSimpson <- subset(melt_plot, variable == "InvSimpson")

# making function (kruskal test, dunn post-hoc test, and plotting with letter)

alpha_diversity_plot_with_letter <- function(melt_dataset, vjust, hjust, index){
  dunn <- dunnTest(value ~ plot_x, data=melt_dataset, method= "hochberg")
  Dunnx <- dunn$res
  cld <- cldList(P.adj ~ Comparison,
                 data = Dunnx,
                 threshold = 0.05)
  data_for_plotting <- group_by(melt_dataset, plot_x) %>%
    summarise(mean=mean(value), quant=quantile(value, probs=0.75))
  data_for_plotting$cld <- toupper(cld$Letter)
  plot <- ggplot(data = melt_dataset, aes(y = value, x = plot_x, fill = plot_x)) +
    geom_boxplot() +
    theme_minimal() + 
    scale_fill_brewer(palette="Dark2") +
    theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank()) +
    geom_text(data= data_for_plotting, aes(x = plot_x, y = quant, label = cld), size=3, vjust=vjust, hjust=hjust) +
    theme(legend.position = "none") +
    labs(x="Samples", y=index)
  return(plot)
}


dunn_observed <- kruskal.test(value ~ Host, data=melt_plot_observed)
Dunnx <- dunn_observed$res
cld <- cldList(P.adj ~ Comparison,
               data = Dunnx,
               threshold = 0.05)
data_for_plotting <- group_by(melt_plot_observed, plot_x) %>%
  summarise(mean=mean(value), quant=quantile(value, probs=0.75))
data_for_plotting$cld <- toupper(cld$Letter)

write.table(dunn_observed$res, "dunn_observed.csv")

dunn_shannon <- dunnTest(value ~ plot_x, data=melt_plot_shannon, method = "hochberg")
write.csv(dunn_shannon$res, "dunn_shannon.csv")


dunn_simpson <- dunnTest(value ~ plot_x, data= melt_plot_InvSimpson, method = "hochberg")
write.csv(dunn_simpson$res, "dunn_simpson.csv")






observed <- alpha_diversity_plot_with_letter(melt_plot_observed, vjust =-2, hjust= -0.5, "Observed")
ggsave("observed.jpg",plot = observed ,width=7.5, height=5, device="jpg", dpi=600)
ggsave("observed.svg",plot = observed ,width=7.5, height=5, device="svg", dpi=600)
ggsave("observed.pdf",plot = observed ,width=7.5, height=5, device="pdf", dpi=600)

shannon <- alpha_diversity_plot_with_letter(melt_plot_shannon, vjust =-2, hjust= -0.5, "Shannon")
ggsave("shannon.jpg",plot = shannon ,width=7.5, height=5, device="jpg", dpi=600)
ggsave("shannon.svg",plot = shannon ,width=7.5, height=5, device="svg", dpi=600)
ggsave("shannon.pdf",plot = shannon ,width=7.5, height=5, device="pdf", dpi=600)

Simpson <- alpha_diversity_plot_with_letter(melt_plot_InvSimpson, vjust =-2, hjust= -0.5, "InvSimpson")
ggsave("Simpson.jpg",plot = Simpson ,width=7.5, height=5, device="jpg", dpi=600)
ggsave("Simpson.svg",plot = Simpson ,width=7.5, height=5, device="svg", dpi=600)
ggsave("Simpson.pdf",plot = Simpson ,width=7.5, height=5, device="pdf", dpi=600)
library(cowplot)
grid_plot_alpha <- plot_grid(observed, shannon, Simpson, ncol = 3, labels =c("A", "B", "C"),
          rel_widths = 10, rel_heights = 5)

ggsave("alpha_plot.jpg",plot = grid_plot_alpha ,width=20, height=5, device="jpg", dpi=600)
ggsave("alpha_plot.svg",plot = grid_plot_alpha ,width=20, height=5, device="svg", dpi=600)
ggsave("alpha_plot.pdf",plot = grid_plot_alpha ,width=20, height=5, device="pdf", dpi=600)


## 2. Beta diversity
ps_16s_bacteria_wo_mitochondria
# Data normalization (CLR - Central log ratio transformation)
# https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
#Following step requires samples on rows and ASVs in columns
otus <- otu_table(ps_16s_bacteria_wo_mitochondria)
taxa_are_rows(otus) # Should be "FALSE" if "TRUE" use this command to transpose matrix "otus <- t(otus)"
#Replace zero values before clr transformation
#Use CZM method to replace zeros and outputs pseudo-counts (1)
require(zCompositions)
otu.n0 <- t(cmultRepl(otus, label =0, method="CZM", output="p-counts"))
otu.n0 <-ifelse(otu.n0 < 0, otu.n0*(-1), otu.n0)
#Convert data to proportions
otu.n0_prop <- apply(otu.n0, 2, function(x) {x/sum(x)})
#CLR transformation
otu.n0.clr<-t(apply(otu.n0_prop, 2, function(x){log(x)-mean(log(x))}))
phyloesq_genus_CLR <- phyloseq(otu_table(otu.n0.clr, taxa_are_rows=F), 
                               tax_table(ps_16s_bacteria_wo_mitochondria), 
                               sample_data(ps_16s_bacteria_wo_mitochondria))
#PCA
pc.clr <- prcomp(otu.n0.clr)
require(compositions)
# Calculate total variance (necessary for calculating %variance explained)
mvar.clr <- mvar(otu.n0.clr)
row <- rownames(otu.n0.clr)
# extract first two PCs
pc_out <- as.data.frame(pc.clr$x[,1:2])
# combine first two PCs with metadata
pc_out_meta <- as.data.frame(cbind(pc_out, sample_data(phyloesq_genus_CLR)))
library(vegan)
# Calculate euclidean distance between sample
dist <- vegdist(otu.n0.clr, method = "euclidean")

# PERMANOVA (Permutational based ANOVA) 
# https://onlinelibrary.wiley.com/doi/full/10.1002/9781118445112.stat07841
permanova <- vegan::adonis2(dist~Type_of_sample * Host * Collection_site, data = pc_out_meta, perm=999)
permanova
permanova_plus <- vegan::adonis2(dist~Type_of_sample + Host + Collection_site, data= pc_out_meta, perm=999)
permanova_plus
permanova_plus_int <- vegan::adonis2(dist~Type_of_sample + Host + Collection_site + Type_of_sample:Host + Host:Collection_site + Type_of_sample:Collection_site + Type_of_sample:Host:Collection_site, data= pc_out_meta, perm = 999)
permanova_plus_int
permanova # Check PERMANOVA results
label = paste0("p-value = ",permanova$`Pr(>F)`[1]," (PERMANOVA)") #prepare label for plot
label # Check p-value

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pc_out_meta$plot_x <- paste0(pc_out_meta[,"Host"], "\n",
                             pc_out_meta[,"Type_of_sample"])
pairwise.adonis2(dist~ plot_x, data= pc_out_meta, perm=999)

#Getting colors based on your comparisons
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category =='qual',]
col_vector <- unlist(mapply(brewer.pal, 8, "Dark2"))
col=sample(col_vector, length(unique(pc_out_meta[,"plot_x"])))

#Making PCA plot
PCA_plot <- ggplot(pc_out_meta, aes(x=PC1, y=PC2, color = plot_x, shape= Collection_site)) +
  geom_point(size=3)+
  theme(legend.position = 'right')+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+ # edit backgroud
  scale_x_continuous(name = paste("PC1: ", round(pc.clr$sdev[1]^2/mvar.clr*100, digits=1), "%", sep="")) + # %variance explained for PC1
  scale_y_continuous(name = paste("PC2: ", round(pc.clr$sdev[2]^2/mvar.clr*100, digits=1), "%", sep="")) + # %variance explained for PC2
  scale_color_brewer(palette="Dark2") +
  theme(legend.title = element_blank())
PCA_plot

ggsave("pca.jpg", plot = PCA_plot, width= 10, height=7.5, units="in", dpi=600)
ggsave("pca.svg", plot = PCA_plot, width= 10, height=7.5, units="in", dpi=600, device="svg")
ggsave("pca.pdf", plot = PCA_plot, width= 10, height=7.5, units="in", dpi=600, device="pdf")

## ALDEx2
## Running ALDEx2 as ASVs, and family level
aldex2_analysis <- function(phyloseq, family = FALSE){
  if(family == FALSE){
    phyloseq <- phyloseq
  } else {
    phyloseq <- tax_glom(phyloseq, taxrank = "family", NArm = FALSE)
  }
  otus <- otu_table(phyloseq)
  library(ALDEx2)
  meta <- as.data.frame(sample_data(phyloseq))
  aldex_16s <- aldex.clr(t(otus), mc.samples= 128, 
                         conds = meta$Type_of_sample)
  aldex_16s_e <- aldex.effect(aldex_16s)
  aldex_16s_t <- aldex.ttest(aldex_16s)
  aldex_16s_all <- data.frame(aldex_16s_e, aldex_16s_t)
  aldex_16s_sig <-which(aldex_16s_all$wi.eBH <=0.05)
  aldex_16s_sig.row <- rownames(aldex_16s_all)[which(aldex_16s_all$wi.eBH <= 0.05)]
  
  asv_16s_melt <- phyloseq %>% 
    transform_sample_counts(function(x) {x * 100} ) %>% 
    #Recalculates the relative abundance 
    psmelt() %>%  #Melts to long format
    arrange(desc(Abundance))
  
  asv_16s_filter_sig <- filter(asv_16s_melt, OTU %in% aldex_16s_sig.row)
  stat_16s <- asv_16s_filter_sig %>%
    group_by(OTU, Type_of_sample, family) %>%
    summarise(Mean=mean(Abundance, na.rm=TRUE))
  stat_16s_mean <- dcast(stat_16s, OTU + family ~ Type_of_sample, 
                         value.var = "Mean")
  stat_16s_mean$log2change<-
    log2(stat_16s_mean$'Pseudoflower'/stat_16s_mean$'Inflorescence')
  stat_16s_mean$Type_of_sample<-ifelse(stat_16s_mean$log2change <0 , 
                                       "Pseudoflower", "Inflorescence")
  
  # Check number of families 
  Count = length(unique(stat_16s_mean$family))
  
  # Getting color code for plotting - randomly assigned color

  Logfold<-ggplot(stat_16s_mean, aes(x=log2change, y=reorder(OTU,log2change), 
                                     fill=family))+
    geom_bar(stat='identity', color='black')+
    theme(axis.title.y = element_blank(), 
          axis.ticks=element_line(color='black'),
          axis.text=element_text(size=10, color='black')) + 
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, 
                               keyheight = 1, ncol=1)) +
    xlab("Log fold change (log2 Negative/Positive)") + 
    theme(panel.background = element_rect(fill=NA, color =NA),
          plot.background = element_rect(fill="white", color =NA),
          panel.border = element_rect(color="black", fill=NA))  +
    scale_fill_brewer(palette = "Dark2")  
  results <- list(stat_16s_mean, Logfold)
  return(results)
}

aldex_16s_asvs <- aldex2_analysis(ps_16s_bacteria_wo_mitochondria)

ggsave("DA.jpg", plot = Logfold, width= 5, height=5, units="in", dpi=600, device="jpg")
ggsave("DA.svg", plot = Logfold, width= 5, height=5, units="in", dpi=600, device="svg")
ggsave("DA.pdf", plot = Logfold, width= 5, height=5, units="in", dpi=600, device="pdf")
write.csv(stat_16s_mean, "DA_results_significant.csv")  

#Core microbiota analysis (ASVs, genus, and family level)
#https://microbiome.github.io/tutorials/Core.html

library(microbiome)
prevalance_phyloseq <- function(phyloseq, taxrank){
  if(taxrank == "ASVs"){
    phyloseq_prev <- phyloseq
  } else {
  phyloseq_prev <- tax_glom(physeq = phyloseq, taxrank = taxrank, NArm=TRUE)
  }
  phyloseq_rel <- microbiome::transform(phyloseq_prev, "compositional")
  phyloseq_rel_influ <- subset_samples(physeq = phyloseq_rel, Type_of_sample == "Inflorescence")
  phyloseq_rel_pseudo <- subset_samples(physeq = phyloseq_rel, Type_of_sample == "Pseudoflower")
  
  if(taxrank == "ASVs"){
    phyloseq_rel_influ_plot <- phyloseq_rel_influ
    phyloseq_rel_pseudo_plot <- phyloseq_rel_pseudo
  } 
  if(taxrank == "genus"){
    phyloseq_rel_influ_plot <- phyloseq_rel_influ
    taxa_names(phyloseq_rel_influ_plot) <- tax_table(phyloseq_rel_influ_plot)[,"genus"]
    phyloseq_rel_pseudo_plot <- phyloseq_rel_pseudo
    taxa_names(phyloseq_rel_pseudo_plot) <- tax_table(phyloseq_rel_pseudo)[,"genus"]
  }
  if(taxrank == "family"){
    phyloseq_rel_influ_plot <- phyloseq_rel_influ
    taxa_names(phyloseq_rel_influ_plot) <- tax_table(phyloseq_rel_influ_plot)[,"family"]
    phyloseq_rel_pseudo_plot <- phyloseq_rel_pseudo
    taxa_names(phyloseq_rel_pseudo_plot) <- tax_table(phyloseq_rel_pseudo_plot)[,"family"]
  }
  
  ##plot core microbiome
  library(RColorBrewer)
  prevalences <- seq(.05, 1, .05)
  detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)
  p1 <- plot_core(phyloseq_rel_influ_plot,
                  plot.type="heatmap",
                  colours = rev(brewer.pal(5, "RdBu")),
                  prevalences = prevalences,
                  detections = detections, min.prevalence = 0.5) +
    xlab("Detection Threshold (Relative Abundance (%))") + theme_bw() + ylab("ASVs") +
    ggtitle("Inflorescence")
  
  p2 <- plot_core(phyloseq_rel_pseudo_plot,
                  plot.type="heatmap",
                  colours = rev(brewer.pal(5, "RdBu")),
                  prevalences = prevalences,
                  detections = detections, min.prevalence = 0.5) +
    xlab("Detection Threshold (Relative Abundance (%))") + theme_bw() + ylab("ASVs") +
    ggtitle("Pseudoflower")
  
  influ_prev <- as.data.frame(prevalence(phyloseq_rel_influ, sort = FALSE, detection = 0))
  pseudo_prev <- as.data.frame(prevalence(phyloseq_rel_pseudo, sort = FALSE, detection = 0))
  tax <- as.data.frame(tax_table(phyloseq_prev))
  tax$influ <- influ_prev[,1]
  tax$pseudo <- pseudo_prev[,1]
  
  results <- list(tax, p1, p2)
  return(results)
}

asvs_prev <- prevalance_phyloseq(ps_16s_bacteria_wo_mitochondria, "ASVs")
write.csv(asvs_prev, "prev_ASVs.csv")
ggsave("Core_ASVs_influ.jpg", plot=asvs_prev[[2]], width=5, height=5, units="in", dpi=600, device="jpg")
ggsave("Core_ASVs_influ.pdf", plot=asvs_prev[[2]], width=5, height=5, units="in", dpi=600, device="pdf")
ggsave("Core_ASVs_influ.svg", plot=asvs_prev[[2]], width=5, height=5, units="in", dpi=600, device="svg")

ggsave("Core_ASVs_psuedo.jpg", plot=asvs_prev[[3]], width=5, height=5, units="in", dpi=600, device="jpg")
ggsave("Core_ASVs_psuedo.pdf", plot=asvs_prev[[3]], width=5, height=5, units="in", dpi=600, device="pdf")
ggsave("Core_ASVs_psuedo.svg", plot=asvs_prev[[3]], width=5, height=5, units="in", dpi=600, device="svg")

genus_prev <- prevalance_phyloseq(ps_16s_bacteria_wo_mitochondria, "genus")
write.csv(genus_prev, "prev_genus.csv")
ggsave("Core_genus_influ.jpg", plot=genus_prev[[2]], width=7, height=7, units="in", dpi=600, device="jpg")
ggsave("Core_genus_influ.pdf", plot=genus_prev[[2]], width=7, height=7, units="in", dpi=600, device="pdf")
ggsave("Core_genus_influ.svg", plot=genus_prev[[2]], width=7, height=7, units="in", dpi=600, device="svg")

ggsave("Core_genus_psuedo.jpg", plot=genus_prev[[3]], width=5, height=5, units="in", dpi=600, device="jpg")
ggsave("Core_genus_psuedo.pdf", plot=genus_prev[[3]], width=5, height=5, units="in", dpi=600, device="pdf")
ggsave("Core_genus_psuedo.svg", plot=genus_prev[[3]], width=5, height=5, units="in", dpi=600, device="svg")


family_prev <- prevalance_phyloseq(ps_16s_bacteria_wo_mitochondria, "family")
write.csv(family_prev, "prev_family.csv")
ggsave("Core_family_influ.jpg", plot=family_prev[[2]], width=5, height=5, units="in", dpi=600, device="jpg")
ggsave("Core_family_influ.pdf", plot=family_prev[[2]], width=5, height=5, units="in", dpi=600, device="pdf")
ggsave("Core_family_influ.svg", plot=family_prev[[2]], width=5, height=5, units="in", dpi=600, device="svg")

ggsave("Core_family_psuedo.jpg", plot=family_prev[[3]], width=5, height=5, units="in", dpi=600, device="jpg")
ggsave("Core_family_psuedo.pdf", plot=family_prev[[3]], width=5, height=5, units="in", dpi=600, device="pdf")
ggsave("Core_family_psuedo.svg", plot=family_prev[[3]], width=5, height=5, units="in", dpi=600, device="svg")


## Core microbiome
install.packages("eulerr")
BiocManager::install("microbiome")

library(eulerr)
library(microbiome)
library(microbiomeutilities)

rownames(alpha_comparison) <- sample_names(ps_16s_bacteria_wo_mitochondria)
alpha_comparison_phy <- sampleData(alpha_comparison)
ps_16s_bacteria_wo_mitochondria_2 <- merge_phyloseq(otu_table(ps_16s_bacteria_wo_mitochondria),
                                                    tax_table(ps_16s_bacteria_wo_mitochondria),
                                                    alpha_comparison_phy)


table(meta(ps_16s_bacteria_wo_mitochondria_2)$plot_x, useNA = "always")
pseq.rel <- microbiome::transform(ps_16s_bacteria_wo_mitochondria_2, "compositional")
pseq.rel_Infloresence <- subset_samples(pseq.rel, 
                                        Type_of_sample == "Inflorescence")

venn_diagram <- function(pseq.rel){
  disease_states <- unique(as.character(meta(pseq.rel_Infloresence)$plot_x))
  list_core <- c() # an empty object to store information
  for (n in disease_states){ # for each variable n in DiseaseState
    #print(paste0("Identifying Core Taxa for ", n))
    
    ps.sub <- subset_samples(pseq.rel_Infloresence, plot_x == n) # Choose sample from DiseaseState by n
    
    core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                           detection = 0, # 0.001 in atleast 90% samples 
                           prevalence = 0)
    print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
    list_core[[n]] <- core_m # add to a list core taxa for each group.
    #print(list_core)
  }
  mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
  venn <- plot(venn(list_core),fills = mycols)
  print(venn)
}

ggsave("Venn_Inflorescence.jpg", venn, dpi=600)
ggsave("Venn_Inflorescence.svg", venn, dpi=600, device="svg")
ggsave("Venn_Inflorescence.pdf", venn, dpi=600, device="pdf")


ggsave("Venn_surinamensis.jpg", venn, dpi=600)
ggsave("Venn_surinamensis.svg", venn, dpi=600, device="svg")
ggsave("Venn_surinamensis.pdf", venn, dpi=600, device="pdf")


ggsave("Venn_all.jpg", venn, dpi=600)
ggsave("Venn_all.svg", venn, dpi=600, device="svg")
ggsave("Venn_all.pdf", venn, dpi=600, device="pdf")

a <- list_core[[1]]
b <- list_core[[2]]
c <- list_core[[3]]

a <- as.data.frame(a)
rownames(a) <- a[,1]
b <- as.data.frame(b)
rownames(b) <- b[,1]
c <- as.data.frame(c)
rownames(c) <- c[,1]

merged <- merge(a,b,by = "row.names", all=TRUE)
rownames(merged) <- merged[,1]
merged[,1] <- NULL
merged <- merge(merged,c, by ="row.names", all=TRUE)
rownames(merged) <- merged[,1]
merged[,1] <- NULL
tax_name <- as.data.frame(tax_table(pseq.rel_Infloresence))
merged <- merge(merged, tax_name, by="row.names", all=TRUE)
write.csv(merged, "Inflorescence_shared.csv")


