# Analisi de dades omiques. PAC 1
# Marina Gay
# Exploratory data analysis of a phosphoproteomics experiment


# Libraries and folders --------------------------------------------------------
# Load libraries
library(tidyverse)
library(ggplot2)
library(xlsx)
library(SummarizedExperiment)
library(mixOmics)
library(ggpubr)
library(plotly)
library(htmlwidgets)

# Getting the dir where this R script is located: mainDir
p_mainDir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Setting mainDir as the current working dir
setwd(p_mainDir)

# Define folders
p_mainDir <-  paste0(p_mainDir,"/")
p_figures <-  paste0(p_mainDir,"figures/")
p_tables <-  paste0(p_mainDir,"tables/")

# Create figures and table folders
dir.create(file.path(p_figures), showWarnings=F)
dir.create(file.path(p_tables), showWarnings=F)

# Load data and prepare for the SummarizedExperiment ---------------------------
# Load data expression and feature data
df.raw <- read.xlsx(file=paste0(p_mainDir,"metaboData/Datasets/2018-Phosphoproteomics/TIO2+PTYR-human-MSS+MSIvsPD.XLSX"), sheetIndex=1)
# df.raw <- read.delim(file=paste0(p_mainDir,"phosphopep.txt"), sep = "\t", stringsAsFactors = FALSE)

# There is a phosphopeptide entry duplicate (GEPNVSYICSR[7] Phospho|[9] Carbamidomethyl):
# One is Y and the other is S/T
# I will add this information to distinguish them
df.raw$SequenceModifications <- with(df.raw, ifelse(SequenceModifications == "GEPNVSYICSR[7] Phospho|[9] Carbamidomethyl",
                                                     paste0(SequenceModifications, "_", PHOSPHO),
                                                     SequenceModifications))

# Select only expression data
df.expr <- df.raw %>% 
  dplyr::select(SequenceModifications, ends_with("_MSS"), ends_with("_PD"))

# Convert to matrix
m.expr <- as.matrix(df.expr %>% dplyr::select(-SequenceModifications))

# Add phosphopeptides as rownames
rownames(m.expr) <- df.expr$SequenceModifications

# Select feature data and add column with number of phospho in the sequence
df.fd <- df.raw %>% 
  dplyr::select(-ends_with("_MSS"), -ends_with("_PD")) %>% 
  mutate(phos = str_count(SequenceModifications, "Phospho") # Count number of phospho
  )

# Load pheno data
df.pheno <- read.xlsx(file=paste0(p_mainDir, "TIO2+PTYR-human-MSS+MSIvsPD.XLSX"), sheetIndex=2)
# df.pheno <- read.delim(file=paste0(p_mainDir, "targets.txt"), sep = "\t", stringsAsFactors = TRUE)

# load metadata
meta <- readLines(paste0(p_mainDir, "description.md"))

# Create the summarized experiment container -----------------------------------
se <- SummarizedExperiment(assays = list(counts = m.expr),
                           colData = df.pheno,
                           rowData = df.fd,
                           metadata = meta)

# Data pre processing ----------------------------------------------------------
# Remove phosphopeptides with no quantification values in any sample
# Check for rows with non-zero values in any sample
non.zero.pep <- rowSums(assay(se, "counts") != 0) > 0

# Subset the SummarizedExperiment object to keep only non-zero rows
se <- se[non.zero.pep, ]

# Replace 0 with NA
assays(se)$counts[assays(se)$counts == 0] <- NA

# log 2 transform 
assays(se)$counts <- log2(assays(se)$counts)

# Tidy data for ggplot
df.t <- df.raw %>% 
  pivot_longer(cols = ends_with("_MSS") | ends_with("_PD"),
               names_to = "replicate",
               values_to = "intensity") %>% 
  mutate_at(c('intensity'), ~na_if(., 0)) %>% # replace 0 with NA
  mutate(log2Int = log2(intensity), # log transform
         phos = str_count(SequenceModifications, "Phospho"), # Count number of phospho
         Phenotype =  str_split_fixed(replicate, "_", 3), # Add pheno data
         sample = Phenotype[,1],
         Individual = Phenotype[,2],
         Phenotype = Phenotype[,3],
         IsNAN = is.na(log2Int)) # Find NA

# Average data by phenotype and calculate FC
df.av <- df.t %>% 
  group_by(SequenceModifications, Accession, Description, Phenotype) %>% 
  summarise(abundance = mean(log2Int)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Phenotype,
              values_from = abundance) %>% 
  mutate(Abu_mean = rowMeans(dplyr::select(., MSS, PD), na.rm = TRUE),
         log2FC = MSS - PD)

# Prepare data form missigness plot
df.miss <- df.t %>%
  mutate(replicate = factor(replicate)) %>%  
  group_by(replicate, Phenotype, IsNAN) %>%  
  summarise(n = n(), .groups = "drop") %>%  
  group_by(replicate, Phenotype) %>%  
  mutate(Percentage = n / sum(n),
         Sample = str_split_fixed(replicate, "_", 3),
         Sample = paste0(Sample[,1], "_", Sample[,2]))
  

# Add technical variables to pheno data
# Filter missing data
df.miss.filt <- df.miss %>% 
  filter(IsNAN == TRUE) %>% 
  rename("miss_count" = "n",
         "miss_per" = "Percentage")

# Extract total intensity
df.totInt <- df.t %>% 
  drop_na(log2Int) %>% 
  group_by(replicate, Phenotype) %>% 
  summarize(totalInt = sum(intensity, na.rm = TRUE),
            phos_count = n()) %>% 
  mutate(Sample = str_split_fixed(replicate, "_", 3),
         Sample = paste0(Sample[,1], "_", Sample[,2]))

# Join data frames
df.pheno <- merge(df.pheno, df.miss.filt, by = c("Sample", "Phenotype"),
                  all.x = TRUE, sort = FALSE)

df.pheno <- merge(df.pheno, df.totInt, by = c("Sample", "Phenotype", "replicate"),
                  all.x = TRUE, sort = FALSE)

# Update pheno data
colData(se)$miss_count <- df.pheno$miss_count
colData(se)$totalInt <- df.pheno$totalInt

# Exploratory analysis ---------------------------------------------------------
df.raw$PHOSPHO <- as.factor(df.raw$PHOSPHO)
summary(df.raw)

# number of phosphopeptides identified
length(unique(df.raw$SequenceModifications))

# number of proteins identified
length(unique(df.raw$Accession))

# Data visualization -----------------------------------------------------------
# Define phenotype color
col <- c("PD" = "darkblue",
         "MSS" = "orange")

# Boxplot of log2 Intensities
plot.boxInt <- ggplot(data = df.t,
       aes(x = replicate, y = log2Int, color = Phenotype)) +
  theme_bw() +
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=90, vjust= 0.3, hjust= 1)) +
  scale_color_manual(values = col)
plot.boxInt

# Number of phosphopeptides identified by sample
plot.barNum <- ggplot(data = df.totInt,
                      aes(x = replicate, y = phos_count, fill = Phenotype)) +
  theme_bw() +
  geom_col() +
  theme(axis.text.x  = element_text(angle=90, vjust= 0.3, hjust= 1)) +
  scale_fill_manual(values = col)
plot.barNum

# Scatter plot of phenotype
plot.pt_groups <- ggplot(data = df.av,
                  aes(x = PD, y = MSS)) +
  theme_bw() +
  geom_point() +
  stat_cor(method="pearson") +
  geom_abline(slope = 1, intercept = 0, linetype = 2)
plot.pt_groups

# MA plot
plot.MA <- ggplot(data = df.av,
               aes(x = Abu_mean, y = log2FC,
                   label = paste(Accession, SequenceModifications))) +
  theme_bw() +
  geom_point() +
  ylab("log2(MSS/PD)") +
  geom_line(y = 0, linetype = 2)
plot.MA

# Barplot with the percentage of missing by sample
plot.miss <- ggplot(df.miss,
                  aes(x=n, y=replicate, alpha=IsNAN, fill=Phenotype, label=scales::percent(Percentage))) +
  theme_bw() +
  geom_col() +
  geom_text(position = position_stack(0.5), size = 2) +
  scale_alpha_manual(values=c(1.0, 0.5)) +
  scale_fill_manual(values = col)
plot.miss

# Extract expression matrix form the se object
m.explog <- assay(se)

# remove NA
m.explog <- m.explog[complete.cases(m.explog), ]

# Computing PCA
pca <- prcomp(t(m.explog))
pcvara <- pca$sdev^2/sum(pca$sdev^2)*100;

# Convert to dataframe
df.pca <- as.data.frame(pca$x)
df.pca$sample <- rownames(df.pca)

plot.pca <- ggplot(data=df.pca,
                      aes(x=PC1, y=PC2, color=df.pheno$Phenotype)) +
  geom_point()+
  geom_text(label=df.pheno$Sample, show.legend = FALSE, vjust = -1) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pcvara[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(pcvara[2], 2), "%)")) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = col)
plot.pca

# PC1 correlation with technical variables
plot(df.pca$PC1, log2(df.pheno$totalInt))
cor(df.pca$PC1, log2(df.pheno$totalInt))

plot(df.pca$PC1, df.pheno$miss_count)
cor(df.pca$PC1, df.pheno$miss_count)


# heatmap sample correlation
manDist <- dist(t(m.explog))
heatmap (as.matrix(manDist), col=heat.colors(16))

# Dendogram
clust.euclid.average <- hclust(dist(t(m.explog)),method="average")
plot(clust.euclid.average, hang=-1)

# MixOmics plots -------
MyResult.pca <- pca(t(m.explog))  # 1 Run the method
plotIndiv(MyResult.pca) # 2 Plot the samples
plotVar(MyResult.pca)   # 3 Plot the variables

tune.pca.multi <- tune.pca(t(m.explog), ncomp = 10, scale = TRUE)
plot(tune.pca.multi)

final.pca.multi <- pca(t(m.explog), ncomp = 3, center = TRUE, scale = TRUE)

# Top variables on the first component only:
head(selectVar(final.pca.multi, comp = 1)$value)

biplot(final.pca.multi)

# Save plots and data ----------------------------------------------------------
# Save plots in png
# List all objects in the environment that start with "plot."
plot_names <- ls(pattern = "^plot\\.")

# Initialize an empty list to store the plots
plot_list <- list()

# Populate the list with plots 
for (plot_name in plot_names) {
  plot_list[[plot_name]] <- get(plot_name)
}

# Save the plots 
lapply(names(plot_list), function(f) {
  ggsave(paste0(p_figures, f, '4x4.png'),
         plot_list[[f]], width=4, height=4, dpi=1200)
})

lapply(names(plot_list), function(f) {
  ggsave(paste0(p_figures, f, '2x4.png'),
         plot_list[[f]], width=4, height=2, dpi=1200)
})

png(filename = paste0(p_figures, "pca_comp.png"), width = 240, height = 240, units = "px")
plot(tune.pca.multi)
dev.off()

png(filename = paste0(p_figures, "heatmap.png"), width = 300, height = 300, units = "px")
heatmap (as.matrix(manDist), col=heat.colors(16))
dev.off()

png(filename = paste0(p_figures, "dendogram.png"), width = 240, height = 240, units = "px")
plot(clust.euclid.average, hang=-1)
dev.off()

png(filename = paste0(p_figures, "PCA_miss.png"), width = 240, height = 240, units = "px")
plot(df.pca$PC1, df.pheno$miss_count)
dev.off()

saveWidget(ggplotly(plot.MA), file = paste0(p_figures, "MAplot.html"))

# Save Rda object
save(se, df.av, plot_list, file=paste0(p_mainDir, "Phosphopeptide_SE", '.Rda'))

