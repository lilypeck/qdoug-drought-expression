library(tidyverse)
library(nlme)

setwd("/Users/lilypeck/Library/CloudStorage/Box-Box/SorkLab/Writing/Steele_etal/github/Files")
## identify shared differentially expressed PFAMs in Gugger et al., 2017 and Mead et al., 2019

## pull in DE files 
gugger.de <- read_delim("Gugger-TableS2.csv")
mead.de <- read_delim("Mead-DE_genes_significant_all_tests.csv")
# Q. douglasii gene annotations
annots <- read_delim("../Output-data/Qdoug_transcriptome_annotations.tsv")
## drought t2 vs t1 all genes
all.genes <- read_delim("../Output-data/salmon.isoform.counts.matrix.QUDO_Drought_T2_vs_QUDO_Drought_T1.DESeq2.DE_results")
# raw expression data from well-watered treatment
ww.genes <- read_delim("../Output-data/salmon.isoform.counts.matrix.QUDO_Control_T2_vs_QUDO_Control_T1.DESeq2.DE_results")
# Q. lobata protein family annotations - downloaded from https://valleyoak.ucla.edu/genomic-resources/
qlob.pfam <- read_delim(col_names = F, "Qlobata.v3.0.PCG.InterProScan5.34-73.0.tsv")

########################
# Figure S1
########################
## tidy and identify pfams DE in both studies
## select all pfams, not just first in list
gugger.de2 <- gugger.de %>%
  dplyr::select(Contig, `TAIR or Pfam ID`) %>%
  rename(Pfam.gugger = `TAIR or Pfam ID`) %>%
  drop_na() %>%
  filter(!grepl("AT",Pfam.gugger)) %>%
  separate_wider_delim(Pfam.gugger, names_sep = "", delim = " ", too_few = "debug") %>%
  dplyr::select(Contig:Pfam.gugger15) %>%
  pivot_longer(cols = -c(Contig), names_to = "nams", values_to = "vals") %>%
  drop_na() %>%
  mutate(Pfam.gugger = gsub(".*:","",vals)) %>%
  dplyr::select(Contig, Pfam.gugger)

## filtered by pval < 0.01 to make more stringent
mead.de2 <- mead.de %>%
  filter(adj.P.Val < 0.01) %>%
  separate_wider_delim(cols = pfam_id, delim = ", ", names_sep = "",too_few = "debug") %>%
  dplyr::select(gene,pfam_id1:pfam_id5) %>%
  pivot_longer(cols = -c(gene), names_to = "nams", values_to = "vals") %>%
  dplyr::select(!nams) %>%
  filter(!vals == "none") %>%
  drop_na() %>%
  mutate(vals = gsub(".*:","",vals),
         vals = gsub("\\..*$","",vals)) %>%
  rename(pfam = vals) %>%
  mutate(both = pfam %in% gugger.de2$Pfam.gugger) %>%
  filter(both == T)
# get Q. lobata pfam names as list
pfam.n <- mead.de2 %>%
  dplyr::select(pfam) %>%
  distinct()
## now check expression of these same pfams in Q. douglasii
annots2 <- annots %>%
  separate_wider_delim(cols = Pfam, delim = "`", names_sep = "",too_few = "debug") %>%
  dplyr::select(transcript_id,Pfam1:Pfam10) %>%
  pivot_longer(cols = Pfam1:Pfam10, names_to = "nams",values_to = "vals") %>%
  drop_na() %>%
  filter(!vals == ".") %>%
  mutate(vals = gsub("\\^.*$","",vals),
         vals = gsub("\\..*$","",vals),
         temp = vals %in% pfam.n$pfam) %>%
  rename(pfam = vals, gene = transcript_id)
## get pfams which were DEG in both mead and gugger
annots.deg <- annots2 %>%
  filter(temp == TRUE) %>%
  dplyr::select(gene,pfam)

## filter by DEG pfam and tidy data for plotting
# drought treatment
dr.genes2 <- all.genes %>%
  left_join(annots.deg) %>%
  mutate(temp = pfam %in% pfam.n$pfam) %>%
  filter(temp == T) %>% 
  pivot_longer(cols = baseMeanA:baseMeanB, names_to = "nams", values_to = "vals") %>%
  group_by(pfam,nams) %>%
  summarise(Expression = mean(vals),gene.n = n()) %>%
  rename(Pfam = pfam, Treatment = nams) %>%
  mutate(Treatment = gsub("baseMeanA","Drought_Time 2",Treatment),
         Treatment = gsub("baseMeanB","Drought_Time 1",Treatment)) %>%
  mutate(Pop = "ONEA",
         Species= "Qdoug")
# well-watered treatment
ww.genes2 <- ww.genes %>%
  left_join(annots.deg) %>%
  mutate(temp = pfam %in% pfam.n$pfam) %>%
  filter(temp == T) %>% 
  pivot_longer(cols = baseMeanA:baseMeanB, names_to = "nams", values_to = "vals") %>%
  group_by(pfam,nams) %>%
  summarise(Expression = mean(vals),gene.n = n()) %>%
  rename(Pfam = pfam, Treatment = nams) %>%
  mutate(Treatment = gsub("baseMeanA","Well-watered_Time 2",Treatment),
         Treatment = gsub("baseMeanB","Well-watered_Time 1",Treatment)) %>%
  mutate(Pop = "ONEA",
         Species= "Qdoug")
## plot constitutive expression for 81 pfams in qdoug - figure S1
qdoug.both <- bind_rows(dr.genes2,ww.genes2)
pal <- c("#ffad73","#26b3ff")
qdoug.both %>%
  separate_wider_delim(Treatment, names = c("Treatment","Time"),delim = "_") %>%
  ggplot(aes(x = Time, y = sqrt(Expression), group = Treatment, color = Treatment)) +
  scale_fill_manual(values= pal) +
  scale_color_manual(values = pal) +
  geom_smooth(aes(color = Treatment, fill = Treatment)) +
  #geom_line(size = 1, alpha = 0.5) +
  facet_wrap(~factor(Treatment, levels = c("Well-watered","Drought"))) +
  theme_minimal() +
  coord_cartesian(ylim = c(0,20)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  ylab("sqrt(Expression)")

ggsave(bg = "white", 
       height = 3, width = 5, units = "in",
       "plots/FigureS1.svg")

########################
# Tables S1-S2 linear models to test effect of family
########################
# up-regulated
drought.deseq2 <- drought.deseq %>%
  filter(padj < 0.01 & log2FoldChange > 0) %>%
  dplyr::select(!sampleA:padj & !edger:ww) %>%
  pivot_longer(cols = QUDO_10_Control_T1_1:QUDO_9_Drought_T2_8,
               names_to = "nams",values_to = "vals") %>%
  separate_wider_delim(cols = nams, delim = "_", names = c("q","site","Treatment","time","indiv")) %>%
  relocate(indiv, .before = Treatment) %>%
  rename(Expression = vals) %>%
  unite(q:site, col = "Family", sep = "_", remove = F) %>%
  unite(site:indiv, col = "Indiv", sep = "_") %>%
  dplyr::select(!q)

# Fit the model - include family in model . expression ~ treatment * time * family, random = ~ 1 | gene
lm.results <- anova(lme(Expression ~ Treatment * time * Family, random = ~ 1 | gene, data = drought.deseq2,
                        control = lmeControl(opt = "optim", msMaxIter = 200, msMaxEval = 200)))

# Tidy the output
results <- lm.results %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(pval = "p-value",
         fval = "F-value") %>%
  mutate(pval.adj = p.adjust(pval, method='BH'),
         key = case_when(pval.adj < 0.05 ~ "*", .default = "ns"))
results
write_delim(results, "plots/TableS1.tsv", delim = "\t")

# then down-regulated
drought.deseq2 <- drought.deseq %>%
  filter(padj < 0.01 & log2FoldChange < 0) %>%
  dplyr::select(!sampleA:padj & !edger:ww) %>%
  pivot_longer(cols = QUDO_10_Control_T1_1:QUDO_9_Drought_T2_8,
               names_to = "nams",values_to = "vals") %>%
  separate_wider_delim(cols = nams, delim = "_", names = c("q","site","Treatment","time","indiv")) %>%
  relocate(indiv, .before = Treatment) %>%
  rename(Expression = vals) %>%
  unite(q:site, col = "Family", sep = "_", remove = F) %>%
  unite(site:indiv, col = "Indiv", sep = "_") %>%
  dplyr::select(!q)

# Fit the model
lm.results <- anova(lme(Expression ~ Treatment * time * Family, random = ~ 1 | gene, data = drought.deseq2,
                        control = lmeControl(opt = "optim", msMaxIter = 200, msMaxEval = 200)))

# Tidy the output
results <- lm.results %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(pval = "p-value",
         fval = "F-value") %>%
  mutate(pval.adj = p.adjust(pval, method='BH'),
         key = case_when(pval.adj < 0.05 ~ "*", .default = "ns"))
results
write_delim(results, "plots/TableS2A.tsv", delim = "\t")

## post-hoc
# select T2
t2 <- drought.deseq2 %>% 
  filter(time == "T2", Treatment == "Drought") %>%
  mutate(Family = as.factor(Family))
# Fit the model - family
model <- aov(Expression ~ Family + (1+gene), data = t2)
# Run Tukey  
values <- summary(glht(model, linfct = mcp(Family = "Tukey")), test = adjusted(type = "bonferroni"))
results <- tidy(values)
write_delim(results, "plots/TableS2B.tsv", delim = "\t")

########################
#Table S3
########################
# table of drought-responsive pfams
# tidy pfams
qlob.pfam2 <- qlob.pfam %>%
  filter(X4 == "Pfam") %>%
  dplyr::select(X1,X5,X6) %>%
  distinct() %>%
  rename(gene = X1, Pfam = X5, Pfam.name = X6)
# make table
tab2 <- annots.deg %>% 
  dplyr::select(pfam) %>% 
  rename(Pfam = pfam) %>%
  left_join(qlob.pfam2) %>%
  dplyr::select(!gene) %>%
  distinct()

# then statistical tests 
## test slope drought
dr.genes.slope <- all.genes %>%
  left_join(annots.deg) %>%
  mutate(temp = pfam %in% pfam.n$pfam) %>%
  filter(temp == T) %>% 
  dplyr::select(gene,baseMeanA,baseMeanB,pfam) %>%
  pivot_longer(cols = baseMeanA:baseMeanB, names_to = "nams", values_to = "vals") %>%
  rename(Pfam = pfam, Time = nams) %>%
  mutate(Time = gsub("baseMeanA","Drought",Time),
         Time = gsub("baseMeanB","Control",Time)) %>%
  mutate(Pop = "ONEA",
         Species= "Qdoug")
## run model by site
## get slope by pfam
mod <- lm(scale(vals) ~ Time * Pfam, data = dr.genes.slope)
tidy(mod)
# Compute the effect (difference) of time within each Pfam
slopes <- emmeans(mod, pairwise ~ Time | Pfam)
tab.dr <- as_tibble(slopes$contrasts)

## test slope well-watered
ww.genes.slope <- ww.genes %>%
  left_join(annots.deg) %>%
  mutate(temp = pfam %in% pfam.n$pfam) %>%
  filter(temp == T) %>% 
  dplyr::select(gene,baseMeanA,baseMeanB,pfam) %>%
  pivot_longer(cols = baseMeanA:baseMeanB, names_to = "nams", values_to = "vals") %>%
  rename(Pfam = pfam, Time = nams) %>%
  mutate(Time = gsub("baseMeanA","Drought",Time),
         Time = gsub("baseMeanB","Control",Time)) %>%
  mutate(Pop = "ONEA",
         Species= "Qdoug")
## run model
mod <- lm(scale(vals) ~ Time * Pfam, data = ww.genes.slope)
# Compute the effect (difference) of time within each Pfam
slopes <- emmeans(mod, pairwise ~ Time | Pfam)
tab.ww <- as_tibble(slopes$contrasts)

## then combine tabs for supp table 3
tab.dr2 <- tab.dr %>%
  dplyr::select(Pfam, estimate) %>%
  rename(Drought = estimate)

tab.ww2 <- tab.ww %>%
  dplyr::select(Pfam, estimate) %>%
  rename(Well_watered = estimate) 

tab2.final <- tab2 %>%
  left_join(tab.dr2) %>%
  left_join(tab.ww2) %>%
  arrange(Pfam)

write_delim(delim = "\t", tab2.final, "plots/TableS3.tsv")

########################
## table s5 - test difference between drought T2 and ww T2
########################
dr.genes.slope2 <- dr.genes.slope %>% mutate(Treatment = "drought")
ww.genes.slope2 <- ww.genes.slope %>% mutate(Treatment = "ww")
both <- bind_rows(dr.genes.slope2, ww.genes.slope2)
mod <- lme(vals ~ Treatment * Time, random = ~ 1 | Pfam, data = both)
summary(mod)
tab.mod <- anova(mod)
write_delim(delim = "\t", tab.mod, "plots/TableS5.tsv")


