# Load library and data
load(paste0(wd,"RData/DO188b_kidney_noprobs.RData"))
data <- read_csv(paste0(wd,"ANOVA/kidney_anova_slope_output.csv"))

# Subset to columns with slopes and pvalues
mcols <- grep("^m.", colnames(data))
sigcols <- grep("^p.", colnames(data))
data <- data[, c(1:8, mcols, sigcols)]

# Using only age realated mRNA and Protein that are significant < 0.05
df_age <- data[data$p.mRNA_Age.Sex <= 0.05,]
df_age <- df_age[df_age$p.Prot_Age.Sex <= 0.05,]

# Count number of genes in each quadrant
# mRNA-Protein
total <- nrow(df_age)
quadI <- nrow(df_age[((df_age$m.mRNA_Age.Sex > 0) & (df_age$m.Prot_Age.Sex > 0)),])
quadIV <- nrow(df_age[((df_age$m.mRNA_Age.Sex > 0) & (df_age$m.Prot_Age.Sex < 0)),])
quadII <- nrow(df_age[((df_age$m.mRNA_Age.Sex < 0) & (df_age$m.Prot_Age.Sex > 0)),])
quadIII <- nrow(df_age[((df_age$m.mRNA_Age.Sex < 0) & (df_age$m.Prot_Age.Sex < 0)),])

# Plot age realted mRNA/Protein slopes
SlopePlot <- ggplot() +
  geom_point(data = data, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour =  "#64676d", alpha = 0.3) + # set colour to light grey
  geom_point(data = df_age, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#66125c") + # set colour of significant genes to eggplant
  scale_x_continuous("mRNA-Age.Sex",
                     breaks = seq(-1, 1, by = 0.1)) +
  scale_y_continuous("Protein-Age.Sex",
                     breaks = seq(-1, 1, by = 0.1)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  labs( title = "Direction of change in gene expression with Age",
        subtitle = paste0("Slope: RNA v. Protein, Total = ", total, ", Sig = ", nrow(df_age))) +
  annotate("text", x = 0.25, y = 0.85, label = paste("A = ", quadI)) +
  annotate("text", x = 0.25, y = -0.85, label = paste("C = ", quadIV)) +
  annotate("text", x = -0.25, y = 0.85, label = paste("D = ", quadII)) +
  annotate("text", x = -0.25, y = -0.85, label = paste("B = ", quadIII)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect( colour = "black", fill = NA))