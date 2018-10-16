# R/3.4.1
library(dplyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE) # args <- "kidney_anova_slope_output.csv"
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("RNAseq_data/DO188b_kidney_noprobs.RData")
output <- list.files(path = "./Anova_output/", pattern = paste0("^",args[[1]]), recursive = TRUE)
data <- read.csv(paste0("./Anova_output/",output[[1]]), header = T)

# For only slopes and sigcols
mcols <- grep("^m.", colnames(data))
sigcols <- grep("^p.", colnames(data))
data <- data[, c(1:8, mcols, sigcols)]

# Age (not interacting with sex)------------------------------------------------
# Using only age realated mRNA and Protein that are significant < 0.05
df_age <- data[data$p.mRNA_Age.Sex <= 0.05,]
df_age <- df_age[df_age$p.Prot_Age.Sex <= 0.05,]
# N = 693

# Count number of genes in each quadrant
# mRNA-Protein
total <- nrow(df_age)
quadI <- nrow(df_age[((df_age$m.mRNA_Age.Sex > 0) & (df_age$m.Prot_Age.Sex > 0)),])
quadIV <- nrow(df_age[((df_age$m.mRNA_Age.Sex > 0) & (df_age$m.Prot_Age.Sex < 0)),])
quadII <- nrow(df_age[((df_age$m.mRNA_Age.Sex < 0) & (df_age$m.Prot_Age.Sex > 0)),])
quadIII <- nrow(df_age[((df_age$m.mRNA_Age.Sex < 0) & (df_age$m.Prot_Age.Sex < 0)),])
# Export gene list in each quadrant
Gene_list <- df_age
Gene_list$quadI <- ((Gene_list$m.mRNA_Age.Sex > 0) & (Gene_list$m.Prot_Age.Sex > 0))
Gene_list$quadII <- ((df_age$m.mRNA_Age.Sex < 0) & (df_age$m.Prot_Age.Sex > 0))
Gene_list$quadIII <- ((df_age$m.mRNA_Age.Sex < 0) & (df_age$m.Prot_Age.Sex < 0))
Gene_list$quadIV <- ((df_age$m.mRNA_Age.Sex > 0) & (df_age$m.Prot_Age.Sex < 0))

Age <- grep("_Age.Sex$", colnames(Gene_list))
Quad <- grep("^quad", colnames(Gene_list))
Gene_list <- Gene_list[, c(1:8, Age, Quad)]
write.csv(Gene_list, "./Anova_output/gene_lists/Quad_Age_Sig.csv", row.names = FALSE, quote = FALSE)

# Plot age realted mRNA/Protein slopes
pdf("./Plot/slope_mRNA_Prot_Age.pdf", width = 6, height =6)
ggplot(df_age, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex)) +
      geom_point() +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.25)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.25)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age",
            subtitle = paste0("Slope: RNA v. Protein, N = ", total)) +
      annotate("text", x = 0.25, y = 0.85, label = paste("I, N = ", quadI)) +
      annotate("text", x = 0.25, y = -0.85, label = paste("IV, N = ", quadIV)) +
      annotate("text", x = -0.25, y = 0.85, label = paste("II, N = ", quadII)) +
      annotate("text", x = -0.25, y = -0.85, label = paste("III, N = ", quadIII)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Age interacting with Sex -----------------------------------------------------
# Using only age realated mRNA and Protein that are significant < 0.05
df_int <- data[data$p.mRNA_Interaction <= 0.05,]
df_int <- df_int[df_int$p.Prot_Interaction <= 0.05,]
# N = 168

# Count number of genes in each quadrant
# mRNA-Protein
total <- nrow(df_int)
quadI <- nrow(df_int[((df_int$m.mRNA_Interaction > 0) & (df_int$m.Prot_Interaction > 0)),])
quadIV <- nrow(df_int[((df_int$m.mRNA_Interaction > 0) & (df_int$m.Prot_Interaction < 0)),])
quadII <- nrow(df_int[((df_int$m.mRNA_Interaction < 0) & (df_int$m.Prot_Interaction > 0)),])
quadIII <- nrow(df_int[((df_int$m.mRNA_Interaction < 0) & (df_int$m.Prot_Interaction < 0)),])
# Export gene list in each quadrant
Gene_list <- df_int
Gene_list$quadI <- ((Gene_list$m.mRNA_Interaction > 0) & (Gene_list$m.Prot_Interaction > 0))
Gene_list$quadII <- ((df_int$m.mRNA_Interaction < 0) & (df_int$m.Prot_Interaction > 0))
Gene_list$quadIII <- ((df_int$m.mRNA_Interaction < 0) & (df_int$m.Prot_Interaction < 0))
Gene_list$quadIV <- ((df_int$m.mRNA_Interaction > 0) & (df_int$m.Prot_Interaction < 0))

Int <- grep("Interaction$", colnames(Gene_list))
Quad <- grep("^quad", colnames(Gene_list))
Gene_list <- Gene_list[, c(1:8, Int, Quad)]
write.csv(Gene_list, "./Anova_output/gene_lists/Quad_AgeIntSex_Sig.csv", row.names = FALSE, quote = FALSE)
# Plot age realted mRNA/Protein slopes
pdf("./Plot/slope_mRNA_Prot_AgeIntSex.pdf", width = 6, heigh =6)
ggplot(df_int, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction)) +
      geom_point() +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age*Sex",
            subtitle = paste0("Slope: RNA v. Protein, N = ", total)) +
      annotate("text", x = 0.06, y = 0.14, label = paste("I, N = ", quadI)) +
      annotate("text", x = 0.06, y = -0.14, label = paste("IV, N = ", quadIV)) +
      annotate("text", x = -0.06, y = 0.14, label = paste("II, N = ", quadII)) +
      annotate("text", x = -0.06, y = -0.14, label = paste("III, N = ", quadIII)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Plot total with subset overlay ----------------------------------------------
# Age
total <- nrow(data)
pdf("./Plot/Total_slope_mRNA_Prot_Age.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex)) +
      geom_point(data = df_age, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#ed0000") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age",
            subtitle = paste0("Slope: RNA v. Protein, Total = ", total, ", Sig = ", nrow(df_age))) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()
# interaction
pdf("./Plot/Total_slope_mRNA_Prot_AgeIntSex.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction)) +
      geom_point(data = df_int, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction), colour = "#ed0000") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age*Sex",
            subtitle = paste0("Slope: RNA v. Protein, Total = ", total, ", Sig = ", nrow(df_int))) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Plot total with subset overlay -----------------------------------------------
cols <- c("m.mRNA_Age.Sex", "m.Prot_Age.Sex", "m.mRNA_Interaction", "m.Prot_Interaction",
          "p.mRNA_Age.Sex", "p.Prot_Age.Sex", "p.mRNA_Interaction", "p.Prot_Interaction")
# Age
# Change in mRNA while protein stays constant
df_deltaAge <- data[data$p.mRNA_Age.Sex <= 0.05,]
df_deltaAge <- df_deltaAge[df_deltaAge$p.Prot_Age.Sex > 0.05,]
sub <- which(colnames(df_deltaAge) %in% cols)
df_deltaAge <- df_deltaAge[, c(1:8,sub)]
write.csv(df_deltaAge, file = "./Anova_output/gene_lists/ChangeRNA_ConstProt_Age.csv", row.names = FALSE, quote = FALSE)

pdf("./Plot/Total_slope_deltaRNA_Age.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex)) +
      geom_point(data = df_deltaAge, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#ed0000") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age",
            subtitle = paste0("Slope: RNA v. Protein (constant), Total = ", nrow(data), ", Sig = ", nrow(df_deltaAge))) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Change in protein while mRNA stays constant
df_deltaProt <- data[data$p.Prot_Age.Sex <= 0.05,]
df_deltaProt <- df_deltaProt[df_deltaProt$p.mRNA_Age.Sex > 0.05,]
df_deltaProt_inc <- data[data$p.Prot_Age.Sex <= 1e-10,]
df_deltaProt_inc <- df_deltaProt_inc[df_deltaProt_inc$p.mRNA_Age.Sex > 0.05,]

sub <- which(colnames(df_deltaProt_inc) %in% cols)
df_deltaProt_inc <- df_deltaProt_inc[, c(1:8,sub)]
write.csv(df_deltaProt_inc, file = "./Anova_output/gene_lists/ChangeProt_ConstRNA_Age.csv", row.names = FALSE, quote = FALSE)

pdf("./Plot/Total_slope_deltaProt_Age.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex)) +
      geom_point(data = df_deltaProt, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#ed0000") +
      geom_point(data = df_deltaProt_inc, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#ffdd00") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age",
            subtitle = paste0("Slope: RNA (constant) v. Protein, Total = ", nrow(data), ", Sig (0.05) = ", nrow(df_deltaProt), ", Sig (1e^-10) = ", nrow(df_deltaProt_inc))) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Interaction
# Change in mRNA while protein stays constant
df_deltaAge <- data[data$p.mRNA_Interaction <= 0.05,]
df_deltaAge <- df_deltaAge[df_deltaAge$p.Prot_Interaction > 0.05,]

sub <- which(colnames(df_deltaAge) %in% cols)
df_deltaAge <- df_deltaAge[, c(1:8,sub)]
write.csv(df_deltaAge, file = "./Anova_output/gene_lists/ChangeRNA_ConstProt_AgeIntSex.csv", row.names = FALSE, quote = FALSE)

pdf("./Plot/Total_slope_deltaRNA_AgeIntSex.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction)) +
      geom_point(data = df_deltaAge, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction), colour = "#ed0000") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age*Sex",
            subtitle = paste0("Slope: RNA v. Protein (constant), Total = ", nrow(data), ", Sig = ", nrow(df_deltaAge))) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Change in protein while mRNA stays constant
df_deltaProt <- data[data$p.Prot_Interaction <= 0.05,]
df_deltaProt <- df_deltaProt[df_deltaProt$p.mRNA_Interaction > 0.05,]
df_deltaProt_inc <- data[data$p.Prot_Interaction <=  0.01,]
df_deltaProt_inc <- df_deltaProt_inc[df_deltaProt_inc$p.mRNA_Interaction > 0.05,]

sub <- which(colnames(df_deltaProt_inc) %in% cols)
df_deltaProt_inc <- df_deltaProt_inc[, c(1:8,sub)]
write.csv(df_deltaProt_inc, file = "./Anova_output/gene_lists/ChangeProt_ConstAge_AgeIntSex.csv", row.names = FALSE, quote = FALSE)

pdf("./Plot/Total_slope_deltaProt_AgeIntSex.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction)) +
      geom_point(data = df_deltaProt, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction), colour = "#ed0000") +
      geom_point(data = df_deltaProt_inc, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction), colour = "#ffdd00") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age*Sex",
            subtitle = paste0("Slope: RNA (constant) v. Protein, Total = ", nrow(data), ", Sig (0.05) = ", nrow(df_deltaProt), ", Sig (0.01) = ", nrow(df_deltaProt_inc))) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()
