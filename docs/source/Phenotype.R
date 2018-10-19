Upheno <- read.delim(paste0(wd,"Phenotype/JAC_CS_urine_chem_v1.txt"), sep = "\t")

pheno <- Upheno %>%
  mutate(
    Mouse.ID = mouse.id,
    duplicated = (duplicated(Upheno$mouse.id) | duplicated(Upheno$mouse.id, fromLast = TRUE)),
    Alb = ma.u,
    Phs = phs.u) %>%
  filter((Mouse.ID %in% annot.samples$Mouse.ID) & duplicated == FALSE) %>%
  select(Mouse.ID, Alb, Phs)


annot.samples <- annot.samples[annot.samples$Mouse.ID %in% pheno$Mouse.ID,]
expr.mrna <- expr.mrna[rownames(expr.mrna) %in% pheno$Mouse.ID,]
expr.protein <- expr.protein[rownames(expr.protein) %in% pheno$Mouse.ID,]

# Identify gene name
Slc34a1_m <- annot.mrna[annot.mrna$symbol == "Slc34a1",]$id
Slc34a1_p <- annot.protein[annot.protein$symbol == "Slc34a1",]$id
Slc34a3_m <- annot.mrna[annot.mrna$symbol == "Slc34a3",]$id
Slc34a3_p <- annot.protein[annot.protein$symbol == "Slc34a3",]$id
Lrp2_m <- annot.mrna[annot.mrna$symbol == "Lrp2",]$id
Lrp2_p <- annot.protein[annot.protein$symbol == "Lrp2",]$id

df <- annot.samples %>%
  mutate(
    Slc34a1_mRNA = expr.mrna[, Slc34a1_m],
    Slc34a1_protein = expr.protein[, Slc34a1_p],
    Slc34a3_mRNA = expr.mrna[, Slc34a3_m],
    Slc34a3_protein = expr.protein[, Slc34a3_p],
    Lrp2_mRNA = expr.mrna[, Lrp2_m],
    Lrp2_protein = expr.protein[, Lrp2_p],
    Sex = as.factor(Sex),
    Generation = factor(Generation, levels = c("G8","G9","G10","G11","G12")),
    Age = factor(as.character(Age), levels = c("6","12","18")),
    Alb = pheno$Alb,
    Phs = pheno$Phs) %>%
  select(Mouse.ID, Sex, Generation, Age, Slc34a1_mRNA, Slc34a1_protein, Slc34a3_mRNA, Slc34a3_protein, Lrp2_mRNA, Lrp2_protein, Alb, Phs)

#Slc34a1 protein ~ RNA
df_r1c3 <- df %>%
  mutate(
    Fitted = fitted(lm(Slc34a1_protein ~ Sex+Age+Slc34a1_mRNA+Sex:Age, data=.)))

df_r1c3_summary <- df_r1c3 %>%
  group_by(Sex, Age) %>%
  summarize(mRNA=mean(Slc34a1_mRNA), sdRNA=sd(Slc34a1_mRNA),
            mProtein=mean(Slc34a1_protein), sdProtein=sd(Slc34a1_protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Slca1_mp <- ggplot(df_r1c3_summary, aes(x=mRNA, y=mProtein, color=Age, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df_r1c3, aes(x=Slc34a1_mRNA, y=Slc34a1_protein), alpha = 0.6) +
  geom_line(data=df_r1c3, aes(x = Slc34a1_mRNA, y = Fitted)) +
  theme_bw() +
  labs(y = "Slc34a1 protein",
       x = "Slc34a1 mRNA") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black")) +
  theme_bw() +
  scale_colour_aaas()

#Slc34a3 protein ~ RNA
df_r2c3 <- df %>% filter(!is.na(df$Slc34a3_protein)) %>%
  mutate( Fitted = fitted(lm(Slc34a3_protein ~ Sex+Age+Slc34a3_mRNA+Sex:Age, data=.)))

df_r2c3_summary <- df_r2c3 %>%
  group_by(Sex, Age) %>%
  summarize(mRNA=mean(Slc34a3_mRNA), sdRNA=sd(Slc34a3_mRNA),
            mProtein=mean(Slc34a3_protein), sdProtein=sd(Slc34a3_protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Slca3_mp <- ggplot(df_r2c3_summary, aes(x=mRNA, y=mProtein, color=Age, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df_r2c3, aes(x=Slc34a3_mRNA, y=Slc34a3_protein), alpha = 0.6) +
  geom_line(data=df_r2c3, aes(x = Slc34a3_mRNA, y = Fitted)) +
  theme_bw() +
  labs(y = "Slc34a3 protein",
       x = "Slc34a3 mRNA") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black")) +
  theme_bw() +
  scale_colour_aaas()

#Lrp2 protein ~ RNA
df_r3c3 <- df %>% filter(!is.na(df$Lrp2_protein)) %>%
  mutate( Fitted = fitted(lm(Lrp2_protein ~ Sex+Age+Lrp2_mRNA+Sex:Age, data=.)))

df_r3c3_summary <- df_r3c3 %>%
  group_by(Sex, Age) %>%
  summarize(mRNA=mean(Lrp2_mRNA), sdRNA=sd(Lrp2_mRNA),
            mProtein=mean(Lrp2_protein), sdProtein=sd(Lrp2_protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Lrp2_mp <- ggplot(df_r3c3_summary, aes(x=mRNA, y=mProtein, color=Age, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df_r3c3, aes(x=Lrp2_mRNA, y=Lrp2_protein), alpha = 0.6) +
  geom_line(data=df_r3c3, aes(x = Lrp2_mRNA, y = Fitted)) +
  theme_bw() +
  labs(y = "Lrp2 protein",
       x = "Lrp2 mRNA") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black")) +
  theme_bw() +
  scale_colour_aaas()


#Slc34a1 Phos ~ protein
df_r1c4 <- df %>% filter(!is.na(Phs)) %>%
  mutate(
    Fitted = fitted(lm(Phs ~ Sex+Age+Slc34a1_protein+Sex:Age, data=.)))

df_r1c4_summary <- df_r1c4 %>%
  group_by(Sex, Age) %>%
  summarize(mRNA=mean(Slc34a1_protein), sdRNA=sd(Slc34a1_protein),
            mProtein=mean(Phs), sdProtein=sd(Phs), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Slca1_pPhos <- ggplot(df_r1c4_summary, aes(x=mRNA, y=mProtein, color=Age, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df_r1c4, aes(x=Slc34a1_protein, y=Phs), alpha = 0.6) +
  geom_line(data=df_r1c4, aes(x = Slc34a1_protein, y = Fitted)) +
  theme_bw() +
  labs(y = "Phosphate",
       x = "Slc34a1 protein") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black")) +
  coord_cartesian(ylim = c(0,800))+
  theme_bw() +
  scale_colour_aaas()

#Slc34a1 Phos ~ protein
df_r2c4 <- df %>% filter(!is.na(Phs) & !is.na(Slc34a3_protein)) %>%
  mutate(
    Fitted = fitted(lm(Phs ~ Sex+Age+Slc34a3_protein+Sex:Age, data=.)))

df_r2c4_summary <- df_r2c4 %>%
  group_by(Sex, Age) %>%
  summarize(mRNA=mean(Slc34a3_protein), sdRNA=sd(Slc34a3_protein),
            mProtein=mean(Phs), sdProtein=sd(Phs), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Slca3_pPhos <- ggplot(df_r2c4_summary, aes(x=mRNA, y=mProtein, color=Age, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df_r2c4, aes(x=Slc34a3_protein, y=Phs), alpha = 0.6) +
  geom_line(data=df_r2c4, aes(x = Slc34a3_protein, y = Fitted)) +
  theme_bw() +
  labs(y = "Phosphate",
       x = "Slc34a3 protein") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black")) +
  coord_cartesian(ylim = c(0,800))+
  theme_bw() +
  scale_colour_aaas()

#Slc34a1 Phos ~ protein
df_r3c4 <- df %>% filter(!is.na(Alb) & !is.na(Lrp2_protein)) %>%
  mutate(
    Fitted = fitted(lm(Alb ~ Sex+Age+Lrp2_protein+Sex:Age, data=.)))

df_r3c4_summary <- df_r3c4 %>%
  group_by(Sex, Age) %>%
  summarize(mRNA=mean(Lrp2_protein), sdRNA=sd(Lrp2_protein),
            mProtein=mean(Alb), sdProtein=sd(Alb), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Lrp2_pAlb <- ggplot(df_r3c4_summary, aes(x=mRNA, y=mProtein, color=Age, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df_r3c4, aes(x=Lrp2_protein, y=Alb), alpha = 0.6) +
  geom_line(data=df_r3c4, aes(x = Lrp2_protein, y = Fitted)) +
  theme_bw() +
  labs(y = "Albumin",
       x = "Lrp2 protein") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black")) +
  #  coord_cartesian(ylim = c(0,20))+
  theme_bw() +
  scale_colour_aaas()

## rest of plot
Slca1_m <- ggplot(df, aes(Age, y = Slc34a1_mRNA, colour = Sex)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
  geom_point(position = position_jitterdodge(0.4)) +
  theme_bw() +
  labs(y = "Slc34a1 mRNA",
       x = "Age") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_colour_aaas()

Slca1_p <- ggplot(df, aes(Age, y = Slc34a1_protein, colour = Sex)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
  geom_point(position = position_jitterdodge(0.4)) +
  theme_bw() +
  labs(y = "Slc34a1 protein",
       x = "Age") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_colour_aaas()

Slca3_m <- ggplot(df, aes(Age, y = Slc34a3_mRNA, colour = Sex)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
  geom_point(position = position_jitterdodge(0.4)) +
  theme_bw() +
  labs(y = "Slc34a3 mRNA",
       x = "Age") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_colour_aaas()

Slca3_p <- ggplot(df, aes(Age, y = Slc34a3_protein, colour = Sex)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
  geom_point(position = position_jitterdodge(0.4)) +
  theme_bw() +
  labs(y = "Slc34a3 protein",
       x = "Age") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_colour_aaas()

Lrp2_m <- ggplot(df, aes(Age, y = Lrp2_mRNA, colour = Sex)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
  geom_point(position = position_jitterdodge(0.4)) +
  theme_bw() +
  labs(y = "Lrp2 mRNA",
       x = "Age") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_colour_aaas()

Lrp2_p <- ggplot(df, aes(Age, y = Lrp2_protein, colour = Sex)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
  geom_point(position = position_jitterdodge(0.4)) +
  theme_bw() +
  labs(y = "Lrp2 protein",
       x = "Age") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_colour_aaas()

