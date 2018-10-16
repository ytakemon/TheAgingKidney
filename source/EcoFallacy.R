library(ggsci)

# Create new fuction to gather gene information
get.gene <- function(gene.symbol){
  
  # pull covariates and convert to factors
  # note Age is numeric and fAge is factor
  Mouse.ID <- rownames(annot.samples)
  Sex <- factor(as.character(model.matrix(~annot.samples[,"Sex"])[,2]+1), levels=c("1","2"), labels=c("F","M"))
  Age <- annot.samples$Age
  fAge <- factor(as.character(Age),levels=c("6","12","18"))
  
  #pull gene identifiers from mrna and protein annotations
  mrna.id <- annot.mrna[annot.mrna$symbol==gene.symbol,"id"]
  protein.id <- annot.protein[annot.protein$symbol==gene.symbol,"id"]
  
  # strat building the tiblle to hold data
  my.data <- data_frame(Mouse.ID=Mouse.ID,
                        Sex=Sex,
                        Age=Age,
                        fAge=fAge)
  # add RNA and Protein columns to the tibble - if it exists
  if(length(mrna.id)>0){
    my.data <- mutate(my.data, RNA=expr.mrna[,mrna.id])
  }
  if(length(protein.id)>0){
    my.data <- mutate(my.data, Protein=expr.protein[,protein.id])
  }
  # return value is a tibble with covariates, RNA and Protein data
  my.data
}

# Slc5a12
df <- get.gene("Slc5a12") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
            mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Slc5a12 <- ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Slc5a12 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "SLC5A12 expression (Rank normalized protein)",
        x = "Slc5a12 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()

# Flot1
df <- get.gene("Flot1") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
            mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Flot1 <- ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Flot1 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "FLOT1 expression (Rank normalized protein)",
        x = "Flot1 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()