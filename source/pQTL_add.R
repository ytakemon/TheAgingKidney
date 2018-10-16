# R/3.4.1
pQTL_best <- read.csv("~/Dropbox/TheAgingKidneyData/QTLscan/protein/BestMarker_add_protein_nobatch.csv")

# Usinb pQTL_best to create plot
# need to reorder "chr" factors
chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
pQTL_best$chr <- factor(pQTL_best$chr, levels = chr_full)
pQTL_best$AdditiveChr <- factor(pQTL_best$AdditiveChr, levels= chr_full)

# Subset out chr 1-19,X from data
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
pQTL_best <- pQTL_best[pQTL_best$chr %in% chr, ]


# check distribution of LOD scores
hist(pQTL_best$AdditiveLOD, breaks =100)
# Set LOD threshold
LODthreshold_diff <- 7.5

# Plot Interactive-Age pQTLs
# Subset Int-Age LOD above total/full and diff
addScan <- pQTL_best[(pQTL_best$AdditiveLOD >= LODthreshold_diff),] # above diff threshold


# Annotate Interactive-Age postion with genes and save file for sharing
save_addScan <- addScan[,c("id", "gene_id","symbol","chr","start","end", "biotype", "AdditiveChr","AdditivePos","AdditiveLOD")]
save_addScan <- arrange(save_addScan, AdditiveChr, AdditivePos)
# save annotated list for sharing
# write.csv(save_addScan, "./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", row.names = FALSE, quote = FALSE)

# Convert transcript and qtl position relative to chromosome positions
# Convert to megabases
chrlen <- sapply(split(addScan$end, addScan$chr), max) * 1e-6
chrsum <- c(0, cumsum(chrlen))
names(chrsum) = names(chrlen)

t.gmb <- addScan$start * 1e-6 # Transcript
q.gmb <- addScan$AdditivePos * 1e-6 # qtl

# Cumulative sum of previosu positions
for(i in c(2:19, "X")) {
  
  wh <- which(addScan$chr == i)
  t.gmb[wh] <- t.gmb[wh] + chrsum[i]
  
  wh <- which(addScan$AdditiveChr == i)
  q.gmb[wh] <- q.gmb[wh] + chrsum[i]
}
addScan$t_gbm <- t.gmb
addScan$q_gbm <- q.gmb

# Custom lablels & lines
# Only display chr1:19,X
chrtick <- chrsum[1:20]
# Shift thick to half way point
max <- max(addScan$q_gbm)
chrtick_half <- NULL
for (i in 1:length(chrtick)){
  if (i == 20){
    
    x <- (chrtick[i] + max)/2
    chrtick_half <- c(chrtick_half, x)
  } else {
    
    x <- (chrtick[i] + chrtick[i + 1])/2
    chrtick_half <- c(chrtick_half, x)
  }
}
#to adjust eQTL plot a bit to the right to match density
chrtick_halfy <- chrtick_half
names(chrtick_halfy)[20] <- "X  "

pPlot <- ggplot(addScan, aes(x= q_gbm, y= t_gbm), color = AdditiveLOD) +
  geom_point(alpha = 0.8) +
  scale_x_continuous("QTL position",
                     breaks = chrtick_half,
                     limits = c(min(addScan$q_gbm), max(addScan$q_gbm)),
                     expand = c(0,0)) +
  scale_y_continuous("Gene position",
                     breaks = chrtick_half,
                     limits = c(min(addScan$t_gbm), max(addScan$t_gbm)),
                     expand = c(0,0)) +
  geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
  geom_hline(yintercept = chrtick[2:20], colour = "grey", size = 0.2) +
  labs( title = "Additive scan pQTLs") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", size = 0.2, fill = NA))+
  scale_colour_gradient2(low = "blue", high = "red", mid = "blue", midpoint = mean(addScan$AdditiveLOD))

interval <- seq(0,max(addScan$q_gbm), by = 10)
interval <- interval + 5
avgLOD <- NULL
for ( i in interval){
  df <- addScan %>% filter((q_gbm > (i-5)) & (q_gbm < (i+5)))
  if(nrow(df) == 0){
    avgLOD <- c(avgLOD,0)
  } else{
    avgLOD <- c(avgLOD, mean(df$AdditiveLOD))
  }
}
df <- data.frame(interval = interval, avgLOD = avgLOD)
z <- c(0,0)
df <- rbind(z, df)

pdensity <- ggplot(addScan, aes(q_gbm, colour = "grey", fill = "grey")) +
  geom_histogram(breaks = seq(0,max(addScan$q_gbm), by = 10)) +
  scale_colour_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
  scale_fill_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
  scale_x_continuous("QTL position",
                     breaks = chrtick_half,
                     limits = c(min(addScan$q_gbm), max(addScan$q_gbm)),
                     expand = c(0,0)) +
  scale_y_continuous(name ="Density",
                     breaks = seq(0,100, by = 10)) +
  geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.2, fill = NA))