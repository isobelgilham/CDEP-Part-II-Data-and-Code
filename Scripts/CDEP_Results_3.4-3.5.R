##### 1. Permanova ####
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

# subset to library and make numeric 
ab_full<- sp_full[10:3038, 8:152]

for (i in 1:ncol(ab_full)){
  x<- as.numeric(ab_full[,i])
  ab_full[,i]<-x 
}

library(vegan)

lib_full<- sp_full[10:3038,]
lib_full[,8:ncol(sp_full)] <- lapply(lib_full[,8:ncol(sp_full)], as.numeric)
library(vegan)

# make community matrix 
merged_abundance<- lib_full[,8:ncol(lib_full)]
# turn abundances into matrix 
merged_matrix<- as.matrix(merged_abundance)
# run and plot nmds 
nmds= metaMDS(t(merged_matrix), distance= 'bray')
nmds
plot(nmds)


dist_matrix<- vegdist(t(merged_matrix), method = "bray")

# is this right? JF says yes 
adonis2(dist_matrix ~ sample_type, data = lib_full[,8:ncol(lib_full)], permutations = 999, strata = sites)
# way more significant when having random effects controlled p<0.01!! beautiful 
# need to figure out a way to visually present this 

##### 2. RDAs ####
library(vegan)
length(sample_type)
dim(t(dist_matrix))

RDA1<- capscale(t(dist_matrix)~sample_type, data = lib_full[,8:ncol(lib_full)]) 
# + sites
summary(RDA1)
RsquareAdj(RDA1)

RDA2<- capscale(t(dist_matrix)~sites+sites:sample_type, data = lib_full[,8:ncol(lib_full)])
summary(RDA2)
RsquareAdj(RDA2)

# can test 
anova.cca(RDA1)
anova.cca(RDA2, by='terms')
# axes instead of terms 

anova(RDA1, permutations = 999)
anova(RDA2, permutations=999)

colour_types<- as.factor(sample_type)
length(levels(colour_types))

col_funct<- c("blue", "green", "orange")

levels(colour_types)<- col_funct

group_labels <- c("Deadwood", "Living wood", "Soil")
group_colors <- c("blue", "green", "orange")


plot(RDA1$CCA$wa[,1], RDA1$CCA$wa[,2], pch= 21, bg=colour_types, cex=2 )

## resp
# Convert sample_type to a factor
colour_types <- as.factor(sample_type)

# Define color mapping
col_funct <- c("blue", "green", "orange") 

# Map the colors correctly
colour_map <- setNames(col_funct, levels(colour_types))  

# Apply the colors to the points
point_colors <- colour_map[colour_types]

# Define legend labels and colors
group_labels <- levels(colour_types)  
group_colors <- col_funct  

# Plot
svg("RDA1.svg", width = 8, height = 6)
plot(RDA1$CCA$wa[,1], RDA1$CCA$wa[,2], pch=21, bg=point_colors, cex=2, xlab= "RDA[,1] 1%", ylab="PCA[,1] 0.7%")

# Add legend
legend("topright", legend = c("Deadwood", "Living wood", "Soil"), pch = 21, pt.bg = group_colors, title = "Microhabitat")
dev.off()


point_colors <- adjustcolor(point_colors, alpha.f = 0.8)
################################################# Plot RDA1
svg("RDA1_new.svg", width = 8, height = 6)

# make bigger margins
par(mar=c(5,5,5,5))

plot(RDA1$CCA$wa[,1], RDA1$CCA$wa[,2], pch = 21, bg = point_colors, cex = 2,
     xlab = "RDA[,1] 0.951%", ylab = "RDA[,2] 0.687%", cex.lab=2, cex.axis=1.5)

# Add legend
legend("topright", legend = c("Deadwood", "Living wood", "Soil"),
       pch = 21, pt.bg = group_colors, title = "Microhabitat", cex=1.5)

dev.off()
###################################################
################################################### Plot RDA2 
svg("RDA2_new.svg", width = 8, height = 6)
par(mar = c(5, 6, 6, 10))
plot(RDA2$CCA$wa[,1], RDA2$CCA$wa[,2], pch = 21, bg = point_colors, cex = 2,
     xlab = "RDA[,1] 2.51%", ylab = "RDA[,2] 0.956%", cex.lab=2, cex.axis=1.5)

# Arrows for biplots
Arrows(x0 = 0, y0 = 0, x1 = (RDA2$CCA$biplot[1,1]) / 3, 
       y1 = (RDA2$CCA$biplot[1,2]) / 5, col = 'black')
Arrows(x0 = 0, y0 = 0, x1 = (RDA2$CCA$biplot[2,1]) / 4, 
       y1 = (RDA2$CCA$biplot[2,2]) / 5, col = 'black')

par(xpd = TRUE)
# Add legend
legend("topright", 
       inset = c(-0.4, 0), 
       legend = c("Deadwood", "Living wood", "Soil"),
       pch = 21, pt.bg = group_colors, title = "Microhabitat", cex=1.5)
par(xpd = FALSE)
dev.off()
################################################

##### 4. RDA Score code ####

# Run RDA
RDA2 <- capscale(t(dist_matrix) ~ sample_type + sites, data = lib_full[, 8:ncol(lib_full)])
summary(RDA2)
RsquareAdj(RDA2)

# Extract site scores
scores_df <- as.data.frame(RDA2$CCA$wa)
scores_df$sample_type <- as.factor(sites) 
colnames(scores_df)[1:2] <- c("RDA1", "RDA2")  
scores_df$sites <- sample_type  

# Define color mapping correctly
group_labels <- unique(sample_type)
group_colors <- setNames(c("orange", "blue", "green"), group_labels)

# Compute mean and confidence intervals

summary_stats <- scores_df %>%
  group_by(sample_type) %>%
  summarise(
    mean_x = mean(RDA1, na.rm = TRUE),
    mean_y = mean(RDA2, na.rm = TRUE),
    ci_x = qt(0.975, df = n() - 1) * sd(RDA1, na.rm = TRUE) / sqrt(n()),
    ci_y = qt(0.975, df = n() - 1) * sd(RDA2, na.rm = TRUE) / sqrt(n())
  )

##### 5. SB Lichen ####
# load in full data 
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[-c(1)]

# subset to libraries only
lib_full<- sp_full[10:nrow(sp_full),]

# set as numeric 
lib_full[,8:ncol(lib_full)]<- lapply(lib_full[, 8:ncol(lib_full)], as.numeric) 

# get relative abundances 
library(vegan)

RA_sp<- decostand(t(lib_full[,8:152]),method = 'total')

RA_sp_cl<- matrix(NA, nrow=nrow(RA_sp), ncol=length(unique(lib_full$O)))

Classes<-unique(lib_full$O)

i<- 1
#which(lib_full$C[DSL]%in% families)

for(i in 1:length(Classes)){
  x<- which(lib_full$O== Classes[i])
  if(length(x)>1) (y<- rowSums(RA_sp[,x]))  
  if(length(x)==1)(y<- RA_sp[,x])
  RA_sp_cl[,i]<- y 
  
}


rownames(RA_sp_cl)<- rownames(RA_sp)
colnames(RA_sp_cl)<- Classes 

groups<-c(rep("D", times=length(deadwood)),
          rep("S", times=length(soil)),
          rep("L", times=length(living)))

# sum by sample type
sample_type<- unlist(lapply(rownames(RA_sp_cl),function(x)(substr(x,start=3, stop=3))))

deadwood<- which(sample_type== 'D')
sum_d<- colMeans(RA_sp_cl[deadwood,])

soil<- which(sample_type== 'S')
sum_s<- colMeans(RA_sp_cl[soil,])

living<- which(sample_type== 'L')
sum_l<- colMeans(RA_sp_cl[living,])

RA_sp_cl_2<- cbind(sum_d, sum_s, sum_l)

top_100_fa_2<- cbind(sum_d, sum_s, sum_l)
# Identify the top 60% of abundance in each sample type
select_top_60 <- function(RA_sp) {
  sorted_abundance <- sort(RA_sp, decreasing = TRUE)
  cumulative_abundance <- cumsum(sorted_abundance) / sum(sorted_abundance)
  top_60_index <- which(cumulative_abundance <= 1)
  return(names(sorted_abundance)[top_60_index])
}

# Apply function to each sample type
top_species_d <- select_top_60(sum_d)
top_species_s <- select_top_60(sum_s)
top_species_l <- select_top_60(sum_l)

# Get unique species across all sample types
top_species <- unique(c(top_species_d, top_species_s, top_species_l))
length(top_species)

# Filter top_100_fa_2 for selected species
top_100_fa_2_filtered <- top_100_fa_2[top_species, ]
col_funct<- rainbow(48)
col_funct<- c( 
'#ff7f00',
'#e41a1c',
'#80cdc1',
'#f781bf',
'#377eb8',
'#4daf4a',
'#ffff33',
'#984ea3',
'#999999',
'#a65628',
'#fff',
'#000')
colour_families<- as.factor(top_species)

colour_map<- setNames(col_funct, levels(colour_families))
colour_families <- colour_map[colour_families]

cbind(colour_families, top_species)

par(mar = c(7, 8, 6, 24))
par(mgp = c(5, 1, 0)) 
bp <- barplot(top_100_fa_2_filtered, col = colour_families, 
              ylab = "Relative Abundance", xlab = "Microhabitat", las = 2, 
              cex.lab=2, cex.axis=1.5)

par(xpd = TRUE)
legend("topright", inset = c(-0.8, 0), legend = names(colour_map), 
       fill = col_funct, title = "Most abundant classes", cex = 1 )

# Reset the plotting behavior
par(xpd = FALSE)


# Barplot with filtered species
barplot(top_100_fa_2_filtered, col = colour_families, las = 2, main = "Top 60% Abundance by Sample Type")


# match with lichen fungi orders 
# List of fungal orders to match
fungal_orders <- c("Acarosporales", "Agaricales", "Arthoniales", "Atheliales", "Baeomycetales",
                   "Caliciales", "Candelariales", "Cantharellales", "Capnodiales", "Chaetothyriales",
                   "Collemopsidiales", "Coniocybales", "Corticiales", "Eremithallales", "Lecanorales",
                   "Lecideales", "Lepidostromatales", "Leprocaulales", "Lichinales", "Monoblastiales",
                   "Odontotrematales", "Ostropales", "Peltigerales", "Pertusariales", "Phaeomoniellales",
                   "Pleosporales", "Pyrenulales", "Rhizocarpales", "Sarrameanales", "Schaereriales",
                   "Strigulales", "Teloschistales", "Thelenellales", "Thelocarpales", "Trypetheliales",
                   "Umbilicariales", "Verrucariales", "Vezdaeales", "Xylariales")

# Filter top species that match the list of fungal orders
matching_species <- top_species[top_species %in% fungal_orders]

# Subset the top_100_fa_2_filtered matrix
matching_fa_filtered <- top_100_fa_2_filtered[matching_species, ]



barplot(matching_fa_filtered, las = 2, main = "Lichenised Fungi Abundance by Sample Type")

# Sum the abundance for each fungal order
order_abundance <- rowSums(matching_fa_filtered)

# Sort orders by total abundance in descending order
sorted_orders <- sort(order_abundance, decreasing = TRUE)

# Select the top 10 most abundant orders
top_10_orders <- names(sorted_orders)[1:10]

# Subset the original abundance matrix to keep only the top 10 orders
top_10_fa_filtered_matching <- matching_fa_filtered[top_10_orders, ]

barplot(top_10_fa_filtered_matching, las = 2, main = "Lichenised Fungi Abundance by Sample Type")
colour_families<- as.factor(top_10_orders)

colour_map<- setNames(col_funct, levels(colour_families))
colour_families <- colour_map[colour_families]

cbind(colour_families, top_10_orders)

svg("mha_lichen_SB2.svg", width = 12, height = 10)
par(mar = c(7, 8, 6, 24))
par(mgp = c(5, 1, 0)) 

# Scale y-values to percentage
top_10_fa_filtered_matching_percent <- top_10_fa_filtered_matching * 100

# ploy
p <- barplot(top_10_fa_filtered_matching_percent, 
              col = colour_families, 
              ylab = "Relative Abundance (%)", 
              xlab = "Microhabitat", 
              las = 1,  
              cex.lab = 2, 
              cex.axis = 1.8, 
              names.arg = c("Deadwood", "Soil", "Living Wood"), 
              cex.names = 1.8) 

# Add the legend with the desired settings
par(xpd = TRUE)
legend("topright", inset = c(-0.6, 0), legend = names(colour_map), 
       fill = col_funct, title = "Top 10 Lichenised Orders", cex = 1)


dev.off()


##### 6. ENS mixed model ####

sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

# subset to library and make numeric 
ab_full<- sp_full[10:3038, 8:152]

for (i in 1:ncol(ab_full)){
  x<- as.numeric(ab_full[,i])
  ab_full[,i]<-x 
}



full_rich<- specnumber(t(ab_full))
full_div<- diversity(t(ab_full), index= 'shannon')
#full_div<- as.matrix(full_div)

length(full_div)
#t_full_div<- t(full_div)

# test difference due to sample type without having site random effect 
sample_type<- unlist(lapply(colnames(ab_full),function(x)(substr(x,start=3, stop=3))))
# so dw is last? 
sample_type[sample_type=='D']<- 'D'
# make a factor 
sample_type<- as.factor(sample_type)
levels(sample_type)


# subset without Rothie and C7 and C8 
nosoilsites<- which(sites =='C7' | sites =='C8' |sites =='RM' | sites =='RO')
boxplot(full_rich[-nosoilsites]~sample_type[-nosoilsites]) 
model1<- summary(lm(full_rich[-nosoilsites]~sample_type[-nosoilsites]))

DSL_rich<-full_rich[-nosoilsites]
DSL_div<- full_div[-nosoilsites]
DSL_st<- sample_type[-nosoilsites]
DSL_expdiv<- exp(DSL_div)

boxplot(exp(DSL_div)~DSL_st) 
boxplot(DSL_div~DSL_st) 

model5 <- lmer(DSL_rich ~ DSL_st + (1 | sites[-nosoilsites]))
summary(model5)


model5 <- lmer(DSL_expdiv ~ DSL_st + (1 | sites[-nosoilsites]))
summary(model5)

##### 7. Pairwise boxplots #####

sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

# subset to library and make numeric 
ab_full<- sp_full[10:3038, 8:152]

for (i in 1:ncol(ab_full)){
  x<- as.numeric(ab_full[,i])
  ab_full[,i]<-x 
}


sites<- unlist(lapply(colnames(ab_full),function(x)(substr(x,start=1, stop=2))))
sample_type<- unlist(lapply(colnames(ab_full),function(x)(substr(x,start=3, stop=3))))

library(vegan)
full_rich<- specnumber(t(ab_full))
full_div<- diversity(t(ab_full), index= 'shannon')
full_expdiv<- exp(full_div)


##### 

dw_s_R<- NULL
dw_l_R<- NULL
s_l_R<- NULL

dw_s_D<- NULL
dw_l_D<- NULL
s_l_D<- NULL

dw_s_E<- NULL
dw_l_E<- NULL
s_l_E<- NULL

#factor_NULL
factor_sites<- unique(sites)
factor_st<- unique(sample_type)

i<-3

for (i in 1:length(factor_sites)){
  x<- which(sites==factor_sites[i])
  y<- which(sample_type[x]=="D")
  z<- which(sample_type[x]=='S')
  q<- which(sample_type[x]=='L')
  dw_s_R1<- full_rich[x][y]-full_rich[x][z]
  dw_l_R1<- full_rich[x][y]-full_rich[x][q]
  s_l_R1<- full_rich[x][z]- full_rich[x][q]
  
  dw_s_D1<- full_div[x][y]-full_div[x][z]
  dw_l_D1<- full_div[x][y]-full_div[x][q]
  s_l_D1<- full_div[x][z]- full_div[x][q]
  
  dw_s_E1<- full_expdiv[x][y]-full_expdiv[x][z]
  dw_l_E1<- full_expdiv[x][y]-full_expdiv[x][q]
  s_l_E1<- full_expdiv[x][z]- full_expdiv[x][q]
  
  dw_s_R<- c(dw_s_R,dw_s_R1)
  dw_l_R<- c(dw_l_R,dw_l_R1)
  s_l_R<- c(s_l_R,s_l_R1)
  
  dw_s_D<- c(dw_s_D,dw_s_D1)
  dw_l_D<- c(dw_l_D,dw_l_D1)
  s_l_D<- c(s_l_D,s_l_D1)
  
  dw_s_E<- c(dw_s_E,dw_s_E1)
  dw_l_E<- c(dw_l_E,dw_l_E1)
  s_l_E<- c(s_l_E,s_l_E1)
  
}

# boxplot dw_l 

mean(dw_l_E)
sd(dw_l_E)

boxplot(dw_s_R, dw_l_R, s_l_R)

boxplot(dw_s_D, dw_l_D, s_l_D)
lines(x=c(-1000,1000), y=c(0,0), lwd=3, col='red')

svg("boxplot_microhabitat.svg", width = 8, height = 6)
par(mar=c(5,5,5,5))
boxplot(dw_s_E, dw_l_E, s_l_E, xlab= 'Microhabitat', ylab= 'Effective number of species', cex.lab=2, cex.axis=1.5)
lines(x=c(-1000,1000), y=c(0,0), lwd=3, col='red')


dev.off()

library(svglite)

library(ggplot2)
library(dplyr)

# Create a data frame for ggplot
df_boxplot <- data.frame(
  Microhabitat = rep(c("dw_s_E", "dw_l_E", "s_l_E"), 
                     times = sapply(list(dw_s_E, dw_l_E, s_l_E), length)),
  ENS = c(dw_s_E, dw_l_E, s_l_E)  # Combine data into one column
)

# Create the boxplot with a horizontal red reference line at y = 0
ggplot(df_boxplot, aes(x = Microhabitat, y = ENS)) +
  geom_boxplot(fill = "lightblue", color = "black") +  # Boxplot styling
  geom_hline(yintercept = 0, linetype = "solid", color = "red", size = 1) +  # ref 
  labs(x = "Microhabitat", y = "Effective Number of Species") +  # lab
  theme_minimal()



boxplot(dw_s_D)
boxplot(dw_l_R)
boxplot(dw_l_D)
boxplot(dw_s_R)
boxplot(dw_s_D)
boxplot(s_l_R)
boxplot(s_l_D)


# Perform one-sample t-tests for richness
t_dw_s_R <- t.test(dw_s_R, mu = 0)
t_dw_l_R <- t.test(dw_l_R, mu = 0)
t_s_l_R <- t.test(s_l_R, mu = 0)

# Perform one-sample t-tests for diversity
t_dw_s_D <- t.test(dw_s_D, mu = 0)
t_dw_l_D <- t.test(dw_l_D, mu = 0)
t_s_l_D <- t.test(s_l_D, mu = 0)

t_dw_s_E <- t.test(dw_s_E, mu = 0)
t_dw_l_E <- t.test(dw_l_E, mu = 0)
t_s_l_E <- t.test(s_l_E, mu = 0)

