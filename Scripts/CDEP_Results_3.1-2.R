# 1. Make dwc variable 
# 2. Decay stage linear models dwc/control
# 3. Line plots decay stage and diversity 
# 4. Histogram decay stage and dwc 
# 5. ANCOMBC 
# 5.1 by decay stage 
# 5.2 by dwc 

# load in data 
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

# subset to library and make numeric 
ab_full<- sp_full[10:3038, 8:152]

for (i in 1:ncol(ab_full)){
  x<- as.numeric(ab_full[,i])
  ab_full[,i]<-x 
}

# establish all factors 
sample_type<- unlist(lapply(colnames(sp_full),function(x)(substr(x,start=3, stop=3))))
deadwood<- which(sample_type =='D')
lib_dw<- sp_full[,deadwood]

sample_type_ab<- unlist(lapply(colnames(ab_full),function(x)(substr(x,start=3, stop=3))))
deadwood<- which(sample_type_ab =='D')
ab_dw<- ab_full[,deadwood]

# add sp back 
lib_dw<- cbind(sp_full[,1:7], lib_dw)


sample_loc<- unlist(lapply(colnames(lib_dw[, 8:99]),function(x)(substr(x,start=7, stop=7))))
sites<- unlist(lapply(colnames(lib_dw[,8:99]),function(x)(substr(x,start=1, stop=2))))

##### 1. Establish factors #####

tree<- unlist(lib_dw[3,8:99]) 
# numeric to make DBH and length continuous 
DBH<- as.numeric(unlist(lib_dw[5, 8:99 ]))
length<- as.numeric(unlist(lib_dw[6,8:99 ]))
lib_dw[7,][lib_dw[7,] == 6] <- 5
decay_stage<- as.numeric(unlist(lib_dw[7,8:99 ])) 
avg_soil_moist<- as.numeric(unlist(lib_dw[8,8:99 ]))
avg_soil_ph<- as.numeric(unlist(lib_dw[9,8:99 ]))


##### 2. Make dwc variable #####
# treatment 
dwc<- rep(NA, length(sites))

dwc[sites %in% c('GO', 'RO','GM', 'RM')] <- "control"
dwc[sites %in% c('C7', 'C8', 'C5', 'C4', 'A5', 'A4', 'A2', 'A6' )] <- "dwcd"
dwc <- factor(dwc)

length(dwc)


# define diversity metrics 
# full lib 
library(vegan)
lib_dead<- lib_dw[10:3038,]
lib_dead[,8:ncol(lib_dead)]<- lapply(lib_dead[, 8:ncol(lib_dead)], as.numeric) 

# richness 
richness<- specnumber(t(lib_dead[,8:ncol(lib_dead)]))

# diversity shannon and ENS 
div<- diversity(t(lib_dead[,8:ncol(lib_dead)]), index= 'shannon')
expdiv<- exp(div)


plot(expdiv~decay_stage)
model<- lm(expdiv~dwc)
summary(model)
View(expdiv)

model<-lm(expdiv~decay_stage)
summary(model)




##### 3. Linear models #####
plot(expdiv~decay_stage)
model<- lm(expdiv~decay_stage)
summary(model)


model<- lm(decay_stage~dwc)
summmary(model)

# intial model 
model<- (lm(expdiv~sample_loc+decay_stage+dwc))
summary(model)
drop1(model)

# drop sample_loc 
model<- (lm(expdiv~decay_stage+dwc))
summary(model)
drop1(model)

# drop DBH 
model2<- (lm(div~tree+treatment+decay_stage))
summary(model2)
drop1(model2)

# drop tree 
model2<- (lm(div~treatment+decay_stage))
summary(model2)
drop1(model2)

# separate stepwise models for control and 
expdiv_cont <- subset(expdiv, dwc == "control")
sl_cont <- subset(sample_loc, dwc == "control")
ds_cont<- subset(decay_stage, dwc == "control")
dbh_cont<- subset(DBH, dwc == "control")
length_cont<- subset(length, dwc == "control")


expdiv_dwcd <- subset(expdiv, dwc == "dwcd")
sl_dwcd <- subset(sample_loc, dwc == "dwcd")
ds_dwcd<- subset(decay_stage, dwc == "dwcd")
dbh_dwcd<- subset(DBH, dwc == "dwcd")
length_dwcd<- subset(length, dwc == "dwcd")

# stepwise dwcd
model<- (lm(expdiv_dwcd~sl_dwcd+ds_dwcd +dbh_dwcd+length_dwcd))
summary(model)
drop1(model)

# drop length 
model<- (lm(expdiv_dwcd~ds_dwcd +dbh_dwcd +sl_dwcd))
summary(model)
drop1(model)

#repeat for control 
model<- (lm(expdiv_cont~sl_cont +ds_cont))
summary(model)
drop1(model)

# drop length 
model<- (lm(expdiv_dwcd~ds_dwcd +dbh_dwcd +sl_dwcd))
summary(model)
drop1(model)

# check ds ens 
model<- lm(expdiv_cont~ds_cont)
summary(model)



##### 4. Plot #####
library(svglite)


# df for plot 
df <- data.frame(
  decay_stage = decay_stage,
  expdiv = expdiv,
  dwc = dwc
)
levels(df$dwc)

dwc_colors <- c("dwcd" = "red", "control" = "cyan")
dwc_colors <- setNames(c("red", "cyan"), levels(df$dwc))


# dwc to factor
df$dwc <- factor(trimws(df$dwc))

# set colours 
dwc_colors <- setNames(c("blue", "red"), levels(df$dwc))

# Plot with ggplot
p <- ggplot(df, aes(x = decay_stage, y = expdiv, color = dwc, linetype = dwc)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(data = df[df$dwc == levels(df$dwc)[1],], method = "lm", se = TRUE, alpha = 0.2, color = dwc_colors[1]) +  
  geom_smooth(data = df[df$dwc == levels(df$dwc)[2],], method = "lm", se = TRUE, alpha = 0.2, color = dwc_colors[2]) +  
  scale_color_manual(values = dwc_colors) +  
  scale_linetype_manual(values = c("solid", "dotted")) +  
  labs(
    title = "Effective number of species vs. Decay Stage",
    x = "Decay Stage",
    y = "ENS (no. of species)",
    color = "Treatment",
    linetype = "Treatment"
  ) +
  theme_bw() + 
  theme(
    axis.text = element_text(size = 14, color = "black"),  
    axis.title = element_text(size = 16),  
    plot.title = element_text(size = 18),  
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),  
    legend.key = element_rect(fill = "white", colour = "black", size = 0.5),  
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12),  
    legend.box.background = element_rect(colour = "black", size = 1)  
  )

# download plot
ggsave("ENSvsDSbyT.svg", plot = q, width = 8, height = 6, dpi = 300)


##### 6. PERMANOVA  #####

lib_dw_notaxa <- as.data.frame(sapply(lib_dw[, 8:99], as.numeric))
# dwc 
Model<- adonis(t(lib_dw_notaxa)~dwc, permutations = 9999)
sum_model<- summary(Model)
Model$aov.tab

##### 7. ANCOMBC dwc  #####
#install.packages('BiocManager')
#BiocManager::install("ANCOMBC")
library(ANCOMBC)

#BiocManager::install("phyloseq")

library(phyloseq)

#BiocManager::install('microbiome')

library('microbiome')


# subset to abundances 
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

ab_full<- sp_full[10:3038, 8:152]

for (i in 1:ncol(ab_full)){
  x<- as.numeric(ab_full[,i])
  ab_full[,i]<-x 
}

# establish all factors 
sample_type<- unlist(lapply(colnames(sp_full),function(x)(substr(x,start=3, stop=3))))
deadwood<- which(sample_type =='D')
lib_dw<- sp_full[,deadwood]

sample_type_ab<- unlist(lapply(colnames(ab_full),function(x)(substr(x,start=3, stop=3))))
deadwood<- which(sample_type_ab =='D')
ab_dw<- ab_full[,deadwood]

# add sp back 
lib_dw<- cbind(sp_full[,1:7], lib_dw)


sample_loc<- unlist(lapply(colnames(lib_dw[, 8:99]),function(x)(substr(x,start=7, stop=7))))
sites<- unlist(lapply(colnames(lib_dw[,8:99]),function(x)(substr(x,start=1, stop=2))))

#1. Establish factors

tree<- unlist(lib_dw[3,8:99]) 
# numeric to make DBH and length continuous 
DBH<- as.numeric(unlist(lib_dw[5, 8:99 ]))
length<- as.numeric(unlist(lib_dw[6,8:99 ]))
decay_stage<- as.numeric(unlist(lib_dw[7,8:99 ])) 
avg_soil_moist<- as.numeric(unlist(lib_dw[8,8:99 ]))
avg_soil_ph<- as.numeric(unlist(lib_dw[9,8:99 ]))
lib_dw[7,][lib_dw[7,] == 6] <- 5

# 1. Make dwc variable
# treatment 
dwc<- rep(NA, length(sites))
dwc[sites %in% c('GO', 'RO','GM', 'RM')] <- "control"
dwc[sites %in% c('C7', 'C8', 'C5', 'C4', 'A5', 'A4', 'A2', 'A6' )] <- "dwcd"
dwc <- factor(dwc)


# otu table 
OTU_TABLE<- otu_table(ab_dw, taxa_are_rows = T)
OTU_TABLE<- OTU_TABLE[-c(which(rowSums(OTU_TABLE)<500)),]

dim(OTU_TABLE)

# samples df 
samples<- as.data.frame(cbind(dwc, sites))
SAMPLE_TABLE<- sample_data(samples)
rownames(SAMPLE_TABLE)<- colnames(OTU_TABLE)
ps1<- phyloseq(OTU_TABLE, SAMPLE_TABLE)
model_ancom<- ancombc(ps1,formula= 'dwc', n_cl = 5, global =T, group = NULL)
#, group= 'sites', global= T )
sign_dif<- which(model_ancom$res$diff_abn[,4]==T)
q_values <- model_ancom$res$q_val

# View q-values
head(q_values)

sign_dif_dwc<- which(model_ancom$res_global$diff_abn==T)
View(sp_full[c(10:3038)[sign_dif_dwc],-c(8:152)])
length(sign_dif_dwc)


q <- RDAqh$CCA$rank

# species responsible for the difference between sample types 

indicators_sample_type3<- sp_full[c(10:3038)[sign_dif3],-c(8:152)]

Working ancombc by dwc  
library(phyloseq)
library(ANCOMBC)

# Create OTU Table (Filter low-abundance taxa)
OTU_TABLE <- otu_table(ab_dw, taxa_are_rows = TRUE)
OTU_TABLE <- OTU_TABLE[-c(which(rowSums(OTU_TABLE) < 500)),]

# Create Sample Data (Ensure Factors Are Correct)
samples <- data.frame(dwc = factor(dwc), sites = sites, ds= decay_stage)
rownames(samples) <- colnames(OTU_TABLE)  
SAMPLE_TABLE <- sample_data(samples)

# Create Phyloseq Object
ps1 <- phyloseq(OTU_TABLE, SAMPLE_TABLE)

# Run ANCOM-BC
model_ancom <- ancombc(ps1, formula = "dwc", 
                       lib_cut = 1000,        
                       n_cl = 5,              
                       group = "dwc",          
                       global = TRUE)         

# Extract diff ab spec
sign_dif_dwc_2<- which(model_ancom$res$diff_abn[,2]==T)
length(sign_dif_dwc_2)

sign_dif_dwcdwcd <- which(model_ancom$res$diff_abn[, "dwcdwcd"] == TRUE)  
q_values <- model_ancom$res$q_val
length(sign_dif_dwcdwcd)

# how are these different? the second one is surely the correct one 
# will use both and see what happens but I think the dwcdwcd is better as im not sure the first is right

# View q-values
head(q_values)

# Extract global sig species
sign_dif_dwc_global <- which(model_ancom$res_global$diff_abn == TRUE)
length(sign_dif_dwc_global)

# View sig spec
View(sp_full[c(10:3038)[sign_dif_dwc_2], -c(8:152)])


# Count no of sig sp
length(sign_dif_dwc)


## save as csv 
indicators_dwc_2<- sp_full[c(10:3038)[sign_dif_dwc_2],-c(8:152)]
write.csv(indicators_dwc_2,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/indicators_dwc_2.csv')

indicators_dwcdwcd<-sp_full[c(10:3038)[sign_dif_dwcdwcd],-c(8:152)]
write.csv(indicators_dwcdwcd,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/indicators_dwcdwcd.csv')


##### 8. stacked bars dwc  #####

# load in data 
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[-c(1)]


# subset to libraries only
lib_full<- sp_full[10:nrow(sp_full),]

# set as numeric 
lib_full[,8:ncol(lib_full)]<- lapply(lib_full[, 8:ncol(lib_full)], as.numeric) 

# get relative abundances 
library(vegan)
RA_sp<- decostand(t(lib_full[,8:152]),method = 'total')


# subset deadwood only 
sample_type<- unlist(lapply(rownames(RA_sp[1:145,]),function(x)(substr(x,start=3, stop=3))))

# deadwood 
deadwood<- which (sample_type == 'D')

# subset deadwood 
RA_dw<- RA_sp[deadwood,]


# deadwood library 
lib_dw<- lib_full[,deadwood +7]
lib_dw<- cbind(lib_full[,1:7], lib_dw)


# subset indicator species 
RA_dwc<- RA_dw[,sign_dif_dwcdwcd]
# RA_dwc= relative abundances, no taxa 

lib_dwc<- lib_dw[sign_dif_dwcdwcd,]
# lib_dwc = raw counts, with taxa but same index


# for loop to get to order level 

dwc_RA_G<- matrix(NA, nrow=nrow(RA_dwc), ncol=length(unique(lib_dwc$C))) 

Genera<-unique(lib_dwc$C)

i<- 1

for(i in 1:length(Genera)){
  x<- which(lib_dwc$C== Genera[i])
  if(length(x)>1) (y<- rowSums(RA_dwc[,x]))  
  if(length(x)==1)(y<- RA_dwc[,x])
  dwc_RA_G[,i]<- y 
  
}

rownames(dwc_RA_G)<- rownames(RA_dwc)
colnames(dwc_RA_G)<- Genera 

# group by dw treatment 
sites<- unlist(lapply(rownames(RA_dwc),function(x)(substr(x,start=1, stop=2))))
dwc<- rep(NA, length(sites))
dwc[sites %in% c('GO', 'RO','GM', 'RM')] <- "control"
dwc[sites %in% c('C7', 'C8', 'C5', 'C4', 'A5', 'A4', 'A2', 'A6' )] <- "dwcd"
dwc <- factor(dwc)

Controls<- which(dwc =='control')
Dwcreation<- which(dwc=='dwcd')

groups<- c(rep("Control", times= length(Controls)), 
           rep("Dwcd", times= length(Dwcreation)))

# sum by group 
sites<- unlist(lapply(rownames(dwc_RA_G),function(x)(substr(x,start=1, stop=2))))
dwc<- rep(NA, length(sites))
dwc[sites %in% c('GO', 'RO','GM', 'RM')] <- "control"
dwc[sites %in% c('C7', 'C8', 'C5', 'C4', 'A5', 'A4', 'A2', 'A6' )] <- "dwcd"
dwc <- factor(dwc)

Controls<- which(dwc =='control')
Dwcreation<- which(dwc=='dwcd')

sum_c<- colMeans(dwc_RA_G[Controls,])
sum_dwcd<- colMeans(dwc_RA_G[Dwcreation,])

dwc_RA_G_2<- cbind(sum_c, sum_dwcd)

# Plot by dwc by class (9 levels) 

# plot and sort out colours 
colour_families<- as.factor(Genera)
length(levels(colour_families))

# select 17 colours to colour by site 
col_funct<- c('#377eb8', '#4daf4a', '#e41a1c', '#984ea3', '#ff7f00', '#999999', '#a65628',  '#f781bf','#ffff33' )



colour_map<- setNames(col_funct, levels(colour_families))
colour_families <- colour_map[colour_families]

cbind(colour_families, Genera)


# final plot dwc 
# Assuming dwc_RA_G_2 and colour_families are already defined

# Create the barplot in the svg output file
svg("DWC_SB_new3.svg", width = 10, height = 10)
par(mar = c(7, 8, 6, 18))  # Set margins
par(mgp = c(5, 1, 0))  # Set margin for axis labels

# Convert relative abundance to percentage
dwc_RA_G_2_percent <- dwc_RA_G_2 * 100

# Barplot with specified settings
bp <- barplot(dwc_RA_G_2_percent, 
              col = colour_families, 
              ylab = "Relative Abundance (%)", 
              xlab = "Treatment", 
              las = 1,  # Rotate x-axis labels to be horizontal
              cex.lab = 2,  # Make axis labels larger
              cex.axis = 1.5,  # Make axis text larger
              names.arg = c("Control", "Deadwood Creation"),  # Modify bar names
              cex.names = 2  # Make bar names larger
)

# Add a legend outside the plot on the right
par(xpd = TRUE)  # Allow drawing outside the plot region
legend("topright", inset = c(-0.6, 0), legend = names(colour_map), 
       fill = col_funct, title = "Class", cex = 1)

# Reset the plotting behavior
par(xpd = FALSE)

# Close the svg device
dev.off()
