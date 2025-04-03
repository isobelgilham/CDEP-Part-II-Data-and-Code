##### 1. Mixed model ENS DS #####

library(vegan)
library(lme4)


# load in full df 
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[-c(1)]

# establish all factors 
sample_type<- unlist(lapply(colnames(sp_full),function(x)(substr(x,start=3, stop=3))))
deadwood<- which(sample_type =='D')
lib_dw<- sp_full[,deadwood]


# add sp back 
lib_dw<- cbind(sp_full[,1:7], lib_dw)


sample_loc<- unlist(lapply(colnames(lib_dw[, 8:99]),function(x)(substr(x,start=7, stop=7))))
sites<- unlist(lapply(colnames(lib_dw[,8:99]),function(x)(substr(x,start=1, stop=2))))

# env factors 
tree<- unlist(lib_dw[3,8:99]) 
# numeric to make DBH and length continuous 
DBH<- as.numeric(unlist(lib_dw[5, 8:99 ]))
length<- as.numeric(unlist(lib_dw[6,8:99 ]))
decay_stage<- as.numeric(unlist(lib_dw[7,8:99 ])) 
avg_soil_moist<- as.numeric(unlist(lib_dw[8,8:99 ]))
avg_soil_ph<- as.numeric(unlist(lib_dw[9,8:99 ]))

# fix where decay stage = 6 
# Replace decay stage 6 with 5
# Replace decay stage 6 with 5
lib_dw[7,][lib_dw[7,] == 6] <- 5
# Replace decay stage 6 with 5


###### Define richness and diversity ######

# full lib 
lib_dead<- lib_dw[10:3038,]
lib_dead[,8:ncol(lib_dead)]<- lapply(lib_dead[, 8:ncol(lib_dead)], as.numeric) 

# richness 
richness<- specnumber(t(lib_dead[,8:ncol(lib_dead)]))


# diversity 
div<- diversity(t(lib_dead[,8:ncol(lib_dead)]), index= 'shannon')

# TEST control for treatment 
model1<- lmer(exp(div)~decay_stage+ (1|sites))
summary(model1)


##### 2. ANCOMBC DS #####
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


str(ab_full)

# define variables 
sample_type<- unlist(lapply(colnames(sp_full[, 8:152]),function(x)(substr(x,start=3, stop=3))))
sites<- unlist(lapply(colnames(sp_full[,8:152]),function(x)(substr(x,start=1, stop=2))))


dwc<- rep(NA, length(sites))
dwc[sites %in% c('GO', 'RO','GM', 'RM')] <- "control"
dwc[sites %in% c('C7', 'C8', 'C5', 'C4', 'A5', 'A4', 'A2', 'A6' )] <- "dwcd"
dwc <- factor(dwc)
length(dwc)

decay_stage<- as.numeric(unlist(lib_dw[7,8:99 ])) 


sample_type[sample_type=='D']<- 'D'

# otu table 
OTU_TABLE<- otu_table(ab_full, taxa_are_rows = T)

OTU_TABLE<- OTU_TABLE[-c(which(rowSums(OTU_TABLE)<500)),]


dim(OTU_TABLE)

samples<- as.data.frame(cbind(sample_type, sites))

SAMPLE_TABLE<- sample_data(samples)

rownames(SAMPLE_TABLE)<- colnames(OTU_TABLE)

ps1<- phyloseq(OTU_TABLE, SAMPLE_TABLE)

# Run ANCOM-BC
model_ancom <- ancombc(ps1, formula = "ds", 
                       lib_cut = 1000,        # Min library size cutoff
                       n_cl = 5,              # Number of cores
                       group = "ds",          # No repeated measures
                       global = TRUE)     # Global significance test

# Extract Global Significance
sign_dif_ds <- which(model_ancom$res_global$diff_abn == TRUE)


library(phyloseq)
library(ANCOMBC)

# Convert decay_stage to a factor
samples$ds <- factor(samples$ds)

# Create OTU Table (Filter low-abundance taxa)
OTU_TABLE <- otu_table(ab_dw, taxa_are_rows = TRUE)
OTU_TABLE <- OTU_TABLE[-c(which(rowSums(OTU_TABLE) < 500)),]  # Remove low-count taxa

# Create Sample Data (Ensure row names match OTU table)
rownames(samples) <- colnames(OTU_TABLE)
SAMPLE_TABLE <- sample_data(samples)

# Create Phyloseq Object
ps1 <- phyloseq(OTU_TABLE, SAMPLE_TABLE)

# Run ANCOM-BC with decay_stage as the grouping variable
model_ancom <- ancombc(ps1, formula = "ds", 
                       lib_cut = 1000,        # Min library size cutoff
                       n_cl = 5,              # Number of cores
                       group = "ds",          # No repeated measures
                       global = TRUE)     # Global significance test

# Extract Differentially Abundant Features
sign_dif_ds <- which(model_ancom$res$diff_abn[, "ds2"] == TRUE)  # Adjust for actual levels
q_values <- model_ancom$res$q_val
length(sign_dif_ds)

# View q-values
head(q_values)

# Extract Global Significance
sign_dif_decay <- which(model_ancom$res_global$diff_abn == TRUE)
length(sign_dif_decay)

# View Significant Features
View(lib_dw[sign_dif_decay,])

# save as csv 
indicators_decay<- sp_full[c(10:3038)[sign_dif_decay],-c(8:152)]
write.csv(indicators_decay,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/indicators_decay.csv')


##### 3. Stacked barplots #####

# subset indicator species 
RA_decay<- RA_dw[,sign_dif_decay]
# RA_decay= relative abundances, no taxa 

# library subset with decay stages 
sample_type<- unlist(lapply(colnames(sp_full),function(x)(substr(x,start=3, stop=3))))
deadwood<- which(sample_type =='D')
lib_dw_taxa<- sp_full[,deadwood]

# decay stage as factor 
lib_dw_taxa[7,][lib_dw_taxa[7,] == 6] <- 5
decay_stage<- as.factor(unlist(lib_dw_taxa[7,1:92 ]))

lib_decay<- lib_dw[sign_dif_decay,]


decay_RA_G<- matrix(NA, nrow=nrow(RA_decay), ncol=length(unique(lib_decay$O))) 

Genera<-unique(lib_decay$O)

i<- 1

for(i in 1:length(Genera)){
  x<- which(lib_decay$O== Genera[i])
  if(length(x)>1) (y<- rowSums(RA_decay[,x]))  
  if(length(x)==1)(y<- RA_decay[,x])
  decay_RA_G[,i]<- y 
  
}

rownames(decay_RA_G)<- rownames(RA_decay)
colnames(decay_RA_G)<- Genera 

# group by decay stage treatment 

ds_1<- which(decay_stage =='1')
ds_2<- which(decay_stage =='2')
ds_3<- which(decay_stage =='3')
ds_4<- which(decay_stage =='4')
ds_5<- which(decay_stage =='5')


groups<- c(rep("ds_1", times= length(ds_1)), 
           rep("ds_2", times= length(ds_2)),
           rep("ds_3", times= length(ds_3)),
           rep("ds_4", times= length(ds_4)),
           rep("ds_5", times= length(ds_5)))

# sum by group 

sum_1<- colMeans(decay_RA_G[ds_1,])
sum_2<- colMeans(decay_RA_G[ds_2,])
sum_3<- colMeans(decay_RA_G[ds_3,])
sum_4<- colMeans(decay_RA_G[ds_4,])
sum_5<- colMeans(decay_RA_G[ds_5,])


decay_RA_G_2<- cbind(sum_1, sum_2, sum_3, sum_4, sum_5)



# plot and sort out colours 
colour_families<- as.factor(Genera)
length(levels(colour_families))

# select 12 colours to colour by site 
col_funct<- c( '#e41a1c','#fff','#377eb8','#4daf4a','#984ea3','#000','#ff7f00','#80cdc1','#ffff33','#a65628','#f781bf','#999999' )


colour_map<- setNames(col_funct, levels(colour_families))
colour_families <- colour_map[colour_families]

cbind(colour_families, Genera)


levels(colour_families)<- col_funct
colour_families<- as.character(colour_families)

cbind(colour_families, Genera)

barplot(decay_RA_G_2, col= colour_families)

# plot
svg("DS_SB5.svg", width = 10, height = 10)
par(mar = c(7, 8, 6, 24))
par(mgp = c(5, 1, 0)) 
bp <- barplot(decay_RA_G_2, col = colour_families, 
              ylab = "Relative Abundance", xlab = "Decay stage", las = 2, 
              cex.lab=2, cex.axis=1.5)

par(xpd = TRUE)

# Add a legend outside the plot on the right
legend("topright", inset = c(-1.1, 0), legend = names(colour_map), 
       fill = col_funct, title = "Orders", cex = 1 )

# Reset the plotting behavior
par(xpd = FALSE)

dev.off()

############################ Plot decay 
svg("DS_SBnew.svg", width = 10, height = 10)
par(mar = c(7, 8, 6, 24))
par(mgp = c(5, 1, 0)) 

# Convert decay_RA_G_2 to percentage
decay_RA_G_2_percent <- decay_RA_G_2 * 100

# Bar plot with percentage y-axis and x-axis labels 1, 2, 3, 4, 5
bp <- barplot(decay_RA_G_2_percent, col = colour_families, 
              ylab = "Relative Abundance (%)", xlab = "Decay stage", 
              las = 1, cex.lab = 2, cex.axis = 1.5, 
              names.arg = c(1, 2, 3, 4, 5), cex.names = 2)

par(xpd = TRUE)

# Add a legend outside the plot on the right
legend("topright", inset = c(-1.1, 0), legend = names(colour_map), 
       fill = col_funct, title = "Orders", cex = 1)

# Reset the plotting behavior
par(xpd = FALSE)

dev.off()

