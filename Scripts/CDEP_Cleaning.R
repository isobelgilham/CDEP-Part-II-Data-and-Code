
rm(list=ls())

# 1./ import libraries, split classification 
# 2./ import barcode ids, including DBH and decay stage data 
# 3.X add environmental variables to the dataframe before cleaning 
# 4./ remove columns and rows with 0 reads, remove erroneous C5D5.T
# 5./cleaning 
      # remove rare species, controls and eukaryote sp., insertae species 
# 6./ run nmds to find outliers 
# 7./ remove those outliers 
# Should be left with the df with which all other analysis will be conducted 
# 8./ Separate this into fungi and non fungi species
# 9./ Add environmental variables
# 10./ Download as csv files
# 12. GHG data 

# load libraries 
library(vegan)
library(tidyr)


##### 1. Import: libraries #####
# load in data 
library_1<-read.csv("Data/Library_1_Table.csv")
library_2<-read.csv("Data/Library_2_Table.csv")

# split classification columns up using separate wider regex 
library(tidyr)

library_1<- separate_wider_regex(data= library_1, col= X, 
                                 patterns = c(K= "k_+\\w+", 
                                              ";", 
                                              P= "p_+\\w+", 
                                              ";",
                                              C="c_+\\w+",
                                              ";", 
                                              O="o_+\\w+",
                                              ";",
                                              Fa="f_+\\w+",
                                              ";",
                                              G="g_+\\w+",
                                              ";", 
                                              S="s_+\\w+"
                                 ),
                                 too_few= 'align_start'
)

library_2<- separate_wider_regex(data= library_2, col= X, 
                                 patterns = c(K= "k_+\\w+", 
                                              ";", 
                                              P= "p_+\\w+", 
                                              ";",
                                              C="c_+\\w+",
                                              ";", 
                                              O="o_+\\w+",
                                              ";",
                                              Fa="f_+\\w+",
                                              ";",
                                              G="g_+\\w+",
                                              ";", 
                                              S="s_+\\w+"
                                 ),
                                 too_few= 'align_start'
)

## removing the underscores and markers from each 
library_1$K<-sub('k__','', library_1$K)
library_1$P<-sub('p__','', library_1$P)
library_1$C<-sub('c__','', library_1$C)
library_1$O<-sub('o__','', library_1$O)
library_1$Fa<-sub('f__','', library_1$Fa)
library_1$G<-sub('g__','', library_1$G)
library_1$S<-sub('s__','', library_1$S)

library_2$K <- sub('k__', '', library_2$K)
library_2$P <- sub('p__', '', library_2$P)
library_2$C <- sub('c__', '', library_2$C)
library_2$O <- sub('o__', '', library_2$O)
library_2$Fa <- sub('f__', '', library_2$Fa)
library_2$G <- sub('g__', '', library_2$G)
library_2$S <- sub('s__', '', library_2$S)


# add missing columns to library_2 
library(tibble)

library_2<- add_column(library_2, barcode8 = NA, .after = 'barcode07.csv')
library_2<- add_column(library_2, barcode35 = NA, .after = 'barcode34.csv')
library_2<- add_column(library_2, barcode36 = NA, .after = 'barcode35')


##### 2. Import: IDs and env. variables ######

lib_1_ids_full<- read.csv('lib_1_ids_full.csv', na.strings = c("", "NA"))
lib_2_ids_full<- read.csv('lib_2_ids_full.csv', na.strings = c("", "NA"))

# make all site names 2 characters
lib_2_ids_full$Site<- gsub("RMT","RM", lib_2_ids_full$Site)
lib_2_ids_full$Site<- gsub("GMT","GM", lib_2_ids_full$Site)
lib_2_ids_full$Site_sample<- gsub("RMT","RM", lib_2_ids_full$Site_sample)
lib_2_ids_full$Site_sample<- gsub("GMT","GM", lib_2_ids_full$Site_sample)


# combine sample ids into one column 
lib_1_ids_full$sample_id<- paste0(lib_1_ids_full$Site_sample, '_', lib_1_ids_full$Tree,'_', lib_1_ids_full$sample_location)
lib_2_ids_full$sample_id<- paste0(lib_2_ids_full$Site_sample, '_', lib_2_ids_full$Tree,'_', lib_2_ids_full$sample_location)

# assign library 1 and 2 ids to headings 
colnames(library_1)[8:103]<-lib_1_ids_full$sample_id
colnames(library_2)[8:103]<- lib_2_ids_full$sample_id


# merge libraries
merged_full<- merge(library_1, library_2, all.x = TRUE, all.y= TRUE)


# merge the ID and variables into one df and maintain correct order 
#install.packages('plyr')

library(plyr)
merged_ids_full<- join(lib_1_ids_full[,4:13],lib_2_ids_full[,4:13], by=NULL, type= 'full')


## Counts for report 
m#erged_full[,8:ncol(merged_full)]<- lapply(merged_full[, 8:ncol(merged_full)], as.numeric) 
total_reads<- sum(colSums(merged_full[,8:ncol(merged_full)]))


##### 4. Remove columns and rows with 0 reads ##### 

# remove columns and rows with zero reads 
merged_full <- merged_full %>% replace(is.na(.), 0)

zero<- which(colSums(merged_full[,8:199])==0)
zero2<-which(rowSums(merged_full[,8:199])==0)

# appears to be a species with no reads ?!

merged_full<- merged_full[-zero2,-c(8:199)[zero]]
  # checked and worked

# remove erroneous C5D5.E (column 95)

merged_full<- merged_full[,-95]


##### 5. Cleaning: remove rare species, controls and species for which the genus is unknown #####

# remove rare species 
  # see the quantile values for the sum of the rows to give a good idea of the cutoff point 
quantile(rowSums(merged_full[,8:195]))

  # Remove species with too few reads 
U50_reads_sp<- which(rowSums(merged_full[,8:195])<=50)

  #count how many times species occurs 
library(vegan)
U10_occur<- which(specnumber(merged_full[,8:195])<10)

low_ab<- unique(c(U50_reads_sp,U10_occur))

no_rare_sp <- merged_full[-low_ab,]

# remove the controls 

  # define treatment 
colnames(no_rare_sp)
treatment<- unlist(lapply(colnames(no_rare_sp),function(x)(substr(x,start=1, stop=1))))

  # select non samples 
treatment<- treatment[8:195]
controls<- which(treatment=='b'| treatment== "P" | treatment=="N" | treatment =="L") 

  # remove non samples 
samples<- c(8:195)[-controls]
sp_cleaner<- no_rare_sp[,samples]

  # add taxonomic information columns back 
no_controls<- cbind(no_rare_sp[,1:7], sp_cleaner)

# remove insertae sedis 
no_is<- no_controls[-grep("Incertae_sedis", no_controls$G),]

##### 6. Run NMDS to find outliers #####

# remove empty rows 
zero<- which(colSums(no_is[,8:162])==0)
zero2<-which(rowSums(no_is[,8:162])==0)

no_is<- no_is[,-c(8:162)[zero]]

# make community matrix 
clean_ab<- no_is[,8:ncol(no_is)]

# turn abundances into matrix 
clean_matrix<- as.matrix(clean_ab)

# run and plot nmds 
nmds= metaMDS(t(clean_matrix), distance= 'bray')
nmds
plot(nmds)


# plot by sample type 
colnames(no_is)
sample_type<- unlist(lapply(colnames(no_is),function(x)(substr(x,start=3, stop=3))))

sample_type<- sample_type[8:160]

colour_type<- as.factor(sample_type)
length(levels(colour_type))

col_funct<- rainbow(3)

levels(colour_type)<- col_funct
plot(nmds$points[,1], nmds$points[,2], pch= 21, bg=colour_type, cex=2)

outliers1<- which(nmds$points[,1]>1.6)
outliers2<- which(nmds$points[,2]>2)

# unqiue in case same point selected by both subsets 
Outliers<-unique(c(outliers1,outliers2))

# plot again to check that all outliers are dealt with properly 
plot(nmds$points[-Outliers,1], nmds$points[-Outliers,2], pch= 21, bg=colour_type[-Outliers], cex=2)

##### 7. Remove those outliers #####

# remove outliers from no_is 
outlierless<- c(8:160)[-Outliers]

outliers_removed<- no_is[,outlierless]

# add taxonomic information columns back 
no_outliers<- cbind(no_is[,1:7], outliers_removed)

# re-plot to check the outliers were removed correctly: all seems to be fine 
  # make community matrix 
    # clean_ab<- no_outliers[,8:ncol(no_outliers)]
  # turn abundances into matrix 
    # clean_matrix<- as.matrix(clean_ab)
  # run and plot nmds 
    # nmds= metaMDS(t(clean_matrix), distance= 'bray')
    # nmds
    # plot(nmds)
  # plot by sample type 
    # colnames(no_outliers)
    # sample_type<- unlist(lapply(colnames(no_outliers),function(x)(substr(x,start=3, stop=3))))
    # sample_type<- sample_type[8:151]
    # colour_type<- as.factor(sample_type)
    # length(levels(colour_type))
    # col_funct<- rainbow(3)
    # levels(colour_type)<- col_funct
    # plot(nmds$points[,1], nmds$points[,2], pch= 21, bg=colour_type, cex=2)


##### 8. Separate this into fungi and non fungi species #####

fungi<- no_outliers[(no_outliers$K %in% "Fungi"),]
nonfungi<- no_outliers[no_outliers$K != "Fungi", ]

##### 9. Add environmental variables ##### 
Names<- paste0(merged_ids_full$Site_sample, "_", merged_ids_full$Tree, "_", merged_ids_full$sample_location)

# add environmental variables to df with all species 
sp_full<- rbind(NA, NA, NA, NA, NA, NA, NA, NA, NA, no_outliers)


for (i in 8:ncol(sp_full)){
  #if (i==15| i==26|i==27|i==28|)(next())
  x<- which(Names==colnames(no_outliers)[i])
  sp_full[1,i]<- merged_ids_full$Site[x]
  sp_full[2,i]<- merged_ids_full$sample_type[x]
  sp_full[3,i]<- merged_ids_full$Tree[x]
  sp_full[4,i]<- merged_ids_full$sample_location[x]
  sp_full[5,i]<- merged_ids_full$DBH[x]
  sp_full[6,i]<- merged_ids_full$Length[x]
  sp_full[7,i]<- merged_ids_full$Decay_stage[x]
  sp_full[8,i]<- merged_ids_full$Avg_soil_moisture[x]
  sp_full[9,i]<- merged_ids_full$Avg_soil_pH[x]
  
  
}

sp_full$K[(1:9)]<- c('site', 'sample_type', 'tree', 'sample_loc', 'DBH', 'length', 'decay_stage', 'avg_soil_moist', 'avg_soil_ph')


# add environmental variables to df with fungi only 
sp_fungi<- rbind(NA, NA, NA, NA, NA, NA, NA, NA, NA, fungi)


for (i in 8:ncol(sp_full)){
  #if (i==15| i==26|i==27|i==28|)(next())
  x<- which(Names==colnames(no_outliers)[i])
  sp_fungi[1,i]<- merged_ids_full$Site[x]
  sp_fungi[2,i]<- merged_ids_full$sample_type[x]
  sp_fungi[3,i]<- merged_ids_full$Tree[x]
  sp_fungi[4,i]<- merged_ids_full$sample_location[x]
  sp_fungi[5,i]<- merged_ids_full$DBH[x]
  sp_fungi[6,i]<- merged_ids_full$Length[x]
  sp_fungi[7,i]<- merged_ids_full$Decay_stage[x]
  sp_fungi[8,i]<- merged_ids_full$Avg_soil_moisture[x]
  sp_fungi[9,i]<- merged_ids_full$Avg_soil_pH[x]
  
  
}

sp_fungi$K[(1:9)]<- c('site', 'sample_type', 'tree', 'sample_loc', 'DBH', 'length', 'decay_stage', 'avg_soil_moist', 'avg_soil_ph')

# add environmental variables to df with nonfungi only 
sp_nonfungi<- rbind(NA, NA, NA, NA, NA, NA, NA, NA, NA, nonfungi)


for (i in 8:ncol(sp_full)){
  #if (i==15| i==26|i==27|i==28|)(next())
  x<- which(Names==colnames(no_outliers)[i])
  sp_nonfungi[1,i]<- merged_ids_full$Site[x]
  sp_nonfungi[2,i]<- merged_ids_full$sample_type[x]
  sp_nonfungi[3,i]<- merged_ids_full$Tree[x]
  sp_nonfungi[4,i]<- merged_ids_full$sample_location[x]
  sp_nonfungi[5,i]<- merged_ids_full$DBH[x]
  sp_nonfungi[6,i]<- merged_ids_full$Length[x]
  sp_nonfungi[7,i]<- merged_ids_full$Decay_stage[x]
  sp_nonfungi[8,i]<- merged_ids_full$Avg_soil_moisture[x]
  sp_nonfungi[9,i]<- merged_ids_full$Avg_soil_pH[x]
  
  
}

sp_nonfungi$K[(1:9)]<- c('site', 'sample_type', 'tree', 'sample_loc', 'DBH', 'length', 'decay_stage', 'avg_soil_moist', 'avg_soil_ph')

##### 10. Download as csv files #####

write.csv(sp_full,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

write.csv(sp_fungi,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_fungi.csv')
sp_fungi<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_fungi.csv')
sp_fungi<- sp_fungi[-c(1)]

write.csv(sp_nonfungi,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_nonfungi.csv')
sp_nonfungi<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_nonfungi.csv')
sp_nonfungi<- sp_nonfungi[-c(1)]



# note: remember to use the below line when you load it in to get rid of the annoying extra column 
df<- df[,-c(1)]



# remove zeros
write.csv(sp_soil,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_soil.csv')
zero<- which(colSums(sp_soil[,8:15])==0)
zero2<-which(rowSums(sp_soil[,8:15])==0)

sp_soil<- sp_soil[-zero2,]
write.csv(sp_soil,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/soil_sp.csv')

##### 12. GHG data #####

GHG<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/GHG_emissions.csv')
GHG$sample_ID<- paste0(GHG$Site, 'D_', GHG$Tree,'_', GHG$sample_location, '_', GHG$Repeat)


##### 2. + and - controls #####

# need unclean data (merged_full)

# subset positive and negative controls 


# clean in the same way up until removing controls 
##### 4. Remove columns and rows with 0 reads ##### 
library(tibble)
# remove columns and rows with zero reads 
merged_full <- merged_full %>% replace(is.na(.), 0)

zero<- which(colSums(merged_full[,8:199])==0)
zero2<-which(rowSums(merged_full[,8:199])==0)

# appears to be a species with no reads ?!

merged_full<- merged_full[-zero2,-c(8:199)[zero]]
# checked and worked

# remove erroneous C5D5.E (column 95)

merged_full<- merged_full[,-95]


##### 5. Cleaning: remove rare species, controls and species for which the genus is unknown #####

# remove rare species 
# see the quantile values for the sum of the rows to give a good idea of the cutoff point 
quantile(rowSums(merged_full[,8:195]))

# Remove species with too few reads 
U50_reads_sp<- which(rowSums(merged_full[,8:195])<=50)

#count how many times species occurs 
library(vegan)
U10_occur<- which(specnumber(merged_full[,8:195])<10)

low_ab<- unique(c(U50_reads_sp,U10_occur))

no_rare_sp <- merged_full[-low_ab,]

no_is<- no_controls[-grep("Incertae_sedis", no_controls$G),]

# positve 

positives <- grep("PC", colnames(no_is))
negatives <- grep("NC", colnames(no_is))

PC_lib<- cbind(merged_full[-low_ab,1:7],merged_full[-low_ab,positives])
NC_lib<- cbind(merged_full[-low_ab,1:7],merged_full[-low_ab,negatives])


# save as csvs
write.csv(PC_lib,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/PC_lib.csv')
write.csv(NC_lib,'C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/NC_lib.csv')

