
##### 1. Base vs non base boxplot #####
# subset deadwood 
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[-c(1)]

# subset to libraries only
lib_full<- sp_full[10:nrow(sp_full),]

# set as numeric 
lib_full[,8:ncol(lib_full)]<- lapply(lib_full[, 8:ncol(lib_full)], as.numeric) 

# sample type 
sample_type<- unlist(lapply(colnames(lib_full),function(x)(substr(x,start=3, stop=3))))

deadwood<- which(sample_type== 'D')

lib_deadwood <- lib_full[,deadwood]

# define sample location and treatment 
sample_loc<- unlist(lapply(colnames(lib_deadwood),function(x)(substr(x,start=7, stop=7))))
sites<- unlist(lapply(colnames(lib_deadwood),function(x)(substr(x,start=1, stop=2))))
sample_loc<- as.factor(sample_loc)

# plot
library(vegan)
dw_div<- diversity(t(lib_deadwood), index= 'shannon')

# model and test model 
model<- summary(lm(dw_div~sample_loc))
boxplot(exp(dw_div)~sample_loc)

#plot and save plot 
svg("boxplot_basenonbase.svg", width = 8, height = 6)

par(mar = c(6, 6, 4, 2))
boxplot(exp(dw_div) ~ sample_loc, 
        data = lib_deadwood, 
        xlab = "Sample Location", 
        ylab = "ENS (no. of species)", 
        main = "Boxplot of ENS by Sample Location", 
        names = c("Base", "Non-Base"), 
        cex.axis = 1.5, 
        cex.lab = 1.8,  
        cex.main = 2,   
        col = "lightblue", 
        border = "black")

dev.off()


##### 2. Species accumulation curves #####
sp_full<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full<- sp_full[,-c(1)]

# subset to library and make numeric 
ab_full<- sp_full[10:3038, 8:152]

for (i in 1:ncol(ab_full)){
  x<- as.numeric(ab_full[,i])
  ab_full[,i]<-x 
}


# plot number of samples on x axis and number species (colsums =/= 0) on y axis 
library(vegan)
t_ab_full<- t(ab_full)

# repeat with subsets 
deadwood<- which(sample_type =='D')
ab_dw<- ab_full[, deadwood]
t_ab_dw<- t(ab_dw)

living<- which(sample_type =='L')
ab_li<- ab_full[, living]
t_ab_li<- t(ab_li)

soil<- which(sample_type =='S')
ab_s<- ab_full[, soil]
t_ab_s<- t(ab_s)

sp1<- specaccum(t_ab_s, method ='random')
plot(sp1, xlab="Number of Samples", ylab="Species Richness",
     main="Soil Species Accumulation Curve", col="blue",
     cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2)

svg("species_accumulation_curves.svg", width = 15, height = 10)
par(mfrow = c(2, 2))

# Soil SAC
sp1 <- specaccum(t_ab_s, method = 'random')
plot(sp1, xlab = "Number of Samples", ylab = "Species Richness",
     main = "Soil Species Accumulation Curve", col = "blue", 
     cex.main = 2.4, cex.lab = 1.6, cex.axis = 1.5)

# Living SAC
sp2 <- specaccum(t_ab_li, method = 'random')
plot(sp2, xlab = "Number of Samples", ylab = "Species Richness",
     main = "Living Species Accumulation Curve", col = "blue", 
     cex.main = 2.4, cex.lab = 1.6, cex.axis = 1.5)

# Deadwood SAC
sp3 <- specaccum(t_ab_dw, method = 'random')
plot(sp3, xlab = "Number of Samples", ylab = "Species Richness",
     main = "Deadwood Species Accumulation Curve", col = "blue", 
     cex.main = 2.4, cex.lab = 1.6, cex.axis = 1.5)
dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))
