# read in data 
# calculate emissions 
# emissions against decay stage 
# emissions against diversity 


# load in 
ghg<- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/GHG_emissions.csv')
summary(ghg)

##### Find emissions #####

# find emissions/day/g 

ghg$CO2_emissions<- ((ghg$CO2_fin - ghg$atmos_CO2_fin)-(ghg$CO2_0 - ghg$atmos_CO2_0))/(2.5*14)
ghg$CH4_emissions<- ((ghg$CH4_fin - ghg$atmos_CH4)-(ghg$CH4_0 - ghg$atmos_CH4))/(2.5*14)

ghg$CO2_emissionsnorm<- ((ghg$CO2_emissions-5.714286)*0.000998859)/0.01
ghg$CH4_emissionsnorm<- ((ghg$CH4_emissions-0.0059142857)*0.000998859)/0.01

# plot 
plot(ghg$CO2_emissions)
plot(ghg$CH4_emissions)
# histogram 
hist(ghg$CO2_emissions, breaks=10)



#plot co2 by site 
boxplot(ghg$CO2_emissions~ghg$Site)
boxplot(ghg$CH4_emissions~ghg$Site)

##### Add decay stage #####
# Create Tree ID for matching
Tree_id <- paste0(ghg$Site, "D_", ghg$Tree, "_B")

# Import full data frame
sp_full <- read.csv('C:/Users/isobe/OneDrive/Documents/Part II/Project/Part II Analysis/sp_full.csv')
sp_full <- sp_full[, -1]  # Remove first column (if it's an index)

# Make a copy of sp_full
sp_ghg <- sp_full

# Add a new column to `ghg` for decay stage (initialize as NA)
ghg$decay_stage <- NA  

# Loop through `ghg` and match `Tree_id` with column names in `sp_ghg`
for (i in 1:nrow(ghg)) {
  match_index <- match(Tree_id[i], colnames(sp_ghg)[8:ncol(sp_ghg)])  # Find matching column index
  
  if (!is.na(match_index)) {  # Ensure a match is found
    actual_col <- match_index + 7  # Adjust index to match format
    ghg$decay_stage[i] <- sp_ghg[7, actual_col]  
  } else {
    warning(paste("No match found for Tree ID:", Tree_id[i]))
  }
}

ghg$decay_stage<- as.numeric(ghg$decay_stage)

##### Plot emissions against decay stage #####

plot(ghg$CO2_emissions~ghg$decay_stage)


model<- lm(ghg$CO2_emissionsnorm~ghg$decay_stage)
summary(model)

# find 


# Scatter plot of CO2 emissions vs. decay stage
plot(ghg$CO2_emissions ~ ghg$decay_stage, 
     main = expression("CO"[2] ~ "Emissions vs. Decay Stage"), 
     xlab = "Decay Stage", 
     ylab = expression("CO"[2] ~ "Emissions (gC day"^{-1}~"gD"^{-1}*")"), 
     pch = 16, col = "black")

# Fit linear model
model <- lm(ghg$CO2_emissions ~ ghg$decay_stage)

# Add regression line in red
abline(model, col = "red", lwd = 2)

### need to add confidence interval ###

# Scatter plot
svg("GHG_2.svg", width = 8, height = 8)
par(mar=c(6,6,5,5))
plot(ghg$CO2_emissionsnorm ~ ghg$decay_stage, 
     main = expression("CO"[2] ~ "Emissions vs. Decay Stage"), 
     xlab = "Decay Stage", 
     ylab = expression("CO"[2] ~ "Emissions (g C day"^{-1}~"g Deadwood"^{-1}*")"), 
     pch = 16, col = "black", cex.lab=2, cex.axis=2)

# Fit linear model
model <- lm(CO2_emissionsnorm ~ decay_stage, data = ghg)

# Add regression line
abline(model, col = "red", lwd = 2)

# Add confidence interval
new_data <- data.frame(decay_stage = seq(1, 
                                         5, 
                                         length.out = 100))

conf_int <- predict(model, newdata = new_data, interval = "confidence")

# Add shaded confidence interval
lines(new_data$decay_stage, conf_int[, "lwr"], col = rgb(0.8, 0.8, 0.8, 0.5), lty = 2)
lines(new_data$decay_stage, conf_int[, "upr"], col = rgb(0.8, 0.8, 0.8, 0.5), lty = 2)
polygon(c(new_data$decay_stage, rev(new_data$decay_stage)),
        c(conf_int[, "lwr"], rev(conf_int[, "upr"])),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

# Add the regression line again for clarity
lines(new_data$decay_stage, conf_int[, "fit"], col = "red", lwd = 2)

dev.off()
# Display model summary
summary(model)

# test 
model<- lm(ghg$CH4_emissions~ghg$decay_stage)
summary(model)


