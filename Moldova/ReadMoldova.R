library(ape)
library(TransPhylo)

file_phylogenetic <-  "C:/Users/Omid/Documents/TransPhyloFork/Moldova/MCCtrees/A4_C1_with_dates.mcc"

setwd("C:/Users/Omid/Documents/TransPhyloFork/Moldova")
A4_C1 <- readRDS("A4_C1.RData")


data <- read.nexus(file_phylogenetic)

# Example usage with a vector of dates
date_vector <- extract_and_convert_date(data$tip.label)  # replace with your actual dates
decimal_years <- convertDatesToDecimalYears(date_vector)
print(decimal_years)
max(decimal_years)

# 1) when you feed the time tinf the data_time first entry should be 
#    earlier than the tinf chosen for the root of tree achievable through 
#    is TP=2 
# 2) DateT is very essential the all of other slicing fo time will be according to
#    DateT so the idea is to pick a floor rounding of the dateLastSample 



create_grid_and_values <- function(start_year, end_year, delta_t, default_value, specific_values) {
  # Create the grid of years with step size delta_t
  grid <- seq(from = start_year, to = end_year, by = delta_t)
  
  # Initialize values with the default value
  values <- rep(default_value, length(grid))
  
  # Set specific values for the years provided in the specific_values list
  for(i in seq_along(grid)) {
    # Find the year part of the grid value
    year_part <- floor(grid[i])
    
    # Check if this year part has a specific value
    if(as.character(year_part) %in% names(specific_values)) {
      values[i] <- specific_values[[as.character(year_part)]]
    }
  }
  
  return(list(grid = grid, values = values))
}

# Define specific values from the table
specific_values <- list(`2018` = 0.8, `2019` = 0.8, `2020` = 0.8, `2021` = 0.8, `2022` = 0.8)
# Example usage:
start_year <- 2000  # Starting year
end_year <- ceiling(max(decimal_years))+5 # Ending year
delta_t <- 0.01        # Step size (1 year)
default_value <- 0.05 # Default value for years outside 2019-2023

# Generate the grid and values
result <- create_grid_and_values(start_year, end_year, delta_t, default_value, specific_values)

# The grid array
grid <- result$grid
# print(grid)

# The values array
values <- result$values
# print(values)

# Plot the values with specified type and color
plot(grid, values, type = 'l', col = 'blue', xlab = 'Grid', ylab = 'Values', main = 'Grid vs. Constant Values', xaxt='n', ylim=c(min(values) - 0.1, max(values) + 0.1))

# Add x-axis ticks for every year
axis(1, at=seq(floor(min(grid)), ceiling(max(grid)), by=1), las=2)

# Add grid lines for y-axis
abline(h=seq(floor(min(values)), ceiling(max(values)), by=0.1), col="gray", lty="dotted")
abline(v=seq(floor(min(grid)), ceiling(max(grid)), by=1), col="gray", lty="dotted")








ptree<-ptreeFromPhylo(data,dateLastSample=ceiling(max(decimal_years)))
plot(ptree)


w.shape=1.5
w.scale=0.1

dateT=2022

set.seed(1)
resTPS<-inferTTree(ptree,mcmcIterations=10000,w.shape=w.shape,w.scale=w.scale
                , isTp =3,
                time_data = result$grid,
                prob_data = result$values,
                delta_t = 0.01,
                updatePi = FALSE,
                dateT=dateT)
# res<-inferTTree(ptree,mcmcIterations=1000,w.shape=w.shape,w.scale=w.scale, isTp = 1,dateT=dateT)
# medTPS=medTTree(resTPS)
# plot(medTPS)



# a=getIncidentCases(resTPS,show.plot = T)
library(coda)
mcmc=convertToCoda(resTPS)
effectiveSize(mcmc)


resTP<-inferTTree(ptree,mcmcIterations=10000,w.shape=w.shape,w.scale=w.scale
                   , isTp =1,
                   time_data = result$grid,
                   prob_data = result$values,
                   delta_t = 0.01,
                   updatePi = FALSE,
                   dateT=dateT)

# a=getIncidentCases(resTPS,show.plot = T)
library(coda)
mcmc=convertToCoda(resTP)
effectiveSize(mcmc)
# medTP =medTTree(resTP)
# plot(medTP)

# a=getIncidentCases(resTP,show.plot = T)

par(mfrow=c(1, 2))

medTP=medTTree(resTP)
medTPS = medTTree(resTPS)
plot(medTP)
plot(medTPS)


par(mfrow=c(1, 2))
a=getIncidentCases(resTP,show.plot = T)
a=getIncidentCases(resTPS,show.plot = T)



# med=medTTree(A4_C1)
# plot(med)
mat=computeMatWIW(A4_C1)
lattice::levelplot(mat,xlab='',ylab='')

mat=computeMatWIW(resTP)
lattice::levelplot(mat,xlab='',ylab='')

mat=computeMatWIW(resTPS)
lattice::levelplot(mat,xlab='',ylab='')


