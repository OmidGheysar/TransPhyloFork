library(TransPhylo)
set.seed(4)
neg=100/365
off.r=1.1
w.shape=10
w.scale=0.1
pi=0.8
set.seed(42)
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2018,dateT=2024)
plot(simu)

ptree<-extractPTree(simu)



res <- inferTTree(ptree, mcmcIterations=1000, w.shape=w.shape, w.scale=w.scale, dateT=2025, 
                  delta_t=0.01, isTp = 2, time_data = c(2,3,66), prob_data = c(20,30,666), updatePi = FALSE, updateNeg =FALSE)



res <- inferTTree(ptree, mcmcIterations=100, w.shape=w.shape, w.scale=w.scale, dateT=2025, 
                  delta_t=0.01, isTp = 3, time_data = c(2,3,66), prob_data = c(20,30,666), updatePi = FALSE, updateNeg =FALSE)






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

# ... (rest of your code)

# Define specific values from the table
specific_values <- list(`2019` = 0.8, `2020` = 0.8, `2021` = 0.75, `2022` = 0.7, `2023` = 0.33)

# Example usage:
start_year <- 2010  # Starting year
end_year <- 2025    # Ending year
delta_t <- 0.01        # Step size (1 year)
default_value <- 0.05 # Default value for years outside 2019-2023

# Generate the grid and values
result <- create_grid_and_values(start_year, end_year, delta_t, default_value, specific_values)

# The grid array
grid <- result$grid
print(grid)

# The values array
values <- result$values
print(values)





# Plot the values with specified type and color
plot(grid, values, type = 'l', col = 'blue', xlab = 'Grid', ylab = 'Values', main = 'Grid vs. Constant Values', xaxt='n', ylim=c(min(values) - 0.1, max(values) + 0.1))

# Add x-axis ticks for every year
axis(1, at=seq(floor(min(grid)), ceiling(max(grid)), by=1), las=2)

# Add grid lines for y-axis
abline(h=seq(floor(min(values)), ceiling(max(values)), by=0.1), col="gray", lty="dotted")
abline(v=seq(floor(min(grid)), ceiling(max(grid)), by=1), col="gray", lty="dotted")


# 1) when you feed the time tinf the data_time first entry should be 
#    earlier than the tinf chosen for the root of tree achievable through 
#    is TP=2 
# 2) DateT is very essential the all of other slicing fo time will be according to
#    DateT so the idea is to pick a floor rounding of the dateLastSample 
# 3) BE mindful about delat_t=0.5

res_TP <- inferTTree(ptree, mcmcIterations=1000, w.shape=w.shape, w.scale=w.scale, dateT=2026, 
                     delta_t=0.01, isTp = 1, time_data = result$grid, prob_data = result$values, updatePi = FALSE, updateNeg =FALSE)

a=getIncidentCases(res_TP,show.plot = T)



res_TPS <- inferTTree(ptree, mcmcIterations=1000, w.shape=w.shape, w.scale=w.scale, dateT=2026, 
                      delta_t=0.01, isTp = 3, time_data = result$grid, prob_data = result$values, updatePi = FALSE, updateNeg =FALSE)


# plot(medTTree(res)) 
# get incident in this example!

a=getIncidentCases(res_TPS,show.plot = T)



































