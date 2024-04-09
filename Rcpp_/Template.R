library(TransPhylo)
set.seed(0)
neg=100/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.25
simu <- simulateOutbreak(neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
                         w.scale=w.scale,dateStartOutbreak=2005,dateT=2008)

ptree<-extractPTree(simu)
value1 <- 0.1
value2 <- 0.9
value3 <- 0.1
threshold1 <-  2000.0;
threshold2 <-  2000 + 0.6;


ptree<-extractPTree(simu)
value1 <- 0.1
value2 <- 0.9
value3 <- 0.1
threshold1 <-  2000.0;
threshold2 <-  2000 + 0.6;


# res<-inferTTree(ptree,mcmcIterations=10000,w.shape=w.shape,w.scale=w.scale,dateT=2020, delta_t=0.01, isTp = 1,
#                 Pi = c(value1,value2,value3,threshold1,threshold2), updatePi = FALSE, updateNeg =FALSE)

res<-inferTTree(ptree,mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=2020, delta_t=0.01, isTp = 1,
                Pi = c(value1,value2,value3,threshold1,threshold2), updatePi = FALSE, updateNeg =FALSE)


# Define your time points and sampling probabilities
time_points = c(1, 2) # replace time1, time2 with actual time values
sampling_probs = c(.5, .3) # replace prob1, prob2 with actual sampling probabilities



# Combine them into a list with named elements
Pi_values <- list(
  time = time_points,
  probs = sampling_probs
)

# Now you can access time and probs with Pi_values$time and Pi_values$probs
# ...

# Assuming you need to pass them as a concatenated vector to the Pi argument
Pi_argument <- c(Pi_values$time, Pi_values$probs, Pi_values$threshold1, Pi_values$threshold2)

# Then use this vector in your function call
res <- inferTTree(ptree, mcmcIterations=100, w.shape=w.shape, w.scale=w.scale, dateT=2011, 
                  delta_t=0.01, isTp = 0, time_data = c(2,3,66), prob_data = c(20,30,666), updatePi = FALSE, updateNeg =FALSE)

# check tinf and date T the time interval will between these two values the grid 
# I put function that shows the value of tinf it can be checked before inference 


res <- inferTTree(ptree, mcmcIterations=100, w.shape=w.shape, w.scale=w.scale, dateT=2011.303, 
                  delta_t=0.01, isTp = 2, time_data = c(2,3,66), prob_data = c(20,30,666), updatePi = FALSE, updateNeg =FALSE)



res <- inferTTree(ptree, mcmcIterations=100, w.shape=w.shape, w.scale=w.scale, dateT=2011, 
                  delta_t=0.01, isTp = 3, time_data = c(2,3,66), prob_data = c(20,30,666), updatePi = FALSE, updateNeg =FALSE)


create_grid_and_values <- function(dateT, tinf, delta_t) {
  n <- round((dateT - tinf) / delta_t)
  grid <- seq(from = tinf, by = delta_t, length.out = n)
  values <- rep(0.7, length(grid))
  return(list(grid = grid, values = values))
}

# Example usage:
dateT <- 2011  # Example current date in numeric format
tinf <- 2000               # Example start date, 10 days before dateT
delta_t <- 0.01                   # Example interval between grid points

# Generate the grid and values
result <- create_grid_and_values(dateT, tinf, delta_t)

# The grid array
grid <- result$grid
print(grid)

# The values array with 0.7
values <- result$values
print(values)



res <- inferTTree(ptree, mcmcIterations=1000, w.shape=w.shape, w.scale=w.scale, dateT=2011, 
                  delta_t=0.01, isTp = 4, time_data = result$grid, prob_data = result$values, updatePi = FALSE, updateNeg =FALSE)
# Plotting the grid against the constant values
plot(grid, values, type = 'l', col = 'blue', xlab = 'Grid', ylab = 'Values', main = 'Grid vs. Constant Values')
# this is the plot for this section
















