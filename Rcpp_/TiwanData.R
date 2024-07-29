#The MCMC procedure to infer the transmission tree given the phylogenetic tree!

set.seed(1)
setwd("C:/Users/Omid/Documents/TransPhyloFork")
library(TransPhylo)
# w.shape = 1.3, w.scale = 3.333, ws.shape = 1.1,
# ws.scale = 2.5

# w.mean= 1.5
# w.std= 3.333
# 
# ws.mean= 1.1
# ws.std= 2.5


w.mean= 1.5
w.std= 1

ws.mean= 1.5
ws.std= 1

pi = 0.7 

library(ape)
tree_c4007<-read.nexus("c4007_n24_original_MCC(1).tree")
tree_c4007$tip.label= substr(tree_c4007$tip.label,1,7)

ptree_c4007=ptreeFromPhylo(tree_c4007,dateLastSample =2023.14 )
plot(ptree_c4007)



# Define specific values from the table
specific_values <- list(`2019` = 0.8, `2020` = 0.8, `2021` = 0.8, `2022` = 0.8, `2023` = 0.8)

# Example usage:
start_year <- 2000  # Starting year: use isTP=2 to get idea, just prior to the first infection
end_year <- 2025    # Ending year
delta_t <- 0.01        # Step size (1 year)
default_value <- 0.05 # Default value for years outside 2019-2023


# Generate the grid and values
result <- create_grid_and_values(start_year, end_year, delta_t, default_value, specific_values)
# The grid & values array
grid <- result$grid
values <- result$values
# Plot the values with specified type and color
plot(grid, values, type = 'l', col = 'blue', xlab = 'Grid', ylab = 'Values', main = 'Grid vs. Constant Values', xaxt='n', ylim=c(min(values) - 0.1, max(values) + 0.1))
axis(1, at=seq(floor(min(grid)), ceiling(max(grid)), by=1), las=2)
abline(h=seq(floor(min(values)), ceiling(max(values)), by=0.1), col="gray", lty="dotted")
abline(v=seq(floor(min(grid)), ceiling(max(grid)), by=1), col="gray", lty="dotted")


res_TP <- inferTTree(ptree_c4007, mcmcIterations=50000,w.mean=w.mean,w.std=w.std,ws.mean=ws.mean,ws.std=ws.std, dateT=2025, 
                     delta_t=0.01, isTp = 1, time_data = c(2019), prob_data = c(0.05, 0.8), updatePi = FALSE, updateNeg =FALSE)

library(coda)
mcmc_TP=convertToCoda(res_TP)
effectiveSize(mcmc_TP)

plot(res_TP)


a_TP=getIncidentCases(res_TP,show.plot = T)


# res_TP <- inferTTree(ptree_c4007, mcmcIterations=1000,w.mean=w.mean,w.std=w.std,ws.mean=ws.mean,ws.std=ws.std, dateT=2026, delta_t=0.01, isTp = 2, time_data = result$grid, prob_data = result$values, updatePi = FALSE, updateNeg =FALSE)
res_TPS <- inferTTree(ptree_c4007, mcmcIterations=50000,w.mean=w.mean,w.std=w.std,ws.mean=ws.mean,
                      ws.std=ws.std, dateT=2025, delta_t=0.01, isTp = 2,time_data = c(2019), prob_data = c(0.05, 0.8), updatePi = FALSE, updateNeg =FALSE)
#*dateT needs to be the same as the last element in time_data, later is fine

mcmc_TPS=convertToCoda(res_TPS)
effectiveSize(mcmc_TPS)

plot(res_TPS)

# plot(medTTree(res_TPS)) 

# get incident in this example
a_TPS=getIncidentCases(res_TPS,show.plot = T)


par(mfrow=c(1, 2))

# med1TP = medTTree(res_TP)
# med1TPS = medTTree(res_TPS)
# plot(med1TP) 
# plot(med1TPS) 

In_TP = getIncidentCases(res_TP,show.plot = T)
In_TPS = getIncidentCases(res_TPS,show.plot = T)



mcmc_TP=convertToCoda(res_TP)
effectiveSize(mcmc_TP)

mcmc_TP=convertToCoda(res_TPS)
effectiveSize(mcmc_TPS)


sum(In_TP$unsampCases)
sum(In_TPS$unsampCases)














