

# res_TP <- inferTTree(ptree_c4007, mcmcIterations=5000,w.mean=w.mean,w.std=w.std,ws.mean=ws.mean,ws.std=ws.std, dateT=2026, 
#                      delta_t=0.01, isTp = 1, time_data = result$grid, prob_data = result$values, updatePi = FALSE, updateNeg =FALSE)
set.seed(42)

record1<-inferTTree(ptree, 
                    w.shape = 1.3, 
                    w.scale = 3.333, 
                    ws.shape = 1.1,
                    ws.scale = 2.5, 
                    mcmcIterations = 10000, thinning = 10,
                    startNeg = 100/365, 
                    startOff.r = 1, 
                    startOff.p = 0.5, 
                    # startPi = 0.8,
                    updateNeg = FALSE, updateOff.r = TRUE, updateOff.p = FALSE, 
                    isTp=1, 
                    delta_t=0.01,
                    time_data = result$grid, prob_data = result$values,
                    updatePi = FALSE,  
                    dateT= dateT, 
                    updateTTree=TRUE,
                    optiStart=2,verbose=F)

library(coda)
mcmc_TP=convertToCoda(record1)
effectiveSize(mcmc_TP)
saveRDS(record1, file = "mcmc_A4_C1_TP_10000_42.rds")
plot(record1)

plot(medTTree(record1))

a_TP=getIncidentCases(record1,show.plot = T)



record2<-inferTTree(ptree, 
                    w.shape = 1.3, 
                    w.scale = 3.333, 
                    ws.shape = 1.1,
                    ws.scale = 2.5, 
                    mcmcIterations = 10000, 
                    thinning = 10,
                    startNeg = 100/365, 
                    startOff.r = 1, startOff.p = 0.5, 
                    # startPi = 0.8,
                    updateNeg = FALSE, 
                    updateOff.r = TRUE, 
                    updateOff.p = FALSE, 
                    isTp=3, 
                    delta_t=0.01,
                    time_data = result$grid, prob_data = result$values,
                    updatePi = FALSE,  
                    dateT=dateT, 
                    updateTTree=TRUE,
                    optiStart=2,verbose=F)

library(coda)
mcmc_TP=convertToCoda(record2)
effectiveSize(mcmc_TP)
# saveRDS(record2, file = "mcmc_A4_C1_TPS_10000_42.rds")
plot(record2)

plot(medTTree(record2))

a_TP=getIncidentCases(record2,show.plot = T)


par(mfrow=c(1, 2))


plot(medTTree(record1)) 
plot(medTTree(record2))



par(mfrow=c(1, 2))


a_TP=getIncidentCases(record1,show.plot = T)
a_TP=getIncidentCases(record2,show.plot = T)

# a_TP=getIncidentCases(A4_C1,show.plot = T)

par(mfrow=c(1, 2))
mat=computeMatWIW(record1)
lattice::levelplot(mat,xlab='',ylab='')

mat=computeMatWIW(record2)
lattice::levelplot(mat,xlab='',ylab='')
