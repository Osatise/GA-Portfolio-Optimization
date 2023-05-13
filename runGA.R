install.packages("quantmod")
install.packages("GA")
install.packages("parallel")
install.packages("doParallel")
 
library(quantmod)
library(GA)
library(parallel)
library(doParallel)


myStocks <- c("AAPL", "GOOG", "IBM", "MSFT", "TSLA", "NFLX", "AMZN", "NVDA", "AMD", "JPM")
#Number of assets in a portfolio
numStocks <- length(myStocks)

#Get the weekly returns for each asset
Stocks <- lapply(myStocks, function(sym) weeklyReturn(na.omit(getSymbols(sym, from="2018-01-01", to="2023-01-01", auto.assign=FALSE))))
#Merge the assets into 1 matrix
Stocks <- do.call(merge, Stocks)

#Set the lower and the upper limit of the generated solutions as well as the chromosome size
lower <- c(rep(0,numStocks))
upper <- c(rep(1,numStocks))

#Scale the chromosomes down to 1
scaleWeights <- function(x){
  return (x/sum(x))
}

#Fitness function
sharpe_ratio <- function(weights) {
  weights <- scaleWeights(weights)
  portfolio_return <- colSums(Stocks * weights)
  portfolio_sd <- sqrt(t(weights) %*% cov(Stocks) %*% weights)
  sharpe <- (portfolio_return - 0.0339) / portfolio_sd
  return(sharpe)
}

maxGenerations <<- 20

monitor <- function(obj){
  # gaMonitor(obj)                      #call the default gaMonitor to print the usual messages during evolution
  iter <- obj@iter                      #get the current iternation/generation number 
  if (iter <= maxGenerations){          #some error checking
    fitness <- obj@fitness              #get the array of all the fitness values in the present population
    #<<- assigns a value to the global variable declared outside the scope of this function.    
    thisRunResults[iter,1] <<- max(fitness)
    thisRunResults[iter,2] <<- mean(fitness)    
    thisRunResults[iter,3] <<- median(fitness)
    cat(paste("\rGA | generation =", obj@iter, "Mean =", thisRunResults[iter,2], "| Best =", thisRunResults[iter,1], "\n"))
    flush.console()
  }  
  else{                               #print error messages
    cat("ERROR: iter = ", iter, "exceeds maxGenerations = ", maxGenerations, ".\n")
    cat("Ensure maxGenerations == nrow(thisRunResults)")
  }
}

runGA <- function(noRuns = 30){
 
  #Set up what stats you wish to note.    
  statnames = c("best", "mean", "median")
  thisRunResults <<- matrix(nrow=maxGenerations, ncol = length(statnames)) #stats of a single run
  resultsMatrix = matrix(1:maxGenerations, ncol = 1)  #stats of all the runs
  
  resultNames = character(length(statnames)*noRuns)
  resultNames[1] = "Generation"
  
  bestFitness <<- -Inf
  bestSolution <<- NULL
  for (i in 1:noRuns){
    cat(paste("Starting Run ", i, "\n"))
  GA <- ga(type="real-valued",
              fitness = sharpe_ratio, maxiter = maxGenerations, popSize = 50,
              lower = lower,
              upper = upper,
              monitor=monitor,seed = 1)
  
  resultsMatrix = cbind(resultsMatrix, thisRunResults)
  
  if (GA@fitnessValue > bestFitness){
    bestFitness <<- GA@fitnessValue
    bestSolution <<- GA@solution
  }
  #Create column names for the resultsMatrix
  for (j in 1:length(statnames)) resultNames[1+(i-1)*length(statnames)+j] = paste(statnames[j],i)
  }
  colnames(resultsMatrix) = resultNames
  return (resultsMatrix)
}

getBestFitness<-function(){
  return(bestFitness)
}

#Best solution scaled to 1
getBestSolution<-function(){
  return(bestSolution/sum(bestSolution))
}
