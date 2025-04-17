########## This code generates data and contour plots for the current 
########## version of the paper

####Set the working directory####

setwd('C:/Users/tymur/OneDrive/Desktop/Work/Staking optimal r and t')

####Load the libraries####

library("nloptr")
library("pracma")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("glue")
library("writexl")
library("openxlsx")

####Utility functions####

tFunc <- function(TT){
  #
  #This function (f()) multiplies the annual staking rate making the
  #staking reward rate paid to be a function of time.
  #
  #INPUTS: - TT: time in years
  #
  #OUTPUT: - f(T) = T
  #
  
  return (TT)
}

objFunc <- function(r,T){
  #
  #This function calculates the value of the objective function optimized.
  #
  #INPUTS: - r: annual staking rewards rate
  #        - T: time in years
  #
  #OUTPUT: - r*f(T)
  #
  
  
  return (r*tFunc(T))}

evalPobaNoDefault <- function(r,T,params){
  #
  # This function calculates the probability of no-default by adopting the 
  # barrier option theory (For details see the appendix of the paper).
  #
  #INPUTS: - r: the annual staking reward rate
  #        - T: time in years
  #        - params: a vector of parameters
  #
  #OUTPUT: - proba(no_default)
  #
  
  #Unpack the parameters
  lam0 <- params[1]
  XX0 <- params[2]
  FF <- params[3]
  SS <- params[4]
  RR_0 <- params[5]
  sigmX <- params[6]
  sigmLambda <- params[7]
  tStar <- params[8]
  rho <- params[11]
  muLambda <- params[9]
  muX <- params[10]
  
  #Set the initial value for U
  U_0 <- lam0*XX0
  
  #Calculate the mu^{U}, simga^{U} and alpha
  mu_U <- muX+muLambda +sigmX*sigmLambda*rho
  sigm_U <- sqrt(sigmX^2+sigmLambda^2+sigmX*sigmLambda*rho)
  alpha <- (1/sigm_U)*(mu_U-0.5*sigm_U^2)
  
  #Calculate the value of the boundary
  m <- (1/sigm_U)*log((r*tFunc(tStar)*SS-RR_0)/(FF*U_0))
  
  #Calculate the constituent elements of the probability and sum them up
  proba1 <- pnorm(-(m+alpha*T)/sqrt(T), 0, 1)
  proba2 <- pnorm((m-alpha*T)/sqrt(T), 0, 1)
  
  A <- exp(2*alpha*m)*proba1 -  exp(2*alpha*m)
  B <- 1 - proba2
  
  #return ((1 - exp(2*alpha*m)+exp(2*alpha*m)*proba1-proba2))
  return (A + B)
  
}

upperR <- function(params){
  #
  #This function computes the upper bound on the annual staking reward rate for the
  #optimization.
  #
  #INPUTS: - params: vector of parameters
  #
  #OUTPUT: - value of the upper bound on the staking rewards rate
  #
  
  #Unpack the parameters
  lam0 <- params[1]
  XX0 <- params[2]
  FF <- params[3]
  SS <- params[4]
  RR_0 <- params[5]
  tStar <- params[8]
  
  #Set the initial value for U
  U_0 <- lam0*XX0
  
  return ((FF*U_0+RR_0)/(tFunc(tStar)*SS))
}

lowerR <- function(params){
  #
  #This function calculates the lower bound on the annual staking rewards rate for
  #the optimization.
  #
  #INPUTS: - params: vector of parameters
  #
  #OUTPUT: - value of the upper bound on the staking rewards rate
  #
  
  #Unpack the parameters
  lam0 <- params[1]
  XX0 <- params[2]
  FF <- params[3]
  SS <- params[4]
  RR_0 <- params[5]
  tStar <- params[8]
  
  return (RR_0/(tFunc(tStar)*SS))
  
  
}

produceContourPlot <- function(params, rGrid, tGrid, FigName, lowerProbaNoDefault){
  #
  #This function produces contour plots for the optimization problem determining 
  #the optimal staking policy. 
  #
  #INPUTS: - params: vector of parameters
  #        - rGrid: a vector containing annual staking reward rate values
  #        - tGrid: a vector containing values of the time when staking rewards
  #          are paid
  #        - FigName: the name assigned to the .jpg file saved (extension is
  #          excluded) 
  #        - lowerProbaNoDefault: the lowest acceptable probability of no default
  #
  #OUTPUT: - contour plots of the objective function of the optimization problem
  #        from the paper
  #        - the solution to the optimization problem in a vector (r,t,tr)
  #
  
  #Check the inputs. The calculation should be vectorized for the grids for r and t
  if (is.vector(tGrid) == 0){ 
    stop("tGrid must be a vector!")
  }else if(is.vector(rGrid) == 0){ 
    stop("rGrid must be a vector!")
  }else if(length(tGrid) != length(rGrid)){ 
    stop("Check the dimensions of rGrid and tGrid: calculations must be
          vectorized for rGrid and tGrid!")
  } 
  
  #Initialize a data frame to store the results
  datab <- data.frame(matrix(NA, nrow = length(tGrid), ncol = 3))
  names(datab)[1] <- "r"
  names(datab)[2] <- "t"
  names(datab)[3] <- "tr"
  
  #Compute the no-default probability
  probNoDefault <- evalPobaNoDefault(rGrid,tGrid,params)
  #Compute the value of the objective function
  o <- objFunc(rGrid,tGrid)
  o[probNoDefault < lowerProbaNoDefault] <- -0.03 # The value of the objective function for
  #the values of t and r that do not satisfy the constraints is set to -0.03
  #Assign the values to the corresponding columns of the data frame
  datab['r'] <- rGrid
  datab['t'] <- tGrid
  datab['tr'] <- o
  
  #Create a ggplot object
  q <- ggplot(datab, aes(x = r, y = t, z = tr, fill = tr)) +
    geom_raster(interpolate = T) +  
    scale_fill_gradientn(colours = rainbow(3))
  #Save the figure
  ggsave(glue('{FigName}.jpg', width = 10, height = 10, dpi = 150, units = "in", plot = q))
  #Return the figure to the screen
  print(q)
  
  #Find the maximum of the objective function st constraints and get the arg maximum
  inxs = which(datab[,3] == max(datab[,3]))
  dd = datab[inxs,]
  
  return(dd)
  
}


####Parameter set up for all the scenarios####

tMax <- 2 # Maximum duration of the contract (Upper bound on t)
tStarr <- tMax # Arg maximum value of f(t) for t in the constraints' set
rhoGrid <- c(-1, -0.5, 0, 0.5, 1) # Grid of values of the correlation coefficient
#between the transactions' intensity and the exchange rate innovations

paramrs <-c()
paramrs[1] <- 500*365*24*60*60 # lam0
paramrs[2] <-  1/0.08 # XX0
paramrs[3] <- 0.01 # FF
paramrs[4] <- 31.9*(10^9) # SS
paramrs[5] <- 250000000 # RR_0
paramrs[6] <- 0.3 # sigmX
paramrs[7] <- 0.1 # sigmLambda
paramrs[8] <- tStarr # tStar 

lb <- c(lowerR(paramrs), 0.00001) # Lower bounds on the variables (r,t)
ub <- c(upperR(paramrs), tMax) # Upper bounds on the variables (r,t)

#Vectorize the calculations over the grids for r and t
rGridd <- as.vector(repmat(linspace(lb[1]+0.0001,ub[1],100),1,400)) # Grid for the 
#staking rate. A small value is added to the lower bound so that the constraints 
#in the optimization problem are satisfied (See the paper for details)
tGridd <- as.vector(repmat(linspace(0.001,tMax,400),100,1)) # Grid for the time 
#when the staking rewards are paid


####Scenario I####
paramrs[9] <- 0 # muLambda
paramrs[10] <- 0 # muX

#Initialize a data frame to store the arg maximum
maximumsI <- data.frame()
#Do the calculations while looping over the rho values
for (i in 1:length(rhoGrid)){
  
  paramrs[11] <- rhoGrid[i]
  mval <- produceContourPlot(paramrs, rGridd, tGridd, glue('ScenarioI{rhoGrid[i]}'),0.99)
  
  maximumsI <- rbind(maximumsI , mval)
  
}

###############SCENARIO II############

paramrs[9] <- 0.1 # muLambda
paramrs[10] <- -0.3 # muX

#Initialize a data frame to store the arg maximum
maximumsII <- data.frame()
#Do the calculations while looping over the rho values
for (i in 1:length(rhoGrid)){
  
  paramrs[11] <- rhoGrid[i]
  mval <- produceContourPlot(paramrs, rGridd, tGridd, glue('ScenarioII{rhoGrid[i]}'),0.99)
  
  maximumsII <- rbind(maximumsII , mval)
  
}

###############SCENARIO III############

paramrs[9] <- 0.1 # muLambda
paramrs[10] <- -0.05 # muX

#Initialize a data frame to store the arg maximum
maximumsIII <- data.frame()
#Do the calculations while looping over the rho values
for (i in 1:length(rhoGrid)){
  
  paramrs[11] <- rhoGrid[i]
  mval <- produceContourPlot(paramrs, rGridd, tGridd, glue('ScenarioIII{rhoGrid[i]}'),0.99)
  
  maximumsIII <- rbind(maximumsIII , mval)
  
}

###############SCENARIO IV############

paramrs[9] <- 0 # muLambda
paramrs[10] <- -0.3 # muX

#Initialize a data frame to store the arg maximum
maximumsIV <- data.frame()
#Do the calculations while looping over the rho values
for (i in 1:length(rhoGrid)){
  
  paramrs[11] <- rhoGrid[i]
  mval <- produceContourPlot(paramrs, rGridd, tGridd, glue('ScenarioIV{rhoGrid[i]}'),0.99)
  
  maximumsIV <- rbind(maximumsIV , mval)
  
}

####Save the data on the optimal staking policy####

write_xlsx(maximumsI, "maximumsIa.xlsx")
write_xlsx(maximumsII, "maximumsIIa.xlsx")
write_xlsx(maximumsIII, "maximumsIIIa.xlsx")
write_xlsx(maximumsIV, "maximumsIVa.xlsx")
