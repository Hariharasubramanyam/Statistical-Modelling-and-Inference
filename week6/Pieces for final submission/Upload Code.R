#Importing the data and forming the matrices

synthetic_regression <- read.csv("synthetic_regression.txt",sep="", nrows=300)

dataset = synthetic_regression[,c(1:31)]

X = as.matrix(cbind(1,dataset[,c(2:31)]))
y = as.vector(dataset[,1])

N <- 300

#Doing the OLS and storing the required values
ols<-lm(y~X)
beta.se <- data.frame(summary(ols)$coefficients[,c(1,2)])
colnames(beta.se) = c("value","se")
#beta.se: The column value stores w and se stores the standard errors
#In the following, res stores the residuals and mledm stores the deviance residuals
res <- as.vector(ols$residuals)
mledm <- (res^2)/var(y)

#Robust Regression and EM
#Setting initial values and then performing 100 iterations
iter <- 100
q <- 1/var(y)
nu <- 10
e <- y - mean(y)
for(i in 1:iter){
  eta <- as.vector((nu + 1)/(nu + q*(e^2) - 2))
  
  W <- solve(t(X) %*% diag(eta) %*% X , t(X) %*% diag(eta) %*% y)
  e <- as.vector(y - X %*% W)
  q <- as.numeric(solve((1/N)*t(e) %*% diag(eta) %*% e))
}

#Using all the stuff from the appendix to calculate the standard errors
ematrix <- matrix(0,31,31)
for(j in 1:N){
  c <- as.numeric((nu + 1)*(nu - 2 - q*(e[j]^2))/((nu + q*(e[j]^2) -2)^2))
  summ <- c*X[j,] %*% t(X[j,])
  ematrix = ematrix + summ
}
var <- solve(q*ematrix)
#The square root of diagonal matrix of variance holds the standard error values
robust.se <- data.frame(value = W,se = sqrt(diag(var)))
row.names(robust.se)[1] <- "Intercept"
#Storing the deviance residuals from the robust regression
robustdm = (e^2)*(q*eta)

##Plotting the answers to question 3, parts 1 and 2 
library(ggplot2)
##I have very limited skill in ggplot2 and hence I am just copy pasting their code which needs to be altered for our submission
#r stores the w +- 1.96 the standard errors for robust regression
#p stores the w +- 1.96 the standard errors for mle regression
#We need to plot r and p beside each other for submission
r <- ggplot() + geom_point(data = robust.se, aes(x = factor(row.names(robust.se),levels = row.names(robust.se)), y = value),colour = 'blue', size = 3) + geom_errorbar(data = robust.se, aes(x = factor(row.names(robust.se),levels = row.names(robust.se)), y = value, ymax = value + 1.96*se, ymin=value - 1.96*se),colour = 'red', width = 0.4) + labs(x= "Predictors",y="Coefficient estimate with error bars",title = "Robust Regression")

p <- ggplot() + geom_point(data = beta.se, aes(x = factor(row.names(beta.se),levels = row.names(beta.se)), y = value),colour = 'blue', size = 3) + geom_errorbar(data = beta.se, aes(x = factor(row.names(beta.se),levels = row.names(beta.se)), y = value, ymax = value + 1.96*se, ymin=value - 1.96*se),colour = 'red', width = 0.4) + labs(x= "Predictors",y="Coefficient estimate with error bars",title = "MLE Regression")

#Code for plotting the deviance residuals
#This part is complete but can be made aesthetically better

layout(matrix(c(1,2), 1, 2))
plot(mledm,ylab = "MLE deviance") + abline(a= quantile(mledm,0.99),b= 0)
plot(robustdm,ylab = "Robust deviance") + abline(a= quantile(robustdm,0.99),b= 0)

##############################################################################
#The solution code for question 3
#Same initial steps as earlier except for loglik and Qtrace
#Loglik compares the values of loglikelihood and is used for breaking the loop once we reach 
#a point where the value barely changes
#Qtrace stores the values of q
iter <- 100
q <- 1/var(y)
nu <- 10
loglik <- rep(0,iter)
Qtrace <- rep(q,iter)
e <- y - mean(y)

for(i in 1:iter){
  eta <- as.vector((nu + 1)/(nu + q*(e^2) - 2))
  
  W <- solve(t(X) %*% diag(eta) %*% X , t(X) %*% diag(eta) %*% y)
  e <- as.vector(y - X %*% W)
  q <- as.numeric(solve((1/N)*t(e) %*% diag(eta) %*% e))

#Coding for the log likelihood based on Q(\theta,\theta^{'})  
  loglik[i] <- (-q/2)*(t(e) %*% diag(eta) %*% e)  + (N/2)*log(q)
  Qtrace[i] <- q
#This part breaks the for loop once the change in value is below 10^{-4}
  if(i > 1){
    if(loglik[i] - loglik[i-1] < 1e-4){
      cat("number of iterations",i)
      break
    }
  }
}

#The graphical output of this solution is the plot of the loglikelihood
#This plot can be made better looking and more informative as well.
plot(loglik[1:i])

## The code for answering part 4 of Question 3
# The code coputes the various values of log likelihood for different values of nu
# It then moves on to the next value of nu once the maximum loglik for the current value is found
#The code breaks when an increase in nu doesn't significantly increase the loglik
#This part is 100% verified yet because the loop break threshold is extremely low.
#I'lltry fixing this but everything else is pretty complete
N <- 300
nu <- seq(2.2,42,0.2)
loglik <- rep(200)

for(j in 1:200){
  W <- beta.se$value
  q <- 1/var(y)
  e <- y - mean(y)
#We set 10 because of what we conclude in question 3  
  for(i in 1:10){
    eta <- as.vector((nu[j] + 1)/(nu[j] + q*(e^2) - 2))
    
    W <- solve(t(X) %*% diag(eta) %*% X , t(X) %*% diag(eta) %*% y)
    e <- as.vector(y - X %*% W)
    q <- as.numeric(solve((1/N)*t(e) %*% diag(eta) %*% e))
    
    loglik[j] <- (-1/2)*(t(e) %*% diag(q*eta) %*% e) + (N/2)*log(q)
  }
  if(j > 1){
    if(loglik[j] - loglik[j-1] < 1e-1){
      cat("Optimum value of nu",nu[j])
      break
    }
  }
}

loglik <- loglik[1:j]
nu = nu[1:j]
layout(matrix(1, 1, 1))
plot(y=loglik,x=nu,ylab = "Log-Likelihood",type="b")

