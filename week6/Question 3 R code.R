#Importing the data
data <- as.matrix(read.table("synthetic_regression.txt", nrow=300, stringsAsFactors = F))[,1:31]
data_i <- cbind(data[,1],1,data[,2:31])
colnames(data_i)[2] <- "X.0"

#Preparing the t and features matrices
output <- data_i[,1]
features <- data_i[,2:32]

###Question 3, Part 1 
##Performing the ols 
ols <- lm( output ~ features)
mle.res <- ols$residuals
#Storing the values of w and standard errors
w <- as.vector(summary(ols)$coefficients[,1])
mle.std <- as.vector(summary(ols)$coefficients[,2])
#Storing the values of deviance residuals for part 2 
mle.dev.res <- (mle.res^2)/var(output)

#Performing the robust regression 
q <- as.numeric(1/var(output))
nu <- 10
iterations <- 100
e <- as.vector(output - mean(output))
N <- nrow(features)
for (i in 1:iterations){
  eta <- ((nu + 1)/(nu + q*(e**2) - 2))*diag(N)
  
  W <-solve(t(features)%*%eta%*%features)%*%t(features)%*%eta%*%output
  q <- solve((1/N)*(t(e)%*%eta%*%e))
  e <- output - features%*%W
}

e %*% t(e)

a <- 5
class(a)
class(q)
