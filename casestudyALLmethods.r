## Read in Data and clean that data
data_backup=data=read.csv("casestudydata.csv")
names(data)

out_sample=which(is.na(data$CKD)==1)
data_out=data[out_sample,]   ## the ones without a disease status
data_in=data[-out_sample,]   ## the ones with a disease status

data_in=na.omit(data_in)
dim(data_in)

CKD=data_in$CKD  # save the dependent variable
data_in=data_in[,-c(1,4,8,34)]

data_backup=data_in
data=data_in

###############
### PCA  ####
#############


## First let's try with the RAW unstandardized data, just in-sample data 

rho=cor(data)  ## we can use a correlation matrix as sigma - it's just like standardizing data
rho

sigma=var(data)
sigma

vars=diag(sigma)
percentvars=vars/sum(vars)


eigenvalues=eigen(sigma)$values
eigenvectors=eigen(sigma)$vectors

eigenvalues
eigenvectors

### Just take the first 5 for plotting reasons later
y1=as.matrix(data)%*%(eigenvectors[,1])
y2=as.matrix(data)%*%(eigenvectors[,2])
y3=as.matrix(data)%*%(eigenvectors[,3])
y4=as.matrix(data)%*%(eigenvectors[,4])
y5=as.matrix(data)%*%(eigenvectors[,5])


y=as.matrix(data)%*%eigenvectors
cor(y)

percentvars  # variance of real data

percentvars_pc=eigenvalues/sum(eigenvalues)
percentvars_pc  ##  Principal Components

names(data)
eigenvectors[,1]  ## name PCA 1  (LDL, Total Chol. - Negative Relationship)
eigenvectors[,2]  ## name PCA 2  (Weight, Waist size - Positive Relationship)

## plot percent of variances
ts.plot(cbind(percentvars,percentvars_pc),col=c("blue","red"),xlab="ith vector",ylab="percent variance")

## plot PCA values
plot(y1,y2, xlab="Low Cholesterol",ylab="Size of Person")
text(y1,y2, labels = rownames(data), pos = 4)

which(rownames(data)=="4042")
data_backup[2802,]   ## the "out of control" patient

## predict CKD with a regression
## all original data, 16% R-Squared
summary(lm(CKD~as.matrix(data_in)))
## all PCs
summary(lm(CKD~as.matrix(y)))

## let's pick the best few
summary(lm(CKD~y1))   ## "LDL and Total Chol"
## not predictive, why?
summary(lm(CKD~y2))   "Waist and Weight"

summary(lm(CKD~y1+y2+y3))   ## first 3, much better (4th doesn't do any better)
eigenvectors[,3]  ## name PCA 2  (Age, Sys. Blood Pressure - Positive Relationship)

#########################
## Next, let's do it with STANDARDIZED data (just like using correlation matrix)
##################

## We do this because we want an EQUAL importance between each variable. 
data=scale(data)
data=as.data.frame(data)



rho=cor(data)  
rho
sigma=var(data)
sigma  ## now this is the same as correlation.  So it doesn't matter which one we use.

vars=diag(sigma)
percentvars=vars/sum(vars)
percentvars  # all same percent  LOL....that's cuz we standardized it

eigenvalues=eigen(sigma)$values
eigenvectors=eigen(sigma)$vectors

eigenvalues  ## see which one is >1
eigenvectors

### Just take the first 5 for plotting reasons later
y1=as.matrix(data)%*%(eigenvectors[,1])
y2=as.matrix(data)%*%(eigenvectors[,2])
y3=as.matrix(data)%*%(eigenvectors[,3])
y4=as.matrix(data)%*%(eigenvectors[,4])
y5=as.matrix(data)%*%(eigenvectors[,5])


y=as.matrix(data)%*%eigenvectors
cor(y)   ## all mostly close to 0 except the last one - it's okay (last one not important)

percentvars  # variance of real data

percentvars_pc=eigenvalues/sum(eigenvalues)
percentvars_pc  ## will plot this later

names(data)
eigenvectors[,1]  ## name PCA 1  (Weight, BMI, Obese, Waist - Positive Relationship) "LARGE PEOPLE"
eigenvectors[,2]  ## name PCA 2  (Age, Hypertension, Blood Pressure  - Positive Relationship) "OLDER PRESSURED People" 

## plot percent of variances, look how it changes...
ts.plot(cbind(percentvars,percentvars_pc),col=c("blue","red"),xlab="ith vector",ylab="percent variance")

## plot PCA values
plot(y1,y2,xlab="LARGE PEOPLE",ylab="OLDER PRESSURED")
text(y1,y2, labels = rownames(data), pos = 4)
which(rownames(data)=="1846")
data_backup[1339,]   ## an "out of control" patient, 30s male, 193 kg weight, very large
which(rownames(data)=="82")
data_backup[54,]   ## 83 year old female, with hypertension 

## predict CKD with a regression

## all original data, 16% R-Squared, nothing changes with standardization except coefficients
summary(lm(CKD~as.matrix(data_in)))
## all PCs
summary(lm(CKD~as.matrix(y)))

## let's pick the best few
summary(lm(CKD~y1))   ## "LARGER PEOPLE"  2.3% R-squared, not great but a little better

summary(lm(CKD~y2))   ## "OLDER PRESSURED"   9.3% R-squared - really good for 1 variable

summary(lm(CKD~y1+y2))   ## first 2, much better than non-standardized, 11.6% R-Squared


##########################
## FACTOR ANALYSIS ##########
#############################


## It starts off the same as PCA. 

# we want to know which basic factors are driving our 30 variables or so.


## Step 1 , Explore the data
rho=cor(data)

## Step 2 - Compute the eigenvalues and eigenvactors of the correlation matrix

eigenvalues=eigen(rho)$values
eigenvalues
(eigenvalues)>1
eigenvectors=eigen(rho)$vectors


## Our Hypothesis , Do we only need 2 factors to represent our 30 variables???
m=2  ## can change
v=30  ## the original number of variables

## Step 3 - Compute Estimated Factor Loadings
L=matrix(nrow=v,ncol=m)
for (j in 1:m){
  L[,j]=sqrt(eigenvalues[j])*eigenvectors[,j]  
}

L  # first column is Factor 1, 2nd column is factor 2, these are ith loadings on the kth variable


## Step 4 - Compute common variance and unique variance. Common is good. Unique is Bad.

common=rowSums(L^2)  ## we want this to be close to 1
unique=1-common  ## this diagonal of error matrix

common  ## desire this to be close to 1
unique   ## desire this to be close to 0, goes down as you use more factors


## Step 5 - Check the model and reproduce correlation matrix using only the 2 loadings

phi=diag(v)*unique

recreate=L%*%t(L)+phi   ## this is our "fake" synthetic correlation matrix
recreate

rho   ## this is our real correlation matrix  - should be similar

## Step 6 - Create Residual Matrix, basically the real one minus the fake one

residual=rho-recreate
residual  ## check to see if off-diagonal elements are "small", actually this is quite small
hist(residual,breaks=50)  ## looks like it does pretty good, tighter with more factors

sum(residual[lower.tri(residual)]^2)  ## sum of square of off-diagonal elements, goes down with more factors

n=length(residual[lower.tri(residual)])  ## how many things you are squaring

sum(sqrt(residual[lower.tri(residual)]^2))/(n-m)  ## RMSE of off-diagonal elements


sum(eigenvalues[(m+1):v]^2)  ## sum of square of non-used eigenvalues, if big, it is bad



## Step 7  - Plot pairs of loadings to interpret factor loadings
## if we can't tell, we may need to do a varimax rotation

plot(L[,1],L[,2],col=1:v,xlab="Loading 1",ylab="Loading 2")
text(L[,1],L[,2],names(data))

## whem m>2, how do you cluster them together?  This is a clustering problem!  
## we will need to learn clustering too.

## variables that are clustered together are from a "similar factor"


## Step 8

install.packages('psych')
library(psych)

## should reproduce our results
fit2 <- principal(data, nfactors=m, rotate="none")

fit2

## more of an issue when you have 4 or 5 variables, separates it more visually
fit <- principal(data, nfactors=m, rotate="varimax")
fit
plot(fit$loadings[,1],fit$loadings[,2],col=1:v)
text(fit$loadings[,1],fit$loadings[,2],names(data))




