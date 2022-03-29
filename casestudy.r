

install.packages("dummies")
library(dummies)
data_backup=data=read.csv("casestudydata.csv")
names(data)

#class(data)
summary(data)

out_sample=which(is.na(data$CKD)==1)
data_out=data[out_sample,]   ## without a disease status
data_in=data[-out_sample,]   ## with a disease status

dim(data_in)

## Step 2  - Missing Data
data_in=na.omit(data_in)
dim(data_in)

CKD=data_in$CKD  # save the dependent variable
data_in=data_in[,-c(1,8, 34)] # omit ID, CareSource, and CKD

## create dummies for racial var
data_in$Activity=as.factor(data_in$Activity)
data_in <- dummy.data.frame(data_in, sep = ".")
names(data_in)

data_backup=data_in
data=data_in

## Step 3 and 4 - Correlation
summary(data)
rho=cor(data)  ## we can use a correlation matrix as sigma - it's just like standardizing data
rho

sigma=var(data)
sigma

vars=diag(sigma)
percentvars=vars/sum(vars)
percentvars

eigenvalues=eigen(sigma)$values
eigenvectors=eigen(sigma)$vectors

eigenvalues
eigenvectors

##########################
## PCA ##########
#############################

### For plotting
y1=as.matrix(data)%*%(eigenvectors[,1])
y2=as.matrix(data)%*%(eigenvectors[,2])
y3=as.matrix(data)%*%(eigenvectors[,3])
y4=as.matrix(data)%*%(eigenvectors[,4])
y5=as.matrix(data)%*%(eigenvectors[,5])
y6=as.matrix(data)%*%(eigenvectors[,6])
y7=as.matrix(data)%*%(eigenvectors[,7])

y=as.matrix(data)%*%eigenvectors
cor(y)

percentvars  # variance of real data

percentvars_pc=eigenvalues/sum(eigenvalues)
percentvars_pc  ##  Principal Components

names(data)
eigenvectors[,1]  ## name PCA 1  (LDL, Total Chol. - Negative Relationship)
eigenvectors[,2]  ## name PCA 2  (Weight, Waist size - Positive Relationship)
eigenvectors[,3]
eigenvectors[,4]
eigenvectors[,5]
eigenvectors[,6]
eigenvectors[,7]

## plot percent of variances
ts.plot(cbind(percentvars,percentvars_pc),col=c("blue","red"),xlab="ith vector",ylab="percent variance")
legend(25, .6, legend=c("Real Data", "Principal Components"), col=c("blue", "red"), lty=c(1,1), cex=.8)

## plot PCA values
plot(y1,y2, xlab="Bad Cholesterol Level",ylab="Size by Weight")
text(y1,y2, labels = rownames(data), pos = 4)

which(rownames(data)=="4042")
data_backup[2802,]   ## the "out of control" patient

plot(y1,y6, xlab="Bad Cholesterol Level",ylab="Systolic Blood Pressure Relative to Age")
text(y1,y3, labels = rownames(data), pos = 4)


## predict CKD with a regression
## all original data, 16% R-Squared
summary(lm(CKD~as.matrix(data_in)))
## all PCs
summary(lm(CKD~as.matrix(y)))

## pick variables with large variances
summary(lm(CKD~data_in$LDL+data_in$Total.Chol+data_in$SBP))
## let's pick the best few
summary(lm(CKD~y1))   ## "LDL and Total Chol"
## not predictive, why?
summary(lm(CKD~y2))   ##"Waist and Weight"

summary(lm(CKD~y1+y2+y3))   ## first 3, much better (4th doesn't do any better)

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

### for plotting later
y1=as.matrix(data)%*%(eigenvectors[,1])
y2=as.matrix(data)%*%(eigenvectors[,2])
y3=as.matrix(data)%*%(eigenvectors[,3])
y4=as.matrix(data)%*%(eigenvectors[,4])


y=as.matrix(data)%*%eigenvectors
cor(y)   ## all mostly close to 0 except the last one - it's okay (last one not important)

percentvars  # variance of real data

percentvars_pc=eigenvalues/sum(eigenvalues)
percentvars_pc  ## will plot this later

names(data)
eigenvectors[,1]  ## name PCA 1  (Weight, BMI, Obese, Waist - Positive Relationship) "LARGE PEOPLE"
eigenvectors[,2]  ## name PCA 2  (Age, Hypertension, Blood Pressure  - Positive Relationship) "OLDER PRESSURED People" 

## plot percent of variances, look how it changes...
ts.plot(cbind(percentvars,percentvars_pc),col=c("blue","red"), xlab="ith vector",ylab="percent variance", xlim=c(1,37))
grid(lty=1, col=gray(.9)) 
legend(25, .1, legend=c("Standardized Data", "Principal Components"), col=c("blue", "red"), lty=c(1,1), cex=.8)

## plot PCA values
plot(y1,y2, xlab="Overweighted People",ylab="Blood Pressure and Age", xlim=c(-6, 11))
text(y1,y2, labels = rownames(data), pos = 4)

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


## Hypothesis: Do we only need 2 factors to represent our 31 variables?
m=2  ## can change
v=37  ## the original number of variables

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

sum(residual[lower.tri(residual)]^2)  ## sum of square of off-diagonal elements, goes down with more factors (sum of squared Σ−(𝐿𝐿^′+𝜓))

n=length(residual[lower.tri(residual)])  ## how many things you are squaring

sum(sqrt(residual[lower.tri(residual)]^2))/(n-m)  ## RMSE of off-diagonal elements


sum(eigenvalues[(m+1):v]^2)  ## sum of square of non-used eigenvalues, if big, it is bad



## Step 7  - Plot pairs of loadings to interpret factor loadings
## if we can't tell, we may need to do a varimax rotation

plot(L[,1],L[,2],col=1:v,xlab="Loading 1",ylab="Loading 2")
text(L[,1],L[,2],names(data), pos=3, cex=.8)

## whem m>2, how do you cluster them together?  This is a clustering problem!  
## we will need to learn clustering too.

## variables that are clustered together are from a "similar factor"


## Step 8

install.packages('psych')
library(psych)

## should reproduce our results
fit2 <- principal(data, nfactors=5, rotate="none")

fit2

## more of an issue when you have 4 or 5 variables, separates it more visually (rotated facto loadings)
fit <- principal(data, nfactors=m, rotate="varimax")
fit
plot(fit$loadings[,1],fit$loadings[,2],col=1:v,xlab="Loading 1",ylab="Loading 2")
text(fit$loadings[,1],fit$loadings[,2],names(data), pos=3, cex=.8)


##########################
## CLUSTERING ANALYSIS ##########
#############################

## default: euclidean distance
D=dist(data)
D

## single linkage
hc<-hclust(D,"single")
summary(hc)

## Plot dendrogram
plot(hc, xlab="Clusters",ylab="Levels")

i = 4
memb<-cutree(hc,k=i)
memb

# install.packages("dendextend")
library(dendextend)

# Plot dendrogram based on the number of clusters
dend <- as.dendrogram(hc)
dend <- color_branches(dend, k = i)
plot(dend, main = "Colored branches")

## get clustser centers
cent<-NULL
for (k in 1:i){cent<-rbind(cent,colMeans(data[memb==k,,drop=FALSE]))}
cent

## calculate sum of total SS . within SS for each cluster
one=sum((data[memb==1,,drop=FALSE]-cent[1,])^2)
two=sum((data[memb==2,,drop=FALSE]-cent[2,])^2)
three=sum((data[memb==3,,drop=FALSE]-cent[3,])^2)
four=sum((data[memb==4,,drop=FALSE]-cent[4,])^2)
one
two
three
four
tss_single=one+two+three+four  ## total sum of squares from cluster centroid
tss_single

## kmeans clustering
library(flexclust)


# Step 1 - Run kmeans and analyze clusters
cl=kmeans(data,i,nstart = 10)  ## let's keep same number of clusters as before
cl

cl$betweenss
cl$totalss

# Step 2 - Plot clusters
plot(data, col = cl$cluster)

#plot(data$Weight, data$SBP, xlab="Weight", ylab="SBP", col = cl$cluster)
#text(data$Weight, data$SBP, labels = rownames(data), cex=.5, col=cl$cluster)
#points(cl$centers[1,], col=cl$cluster, pch = 20)
#text(data,rownames(data),col=cl$cluster)

# Step 3 - Choose k  (plot total sum of squares)
tss<-rep(0,8)
for (k in 1:8){tss[k]=kmeans(data,k)$tot.withinss}
plot(1:8,tss, xlab="# of clusters", ylab="Total Within Sum of Squares")


# Step 4 - Interpret clusters
cl$centers

cl$cluster


# Step 5  - Plot clusters with pca (this makes your intreptation even better)

## here we are putting food names to the clusters...it helps people understand.
install.packages("pca3d")
install.packages("pca2d")

install.packages("rgl")
library(rgl)
library(pca3d)
library(pca2d)

?pca3d

pca <- prcomp(data, scale.= TRUE )
pca2d( pca, group=cl$cluster )
pca2d( pca, group=cl$cluster,show.labels=TRUE)
pca3d( pca, group=cl$cluster,show.labels=TRUE )


## Instead of pca above

install.packages(c("factoextra", "FactoMineR"))
library("factoextra")
library("FactoMineR")
res.pca <- prcomp(data, scale = TRUE)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 12))

var <- get_pca_var(res.pca)
var

library("corrplot")
corrplot(var$cos2, is.corr=FALSE)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PCs
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
)

# Contributions of variables to PC1&PC2
fviz_contrib(res.pca, choice = "var", axes = 1, top = 5, fill="blue", col="blue")
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 5, fill="orange", col="orange")

# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1:5, fill="blue", col="blue")

## individuals contributing the most to the PCs
fviz_pca_ind(res.pca, select.ind = list(contrib = 20))
data_in[c(240,1846,3488,5118,5153)]

## Another way of K-Means

library(tidyverse)  # data manipulation
install.packages("cluster")
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

K2 <- kmeans(data, centers = 2, nstart = 25)
K3 <- kmeans(data, centers = 3, nstart = 25)
K4 <- kmeans(data, centers = 4, nstart = 25)
K5 <- kmeans(data, centers = 5, nstart = 25)

p1 <- fviz_cluster(K2, geom = "point", data = data) + ggtitle(" K = 2")
p2 <- fviz_cluster(K3, geom = "point", data = data) + ggtitle(" K = 3")
p3 <- fviz_cluster(K4, geom = "point", data = data) + ggtitle(" K = 4")
p4 <- fviz_cluster(K5, geom = "point", data = data) + ggtitle(" K = 5")

install.packages("dplyr")
library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

K2
K3

data %>% 
  mutate(Cluster = K3$cluster) %>%
  group_by(Cluster) %>%
  summarize_at(c("Obese","Weight","BMI","Waist","Age","SBP","Hypertension"), median, na.rm = TRUE)



i = 2
k2=kmeans(data, centers = i, nstart = 25)
str(k2)
k2
fviz_cluster(k2, data = data)
distance=get_dist(data)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


set.seed(12345)

# function to compute total within-cluster sum of square 
wss=function(k) {
  kmeans(data, k, nstart = 10 )$tot.withinss
}

wss

## optimal number of clusters
fviz_nbclust(data, kmeans, method = "wss")
fviz_nbclust(data, kmeans, method = "silhouette")

gap_stat=clusGap(data, FUN = kmeans, nstart = 25,
                 K.max = 8, B = 50)
fviz_gap_stat(gap_stat)


##########################
## LOGISTIC REGRESSION ##########
#############################

data=read.csv("casestudydata.csv")
names(data)
data=data[,-1]  ## remove ID
data=data[,-c(7)] # omit CareSource
out_sample=which(is.na(data$CKD)==1)
data_out=data[out_sample,]   ## the ones without a disease status
data_in=data[-out_sample,]   ## the ones with a disease status
summary(data_in)

## Step 1 - Explore and relabel Data
class(data)
summary(data)

#data_in=na.omit(data_in)
## Training and testing sets
smp_siz = floor(0.75*nrow(data_in))
smp_siz

set.seed(123)
train_ind=sample(seq_len(nrow(data_in)),size = smp_siz) 
train=data_in[train_ind,]
test=data_in[-train_ind,]

## Step 2  - Run the Logistic Regression with one variable
## Logistic regression on selected vars
model_tr=glm(CKD~Age+Hypertension+Diabetes+Obese+CVD+SBP+LDL+Weight+BMI+Waist+Female, family="binomial",train)

formula(model_tr)
summary(model_tr)

confint.default(model_tr)
confint(model_tr)

with(model_tr, null.deviance - deviance)
with(model_tr, df.null - df.residual)
with(model_tr, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))

-2*logLik(model_tr)
with(model_tr, pchisq(deviance, df.residual, lower.tail = FALSE))

## Predict probabilities and Odds Ratios of New Data
newdata1=data_out

summary(newdata1)
#newdata1$Activity=as.factor(newdata1$Activity)
#newdata1 <- dummy.data.frame(newdata1, sep = ".")
phatnew=predict(model_tr, newdata = newdata1, type = "response")
phatnew

## odds ratios
phatnew/(1-phatnew)

phat3=predict(model_tr,type="response")  # predicts for ALL in sample data
summary(phat3)

classify=ifelse(phat3>.3,1,0)  # this is a threshold, we say if probability >x% , then say "yes"
summary(classify)  # notice that not many are "yes"  - is this desirable?

c_accuracy(train$CKD,classify)


# let's compare that to a "simple" model with only age
model=glm(CKD~Age,family="binomial",data=data_in)
Age=seq(0,100,by=.1)
phat=predict(model,list(Age = Age),type="response")
summary(phat)  # range from 0.02% to 76.2%, it's similar but doesn't fit as well
plot(data_in$Age, data_in$CKD, pch = 16, xlab = "Age", ylab = "CKD")
lines(Age, phat)

## plot the actual probabilities for the "complicated" model
plot(train$Age, train$CKD, pch = 16, xlab = "Age", ylab = "CKD")
points(train$Age, phat3,col="blue")  # this plots all phat3's according to age

## Step 8 - Classification
summary(phat3)
classify=ifelse(phat3>.5,1,0)  # this is a threshold, we say if probability >50% , then say "yes"
summary(classify)  # notice that not many are "yes"  - is this desirable?

c_accuracy(data_in$CKD,classify)  # to run this you must run my code below first.

# notice these are the accuracy results when you use 50% as the threshold. they are kind of extreme
#    the false positive rate is almost 0% BECAUSE you almost NEVER say "yes"
#         true positive rate is 13.9% which isn't that bad because you have almost no false positives

## Step 9 - Caclculate Costs
acc=c_accuracy(train$CKD,classify)
c1=100   # penalize me  $100 for a false positive
c2=200  #  penalize me $200 for a false negatives
cost=acc[9]*c1+acc[10]*c2
acc
cost   ## my costs are $48,800  because I got hit with a ton of false negative costs 
# i said no too much!  many FN.  hmmm, 

## YOu must change the threshold of 50% above , to lower your costs.
##    How to do this?  One way is to search over 101 thresholds (each % from 0 to 100)

##  you may realize at some point, that plotting an ROC curve with roc()  gives you all possibilities
##   that's a high level understanding

######################### model ##########################
## Logistic regression on selected vars
## The training  model should not use imputed data
data_in=na.omit(data_in)
model_semi=glm(CKD~Age+Female+Hypertension+Diabetes+Obese+CVD+SBP+LDL+Weight+BMI+Waist, family="binomial", data_in)

formula(model_semi)
summary(model_semi)

model_final=glm(CKD~Age+Female+Hypertension+Diabetes+CVD+Weight+BMI, family="binomial", data_in)

formula(model_final)
summary(model_final)

## kNN Imputation
install.packages("purrr")
install.packages("RANN")
library(caret)
library(RANN)

data$Activity=as.factor(data$Activity)

##values will be normalized
preProcValues <- preProcess(data %>% 
                              dplyr::select(Educ,Unmarried,Income,Insured,Weight,Height,BMI,Obese,Waist,SBP,DBP,HDL,
                                            LDL,Total.Chol,Activity,PoorVision,Diabetes,Hypertension,Stroke,CVD,Fam.CVD,CHF,Anemia),
                            method = c("knnImpute"),
                            k = 25,
                            knnSummary = mean)
impute_data <- predict(preProcValues, data, na.action = na.pass)

procNames <- data.frame(col = names(preProcValues$mean), mean = preProcValues$mean, sd = preProcValues$std)

## Denormalize
for(i in procNames$col){
  impute_data[i] <- impute_data[i]*preProcValues$std[i]+preProcValues$mean[i] 
}

summary(impute_data)

out_sample=which(is.na(data$CKD)==1)
data_out=impute_data[out_sample,]   ## the ones without a disease status
data_in=impute_data[-out_sample,]   ## the ones with a disease status
summary(data_in)
summary(data_out)


model_final=glm(CKD~Age+Female+Hypertension+Diabetes+CVD+Weight+BMI, family="binomial", data_in)

formula(model_final)
summary(model_final)

confint.default(model_final)
confint(model_final)

with(model_final, null.deviance - deviance)
with(model_final, df.null - df.residual)
with(model_final, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))

-2*logLik(model_final)
with(model_final, pchisq(deviance, df.residual, lower.tail = FALSE))

## Predict probabilities and Odds Ratios of New Data
newdata1=data_out
summary(newdata1)
#newdata1$Activity=as.factor(newdata1$Activity)
#newdata1 <- dummy.data.frame(newdata1, sep = ".")
phatnew=predict(model_final, newdata = newdata1, type = "response")
phatnew
write.csv(phatnew, file = "prob.csv")

## odds ratios
phatnew/(1-phatnew)

## Test on the data_in (with CKD) sample to see how the model predicts (prob. predicted)
phat=predict(model_final,type="response")
summary(phat)

yes_ckd=which(data$CKD==1)
no_ckd=which(data$CKD==0)
## randomly inspect ckd and no ckd individuals
data_ckd=data_in[yes_ckd,]
data_nockd=data_in[no_ckd,]

data_ckd=data_ckd[sample(nrow(data_ckd), 5), ]
data_nockd=data_nockd[sample(nrow(data_nockd), 5), ]
write.csv(rbind(data_ckd, data_nockd), file = "sample.csv")

## odds ratios
#avg odd ratio of ckd sample
(phat/(1-phat))[yes_ckd]
#avg odd ratio of no ckd sample
(phat/(1-phat))[no_ckd]

phat[c(167,842,1565,2499,3869,988,2095,3368,3703,5229)]

classify1=ifelse(phat>.14,1,0)  # this is a threshold, we say if probability >14% , then say "yes"
summary(classify1)

c_accuracy(data_in$CKD,classify1)

classify2=ifelse(phat>.20,1,0)  # this is a threshold, we say if probability >14% , then say "yes"
summary(classify2)

c_accuracy(data_in$CKD,classify2)

# of the out sample data (2819 patients)
classify=ifelse(phatnew>.14,1,0)
summary(classify)
write.csv(classify, file = "class.csv")

c_accuracy(data_out$CKD,classify)


install.packages("pROC")
install.packages("ROCR")
library(pROC)
library(ROCR)
rocobj <- roc(data_in$CKD, phat)
auc <- round(auc(data_in$CKD, phat),4)
ggroc(rocobj, colour = 'steelblue', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))


## Step 8 - Classification
summary(phat3)
classify=ifelse(phat3>.5,1,0)  # this is a threshold, we say if probability >50% , then say "yes"
summary(classify)  # notice that not many are "yes"  - is this desirable?

c_accuracy(data_in$CKD,classify)  # to run this you must run my code below first.

# notice these are the accuracy results when you use 50% as the threshold. they are kind of extreme
#    the false positive rate is almost 0% BECAUSE you almost NEVER say "yes"
#         true positive rate is 13.9% which isn't that bad because you have almost no false positives

## Step 9 - Caclculate Costs
acc=c_accuracy(data_out$CKD,classify)
c1=100   # penalize me  $100 for a false positive
c2=200  #  penalize me $200 for a false negatives
cost=acc[9]*c1+acc[10]*c2
acc
cost

## CKD by age
h <- hist(data_in$Age[(data_in$CKD==1)], breaks=c(20,30,40,50,60,70,80,90), plot=FALSE)
h$density = h$counts/sum(h$counts) * 100
plot(h, main="Distribution of Age Groups Among CKD Population",
     xlab="Age Group",
     ylab="Percent",
     col="orange",
     freq=FALSE)
abline(h = c(10,20,30), col = "grey", lty = "dotted")
plot(h, main="Distribution of Age Groups",
     xlab="Age Group",
     ylab="Percent",
     col="orange",
     freq=FALSE,
     add=TRUE)

##CKD by hypertension
df=ifelse(data_in$Hypertension[data_in$CKD==1]==1,"Yes","No")
df2=ifelse(data_in$Hypertension[data_in$CKD==0]==1,"Yes","No")
data_perc <- t(prop.table(table(df))) * 100 
data_perc2 <- t(prop.table(table(df2))) * 100 

a <- rbind(data_perc,data_perc2)
colnames(a) <- c("NO", "YES")
rownames(a) <- c("CKD", "NO CKD")

barplot(height=a, beside = TRUE, legend.text = TRUE, main="Distribution of Hypertension", col=c("blue", "orange"), ylab = "Percent", , ylim=c(0,100), xlab="Hypertension")
abline(h = c(20,40,60,80), col = "grey", lty = "dotted")
barplot(height=a, beside = TRUE, legend.text = TRUE, main="Distribution of Hypertension", col=c("blue", "orange"), ylab = "Percent", , ylim=c(0,100), xlab="Hypertension", add=TRUE)


##CKD by diabetes
df=ifelse(data_in$Diabetes[data_in$CKD==1]==1,"Yes","No")
df2=ifelse(data_in$Diabetes[data_in$CKD==0]==1,"Yes","No")
data_perc <- t(prop.table(table(df))) * 100 
data_perc2 <- t(prop.table(table(df2))) * 100 

a <- rbind(data_perc,data_perc2)
colnames(a) <- c("NO", "YES")
rownames(a) <- c("CKD", "NO CKD")

barplot(height=a, beside = TRUE, legend.text = TRUE, main="Distribution of Diabetes", col=c("blue", "orange"), ylab = "Percent", , ylim=c(0,100), xlab="Diabetes")
abline(h = c(20,40,60,80), col = "grey", lty = "dotted")
barplot(height=a, beside = TRUE, legend.text = TRUE, main="Distribution of Diabetes", col=c("blue", "orange"), ylab = "Percent", , ylim=c(0,100), xlab="Diabetes", add=TRUE)



## Function Below, RUN THIS FIRST
## make sure actuals and classifications are 0 (no) or 1 (yes) only 
##  Built by Matthew J. Schneider
c_accuracy=function(actuals,classifications){
  df=data.frame(actuals,classifications);
  
  
  tp=nrow(df[df$classifications==1 & df$actuals==1,]);        
  fp=nrow(df[df$classifications==1 & df$actuals==0,]);
  fn=nrow(df[df$classifications==0 & df$actuals==1,]);
  tn=nrow(df[df$classifications==0 & df$actuals==0,]); 
  
  
  recall=tp/(tp+fn)
  precision=tp/(tp+fp)
  accuracy=(tp+tn)/(tp+fn+fp+tn)
  tpr=recall
  fpr=fp/(fp+tn)
  fmeasure=2*precision*recall/(precision+recall)
  scores=c(recall,precision,accuracy,tpr,fpr,fmeasure,tp,tn,fp,fn)
  names(scores)=c("recall","precision","accuracy","tpr","fpr","fmeasure","tp","tn","fp","fn")
  
  #print(scores)
  return(scores);
}


## produce column 3 for project
out_sample=which(is.na(data$CKD)==1)
data2=data[out_sample,]
summary(data2)

data2$Points<-ifelse(data2$Female==1,2,1)

tf <- which(data2$Age<45)
data2$Points[tf] <- data2$Points[tf] + 1

tf <- which((data2$Age>44)&(data2$Age<65))
data2$Points[tf] <- data2$Points[tf] + 2

tf <- which(data2$Age>64)
data2$Points[tf] <- data2$Points[tf] + 3

tf <- which((data2$SBP>129)&(data2$SBP<140))
data2$Points[tf] <- data2$Points[tf] + 1

tf <- which((data2$SBP>139)&(data2$SBP<150))
data2$Points[tf] <- data2$Points[tf] + 2

tf <- which(data2$SBP>150)
data2$Points[tf] <- data2$Points[tf] + 3

tf <- which(data2$Diabetes==1)
data2$Points[tf] <- data2$Points[tf] + 2

tf <- which(data2$BMI<24.9)
data2$Points[tf] <- data2$Points[tf] + 1

tf <- which((data2$BMI>24.9)&(data2$BMI<30))
data2$Points[tf] <- data2$Points[tf] + 2

tf <- which(data2$BMI>29.9)
data2$Points[tf] <- data2$Points[tf] + 3

tf <- which(data2$CVD==1)
data2$Points[tf] <- data2$Points[tf] + 2

## Binary classification
data2$CKD_tr<-ifelse(data2$Points>4,1,0)

write.csv(data2$Points, file = "here.csv")






