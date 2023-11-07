library(tidyverse)


#reading in data and extracting the considered variables for our models
raw <- read_csv("Datasets/MergedData/MergedData-ChangeInPop.csv", na=c("N/A","NA","(X)",".","N","-"),skip=2)
exclude_vars <- c("PopChg", "Year","GEO_ID","NAME")

# Use lapply to convert all variables (columns) except the excluded ones to numeric
raw[, !(names(raw) %in% exclude_vars)] <- lapply(raw[, !(names(raw) %in% exclude_vars)], as.numeric)

raw[, (names(raw) %in% exclude_vars)] <- lapply(raw[, (names(raw) %in% exclude_vars)], as.factor)


#character vector of model building vars
vars <- read_csv("Datasets/MergedData/Considered Variables.csv", skip=2) 
vars <- names(vars) 

###creating dataframe with only the selected vars
pop <- raw[,vars]
pop <- pop[pop$Year != 2011, ] #2011 has N/A for target variable

#variable plots
library(ggplot2)
library(gridExtra)



# Identify numerical variables in dataframe
numerical_vars <- names(pop)[sapply(pop, is.numeric)]

# Create histograms for the numerical variables
hist_list <- lapply(numerical_vars, function(var) {
  ggplot(pop, aes_string(x = var)) +
    geom_histogram(bins = 20, fill = "blue", color = "blue") +
    labs(title = var) +
    theme_minimal()
})


# Create sublists of 9 plots each
sublists <- split(hist_list, ceiling(seq_along(hist_list) / 9))

# Display the histograms in grids of 9 at a time
for (i in seq_along(sublists)) {
  grid.arrange(grobs = sublists[[i]], ncol = 3)
}




library(caTools)
set.seed(8480)
split <- sample.split(pop$PopChg,SplitRatio=0.5)

#check number of NAs per column
sapply(pop, function(x) sum(is.na(x)))

##imputing by mean##
pop.imp <- pop
vars.na <- pop.imp %>% select_if(colSums(is.na(.))>0) %>% names
mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

pop.imp <-pop.imp %>% 
  mutate(across(where(is.factor) & !PopChg, ~replace_na(.,mode(.[split])))) %>%  # Nominal Input: By Mode #
  mutate(across(where(is.numeric), ~replace_na(.,mean(.[split],na.rm = TRUE))))

# Create Missing Value Flag #
pop.imp[paste(vars.na, "NA", sep=".")] <- ifelse(is.na(pop[vars.na]), 1, 0) 

pop.imp <- pop.imp %>% mutate(across(ends_with(".NA"), as.factor))


##class imbalance
summary(pop$PopChg)/length(pop$PopChg)
summary(pop[split,]$PopChg)/length(pop[split,]$PopChg)



### Variable selection- Logistic Regression ###
pop.mdl <- pop.imp # pass on data

full = glm(PopChg ~., family=binomial, data=pop.mdl[split,])
null = glm(PopChg ~1,family=binomial, data=pop.mdl[split,])
n <- sum(split)

reg.bwd.BIC <- step(full, direction="backward",k=log(n), trace = FALSE)
summary(reg.bwd.BIC)

reg.bwd.AIC <- step(full, direction="backward", trace = FALSE)
summary(reg.bwd.AIC)

reg.step.BIC <- step(full, direction="both",k=log(n), trace=F)

reg.step.AIC <- step(full, direction="both", trace=F)

reg.bwd.BIC.prob<-predict(reg.bwd.BIC,pop.mdl[!split, ], type = "response")
reg.bwd.BIC.class <- ifelse(reg.bwd.BIC.prob > 0.5, 1, 0)
reg.bwd.BIC.misc<- 1-mean(reg.bwd.BIC.class == pop.mdl$PopChg[!split])

reg.bwd.AIC.prob<-predict(reg.bwd.AIC,pop.mdl[!split, ], type = "response")
reg.bwd.AIC.class <- ifelse(reg.bwd.AIC.prob > 0.5, 1, 0)
reg.bwd.AIC.misc<- 1-mean(reg.bwd.AIC.class == pop.mdl$PopChg[!split])

reg.step.BIC.prob<-predict(reg.step.BIC,pop.mdl[!split, ], type = "response")
reg.step.BIC.class <- ifelse(reg.step.BIC.prob > 0.5, 1, 0)
reg.step.BIC.misc<- 1-mean(reg.step.BIC.class == pop.mdl$PopChg[!split])

reg.step.AIC.prob<-predict(reg.step.AIC,pop.mdl[!split, ], type = "response")
reg.step.AIC.class <- ifelse(reg.step.AIC.prob > 0.5, 1, 0)
reg.step.AIC.misc<- 1-mean(reg.step.AIC.class == pop.mdl$PopChg[!split])

c(reg.bwd.BIC.misc,reg.bwd.AIC.misc,reg.step.BIC.misc,reg.step.AIC.misc) #Backwards and stepwise models are identical, AIC model is best by misclassification rate
#reg.bwd.AIC.misc = .4396

# Variable Importance based on Optimal Logistic Regression Model
library(caret)
varImp(reg.bwd.AIC) %>% arrange(desc(Overall))
varImp(reg.bwd.BIC) %>% arrange(desc(Overall))


## Lasso Regression
library(tidyverse)
library(olsrr)
library(data.table)
library(lmtest)
library(car)
library(glmnet)
#Taking certain variables out of the dataset and defining y vector and x matrix
xmat2b <- pop.mdl
yb <- xmat2b$PopChg
xmat2c <- subset(xmat2b, select = -c(PopChg))

# Convert factors to dummy variables
xmat2c <- model.matrix(~ . -1, data = xmat2c)
xmat3b <- as.matrix(xmat2c)

#Finding best lambda value and doing lasso regression
yc <- as.numeric(yb)
cv.lasso1 = cv.glmnet(xmat3b,yc,alpha=1, type.measure="class")
plot(cv.lasso1)
bestlamlasso=cv.lasso1$lambda.min
lasso2=glmnet(xmat3b,yc,alpha=1,lambda=bestlamlasso, family="binomial")
coef(lasso2) #Looking at the parameter estimates.
sort(coef(lasso2))
plot(lasso2, label=TRUE) 

# Extract coefficients as a numeric vector
coefficients <- as.numeric(coef(lasso2))

# Extract variable names
variable_names <- rownames(coef(lasso2))

# Create a data frame with coefficients and variable names
coefficients_df <- data.frame(Coefficient = coefficients, Variable = variable_names)

# Sort the data frame by the absolute values of the coefficients
sorted_coefficients_df <- coefficients_df[order(-abs(coefficients_df$Coefficient)), ]

# Access the sorted variable names and coefficients
sorted_variable_names <- sorted_coefficients_df$Variable
sorted_coefficients <- sorted_coefficients_df$Coefficient




##transformations
library(caret)
pop.xf <- pop.mdl

#removing 0's NAs
# pop.xf <- pop.xf %>% mutate(across(where(is.numeric),~replace(.,.==0,NA)))
# pop.xf <-pop.xf %>% 
#   mutate(across(where(is.factor) & !PopChg, ~replace_na(.,mode(.[split])))) %>%  # Nominal Input: By Mode #
#   mutate(across(where(is.numeric), ~replace_na(.,mean(.[split],na.rm = TRUE))))

#log+1 transform
#pop.xf <- pop.xf %>% 
  #mutate(across(!CL_UP & where(is.numeric), ~ log(.+1)))

#log+1 transform across non-mean/median metrics
pop.xf <- pop.xf %>% 
   mutate(across((!MEDIAN_INCOME | MeanTravel | MEAN_INCOME | MEDIAN_MALE_EARNINGS | MEDIAN_FEM_EARNINGS) & where(is.numeric), ~ log(.+1)))

### Variable selection- Logistic Regression TRANSFORMED VALUES ###
pop.mdl <- pop.xf # pass on data

full = glm(PopChg ~., family=binomial, data=pop.mdl[split,])
null = glm(PopChg ~1,family=binomial, data=pop.mdl[split,])
n <- sum(split)

reg.bwd.BIC <- step(full, direction="backward",k=log(n), trace = FALSE)
summary(reg.bwd.BIC)

reg.bwd.AIC <- step(full, direction="backward", trace = FALSE)
summary(reg.bwd.AIC)

reg.step.BIC <- step(full, direction="both",k=log(n), trace=F)

reg.step.AIC <- step(full, direction="both", trace=F)

reg.bwd.BIC.prob<-predict(reg.bwd.BIC,pop.mdl[!split, ], type = "response")
reg.bwd.BIC.class <- ifelse(reg.bwd.BIC.prob > 0.5, 1, 0)
reg.bwd.BIC.misc<- 1-mean(reg.bwd.BIC.class == pop.mdl$PopChg[!split])

reg.bwd.AIC.prob<-predict(reg.bwd.AIC,pop.mdl[!split, ], type = "response")
reg.bwd.AIC.class <- ifelse(reg.bwd.AIC.prob > 0.5, 1, 0)
reg.bwd.AIC.misc<- 1-mean(reg.bwd.AIC.class == pop.mdl$PopChg[!split])

reg.step.BIC.prob<-predict(reg.step.BIC,pop.mdl[!split, ], type = "response")
reg.step.BIC.class <- ifelse(reg.step.BIC.prob > 0.5, 1, 0)
reg.step.BIC.misc<- 1-mean(reg.step.BIC.class == pop.mdl$PopChg[!split])

reg.step.AIC.prob<-predict(reg.step.AIC,pop.mdl[!split, ], type = "response")
reg.step.AIC.class <- ifelse(reg.step.AIC.prob > 0.5, 1, 0)
reg.step.AIC.misc<- 1-mean(reg.step.AIC.class == pop.mdl$PopChg[!split])

c(reg.bwd.BIC.misc,reg.bwd.AIC.misc,reg.step.BIC.misc,reg.step.AIC.misc) #Stepwise/backwards models are identical again for AIC/BIC
#best model is misc rate of 0.4432726 for AIC.  Worse than pre-transformation

###Decision Tree###
library(rpart)
tree <- rpart(PopChg ~ ., data=pop.imp[split,],control=rpart.control(cp=0))

library(pROC)
cp.seq=tree$cptable[,1]
fscore<-numeric()
fscore[1]<-0  # Set root node F-score zero
for (i in 2:length(cp.seq)) {
  tree.prob = predict(prune(tree, cp=cp.seq[i]), pop.imp[!split,],type="prob")[,2] 
  rocCurve.tree <- roc(pop.imp$PopChg[!split], tree.prob, quiet=TRUE)
  treeThresh <-  coords(rocCurve.tree, x = "best", best.method = "closest.topleft")
  tree.class <- as.factor(ifelse(tree.prob >= treeThresh$threshold, 1,0))
  fscore[i]<-confusionMatrix(tree.class,pop.imp$PopChg[!split],positive = "1")$byClass["F1"]
}

plot(tree$cptable[,'nsplit']+1,fscore,
     type="o", xlab="Number of Leaves", ylab="F-score")

tree.final=prune(tree,cp=cp.seq[fscore==max(fscore)])  #max f-score is 0.5681197
library(partykit)
plot(as.party(tree.final))
tree.pred <- predict(tree.final, newdata = pop.imp[!split,], type = "class")
unname(1-confusionMatrix(tree.pred,pop.imp$PopChg[!split])$overall["Accuracy"]) #misclassification rate is .4414; worse than stepwise

#####Random Forest#####
pop.rf <- pop.imp
library(randomForest)

#internal downsampling
minor<-unname(summary(pop.rf$PopChg[split])[2])
set.seed(1234)
RF <- randomForest(PopChg ~., data=pop.rf[split,],
                   ntree = 1000,
                   strata= pop.rf$PopChg[split],
                   sampsize=c(minor,minor),
                   importance = TRUE)

# Make predictions #
library(caret)
RF.class<- predict(RF, newdata=pop.rf[!split,], type="response")
unname(1-confusionMatrix(RF.class, pop.rf$PopChg[!split],
                        positive = "1")$overall["Accuracy"])  
#misc rate is 0.437722


# Variable importance #
RF$importance
varImpPlot(RF)  
