library(tidyverse)
library(hmisc)
library(corrplot)
library(car)
library(caret)
library(klaR)
library(MASS)
library(psd)


rm(list = ls())
##gettubg labels
##########################################################################################
act_labels = read_delim("activity_labels.txt", " ", col_names=F, trim_ws = T) 
act_labels = act_labels %>% dplyr::select(X1,X2)
act_labels


labels <- read_delim("RawData/Train/labels_train.txt", " ", col_names = F)
colnames(labels) <- c('trial', 'userid', 'activity', 'start', 'end')
labels

train_labels = labels
##########################################################################################


##function to calculate entropy 
##########################################################################################
entropy  <- function(x, nbreaks = nclass.Sturges(x)) {
  r = range(x)
  x_binned = findInterval(x, seq(r[1], r[2], len= nbreaks))
  h = tabulate(x_binned, nbins = nbreaks) # fast histogram
  p = h/sum(h)
  -sum(p[p>0] * log(p[p>0]))
}
##########################################################################################


##DEFNING THE FUNCTION to extract features
##########################################################################################
extractTimeDomainFeatures <- 
  function(filename, labels = train_labels) {
    # extract user and experimental run ID's from file name
    username = gsub(".+user(\\d+).+", "\\1", filename)
    expname =  gsub(".+exp(\\d+).+", "\\1", filename)
    
    # import the data from the file
    user01 <- read_delim(filename, " ", col_names = F, progress = TRUE, col_types = "ddd")
    
    # select this users activity labels from the `labels` frame for this experimental run
    user_labels <- 
      labels %>% 
      dplyr::filter(userid==as.numeric(username) & trial == as.numeric(expname)) %>% # select rows pertaining to current signals
      mutate(segment = row_number()) %>%                                             # segment identifies different WALKING_UPSTAIRS etc
      gather(start_end, vec, -trial, -userid, -activity, -segment) %>%               # stack start and end on top of each other
      arrange(vec) %>%                                                               # arrange rows in natural order
      mutate(activity = ifelse(start_end == 'end', NA, activity), activity_id = row_number()) # remove activity label `end` time rows
    
    # add activity labels to each sample
    user <- 
      user01 %>% 
      mutate(sample = row_number()-1) %>%
      mutate(activity_id = findInterval(sample, user_labels$vec)) %>%
      left_join(user_labels, by = "activity_id") 
    
    # split in epochs of 128 samples and compute features per epoch
    usertimedom <- 
      user %>%
      mutate(epoch = sample %/% 128) %>% # epoch = 2.56 sec
      group_by(epoch) %>%
      summarise(
        user_id = username, # added to identify user in data frame rows
        exp_id = expname,   # added to identify experimental run in data frame rows
        activity = names(which.max(table(c("-", activity)))),
        sample = sample[1],
        m1 = mean(X1),  ###Means
        m2 = mean(X2),
        m3 = mean(X3),
        sd1 = sd(X1), ##Standard deviations
        sd2 = sd(X2),
        sd3 = sd(X3),
        q1_25 = quantile(X1, .25), ##Quartiles
        skew1 = e1071::skewness(X1), ##Skewness
        skew2 = e1071::skewness(X2),
        skew3 = e1071::skewness(X3),
        AR1.2 = cor(X1, lag(X1, n = 3), use = "pairwise"), ##Autocorrelations
        AR2.2 = cor(X2, lag(X2, n = 3), use = "pairwise"),
        AR3.2 = cor(X3, lag(X3, n = 3), use = "pairwise"),
        AR12.1 = cor(X1, lag(X2, n= 2),use = "pairwise"),
        AR23.1 = cor(X2, lag(X3, n=2),use = "pairwise"),
        AR13.1 = cor(X1, lag(X3,n=2), use = "pairwise"),
        max.1 = max(abs(X1)), ##Max and min values
        max.2 = max(abs(X2)),
        max.3 = max(abs(X3)),
        min.1 = min(abs(X1)),
        min.2 = min(abs(X2)),
        min.3 = min(abs(X3)),
        bim.1 = bimodality_coefficient(X1), ##Bimodality coefficients
        bim.2 = bimodality_coefficient(X2),
        bim.3 = bimodality_coefficient(X3),
        entropy.1 = entropy(X1),  #Entropy
        entropy.2 = entropy(X2), 
        entropy.3 = entropy(X3), 
        reg1 = ar(X1, aic=FALSE, order.max=1)$ar, ##Auto regression
        reg2 = ar(X2, aic=FALSE, order.max=1)$ar,
        reg3 = ar(X3, aic=FALSE, order.max=1)$ar,
        fft.1 = mean(abs(fft(X1))), ##Means of the fast Fourier transform (FFT), Fourier analysis converts a signal 
        fft.2 = mean(abs(fft(X2))), #from its original domain (often time or space) to a representation in the
        fft.3 = mean(abs(fft(X3))),  #frequency domain.
        psd.1 = mean(spectrum(X1, plot = F)$spec), #Mean of the power spectral density (A Power Spectral Density 
        psd.2 = mean(spectrum(X2, plot = F)$spec), #(PSD) is the measure of signal's power content versus frequency        
        psd.3 = mean(spectrum(X3, plot = F)$spec),
        n=n()
      ) 
    
    usertimedom 
  }
##########################################################################################


###running the function on accelartion and gyro and combining them
##########################################################################################
filenames_acc <- dir("RawData/Train/", "^acc", full.names = TRUE) # for demo only first 5 files

myData_Acc =  filenames_acc %>%
  map_dfr(extractTimeDomainFeatures) # map_dfr runs `extractTimeDomainFeatures` on all elements in filenames and binds results row wise


filenames_gyro <- dir("RawData/Train/", "^gyr", full.names = TRUE) 

myData_Gyro =  filenames_gyro %>%
  map_dfr(extractTimeDomainFeatures) # map_dfr runs `extractTimeDomainFeatures` on all elements in filenames and binds results row wise


myData_Full <- left_join(myData_Acc, myData_Gyro, by = c("epoch", "user_id", "exp_id")) %>% 
  dplyr::select(-activity.y, -sample.y, -n.y) 

myData_Full$activity.x <-  as.integer(myData_Full$activity.x)

myData_Full$activity.x <- plyr::mapvalues(myData_Full$activity.x, 1:12, act_labels$X2)
data_to_work_with <- myData_Full  %>% dplyr::filter(n.x >30) %>% dplyr::select( -c(1,3, sample.x, n.x)) %>% drop_na() %>% mutate(activity.x = factor(activity.x))

##########################################################################################


##make a function that splits the data into train and test, but based on the user_id and
#returns the cross validated classification accuracy using the model
##########################################################################################
cv_per_user <- function(dataset, model){

ids <- group_indices(dataset, user_id)

data_to_work_with_cv <- cbind(ids, dataset)

data_to_work_with_cv <- as_tibble(data_to_work_with_cv)

means <- c()
for (i in 1:10) {
train_index <- sample(x = ids, size = 12, replace = F)
train_set <- data_to_work_with_cv %>% dplyr::filter(ids %in% train_index)
'%ni%' <- Negate('%in%')
test_set <- data_to_work_with_cv %>% dplyr::filter(ids %ni% train_index)

lda.pred=predict (model, test_set[,-c(1,2)])
lda.class=lda.pred$class

means[i] <- mean(lda.class == test_set$activity.x)
      }
mean(means)

  }
##########################################################################################

##fit a 10/15/20 variable stepwise qda, 
##first starting variable
qda.fit.15.1 <- train(activity.x ~ ., data = dplyr::select(data_to_work_with, -user_id),
              method = "stepQDA",
              trControl = trainControl(method = "cv"),
              tuneGrid = expand.grid(maxvar = 20, direction = "forward"))

qda.fit.15.1 <-qda(activity.x ~ max.1.x + m2.x + AR2.2.x + AR12.1.x + q1_25.y + q1_25.x + m3.y + fft.3.x +
                     fft.2.x + AR12.1.y + reg1.x + reg1.y + entropy.2.x + reg3.x + min.2.x, 
                   data_to_work_with)

##model 1
cv_per_user(data_to_work_with,qda.fit.15.1)



##Our predictions weren't very good so we tried a few things to improve it
#The first was removing values high in VIf (this was done in two ways)
#The second was removing highly correlated variables
#The third was removing extreme outliers
#The forth was centering and scaling the data

##########################################################################################
##REMOVING VALUES HIGH IN VIF
##There are two ways you can do it:
##1 - Simply remove variables with highest VIF
##2 - Remove variables with higherst VIF but also lowest effect on the outcome variables, in hopes that 
##the total collinearity decreases but you still keep the best predicting variables.

##Way 1
vif_model <- lm(as.numeric(activity.x) ~., data_to_work_with[,-1])

vif(vif_model)

vif_model <- lm(as.numeric(activity.x) ~., dplyr::select(data_to_work_with,-user_id, -m1.x,
                                                        -sd1.x, -max.1.x, -sd1.y, -fft.2.y,
                                                        -fft.1.y, -q1_25.x, -max.3.y, 
                                                        -fft.2.x, -fft.1.x, -fft.3.x, 
                                                        -q1_25.y, -min.2.x, -fft.3.y, 
                                                        -psd.3.x, -reg3.x, -psd.2.x, -sd3.y, -max.2.y, 
                                                        -sd3.x, -psd.2.y))

vif(vif_model)

data_no_cor_1 <- dplyr::select(data_to_work_with, -m1.x,
                               -sd1.x, -max.1.x, -sd1.y, -fft.2.y,
                               -fft.1.y, -q1_25.x, -max.3.y, 
                               -fft.2.x, -fft.1.x, -fft.3.x, 
                               -q1_25.y, -min.2.x, -fft.3.y, 
                               -psd.3.x, -reg3.x, -psd.2.x, -sd3.y, -max.2.y, 
                               -sd3.x, -psd.2.y)


qda.fit.15_no_cor_1 <- train(activity.x ~ ., data = data_no_cor_1[,-1],
                    method = "stepQDA",
                    trControl = trainControl(method = "cv"),
                    tuneGrid = expand.grid(maxvar = 15, direction = "forward"))


qda.fit.15_no_cor_1 <- qda(activity.x ~ min.1.x + m2.x + AR12.1.x + AR2.2.x + psd.1.x + m3.y + entropy.1.x +
                         AR1.2.y + AR1.2.x +AR2.2.y + AR12.1.y + entropy.3.y + m3.x + skew2.x + AR13.1.y, 
                       data = data_no_cor_1)


cv_per_user(data_no_cor_1, qda.fit.15_no_cor_1)


#Way 2 ##LEFT TO FIX THIS VIF MODEL, MAYBE MAKE A FUNCTION FOR IT BASED ON THE DISTANCE
##Now try removing correlation by taking the effect of the variables into account
##first we can examine effects of individual predictors 
#make a loop which fits single variables one by one
single_effects <- c()
for (i in 3:length(colnames(data_to_work_with))) {
  
  data = data_to_work_with[,c(2,i)]
  
  qda.fit <- qda(activity.x ~., data = data)
  
  acc <- cv_per_user(data_to_work_with, qda.fit)
                       
single_effects[i] <- paste(colnames(data_to_work_with[,i]), round(acc,3))
}
single_effects


# NA                  NA                  "m1.x 0.39"        
# "m2.x 0.436"        "m3.x 0.382"        "sd1.x 0.411"      
# "sd2.x 0.335"       "sd3.x 0.328"       "q1_25.x 0.55"     
# "skew1.x 0.28"      "skew2.x 0.271"     "skew3.x 0.203"    
# "AR1.2.x 0.362"     "AR2.2.x 0.293"     "AR3.2.x 0.286"    
# "AR12.1.x 0.36"     "AR23.1.x 0.247"    "AR13.1.x 0.274"   
# "max.1.x 0.593"     "max.2.x 0.357"     "max.3.x 0.328"    
# "min.1.x 0.51"      "min.2.x 0.421"     "min.3.x 0.36"     
# "bim.1.x 0.315"     "bim.2.x 0.219"     "bim.3.x 0.214"    
# "entropy.1.x 0.224" "entropy.2.x 0.218" "entropy.3.x 0.237"
# "reg1.x 0.371"      "reg2.x 0.338"      "reg3.x 0.358"     
# "fft.1.x 0.562"     "fft.2.x 0.411"     "fft.3.x 0.374"    
# "psd.1.x 0.406"     "psd.2.x 0.339"     "psd.3.x 0.311"    
# "m1.y 0.293"        "m2.y 0.294"        "m3.y 0.241"       
# "sd1.y 0.338"       "sd2.y 0.31"        "sd3.y 0.324"      
# "q1_25.y 0.363"     "skew1.y 0.23"      "skew2.y 0.219"    
# "skew3.y 0.191"     "AR1.2.y 0.233"     "AR2.2.y 0.238"    
# "AR3.2.y 0.245"     "AR12.1.y 0.245"    "AR23.1.y 0.278"   
# "AR13.1.y 0.209"    "max.1.y 0.35"      "max.2.y 0.314"    
# "max.3.y 0.318"     "min.1.y 0.241"     "min.2.y 0.249"    
# "min.3.y 0.249"     "bim.1.y 0.219"     "bim.2.y 0.197"    
# "bim.3.y 0.203"     "entropy.1.y 0.217" "entropy.2.y 0.205"
# "entropy.3.y 0.215" "reg1.y 0.264"      "reg2.y 0.249"     
# "reg3.y 0.313"      "fft.1.y 0.361"     "fft.2.y 0.333"    
# "fft.3.y 0.346"     "psd.1.y 0.278"     "psd.2.y 0.255"    
# "psd.3.y 0.303"  

##examing vif once again
vif_model <- lm(as.numeric(activity.x) ~., data_to_work_with[,-1])

vif(vif_model)

vif_model <- lm(as.numeric(activity.x) ~., dplyr::select(data_to_work_with,-user_id, -m1.x, -fft.2.y, 
                                                         -sd1.y, -fft.3.y,-max.1.y, -fft.3.x, -psd.2.y, 
                                                         -max.3.y, -m1.y, -psd.2.x, -sd1.x, -max.1.x, 
                                                         -fft.2.x, -fft.1.y, -AR3.2.x, -min.3.x, 
                                                          -fft.1.x, -min.1.x, -sd3.y, -sd2.y, -max.2.x))
vif(vif_model)

data_no_cor_2 <- dplyr::select(data_to_work_with, -m1.x, -fft.2.y, 
                               -sd1.y, -fft.3.y,-max.1.y, -fft.3.x, -psd.2.y, 
                               -max.3.y, -m1.y, -psd.2.x, -sd1.x, -max.1.x, 
                               -fft.2.x, -fft.1.y, -AR3.2.x, -min.3.x, 
                               -fft.1.x, -min.1.x, -sd3.y, -sd2.y, -max.2.x)

qda.fit.15_no_cor_2 <- train(activity.x ~ ., data = data_no_cor_2[,-1],
                         method = "stepQDA",
                         trControl = trainControl(method = "cv"),
                         tuneGrid = expand.grid(maxvar = 15, direction = "forward"))


qda.fit.15_no_cor_2 <- qda(activity.x ~ q1_25.x + m2.x + AR12.1.x + q1_25.y + reg2.x + sd2.x + m3.y + bim.3.y +
                         entropy.1.x + m3.x + AR13.1.y + skew1.x + AR13.1.x + AR23.1.y + reg1.y, 
                       data = data_no_cor_2)


cv_per_user(data_no_cor_2, qda.fit.15_no_cor_2)
##########################################################################################

##REMOVING HIGHLY CORRELATED VALUES
##########################################################################################
###Can also just remove varables with correlations over .80
CUTOFF <- 0.80
cor_matrix <- cor(data_to_work_with[,-c(1,2)]) 
cor_high <- findCorrelation(cor_matrix, CUTOFF)
high_cor_remove <- row.names(cor_matrix)[cor_high] 


rawTrain <- data_to_work_with[,-c(1,2)][,  -cor_high] 


descrCor <- cor(rawTrain)
summary(descrCor[upper.tri(descrCor)])

data_no_correlation <- as_tibble(cbind(data_to_work_with[,c(1,2)], rawTrain))


##modelling 
qda.fit.15_no_cor <- train(activity.x ~ ., data = data_no_correlation[,-1],
                           method = "stepQDA",
                           trControl = trainControl(method = "cv"),
                           tuneGrid = expand.grid(maxvar = 15, direction = "forward"))

qda.fit.15_no_cor_3 <- qda(activity.x ~ q1_25.x + AR12.1.x + psd.1.x + AR2.2.x + reg1.y + m3.y + reg2.y + m1.y +
                           m3.x + skew1.x + AR1.2.x + entropy.2.y + AR13.1.y + max.2.x + entropy.3.x, 
                         data = data_no_correlation)



cv_per_user(data_no_correlation, qda.fit.15_no_cor_3)
##########################################################################################



##REMOVING OUTLIERS
##########################################################################################
##Now trying to remove the outliers and fit a qda, only removing the most extreme of outliers
##this function we found online
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 18 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

data_no_outliers <- data_to_work_with  %>%  mutate_at(vars(m1.x:psd.3.y), remove_outliers) %>% drop_na()


qda.fit.10_no_outliers <- train(activity.x ~ ., data = data_no_outliers[,-1],
                         method = "stepQDA",
                         trControl = trainControl(method = "cv"),
                         tuneGrid = expand.grid(maxvar = 10, direction = "forward"))

qda.fit.10_no_outliers <- qda(activity.x ~ max.1.x + m2.x + q1_25.y + AR12.1.x + AR2.2.x + sd1.x + AR12.1.y +
                                fft.3.x + m2.y + entropy.1.x, data = data_no_outliers)


cv_per_user(data_no_outliers, qda.fit.10_no_outliers)
##########################################################################################


##SCALING DATA 
##########################################################################################
preProcValues <- preProcess(data_to_work_with[,-c(1,2)], method = c("center", "scale"))
trainTransformed <- predict(preProcValues, data_to_work_with)


qda.fit.15_scaled <- train(activity.x ~ ., data = trainTransformed[,-1],
                                method = "stepQDA",
                                trControl = trainControl(method = "cv"),
                                tuneGrid = expand.grid(maxvar = 30, direction = "forward"))


##THIS WAS THE BEST MODEL IN THE END
qda.fit.15_scaled <- qda(activity.x ~ max.1.x + m2.x + AR2.2.x +AR12.1.x + m1.y + q1_25.x + m3.y + reg3.x +
                           reg1.y + reg2.y + psd.3.y + reg1.x + entropy.3.y + skew2.x +AR13.1.y, data = trainTransformed)


cv_per_user(trainTransformed, qda.fit.15_scaled)

##The qda model using centeredand scaled data gave us the best improvement in our predicted variable

##########################################################################################





##Making the test datast
##########################################################################################
filenames_acc_test <- dir("RawData/Test/", "^acc", full.names = TRUE) 
filenames_gyro_test <- dir("RawData/Test/", "^gyr", full.names = TRUE) 


myData_Acc_test =  filenames_acc_test %>%
  map_dfr(extractTimeDomainFeatures)

myData_Gryo_test <- filenames_gyro_test %>%
  map_dfr(extractTimeDomainFeatures)

myData_Full_test <- left_join(myData_Acc_test, myData_Gryo_test, by = c("epoch", "user_id", "exp_id")) %>% 
  dplyr::select(-sample.y, -n.y, -n.x, -activity.x, -activity.y) 


##########################################################################################
##Predicting on scaled and centered data 
preProcValues <- preProcess(myData_Full_test[,-c(1,2)], method = c("center", "scale"))
testTransformed <- predict(preProcValues, myData_Full_test)

lda.pred=predict(qda.fit.15_scaled, testTransformed)
preds <- lda.pred$class
myData_Full_test$activity <- preds


myData_Full_test %>%
  mutate(user_id = paste("user", user_id, sep=""), exp_id = paste("exp", exp_id, sep="")) %>%
  unite(Id, user_id, exp_id, sample.x) %>%
  dplyr::select(Id, Predicted = activity) %>%
  write_csv("test_set_predictions10.csv")

file.show("test_set_predictions10.csv")
##########################################################################################






