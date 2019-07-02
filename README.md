# Classifying-Physical-Activity-Kaggle-
This repository contains R script to the in-class Kaggle activity on recognising physical activity based on phone signal data
This activity was part of the Behavioural Data Science course at the University of Amsterdam. 

The description of the competition on Kaggle is as follows:
This competition involves building a classifier that recognizes different types of physical activity from signals measured by the accelerometer and gyroscope in your smartphone, which both measure aspects of movement and orientation. The data for this competition were collected in a lab using a basic smartphone in experiments with human participants carrying out various daily activities in set order.
Activities classes were three static postures (standing, sitting, lying),  three dynamic activities (walking, walking downstairs and walking upstairs) and postural transitions that occurred between the static postures (stand-to-sit, sit-to-stand, sit-to-lie, lie-to-sit, stand-to-lie, and lie-to-stand)

For feature extraction the signal data was split into time-frames and from then each time period was treated as a histogram of signal readings. The features extracted included the mean, skewness, autocorrelation, entropy and spectral features. 

Due to many variables the variables were checked for multicolinearity using the correlation matrix and VIF index. 
The restrictions set by the course confined the possible methods for classifying the signals into activities to logistic regression and quadratic discriminant analysis (QDA). 

In the end the best model was using the QDA, selected through cross validated forward selection. Despite using a quite simple model, the activities were classified with around a 71% accuracy when tested on a test dataset. 
