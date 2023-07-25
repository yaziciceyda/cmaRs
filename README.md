# cmaRs


Welcome to cmaRs, an R package for Conic Multivariate Adaptive Regression (CMARS).

### Description

An implementation of 'Conic Multivariate Adaptive Regression Splines (CMARS)' in R.
Please check Weber, et al. (2011) CMARS: a new contribution to nonparametric regression with multivariate adaptive regression splines supported by continuous optimization. Inverse Problems in Science and Engineering, 2012, 20.3: 371-400. doi:10.1080/17415977.2011.62477, for more details.
It constructs models by using the terms obtained from the forward step of MARS and then estimates parameters by using  'Tikhonov' regularization and conic quadratic optimization. It is possible to 
construct models for prediction and binary classification. It provides performance 
measures for the model developed. 

## Installation

In order to construct CMARS models, both MOSEK software and Rmosek package needed. Please follow carefully the steps available in https://docs.mosek.com/latest/rmosek/install-interface.html for successful installation.
``` r
library(cmaRs)
```

## Inputs

* **formula**: description of the model
* **degree**: maximum degree of interaction 
* **nk**: maximum number of model terms before pruning 
* **data**: data frame 
* **classification**: binary variable indicating whether the model is for prediction or binary classification 
* **threshold.class**: the threshold value that is used to convert probabilities to classes in binary classification models.

## Prediction Example
cmaRs can construct prediction models which is exemplified below.

``` {r echo=TRUE,  results='hide', eval=FALSE}
prediction.model <- cmaRs(Volume ~ ., degree = 2, nk = 20, data = trees)
summary(prediction.model)
```

Some performance measures are also printed at the end of the output. For instance, the $R^2$ , r and RSS values are given for the prediction models. It is also possible to construct several graphs of a prediction model with the **plot** function. 

## Classification Example
In addition to prediction modeling, it is also possible to construct binary classification models.

``` {r echo=TRUE,  results='hide', eval=FALSE}
library(earth) 
classification.model <- cmaRs(survived ~ age, nk = 35, classification = TRUE, data = etitanic)
summary(classification.model)
```

Note that, the classification argument is set as TRUE here which indicates a binary classification model. Moreover, the degree argument is used as 1 which is its default value indicating a main effect model. This model is constructed by using only continuous variable in the data set, age. Similar to the previous example, the **summary** function presents the details of the model.






