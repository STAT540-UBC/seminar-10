STAT540 - Seminar 10: Supervised learning, classification, cross
validation, variable selection
================

## Attributions

Contributors: Gabriela Cohen Freue, W. Evan Durno, Jasleen Grewal, and
Keegan Korthauer

## Learning objectives

By the end of this tutorial, you should be able to  
- Filter gene expression data to remove uninformative features.  
- Further subset variables using a fixed information criteria, and know
when it may not be required/possible to do so.  
- Understand the concept of held-out test set, training set, and
validation set.  
- Select a suitable classifier (based on size of training dataset) and
train a model with the input training data.  
- Get predicted labels for a validation set, from the trained classifier
model.  
- Understand the utility of cross validation in selecting an optimal
model  
- Be able to identify at what part of the supervised learning process CV
should be used. - Explain how cross validation is different from
bootstrapping.  
- Assess the best selected model’s performance on the held out test
set.  
- Recognize why we can’t go back and reselect models after assessment on
the held out test set.  
- Be able to distinguish between the following scenarios for using
supervised machine learning, and enumerate metrics for assessing
performance in each case:  
- 2 class classification problem  
- Multi-class classification problem  
- Be able to define accuracy, precision, sensitivity, specificity,
F1-score, Kappa score.  
- Replicate analysis in either MLInterface, GLMnet, CMA, or Caret (take
home exercise, optional).

# Introduction

In this Seminar we go over packages and codes for performing supervised
learning and evaluation in R. In supervised learning, one is given a
training data set with a known response, or such as a class, and a set
of covariates or variables. The goal is to use this set to generate a
model that predicts values for a response given covariates in a new
dataset. A supervised learning process can be decomposed into the
following steps:

*Step 1*: Data preprocessing. Before getting started with feature
selection and training our model, we must make sure we have adjusted the
data to account for any batch effects (technical sources of variation
like measuring samples in groups) or biases from outliers (individuals
that look distinctly “different” from the rest of our sample).

*Step 2*: Select Features. Before training a model, in many
applications, it is usually important to perform a pre-filtering step in
which one retains only the most informative features (e.g. genes) as
candidate “biomarkers”. The amount of features retained to select and
train a model is up to the analyst and the methods used in the next
steps. For example,running some methods may be unfeasible or have an
outrageously slow runtime with a large number of features.

*Step 3*: Select and train a classifier. Once the set of candidate
markers have been selected, the next step is to select and train a model
to predict the labels of a test data. We will also tune parameters for
each classifier using cross-validation.

*Step 4*: Test. Finally, a model is chosen and used to predict labels of
a test data. From there we evaluate the model’s performance.

## R packages

There are many packages for performing supervised learning in R, each of
which may implement one or more algorithms. There have also been at
least two major efforts to unify these libraries under a common
framework to make them easier to use: `MLInterfaces` and `CMA`. Although
these may be useful and save you a lot of time in your analysis, it is
important that you understand what these packages are doing and what
they are *not* doing. Thus, I will not use these packages in this
Seminar but I encourage you to reproduce the analysis using at least one
of them! (I recommend `CMA`).

Install the following packages from Bioconductor: `CMA` and `GEOquery`,
and from CRAN: `ROCR`, `car`, `e1071` (for SVM), and `glmnet` along with
their dependencies.

``` r
library(BiocManager)
install('GEOquery')
install('CMA')
install('ROCR')
install(c('e1071','glmnet','mlbench','gbm', 'car','dimRed'))
install('caret')
install('kernlab')
install('gbm')
```

``` r
library(MASS)
library(tidyverse)
theme_set(theme_bw())
library(car)
library(limma)
library(e1071)
library(glmnet)
library(ROCR)
library(CMA)
library(class)
library(GEOquery)
```

## Data Set

This seminar is based on a dataset that comes from a paper by [Smeets et
al. 2010](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23177),
who studied Affymetrix expression profiles from primary breast tumors.
Smeets group was interested in whether tumors which had spread to lymph
nodes (LN positive, generally a bad sign) have different gene expression
profiles than LN negative tumors. If so, a gene expression signature can
be use to predict tumor class.

Their data set contains 24236 genes on 116 samples. The status of the
lymph node is known for each sample, with 59 LN positive and 57 LN
negative. Samples were divided into two parts: 96 samples (48 LN
positive and 48 LN negative) were used as a “training” set and 20
samples (11 LN positive and 9 LN negative) were used as a “test” set.
There is also a quantitative measure, “LnRatio”, the fraction of
affected lymph nodes, presumably reflecting “how bad” the LnStatus is.
Thus, we can use this dataset to illustrate classification and
regularization methods! This seminar will focus on the first task, i.e.,
classification. In the past, Paul selected this dataset to illustrate
the challenges of supervised learning tasks!

In the paper, the authors trained a support vector machine classifier to
distinguish between LN positive and LN negative tumors (i.e.,
classification), and evaluated the results using ROC curves. After some
optimization, they got an area under the ROC curve (AUC) of 0.66 on the
training set and 0.65 on the test data. This is better than chance, but
still not very convincing results about the relevance of the derived
molecular signature (random chance would give an AUC of 0.5; perfect
classification would give an AUC of 1.0).

### Data Preparation

First, let’s retrieve our datasets from GEO with `getGEO` from
`GEOquery` package. Warning: this may take several minutes! So to avoid
re-downloading in the future, save it as an RData object.

``` r
# Returns a list of expressionsets
datgeo <- getGEO('GSE23177', GSEMatrix = TRUE, AnnotGPL = TRUE) 
dat <- datgeo[[1]]   #Note that dat is an ExpressionSet

dat
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 24236 features, 116 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM570498 GSM570499 ... GSM570613 (116 total)
    ##   varLabels: title geo_accession ... patient type:ch1 (49 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1007_s_at 1053_at ... 91952_at (24236 total)
    ##   fvarLabels: ID Gene title ... GO:Component ID (21 total)
    ##   fvarMetadata: Column Description labelDescription
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 21116709 
    ## Annotation: GPL570

Now we’ll do some data wrangling to pull out the metadata variables we
are interested in, and recode some of them.

``` r
str(pData(dat), max.level = 0)
```

    ## 'data.frame':    116 obs. of  49 variables:

``` r
# extract only those variables of interest 
pData(dat) <- pData(dat) %>%
  rename(sample_id = geo_accession,
         LnStatus = characteristics_ch1.2,
         LnRatio = characteristics_ch1.3,
         Set = characteristics_ch1) %>%
  select(sample_id, LnStatus, LnRatio, Set) %>%
  mutate(LnStatus = factor(gsub("ln: ", "", LnStatus))) %>%
  mutate(LnRatio = as.numeric(gsub("lnratio: ", "", LnRatio))) %>%
  mutate(Set = ifelse(Set == "patient type: training set", "training", "test"))

str(pData(dat))
```

    ## 'data.frame':    116 obs. of  4 variables:
    ##  $ sample_id: chr  "GSM570498" "GSM570499" "GSM570500" "GSM570501" ...
    ##  $ LnStatus : Factor w/ 2 levels "neg","pos": 1 1 1 1 1 1 1 1 1 2 ...
    ##  $ LnRatio  : num  0 0 0 0 0 0 0 0 0 0.5 ...
    ##  $ Set      : chr  "test" "test" "test" "test" ...

``` r
#Note: LNRatio will not be used in this Seminar. However, you can use it to try some of the regularization techniques learned in class
```

Great! Next, let’s split the `ExpressionSet` object into two different
parts - one for the training and one for the test set.

``` r
# split the ExpressionSet into training and test sets. 
table(pData(dat)$Set)
```

    ## 
    ##     test training 
    ##       20       96

``` r
train.es <- dat[, dat$Set == "training"]
test.es <- dat[ , dat$Set == "test"]
```

Now, we can do some exploratory analysis of the data before trying some
classification methods.

``` r
# understand your data for classification
table(pData(train.es)$LnStatus)
```

    ## 
    ## neg pos 
    ##  48  48

``` r
table(pData(test.es)$LnStatus)
```

    ## 
    ## neg pos 
    ##   9  11

``` r
# understand the continuous response
tapply(pData(train.es)$LnRatio,pData(train.es)$LnStatus,summary)
```

    ## $neg
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0       0       0       0       0       0 
    ## 
    ## $pos
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0400  0.0700  0.1100  0.1935  0.2275  0.9600

``` r
tapply(pData(test.es)$LnRatio,pData(test.es)$LnStatus,summary)
```

    ## $neg
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0       0       0       0       0       0 
    ## 
    ## $pos
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.0500  0.1800  0.5000  0.4573  0.6100  0.9400

``` r
# look at the expression of 3 randomly picked genes in both training and test sets
set.seed(1234)
rangenes <- sample(1:nrow(dat), size = 3) 

# function to create tidy data table of expression and metadata
toLonger <- function(expset) {
    stopifnot(class(expset) == "ExpressionSet")
    
    expressionMatrix <- longExpressionMatrix <- exprs(expset) %>% 
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(cols = !gene, 
                   values_to = "expression",
                   names_to = "sample_id") %>%
      left_join(pData(expset), by = "sample_id")
  return(expressionMatrix)
}

toLonger(dat[rangenes,]) %>%
  ggplot(aes(y = expression, x = LnStatus)) +
    facet_wrap(Set ~ gene) +
    geom_jitter(width = 0.2, alpha = 0.5)
```

![](sm10_supervisedLearning_files/figure-gfm/eda-1.png)<!-- -->

# Classification

The prediction of a discrete response is usually refer to as
*classification*. A response taking values over a finite set of labels
is essentially the same thing as a factor. We will use the dataset from
Smeets et al. to find the *best-trained* classifier and use it to
predict the `LnStatus` of the 20 samples in the test set, i.e., classify
those as “lymph node positive” or “negative”.

## Data preprocessing

We should check to ensure there are no missing values in our data.

``` r
sum(is.na(exprs(train.es)))
```

    ## [1] 0

``` r
sum(is.na(exprs(test.es)))
```

    ## [1] 0

Here we see there are no missing values in our dataset, so we don’t have
to worry about that.

Other pre-processing operations that can be done are:  
- centering, scaling, normalizing  
- imputing missing data  
- transforming individual features (like boolean measurements)

When you are working with expression data, you may also need to use some
normalization methods to ensure your variables are comparable across all
your samples. For example, when working with count data, it is advised
to log transformed and normalize (e.g. quantile or library size
normalization) your data, so that your samples have similar
distributions.

## Feature and Model Selection

We will now use cross-validation to find the best set of features to
predict our response. Thus, I will divide the training set into 6 folds
(the authors used 10 folds). We also want the proportion of positive and
negative examples in each split to be approximately the same as for the
full data set (i.e., stratified 6-fold CV with 8 positive and 8 negative
samples within each fold). For each round of cross-validation, we use
one fold as the test data and the rest of the data as training to select
features and train different classifier.

### Cross-validation

Although it makes sense to proceed as described above, many methods are
available and many constants within methods need to be selected in these
steps. Thus, *cross-validation* is usually required to *evaluate* how
well different trained models work and select the *best* model to
proceed. Note that although you may only want to select among different
choices available in Step 2, the cross-validation needs to start in
Step 1. Why? The results of the cross-validation will be over-optimistic
and biased if the samples in the test sets of the cross-validation
(i.e., left-out folds) were used to *select* the most promising features
in Step 1!! For example, if the performance of a complex model is
(artificially) good, you may not penalize regression coefficients enough
in Step 2, and may yield to a poor performance in Step 3.

This is why, as a general rule, we **never evaluate a model on the same
data we used to train that model.**

In many studies, in the absence of a test set, cross-validation is used
to estimate performance. In those cases, *nested cross-validation* is
required! The inner cross-validation will be used to select features and
tune parameters, the outer cross-validation will be used to test each
selected model.

In this seminar, we are working with a dataset that has both a training
and a test set. Thus, we will not do a nested cross-validation. However,
keep it in mind for your project or future work, especially if you have
a small number of samples/subjects!

#### Making cross validation splits

This is not the only way to create splits of the training data to run a
cross-validation. Note that if the samples can not be evenly divided
into the nfolds you specified, then you need to complete the matrices
below with NAs and call for entries different from NA at those folds.

``` r
nfold <- 6
splitData <- function(train.es, nfold, seed = NA){
  tabTrain <- table(train.es$LnStatus)

  indlist <- sapply(names(tabTrain), function(z) which(train.es$LnStatus == z), simplify = FALSE)
  if(!is.na(seed)){
    set.seed(seed)
  }
  #Each row contains 8 pos and 8 negative samples for 6 folds
  res <- list()
  res$fold.pos <- matrix(sample(indlist[["pos"]]),nrow=nfold)
  res$fold.neg <- matrix(sample(indlist[["neg"]]),nrow=nfold)
  return(res)
}
result <- splitData(train.es, nfold, seed = 1234)
fold.pos <- result$fold.pos
fold.neg <- result$fold.neg
```

*Note*: with `CMA` you can use the command `GenerateLearningsets` to
split the training data into folds. However, it does not show you how
the data was split. Thus, you either use CMA for all or you write your
own script. Below is how we would generate the splits above using CMA:

``` r
splits <- GenerateLearningsets(y = train.es$LnStatus, method="CV", fold=6, strat= TRUE)
```

### Loop for feature selection and modeling

To illustrate how to select a model, I will use the top-50 genes
selected by `limma` (within each fold). Note that this number is very
arbitrary and other options may make more sense like using a p-value
threshold or testing different options with this CV. For this example,
I’m using only the top-50 genes as methods like LDA and Logit can not be
run on more features than samples. However, other methods like KNN or
SVM will do well with more features.

In this example, I will compare 7 different models: KNN for
k={1,5,10,15}, LDA, Logit, SVM. Feel free to add other methods to the
list!

``` r
#Define here the constants that you will not evaluate. For example, I will use the top-50 limma genes
set.seed(495)
ngenes <- 50
nmethod <- 7 #number of methods you plan to compare. 

#Define an output object here to store results
pr.err <- matrix(-1, nfold,nmethod, dimnames=list(paste0("Fold",1:nfold),c("1NN","5NN","10NN", "15NN","LDA","Logit","SVM")))

 for(i in 1:nfold){

  #Test Fold for the i-th step
  testDat.fold<-exprs(train.es)[,c(fold.pos[i,],fold.neg[i,])]
  #I will create a factor of classes for the test set of the i_th fold
  testclass.fold<-train.es$LnStatus[c(fold.pos[i,],fold.neg[i,])]
  
    
  #The rest of the samples are the training set for the i-th step
  trainDat.fold<-exprs(train.es)[,-c(fold.pos[i,],fold.neg[i,])]
  trainclass.fold<-train.es$LnStatus[-c(fold.pos[i,],fold.neg[i,])]

  #Step 1: feature selection (do you remember limma?). 

  # Note that a different set of genes will be selected for each fold! you can then compare how consistent these sets were.

  limma.dat<-as.data.frame(trainDat.fold)
  desMat <- model.matrix(~ trainclass.fold, limma.dat) #design matrix
  trainFit <- lmFit(limma.dat, desMat)
  eBtrainFit <- eBayes(trainFit)
  
  # top-50 limma genes
  top.fold <- topTable(eBtrainFit, coef = which(colnames(coef(trainFit)) != "(Intercept)"),
                       n = ngenes,sort.by="P")
  
  #Retain the top-50 limma genes from the train and test sets
  trainDat.fold <- trainDat.fold[rownames(top.fold),]
  testDat.fold <-  testDat.fold[rownames(top.fold),]

  
  #STEP 2: select a classifier
  #Set a counter for the method tested
  l <- 0

  #kNN classifiers
  for(kk in c(1,5,10,15)) {
    #every time you get inside this loop, the l counter gets redefined (i.e., 1, 2, etc for         method 1, method 2, etc)
    l <- l+1

    #knn needs samples in rows
    yhat.knn <- knn(train=t(trainDat.fold), test=t(testDat.fold), cl=trainclass.fold,
                    k = kk)
    #Store the prediction error for each kk within this fold
    pr.err[i,l]<- mean(testclass.fold != yhat.knn)
                          } #end of kNN loop

  #LDA method. Note that you can change the prior parameter to reflect a different proportion of case and control samples. The default is to use the class proportions from the training set.
  
  m.lda <- lda(x=t(trainDat.fold), group=trainclass.fold, prior=c(.5, .5))
  yhat.lda <- predict(m.lda, newdata=t(testDat.fold))$class
  pr.err[i,"LDA"] <-mean(testclass.fold != yhat.lda)
   
  #Logit
  glm.dat <- data.frame(t(trainDat.fold), group=trainclass.fold)
  
  # 50 factors still will cause optimization warnings  
  # Try without warning suppression to see 
  # To further reduce parameters, regularized regression can be used 
  # To use regularized regression uncomment lines followed by "uncomment for regularized regression" 
  suppressWarnings( m.log <- glm(group ~ ., data=glm.dat,family=binomial) ) 
  
  # uncomment for regularized regression 
  # m.log <- glmnet( t(trainDat.fold) , trainclass.fold ,family="binomial") 

  pr.log <- predict(m.log,newdata=data.frame(t(testDat.fold)),type="response")
  
  # uncomment for regularized regression 
  # pr.log <- predict(m.log,newdata=data.frame(t(testDat.fold)),type="response",newx=t(testDat.fold)) 
  
  pr.cl <- rep(0,length(testclass.fold))
  pr.cl[pr.log > 1/2] <- "pos"
  pr.cl[pr.log <= 1/2] <- "neg"

  pr.cl <- factor(pr.cl)
  pr.err[i,"Logit"] <- mean( pr.cl != testclass.fold )

  #SVM
  m.svm <- svm(x=t(trainDat.fold), y=trainclass.fold, cost=1, type="C-classification", 
               kernel="linear")
  pr.svm <- predict(m.svm,newdata=t(testDat.fold)) 
   
  pr.err[i,"SVM"] <- mean( pr.svm != testclass.fold )
  } #end of CV loop
```

### Error Rates

Now you can get the average prediction error for all methods. Note that
the prediction errors are high! not too much hope for the real test run!

``` r
cv.err <- colMeans(pr.err)

# mean - 1 sd (sd of the 6 error rates)
ls <- cv.err - apply(pr.err, 2, sd)

# mean + 1 sd (sd of the 6 error rates)
us <- cv.err + apply(pr.err, 2, sd)

# plot the results
plot(1:nmethod, cv.err, ylim=c(0, 1), xlim=c(1, (nmethod+.5)),type='n', 
axes=FALSE, xlab='Classifier', ylab='Error rate',main="6-fold CV Error")

for(j in 1:ncol(pr.err)) 
   points(jitter(rep(j, 6), factor=2), jitter(pr.err[,j]), cex=0.8, pch='X', col='gray')

for(i in 1:nmethod)
   lines(c(i, i), c(ls[i], us[i]), lwd=2, col='gray')
points(1:nmethod, ls, pch=19, col='red')
points(1:nmethod, us, pch=19, col='green')
points(1:nmethod, cv.err, pch=19, cex=1.5, col='black')
axis(2, ylab='Error rate')
axis(1, 1:nmethod, colnames(pr.err))

box()
```

![](sm10_supervisedLearning_files/figure-gfm/err-1.png)<!-- -->

### Results of the CV

According to these results, LDA and 10NN may be the better classifier to
try in the test data. However, remember that this CV results depend on
the first split of the data we did. Thus, we need to repeat this CV

**Exercise 1**: perform 100 runs of this CV before selecting a model to
test! We can do this by running different random splits using the
`splitData` function we created without a seed, and selecting the best
average error across each run of CV. Add at least one rule to select
data for use in the underlying models, such as a P value threshold for
genes selected in limma rather than just the top 50 genes by P value, or
a different cost for the SVM model.

``` r
# your code here
```

**Exercise 2**: Use AUC as a criteria to select a model based on the
training data! Tip: extract the predicted probabilities from each method
and use the calculateAUC function defined below:

``` r
require(ROCR)
calculateAUC <- function(predictions,labels){
  tmp_predictions <- ROCR::prediction(
    as.numeric(predictions),
    as.numeric(labels)
  )
  perf <- performance(tmp_predictions, measure = "auc")
  return(unlist(slot(perf,"y.values")))
}

# Example: Calculate the AUC for our SVM predictions
calculateAUC(pr.svm,testclass.fold)
```

    ## [1] 0.5625

``` r
# your code here
```

If you have never heard of AUC, please work through the rest of this
seminar, and the [Additional metrics for
evaluation](#additional-metrics-for-evaluation) section in particular,
before completing these two tasks.

## Testing the selected model

Now that we decided on which method we are going to use to classify
samples in the test set, we need to train the model using the *FULL*
training set and then classify samples of the test set. I will use the
10NN model.

``` r
yhat.knn <- knn(train=t(exprs(train.es)), test=t(exprs(test.es)), cl=train.es$LnStatus,
                     k = 10)
#Store the prediction error for each kk within this fold
pr.errTest<- mean(test.es$LnStatus != yhat.knn)
pr.errTest
```

    ## [1] 0.45

What does the prediction error mean?  
In this instance, we have evaluated how often the prediction matched the
actual lymph node status, against the total number of cases. This is the
**accuracy** metric.

Not good! In real practice, you should not keep trying until we get a
good result! In fact, you must use cross-validation on the training
dataset to evaluate different parameters and classifiers, and only
evaluate the generalizability of the *best* model on the test set.  
However, in this seminar, I encourage you to try different options **as
an exercise** and to see how much the results can change.

## CMA

Many steps of the CV defined above can be easily done with CMA. For
example, Step 1 in the loop above can also be done using ‘CMA’ with the
function ‘GeneSelection’, which selects the most informative features
(e.g., gene) to build a classifier within each of the splits generated
by ‘GenerateLearningsets’. Some learning algorithms do better if you
only give them “useful” features.

``` r
featureScores<-GeneSelection(X=t(exprs(train.es)), y=train.es$LnStatus, learningsets=splits, method="limma")
```

    ## GeneSelection: iteration 1 
    ## GeneSelection: iteration 2 
    ## GeneSelection: iteration 3 
    ## GeneSelection: iteration 4 
    ## GeneSelection: iteration 5 
    ## GeneSelection: iteration 6

``` r
#Compare list of selected genes using:
toplist(featureScores)
```

    ## top  10  genes for iteration  1 
    ##  
    ##    index importance
    ## 1   9265   28.21333
    ## 2   1702   27.11136
    ## 3  20571   25.15147
    ## 4  10254   24.12519
    ## 5  23504   20.09023
    ## 6  23567   19.72051
    ## 7  11052   18.46525
    ## 8   6936   18.45985
    ## 9  18958   18.14084
    ## 10 19526   17.92896

``` r
#We can aggregate the results across the 6 splits.

seliter<-numeric()
for(i in 1:nfold) seliter<-c(seliter, toplist(featureScores, iter=i, top = 10, show=FALSE)$index)
(sort(table(seliter), dec=T)) # summarize
```

    ## seliter
    ##  9265  1702 10254  6936 19932 23567   808  6938 14544 18958 21592 21919 23504 
    ##     6     5     4     3     3     3     2     2     2     2     2     2     2 
    ##  1478  2690  4580  5553  7183 10171 11052 13581 14249 14378 17804 18127 18918 
    ##     1     1     1     1     1     1     1     1     1     1     1     1     1 
    ## 19410 19526 19697 19821 20571 21064 21679 22343 22524 
    ##     1     1     1     1     1     1     1     1     1

``` r
# Choose the 20 probes which are chosen most commonly in the 6 splits
bestprobes<-as.numeric(names(sort(table(seliter), dec=T)))[1:20]

# examine the feature data for the best probes
fData(dat)[bestprobes, c("ID", "Gene symbol", "Gene title", "Chromosome location")]
```

    ##                        ID                        Gene symbol
    ## 212384_at       212384_at ATP6V1G2-DDX39B///SNORD84///DDX39B
    ## 1569472_s_at 1569472_s_at                      TTC3P1///TTC3
    ## 213593_s_at   213593_s_at                              TRA2A
    ## 208661_s_at   208661_s_at                      TTC3P1///TTC3
    ## 230609_at       230609_at                             CLINT1
    ## 243751_at       243751_at                               CHD2
    ## 1556088_at     1556088_at                              RPAIN
    ## 208663_s_at   208663_s_at                      TTC3P1///TTC3
    ## 222439_s_at   222439_s_at                             THRAP3
    ## 228510_at       228510_at                              ATAT1
    ## 236196_at       236196_at                             ZNF326
    ## 237746_at       237746_at                                   
    ## 243495_s_at   243495_s_at                             ZNF652
    ## 1563475_s_at 1563475_s_at                            ETFBKMT
    ## 201440_at       201440_at                              DDX23
    ## 203537_at       203537_at                            PRPSAP2
    ## 204841_s_at   204841_s_at                               EEA1
    ## 208921_s_at   208921_s_at                                SRI
    ## 213472_at       213472_at                            HNRNPH1
    ## 215220_s_at   215220_s_at                                TPR
    ##                                                                                                         Gene title
    ## 212384_at    ATP6V1G2-DDX39B readthrough (NMD candidate)///small nucleolar RNA, C/D box 84///DEAD-box helicase 39B
    ## 1569472_s_at                    tetratricopeptide repeat domain 3 pseudogene 1///tetratricopeptide repeat domain 3
    ## 213593_s_at                                                                            transformer 2 alpha homolog
    ## 208661_s_at                     tetratricopeptide repeat domain 3 pseudogene 1///tetratricopeptide repeat domain 3
    ## 230609_at                                                                                    clathrin interactor 1
    ## 243751_at                                                              chromodomain helicase DNA binding protein 2
    ## 1556088_at                                                                                 RPA interacting protein
    ## 208663_s_at                     tetratricopeptide repeat domain 3 pseudogene 1///tetratricopeptide repeat domain 3
    ## 222439_s_at                                                          thyroid hormone receptor associated protein 3
    ## 228510_at                                                                        alpha tubulin acetyltransferase 1
    ## 236196_at                                                                                  zinc finger protein 326
    ## 237746_at                                                                                                         
    ## 243495_s_at                                                                                zinc finger protein 652
    ## 1563475_s_at                                  electron transfer flavoprotein beta subunit lysine methyltransferase
    ## 201440_at                                                                                     DEAD-box helicase 23
    ## 203537_at                                             phosphoribosyl pyrophosphate synthetase associated protein 2
    ## 204841_s_at                                                                               early endosome antigen 1
    ## 208921_s_at                                                                                                 sorcin
    ## 213472_at                                                           heterogeneous nuclear ribonucleoprotein H1 (H)
    ## 215220_s_at                                                   translocated promoter region, nuclear basket protein
    ##                Chromosome location
    ## 212384_at    6p///6p21.33///6p21.3
    ## 1569472_s_at      Xq13.3///21q22.2
    ## 213593_s_at                 7p15.3
    ## 208661_s_at       Xq13.3///21q22.2
    ## 230609_at                   5q33.3
    ## 243751_at                    15q26
    ## 1556088_at                 17p13.2
    ## 208663_s_at       Xq13.3///21q22.2
    ## 222439_s_at                 1p34.3
    ## 228510_at                  6p21.33
    ## 236196_at                   1p22.2
    ## 237746_at                         
    ## 243495_s_at               17q21.32
    ## 1563475_s_at              12p11.21
    ## 201440_at                 12q13.12
    ## 203537_at              17p11.2-p12
    ## 204841_s_at                  12q22
    ## 208921_s_at                 7q21.1
    ## 213472_at                   5q35.3
    ## 215220_s_at                   1q25

This looks promising since I get TTC3 and at least a couple of other
genes that show up on Table 3 of the paper.

Similarly, you can use CMA to train and test a classifier within each CV
fold (learningsets). However, there are things you can not do within CMA
or that CMA is not doing right. For example, CMA can not do a full
nested cross-validation. Additionally, it is not trivial to train the
selected in the full dataset and then test it in the test set. CMA is
more designed for CV. Thus, it is good to know how to do this things by
hand as well.

Paul solved this problem in the following way: he made a `learningsets`
object that has just one “split” defined by the samples in the training
set.

``` r
m<-matrix(which(dat$Set == "training"), 1)

full.learningset<-new("learningsets", learnmatrix=m, method="my own", ntrain=96, iter=1)

fullFeatureScores<-GeneSelection(X=t(exprs(dat)), learningsets= full.learningset, y=dat$LnStatus, method="t.test")
```

    ## GeneSelection: iteration 1

``` r
testclassif<-classification(X=t(exprs(dat)), y=dat$LnStatus, learningsets= full.learningset, genesel=fullFeatureScores, nbgene = 100, classifier =pknnCMA, k=5)
```

    ## iteration 1

``` r
#Evaluation:
tres<-testclassif[[1]]
ftable(tres)
```

    ## number of missclassifications:  11 
    ## missclassification rate:  0.55 
    ## sensitivity: 0.545 
    ## specificity: 0.333 
    ##     predicted
    ## true 0 1
    ##    0 3 6
    ##    1 5 6

``` r
roc(tres)
```

![](sm10_supervisedLearning_files/figure-gfm/byhand-1.png)<!-- -->

Note: his optimized classifier did terribly as well.

## Multiclass learning

You won’t always have a binary learning problem, where you are
classifying a sample into 1 of 2 classes. Sometimes we might want to
train a classifier a classifier with more than two classes.  
Here we will use the *caret* and *mlbench* packages for classification
on a multi-class problem.

``` r
library(caret)
library(mlbench)
```

    ## Warning: package 'mlbench' was built under R version 4.4.3

We will be using the Soybean dataset, where our prediction is for the
different problems associated with soybean crops. Our dataset has 683
samples, and 35 features being measured for each sample.  
Our categories are 19.

``` r
cv_folds <- trainControl(method="cv", number=5)
data(Soybean)
unique(Soybean$Class)
```

    ##  [1] diaporthe-stem-canker       charcoal-rot               
    ##  [3] rhizoctonia-root-rot        phytophthora-rot           
    ##  [5] brown-stem-rot              powdery-mildew             
    ##  [7] downy-mildew                brown-spot                 
    ##  [9] bacterial-blight            bacterial-pustule          
    ## [11] purple-seed-stain           anthracnose                
    ## [13] phyllosticta-leaf-spot      alternarialeaf-spot        
    ## [15] frog-eye-leaf-spot          diaporthe-pod-&-stem-blight
    ## [17] cyst-nematode               2-4-d-injury               
    ## [19] herbicide-injury           
    ## 19 Levels: 2-4-d-injury alternarialeaf-spot anthracnose ... rhizoctonia-root-rot

``` r
summary(Soybean)
```

    ##                  Class          date     plant.stand  precip      temp    
    ##  brown-spot         : 92   5      :149   0   :354    0   : 74   0   : 80  
    ##  alternarialeaf-spot: 91   4      :131   1   :293    1   :112   1   :374  
    ##  frog-eye-leaf-spot : 91   3      :118   NA's: 36    2   :459   2   :199  
    ##  phytophthora-rot   : 88   2      : 93               NA's: 38   NA's: 30  
    ##  anthracnose        : 44   6      : 90                                    
    ##  brown-stem-rot     : 44   (Other):101                                    
    ##  (Other)            :233   NA's   :  1                                    
    ##    hail     crop.hist  area.dam    sever     seed.tmt     germ     plant.growth
    ##  0   :435   0   : 65   0   :123   0   :195   0   :305   0   :165   0   :441    
    ##  1   :127   1   :165   1   :227   1   :322   1   :222   1   :213   1   :226    
    ##  NA's:121   2   :219   2   :145   2   : 45   2   : 35   2   :193   NA's: 16    
    ##             3   :218   3   :187   NA's:121   NA's:121   NA's:112               
    ##             NA's: 16   NA's:  1                                                
    ##                                                                                
    ##                                                                                
    ##  leaves  leaf.halo  leaf.marg  leaf.size  leaf.shread leaf.malf  leaf.mild 
    ##  0: 77   0   :221   0   :357   0   : 51   0   :487    0   :554   0   :535  
    ##  1:606   1   : 36   1   : 21   1   :327   1   : 96    1   : 45   1   : 20  
    ##          2   :342   2   :221   2   :221   NA's:100    NA's: 84   2   : 20  
    ##          NA's: 84   NA's: 84   NA's: 84                          NA's:108  
    ##                                                                            
    ##                                                                            
    ##                                                                            
    ##    stem     lodging    stem.cankers canker.lesion fruiting.bodies ext.decay 
    ##  0   :296   0   :520   0   :379     0   :320      0   :473        0   :497  
    ##  1   :371   1   : 42   1   : 39     1   : 83      1   :104        1   :135  
    ##  NA's: 16   NA's:121   2   : 36     2   :177      NA's:106        2   : 13  
    ##                        3   :191     3   : 65                      NA's: 38  
    ##                        NA's: 38     NA's: 38                                
    ##                                                                             
    ##                                                                             
    ##  mycelium   int.discolor sclerotia  fruit.pods fruit.spots   seed    
    ##  0   :639   0   :581     0   :625   0   :407   0   :345    0   :476  
    ##  1   :  6   1   : 44     1   : 20   1   :130   1   : 75    1   :115  
    ##  NA's: 38   2   : 20     NA's: 38   2   : 14   2   : 57    NA's: 92  
    ##             NA's: 38                3   : 48   4   :100              
    ##                                     NA's: 84   NA's:106              
    ##                                                                      
    ##                                                                      
    ##  mold.growth seed.discolor seed.size  shriveling  roots    
    ##  0   :524    0   :513      0   :532   0   :539   0   :551  
    ##  1   : 67    1   : 64      1   : 59   1   : 38   1   : 86  
    ##  NA's: 92    NA's:106      NA's: 92   NA's:106   2   : 15  
    ##                                                  NA's: 31  
    ##                                                            
    ##                                                            
    ## 

Let us pre-process our data.

``` r
#Remove rows (samples) with missing values  
soybean_x = Soybean[rowSums(is.na(Soybean)) == 0,]
#Then remove columns (features) with missing values  
soybean_x = soybean_x[,colSums(is.na(soybean_x)) == 0]
dim(soybean_x)
```

    ## [1] 562  36

We are left with 562 samples and 35 attributes. The first column,
`Class`, describes the categories. We will refactor this column since we
have removed certain columns (and possibly some of the 19 classes)

``` r
soybean_x$Class = (as.factor(as.character(soybean_x$Class)))
unique(soybean_x$Class)
```

    ##  [1] diaporthe-stem-canker  charcoal-rot           rhizoctonia-root-rot  
    ##  [4] phytophthora-rot       brown-stem-rot         powdery-mildew        
    ##  [7] downy-mildew           brown-spot             bacterial-blight      
    ## [10] bacterial-pustule      purple-seed-stain      anthracnose           
    ## [13] phyllosticta-leaf-spot alternarialeaf-spot    frog-eye-leaf-spot    
    ## 15 Levels: alternarialeaf-spot anthracnose ... rhizoctonia-root-rot

Now we have 15 classes!

In this instance, we don’t have an external test set for evaluation, so
we will be assessing the performance on a held-out test set. First, we
create the held-out test set.  
*Note* that we are holding out this test set from our training data, and
then performing data-splitting for validation within our training
subset. In practise, like we did earlier, you would want to loop this
multiple times with holding out a test set and training on the remainder
dataset (using cross validation or bootstrapping to estimate your model
accuracy on the test set).

``` r
trainIndex <- createDataPartition(soybean_x$Class, p = .8,  list = FALSE,  times = 1)
soyTrain <- soybean_x[ trainIndex,]
soyTest  <- soybean_x[-trainIndex,]
```

**Cross validation results**  
We set up our cross-validation folds. Note we can also choose the option
‘cv’ or ‘LOOCV’ instead of ‘repeatedcv’.  
With ‘repeatedcv’, we don’t have to manually set up the multiple loops
we did when we were using CMA.

``` r
#Prepare resampling method for Cross-Validation folds  
set.seed(7)
cv_control = trainControl(method="repeatedcv", number=5, repeats=10, classProbs=FALSE) #summaryFunction=mnLogLoss)
modelSvm_cv <- train(Class~., data=soyTrain, method="svmRadial", metric="Accuracy", trControl=cv_control, trace=FALSE)
# display results
print(modelSvm_cv)
```

    ## Support Vector Machines with Radial Basis Function Kernel 
    ## 
    ## 452 samples
    ##  35 predictor
    ##  15 classes: 'alternarialeaf-spot', 'anthracnose', 'bacterial-blight', 'bacterial-pustule', 'brown-spot', 'brown-stem-rot', 'charcoal-rot', 'diaporthe-stem-canker', 'downy-mildew', 'frog-eye-leaf-spot', 'phyllosticta-leaf-spot', 'phytophthora-rot', 'powdery-mildew', 'purple-seed-stain', 'rhizoctonia-root-rot' 
    ## 
    ## No pre-processing
    ## Resampling: Cross-Validated (5 fold, repeated 10 times) 
    ## Summary of sample sizes: 361, 362, 362, 361, 362, 361, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   C     Accuracy   Kappa    
    ##   0.25  0.6945666  0.6501528
    ##   0.50  0.8793964  0.8645924
    ##   1.00  0.9167649  0.9068591
    ## 
    ## Tuning parameter 'sigma' was held constant at a value of 0.06506239
    ## Accuracy was used to select the optimal model using the largest value.
    ## The final values used for the model were sigma = 0.06506239 and C = 1.

``` r
y_pred_cv = predict(modelSvm_cv, soyTest[,-1])
```

We can also use bootstrap re-sampling instead of cross-validation folds.
This means we take random samples from the dataset (with re-selection),
against which to evaluate our model. This gives us an idea about the
variance of the model itself.

Bootstrapping is different from cross validation in that the latter
splits the entire dataset into folds *without re-selection*.It is a
robust method to estimate the accuracy of our model.

**Bootstrapped results**

This part will take some time.

``` r
#Prepare resampling method for Bootstrapping

set.seed(7)
cv_boot = trainControl(method="boot", number=20, classProbs=FALSE) 
modelSvm_boot <- train(Class~., data=soyTrain, method="svmRadial",
                       metric="Accuracy", trControl=cv_boot, trace=FALSE)
# display results
print(modelSvm_boot)
```

    ## Support Vector Machines with Radial Basis Function Kernel 
    ## 
    ## 452 samples
    ##  35 predictor
    ##  15 classes: 'alternarialeaf-spot', 'anthracnose', 'bacterial-blight', 'bacterial-pustule', 'brown-spot', 'brown-stem-rot', 'charcoal-rot', 'diaporthe-stem-canker', 'downy-mildew', 'frog-eye-leaf-spot', 'phyllosticta-leaf-spot', 'phytophthora-rot', 'powdery-mildew', 'purple-seed-stain', 'rhizoctonia-root-rot' 
    ## 
    ## No pre-processing
    ## Resampling: Bootstrapped (20 reps) 
    ## Summary of sample sizes: 452, 452, 452, 452, 452, 452, ... 
    ## Resampling results across tuning parameters:
    ## 
    ##   C     Accuracy   Kappa    
    ##   0.25  0.7421922  0.7060114
    ##   0.50  0.8776802  0.8617766
    ##   1.00  0.9005277  0.8878170
    ## 
    ## Tuning parameter 'sigma' was held constant at a value of 0.06506239
    ## Accuracy was used to select the optimal model using the largest value.
    ## The final values used for the model were sigma = 0.06506239 and C = 1.

``` r
y_pred_boot = predict(modelSvm_boot, soyTest[,-1])
```

**Fitting other models**

This part will take a while. You can consider other algorithmic
approaches for classifiers too, as follows:

``` r
# trace = FALSE or vocal = FALSE silences these functions depending on the method. it's stupid but blame the train function

# train the GBM model
set.seed(7)
modelGbm <- train(Class~., data=soyTrain, method="gbm", metric="Accuracy",
                  trControl=cv_boot, verbose = FALSE, preProcess = "zv")
# train the SVM model
set.seed(7)
modelSvm <- train(Class~., data=soyTrain, method="svmRadial",
                  metric="Accuracy", trControl=cv_boot, trace=FALSE)

# collect resamples
results <- resamples(list(GBM=modelGbm, SVM=modelSvm))
# summarize the distributions
summary(results)
```

    ## 
    ## Call:
    ## summary.resamples(object = results)
    ## 
    ## Models: GBM, SVM 
    ## Number of resamples: 20 
    ## 
    ## Accuracy 
    ##          Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
    ## GBM 0.8544304 0.8819049 0.8929854 0.8951865 0.9134987 0.9512195    0
    ## SVM 0.8711656 0.8876165 0.8960235 0.9005277 0.9125946 0.9393939    0
    ## 
    ## Kappa 
    ##          Min.   1st Qu.    Median      Mean   3rd Qu.      Max. NA's
    ## GBM 0.8373613 0.8667342 0.8810342 0.8818997 0.9023763 0.9450126    0
    ## SVM 0.8549330 0.8716083 0.8827638 0.8878170 0.9021221 0.9331632    0

``` r
# boxplots of results
bwplot(results)
```

![](sm10_supervisedLearning_files/figure-gfm/othermodels-1.png)<!-- -->

> You can finetune the parameters for different models in Caret using
> tuneGrid.

### Additional metrics for evaluation

We can also consider other metrics to assess the performance of a model.
This becomes particularly important when, for example, your your test
set is imbalanced. In that scenario, evaluating the accuracy of a model
might not be the best indication that your classifier works well. In
fact, it may be biased by the over-represented class in your test set.

#### Kappa-score

Overall accuracy is a misleading statistic in case of unbalanced
datasets. The kappa statistic overcomes this by taking the expected
error into account.

#### ROC, Precision & Recall

Receiver-Operator Characteristic (ROC) Curves can be used to
characterize the performance of our model in a binary classification
setting.  
For binary classification, we might also be interested in precision and
recall, i.e. a metric of our True Positive and True Negative rates.

``` r
binary_preds = as.integer(y_pred_cv)-1
binary_true = as.integer(soyTest$Class)-1
precision <- sum(binary_preds & binary_true) / sum(binary_preds)
recall <- sum(binary_preds & binary_true) / sum(binary_true)
```

Precision tells us how often, when we predict a class ‘y’, do we get it
correct.  
Recall tells us how often, out of all our ‘y’ instances, do we predict
them correctly.

``` r
#Precision    
print(precision)
```

    ## [1] 0.1528998

``` r
#Recall  
print(recall)
```

    ## [1] 0.1389776

In case of multi-class classification, another metric that comes in
handy is the F1-score. This is the harmonic mean of precision and
recall. A high F1-score will tell you that you are quite precise and
sensitive in your prediction of *all* classes.

``` r
Fmeasure <- 2 * precision * recall / (precision + recall)
Fmeasure
```

    ## [1] 0.1456067

> In the models we fit for our multiclass problem, what was the Kappa
> statistic on our test set? What was the Accuracy?  
> Which model selection appraoch worked better, the bootstrapped SVM or
> the repeated 5-CV approach? How do we assess this (on the test set, or
> from the cross-validation results?)

## Discussion points

How good do you think a classifier will have to be to be clinically
relevant? What level of specificity or sensitivity do you think is “good
enough”? The author’s data set has half lymph-node positive and half
negative. The incidence of lymph-node-positive breast cancer is about
33% at the time of diagnosis (according to
\[<http://seer.cancer.gov/statfacts/html/breast.html>\]). How does that
affect your opinion about the classifier?

# References

1: Varma S, Simon R: Bias in error estimation when using
cross-validation for model selection. BMC Bioinformatics 2006, 7:91.

2: Statnikov A, Aliferis CF, Tsamardinos I, Hardin D, Levy S: A
comprehensive evaluation of multicategory classification methods for
microarray gene expression cancer diagnosis. Bioinformatics 2005,
21:631-643.
