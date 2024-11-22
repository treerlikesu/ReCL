library(reda)
library(purrr)
library(dplyr)
library(reReg)
library(nleqslv)
library(frailtypack)
library(mvtnorm)
library(numDeriv)
library(Rsolnp)
library(optimx)
library(BB)
library(rootSolve)
library(survival)
library(MASS)
#library(dummies)
library(ggplot2)
library(ggpubr)
library(parallel)
library(SuperLearner)
library(gam)
library(randomForest)
library(rpart)
library(tidyr)
library(sft)

expit <- function(x) {exp(x)/(1 + exp(x))}

data.generation.bin <- function(N, tau){
  
  x1 <- rnorm(N, mean = 0, sd = 1) 
  x2 <- rnorm(N, mean = 0, sd = 1) 
  x3 <- rnorm(N, mean = 0, sd = 1)
  x <- data.frame(x1, x2, x3)
    
  prob <- expit(0.3*x1 - 0.5*x2)
  A <- rbinom(N, 1, prob)
  g.opt <- (x1 > - 1)*(x2 > -0.5)
  funR <- x2 + abs(1.5*x1 - 0.5)*(A - g.opt)^2 - 0.8
  rr <- simEventData(z = cbind(1, funR), zCoef = c(0, 1), rho = 0.5, endTime = tau, nProcess = N)
  cens <- runif(N, tau - 1, tau)

  rr <- cbind(rr, x[rr$ID,], A = A[rr$ID], cens = cens[rr$ID])
  data <- rr
  data$time <- apply(cbind(data$time, data$cens), 1, min)
  event <- ifelse(data$time < data$cens, 1, 0)
  data$event <- event
  data <- data[!duplicated(data$time, fromLast = T), ]
  data$t.stop <- data$time
  data$t.start <- as.vector(unlist(lapply(split(data$time, data$ID), function(x) c(0, x[-length(x)]))))
  data$id <- data$ID
  num <- as.vector(table(data$id))
  data <- data[,-1]
  data.sur <- data %>% filter(t.start == 0)
  
  return(list(data = data, data.sur = data.sur, A = A, x = x, cens = cens, g.opt = g.opt, prob = prob, num = num))
}

mod = '~ x1 + x2 + A*x1 + A*x2'
SMR.bin <- function(data.all, t.fix, mod = '~ x1 + x2 + A*x1 + A*x2'){
  
  # Data summary
  Data <- data.all$data
  opt.treat <- data.all$g.opt
  cens <- data.all$cens
  A <- data.all$A
  x <- data.all$x
  x1 <- x$x1
  x2 <- x$x2
  x3 <- x$x3
  Ax1 <- x1*A
  Ax2 <- x2*A
  Ax3 <- x3*A
  prob <- data.all$prob
  
  # Estimation
  formula <- as.formula(paste('Recur(t.start %to% t.stop, id, event)', mod))
  fm <- with(Data, formula)
  fit <- reReg(fm, data = Data, model = "cox.LWYY", B = 0)
  coef <- fit$par1
  return(coef)
}

h_function=function(x){
  exp(x)
}

h_deri_function=function(x){
  exp(x)
}

Baseline=function(t, Data, coef, x1, x2, AA){
  
  Obs_time_vec=Data[,"t.stop"]
  Delta_vec=Data[,"event"]
  Cen_time_vec=c()
  for(i in 1:length(table(Data[,"id"]))){
    location=which(Data[,"id"]==i & Data[,"event"]==0)
    Cen_time_vec=c(Cen_time_vec,Data[location,"t.stop"])
  }
  number=as.vector(table(Data[,"id"]))
  if(t==0){
    total=0
  }else{
    V1=unique(sort(c(0,Obs_time_vec[which(Obs_time_vec<=t & Delta_vec==1   )])))####& Delta_vec==1
    total=sum(apply(matrix(1:length(V1),1,length(V1)),2,FUN=function(ind){
      df=data.frame(x=ifelse(Obs_time_vec==V1[ind] & Delta_vec==1,1,0),id=rep(1:length(Cen_time_vec),number))
      denominator=sum(aggregate(df$x, by=list(df$id), FUN=sum)$x)
      numerator=sum(as.numeric(V1[ind]<=Cen_time_vec)*exp((coef[1]*x1)+(coef[2]*x2)+(coef[3]*AA)+(coef[4]*AA*x1)+(coef[5]*AA*x2)))
      ans=(denominator/numerator)
    }## for function
    ))
  }
  return(total)
}

PS.bin <- function(data.all, PS = 'SL'){
  
  # Data summary
  A <- data.all$A
  x <- data.all$x
  x1 <- x$x1
  x2 <- x$x2
  x3 <- x$x3
  prob <- data.all$prob
  
  # Fit the propensity score function
  if(PS == 'SL'){
    PS.library <- c("SL.knn", "SL.randomForest","SL.glm")
    #xx <- data.frame(x1, x2, x3)
    PS.SL <- SuperLearner(A, x, newX=x, family=binomial(), 
                          SL.library=PS.library, method = "method.NNLS",control = list(), cvControl = list())
    ps <- PS.SL$SL.predict
  }
  if(PS == 'GLM'){
    xx = data.frame(A, x1, x2, x3)
    PS.mod <- glm(A ~ x1 + x2, family = "binomial", data = xx) ## fit logistic regression
    ps <- PS.mod$fitted.values
  }
  if(PS == 'TRUE'){
    ps <- prob
  }
  if(PS == 'WRONG'){
    expx3 <- exp(x3)
    xx = data.frame(A, x1, x2, expx3)
    PS.mod <- glm(A ~ x1 + expx3, family = "binomial", data = xx) ## fit logistic regression
    ps <- PS.mod$fitted.values
  }
  
  return(ps)
}

CLrec.bin <- function(data.all, ps, t.fix, mod.class = 'label ~ x1 + x2', ICW = 'PO', method = 'tree') {
  
  # Data summary
  Data <- data.all$data
  Data.sur <- data.all$data.sur
  opt.treat <- data.all$g.opt
  cens <- data.all$cens
  A <- data.all$A
  x <- data.all$x
  x1 <- x$x1
  x2 <- x$x2
  x3 <- x$x3
  Ax1 <- x1*A
  Ax2 <- x2*A
  Ax3 <- x3*A
  prob <- data.all$prob
  
  # Estimate contrast/weight for data space expansion 
  class <- sort(unique(A))
  num <- length(class)
  mus <- matrix(NA, length(A), num)
  if (ICW == 'PO'){
    # IPTW: Pseudo outcome
    Pseudo_value_IPW=apply(matrix(c(1:N),1,N),2,function(t){ 
      location=which(Data[,"id"]==t)
      delete_data=Data[-location,]
      Output1_MCF=mcf(Recur(t.stop,id,  event) ~ 1, data = Data)@MCF
      Output2_MCF=mcf(Recur(t.stop,id,  event) ~ 1, data = delete_data)@MCF
      CMF1<-stepfun(Output1_MCF$time, c(0,Output1_MCF$MCF))
      CMF2<-stepfun(Output2_MCF$time, c(0,Output2_MCF$MCF))
      difference=(N*(CMF1(t.fix)))-((N-1)*(CMF2(t.fix)))
    }
    )
    mus[,1] <- (1 - A)*Pseudo_value_IPW/(1-ps)
    mus[,2] <- A*Pseudo_value_IPW/ps
  }
  
  if (ICW == 'DR'){
    # IPTW: Pseudo outcome
    Pseudo_value_IPW=apply(matrix(c(1:N),1,N),2,function(t){ 
      location=which(Data[,"id"]==t)
      delete_data=Data[-location,]
      Output1_MCF=mcf(Recur(t.stop,id,  event) ~ 1, data = Data)@MCF
      Output2_MCF=mcf(Recur(t.stop,id,  event) ~ 1, data = delete_data)@MCF
      CMF1<-stepfun(Output1_MCF$time, c(0,Output1_MCF$MCF))
      CMF2<-stepfun(Output2_MCF$time, c(0,Output2_MCF$MCF))
      difference=(N*(CMF1(t.fix)))-((N-1)*(CMF2(t.fix)))
    }
    )
    
    # Estimation using SMR model
    formula <- as.formula(paste('Recur(t.start %to% t.stop, id, event)', mod))
    fm <- with(Data, formula)
    fit <- reReg(fm, data = Data, model = "cox.LWYY", B = 0)
    coef <- fit$par1
    
    mus[,1] <- (1 - A)*Pseudo_value_IPW/(1-ps) + 
      (A-ps)*(1/(1-ps))*Baseline(t.fix, Data, coef, x1, x2, AA = 0)*exp(coef[1]*x1 + coef[2]*x2)
    mus[,2] <- A*Pseudo_value_IPW/ps + 
      (1-A/ps)*Baseline(t.fix, Data, coef, x1, x2, AA = 1)*exp(coef[1]*x1 + coef[2]*x2 + coef[3]*1 + coef[4]*x1 + coef[5]*x2)
  }
  
  if (ICW == 'Sur'){
    # IPTW: Pseudo outcome
    Pseudo_value_IPW=apply(matrix(c(1:N),1,N),2,function(t){ 
      location=which(Data.sur[,"id"]==t)
      delete_data=Data.sur[-location,]
      CMF1 <- estimateNAH(Data.sur$time, Data.sur$event)$H
      CMF2 <- estimateNAH(delete_data$time, delete_data$event)$H
      difference=(N*(CMF1(t.fix)))-((N-1)*(CMF2(t.fix)))
    }
    )
    mus[,1] <- (1 - A)*Pseudo_value_IPW/(1-ps)
    mus[,2] <- A*Pseudo_value_IPW/ps
  }
  # Estimate contrast_IPTW
  Contrast <- abs(mus[,2] - mus[,1]) 
  
  label <- ifelse(mus[,2] > mus[,1], 0, 1)
  mod.class <- as.formula(mod.class)
  rpData <- data.frame(x, label)
  
  if(method == 'tree'){
    # fit CART 
    tree <- rpart(mod.class, weights = Contrast, method = "class", data = , na.action = na.omit)
    
    # estimate the OTR
    est.treat <- as.vector(as.numeric(predict(tree, type = "class", newdata = data.frame(x))) - 1)
  }
  
  return(list(tree = tree, est.treat = est.treat))
}

CLrec.bin.test <- function(res, mod.class = '~ x1 + x2', data.test){
  
  Hclassify <- cbind(get_all_vars(mod.class, data.test))
  est.treat <- as.numeric(predict(res$tree, type = 'class', newdata = Hclassify)) - 1
  re <- test.fun(est.treat, data.test)
  return(re)
}

SMR.bin.test <- function(res, data.test){
  
  number <- length(unique(data.test$A))
  mus <- matrix(NA, length(data.test$A), number)
  mus[,1] <- exp(res[1]*data.test$x1 + res[2]*data.test$x2)
  mus[,2] <- exp(res[1]*data.test$x1 + res[2]*data.test$x2 + res[3]*1 + res[4]*data.test$x1 + res[5]*data.test$x2)
  est.treat <- ifelse(mus[,2] < mus[,1], 1, 0)
  re <- test.fun(est.treat, data.test)
  return(re)
}

test.fun <- function(est.treat, data.test, rho = 0.8){
  
  est.acc <- mean(est.treat == data.test$opt.treat)
  est.con <- data.test$num[est.treat == data.test$A]
  est.non <- data.test$num[est.treat != data.test$A]
  freq.con <- c(table(est.con)[1:9],sum(table(est.con)[-c(1:9)]))
  freq.non <- c(table(est.non)[1:9],sum(table(est.non)[-c(1:9)]))
  names(freq.con) <- c(1:10)
  names(freq.non) <- c(1:10)
  funR <- exp(data.test$x2 + abs(1.5*data.test$x1 - 0.5)*(est.treat - data.test$opt.treat)^2 - 0.8)
  est.num = mean(rho * t.fix * funR)
  return(list(est.acc = est.acc, est.num = est.num, est.con = est.con, est.non = est.non,
              freq.con = freq.con, freq.non = freq.non))
}

