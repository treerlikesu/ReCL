rm(list = ls(all = TRUE))

Simu <- function(N, NS, t.fix){
  
  # testing data
  N.test = 5000
  data.all.test <- data.generation.bin(N = N.test, tau = 4, gen = gen)
  data.test <- data.frame(opt.treat = data.all.test$g.opt,
                          x1 = data.all.test$x[,1], x2 = data.all.test$x[,2],
                          x3 = data.all.test$x[,3], A = data.all.test$A, 
                          num = data.all.test$num)
  
  acc.test <- matrix(NA, NS, 7)
  num.test <- matrix(NA, NS, 7)
  
  set.seed(777)
  for(iter in 1:NS){
    print(iter)
    
    pst <- PS.bin(data.all, PS = 'SL')
    psf <- PS.bin(data.all, PS = 'WRONG')
    
    # CL
    ## IPW
    res.sl <- CLrec.bin(data.all, pst, t.fix = t.fix, mod.class = 'label ~ x1 + x2', ICW = 'IPW')
    test.sl <- CLrec.bin.test(res.sl, data.test = data.test)

    ## AIPW-PST
    res.drt <- CLrec.bin(data.all, pst, t.fix = t.fix, mod.class = 'label ~ x1 + x2', ICW = 'AIPW')
    test.drt <- CLrec.bin.test(res.drt, data.test = data.test)

    ## AIPW-PSF
    res.drf <- CLrec.bin(data.all, psf, t.fix = t.fix, mod.class = 'label ~ x1 + x2', ICW = 'AIPW')
    test.drf <- CLrec.bin.test(res.drf, data.test = data.test)

    ## SMR
    res.smr <- SMR.bin(data.all, t.fix = t.fix)
    test.smr <-SMR.bin.test(res.smr, data.test)

    # First
    res.sur <- CLrec.bin(data.all, pst, t.fix = t.fix, mod.class = 'label ~ x1 + x2', ICW = 'First')
    test.sur <- CLrec.bin.test(res.sur, data.test = data.test)
    
    # random
    est.rd <- rbinom(N.test, 1, 0.5)
    test.rd <- test.fun(est.rd, data.test)
    
    # Optimal
    test.opt <- test.fun(data.test$opt.treat, data.test)
    
    acc.test[iter, ] <- c(test.sl$est.acc, test.drt$est.acc, test.drf$est.acc, test.smr$est.acc,
                          test.sur$est.acc, test.rd$est.acc, test.opt$est.acc)
  
    num.test[iter,] <- c(mean(test.sl$est.num), mean(test.drt$est.num), mean(test.drf$est.num),
                         mean(test.smr$est.num), mean(test.sur$est.num),
                         mean(test.rd$est.num), mean(test.opt$est.num))
   
  }
  colnames(acc.test) <- c('CL-IPW', 'CL-AIPW-T', 'CL-AIPW-F', 'CL-SMR', 'First', 'Random', 'Optimal')
  colnames(num.test) <- c('CL-IPW', 'CL-AIPW-T', 'CL-AIPW-F', 'CL-SMR', 'First', 'Random', 'Optimal')
  
  return(list(acc.test = acc.test, num.test = num.test))
}

N = 400; NS = 100; t.fix = 3
res1 <- Simu(N, NS, t.fix = t.fix)
N = 600; NS = 100; t.fix = 3
res2 <- Simu(N, NS, t.fix = t.fix)

