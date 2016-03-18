## Fonction de calcul ds l'AR , theta1, etc...
## 10/05/2015  : New version  : f(x, Z) , g(x)  x= factors,  z = weather factors...
##             :  Correction computing theta (denominator)
## 15/05/2015  : Adapted to parallel tratment
##             : Correction : computing g with abs(ei) on same regressors. 

## 03/06/2015  :Version 4 with computation of bandwidth also for g*
## 02/11/2015 : Updated set of regressors (U = region) 
## 05/11/2015 : Version with direct introduction of  regressors  
## 12/02/2016 : NpRiskAversion5 was computing G() qwith residuals, not residuals squared
## 12/02/2016 : New function NpRiskAversion6 wityh new computation of g based on squared residuals !
## 15/02/2016 : Correction of the theta2 coef (Milk price should be divided per 1000: implicitely done in theta1 (grain should also so it cancelels in theta1) )
#  17/02/2016 : Correction of the computation of sigma, estimated NP

 ##
NpAversion6 <- function(year, 
                        Xvar = x,
                        bwmethod = "CV",
                        bwcompute = "FALSE",
                        samplereduc = 1 ,
                        data = data.all,
                        ampli = 1,
                        ...
) 
{
  options(np.messages=FALSE)
  # We take observations for the considered year 
  data.work <- subset(data,  annee == year )
  N <- nrow(data.work)
  print(N)
  
  # we need to restrict the dataframe to complete observations
  
  data.work <- data.work[complete.cases(data.work),]
  
  ####       ############ For testing  : subset of the data 
  if(samplereduc < 1) {
                set.seed(0924)
                S <-  samplereduc * nrow(data.work)      # Sample size of the training file
                ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
                data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
                  
  }    
  # #       #############
  attach(data.work)
  print(" --- Sample Size for that year ---")
  print(nrow(data.work))
  print(" ---------------")
  # NEW : All the variable are in here 
  W <- as.data.frame(data.work[,Xvar])           
  
  
  ##################   Estimation ###################  
  # New model without  "prixlaiteriemoyen" 
  ## Step 1 -  bandwidth for f 
  
  tic  <- proc.time()
  
  
  if(bwmethod == "CV"){
    if(bwcompute != "TRUE") {
      # 11/05/2015  we load the cross-validated bandwidths
      #load (paste("Results/CV/bw.f.",year,".RData", sep=""))
      load(paste(Myroot,"Results/bw.f.",bwmethod,".",year,".RData", sep=""))
      
    }
    if(bwcompute == "TRUE") {
      bw.f <- npregbw(ydat = data.work$laitproduit, 
                      xdat = W, 
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      ckertype = "epanechnikov",
                      ukertype = "liracine",
                      okertype = "liracine",
                      #  Options to change below 
                      #                       nmulti = 5, 
                      #                       tol = 0.1,
                      #                       ftol = 0.1,
                      #                       itmax = 20, 
                      data = data.work )
      
    }
  }
  
  if(bwmethod == "FIXED") {
    # A déterminer 
    bw.fixed.f <- c(3.766250e+06, 4.150987e+05, 2.359682e+06, 4.393298e+06,
                    1.204398e+06, 8.197868e+07, 1.394551e+09, 1.174233e+00,
                    6.407930e+04, 1.398660e+05)
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = W,                       
                    regtype = "ll",
                    bws = bw.fixed.f,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  if(bwmethod == "SILVERMAN") {
    bw.silver.f <-  apply(W, 2, silverman.bw)
    
    # A déterminer si on modifie la fenêtre de Silverman ... (ici x 2) 
    bw.silver.f <- ampli*bw.silver.f  ### <<<<<<<------- amplificateur  #####
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = W,                       
                    regtype = "ll",
                    bws = bw.silver.f,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  # bandwidth information   
  print(" --- bandwith used for F ---")
  print(" ---------------")
  print(summary(bw.f))
  
  tac  <- proc.time()
  duree.bw <-tac-tic
  
  save(bw.f, duree.bw,  file = paste(Myroot,"Results/bw.f.",bwmethod,".",year,".RData", sep=""))
  
  # npplot(bw.f,  plot.errors.method="bootstrap")
  ##  Np estimation of f  (need to integrate to data.work)  
  
  f.np <- npreg( bws = bw.f,
                 gradients = TRUE,
                 residuals = TRUE)
  
  
  ## Step 2 -- Computing g
  # Residuals Squared !! (12/02/16 - Discussion with M. Simioni)
  
  e.2 <-(f.np$resid)^2
  
  
  # -- Computing bw for g  (same variable than for f )
  if(bwmethod == "CV") {
    if(bwcompute != "TRUE") {
      # 11/05/2015  we load the cross-validated bandwidths
      load (paste("Results/bw.g.",bwmethod,".",year,".RData", sep=""))
    }
    if(bwcompute == "TRUE") {
      
      bw.g <- npregbw(ydat = e.2, 
                      xdat = W,  
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      ckertype = "epanechnikov",
                      ukertype = "liracine",
                      okertype = "liracine",
                      #                       nmulti = 5, 
                      #                       tol = 0.1,
                      #                       ftol = 0.1,
                      #                       itmax = 20, 
                      data = data.work )
      
    }
  }
  
  if(bwmethod == "FIXED") {
    # A déterminer 
    bw.fixed.g <- c(1.522843e+07, 2.113754e+06 ,1.407951e+07, 5.513004e+06, 9.820912e+05 ,
                    1.740240e+08 , 2.020195e+08, 1.940852e-02 , 6.899923e+00, 7.041676e+00)
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = e.2, 
                    xdat = W,    
                    regtype = "ll",
                    bw <- bw.fixed.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  if(bwmethod == "SILVERMAN") {
    
    bw.silver.g <-  apply(W, 2, silverman.bw)
    # A déterminer si on modifie la fenêtre de Silverman .<<<---- ACHTUNG A MODIFIER ? #####
    bw.silver.g <- ampli*bw.silver.g 
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = e.2, 
                    xdat = W,         
                    regtype = "ll",
                    bws = bw.silver.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  
  # bandwidth information   
  print(" --- bandwith used for G  ---")
  print(" ---------------")
  print(summary(bw.g))
  tac  <- proc.time()
  duree.bw <-tac-tic
  save(bw.g, duree.bw, file = paste(Myroot,"Results/bw.g.",bwmethod,".",year,".RData", sep=""))
  
  
  ##  Np estimation of g 
  g.np <- npreg( bws = bw.g,
                 data = data.work,
                 gradients = TRUE,
                 residuals = TRUE )
  # summary(g.np)
  
  save(bw.f, bw.g, f.np, g.np,  file = paste(Myroot,"Results/np.",bwmethod,".",year,".RData", sep=""))
  
  
  #### computing  theta's  ###
  
  # step 4 --  derivative of f and  g 
  # first retrive the names of the variables
  nom <- f.np$xnames
  idx.grain <- match("concvlr", nom) # <<<- correspond to  concvlr
  idx.irri <- match("hasfpirri", nom) # <<<- correspond to  hasfirri
  
  # gradient for  f (grain and irri)
  f1.np <- f.np$grad[,  idx.grain]  
  f2.np <- f.np$grad[, idx.irri]  
  
  # Step 4 -- computing theta1 
  w1 <- prixconcvlr
  p <-  prixlaiteriemoyen
  
  #  12/02/2016: Def de Theta1  sans  signe NEGATIF 
  # gradient for  g (grain)  new computation  is required s ince regression was made on e.2, 
  
  theta1 <- (f1.np - (w1/p)) *  2 * sqrt(g.np$grad[,idx.grain])
  
  #       --- computing theta2
  # ordre de grandeur pour le coût marginal de l'irrigation : environ 15 euros par hectare 
  # c'est une info moyenne pour la campagne 2010 pour la région du sud-ouest.
  
  # 12/02/2016: Def de Theta2  sans  signe NEGATIF  
  # gradient for  g (irri)  new computation  is required s ince regression was made on e.2, 
  
  w2 <- 15    #  Le prix est en euros par litre doit être ramené aux 1000 l comme pour grain ! (SC et CB  16/02/2016 )
  theta2 <-  (f2.np - (w2/(p/1000))) * 2 * sqrt(g.np$grad[, idx.irri])
  
  # computing theta 
  theta <- (theta1  + theta2) /2
  
  # Computing PROFIT  and sigma  (CB 27/11/2015)
  # Profit exist already  !     

  # Computing sigma as in the paper (CB 17/02/2016)
  
  Sigma2 <-  p^2 * g.np$mean  #  estimates of the regression function at the evaluation points
  SigmaProf <- sqrt(Sigma2)
  # Conforme au papier
  AR <- - theta/SigmaProf
  
  # Computing RP ( Checked with Stephane : 12/02/2016)
  
  RP <-  0.5 * AR * Sigma2 * mean(Profit) 
  RP.pc <- RP/Profit
  
  
  ####  RESULTS ############
  summary(theta)
  summary(AR)
  
  ###returned results ####
  
  #     risk.results <- list(sample.used =  data.work, 
  #                        theta = theta, 
  #                        AR = AR )  
  #     
  risk.results <- as.data.frame(cbind(ident, 
                                      annee,
                                     # region,     Since already in data, would cause region.x pb
                                      f.np$mean, 
                                      g.np$mean, 
                                      f1.np,
                                      f2.np,
                                      theta1,
                                      theta2,
                                      theta,
                                     # Profit,    Since already in data, would cause Profit.x pb
                                      SigmaProf,
                                      AR,
                                      RP, 
                                      RP.pc)) 
  save( risk.results, file = paste(Myroot,"Results/Risknp.",bwmethod,".",year,".RData", sep=""))
  
  return(risk.results)
  
}








# ======================OLD VERSION ==========
# This function was computing g based on the residuals (not the residuals squared !!)
NpAversion5 <- function(year, 
                        Xvar = x,
                        bwmethod = "CV",
                        bwcompute = "FALSE",
                        data = data.all,
                        ampli = 1,
                        ...
) 
{
  options(np.messages=FALSE)
   # We take observations for the considered year 
  data.work <- subset(data,  annee == year )
  N <- nrow(data.work)
  print(N)

  # we need to restrict the dataframe to complete observations
  
  data.work <- data.work[complete.cases(data.work),]
  
# # #   #       ############ For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
#               set.seed(0924)
#               S <-  100   # Sample size of the training file
#               ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
#               data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
#                 
# # #   #       
  # #       #############
  attach(data.work)
              
  # NEW : All the variable are in here 
  W <- as.data.frame(data.work[,Xvar])           
  
 
  ##################   Estimation ###################  
  # New model without  "prixlaiteriemoyen" 
  ## Step 1 -  bandwidth for f 
  
  tic  <- proc.time()
  
  
  if(bwmethod == "CV"){
    if(bwcompute != "TRUE") {
      # 11/05/2015  we load the cross-validated bandwidths
      #load (paste("Results/CV/bw.f.",year,".RData", sep=""))
      load(paste(Myroot,"Results/bw.f.",bwmethod,".",year,".RData", sep=""))
      
    }
    if(bwcompute == "TRUE") {
      bw.f <- npregbw(ydat = data.work$laitproduit, 
                      xdat = W, 
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      ckertype = "epanechnikov",
                      ukertype = "liracine",
                      okertype = "liracine",
                      #  Options to change below 
                      #                       nmulti = 5, 
                      #                       tol = 0.1,
                      #                       ftol = 0.1,
                      #                       itmax = 20, 
                      data = data.work )
      
    }
  }
  
  if(bwmethod == "FIXED") {
    # A déterminer 
    bw.fixed.f <- c(3.766250e+06, 4.150987e+05, 2.359682e+06, 4.393298e+06,
                    1.204398e+06, 8.197868e+07, 1.394551e+09, 1.174233e+00,
                    6.407930e+04, 1.398660e+05)
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = W,                       
                    regtype = "ll",
                    bws = bw.fixed.f,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  if(bwmethod == "SILVERMAN") {
    bw.silver.f <-  apply(W, 2, silverman.bw)
    
    # A déterminer si on modifie la fenêtre de Silverman ... (ici x 2) 
    bw.silver.f <- ampli*bw.silver.f  ### <<<<<<<------- amplificateur  #####
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = W,                       
                    regtype = "ll",
                    bws = bw.silver.f,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  # bandwidth information   
  print(" --- bandwith used for F ---")
  print(" ---------------")
  print(summary(bw.f))
  
  tac  <- proc.time()
  duree.bw <-tac-tic
  
  save(bw.f, duree.bw,  file = paste(Myroot,"Results/bw.f.",bwmethod,".",year,".RData", sep=""))
  
  # npplot(bw.f,  plot.errors.method="bootstrap")
  ##  Np estimation of f  (need to integrate to data.work)  
  
  f.np <- npreg( bws = bw.f,
                 gradients = TRUE,
                 residuals = TRUE)
  
  
  ## Step 2 -- Computing g
  # Residuals Squared !! (12/02/16 - Discussion with M. Simioni)
  # e.f <- abs(f.np$resid)
  e.f <-(f.np$resid)^2
  
  
  # -- Computing bw for g  (same variable than for f )
  if(bwmethod == "CV") {
    if(bwcompute != "TRUE") {
      # 11/05/2015  we load the cross-validated bandwidths
      load (paste("Results/bw.g.",bwmethod,".",year,".RData", sep=""))
    }
    if(bwcompute == "TRUE") {
      
      bw.g <- npregbw(ydat = e.f, 
                      xdat = W,  
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      ckertype = "epanechnikov",
                      ukertype = "liracine",
                      okertype = "liracine",
                      #                       nmulti = 5, 
                      #                       tol = 0.1,
                      #                       ftol = 0.1,
                      #                       itmax = 20, 
                      data = data.work )
      
    }
  }
  
  if(bwmethod == "FIXED") {
    # A déterminer 
    bw.fixed.g <- c(1.522843e+07, 2.113754e+06 ,1.407951e+07, 5.513004e+06, 9.820912e+05 ,
                    1.740240e+08 , 2.020195e+08, 1.940852e-02 , 6.899923e+00, 7.041676e+00)
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = e.f, 
                    xdat = W,    
                    regtype = "ll",
                    bw <- bw.fixed.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  if(bwmethod == "SILVERMAN") {
    
    bw.silver.g <-  apply(W, 2, silverman.bw)
    # A déterminer si on modifie la fenêtre de Silverman .<<<---- ACHTUNG A MODIFIER ? #####
    bw.silver.g <- ampli*bw.silver.g 
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = e.f, 
                    xdat = W,         
                    regtype = "ll",
                    bws = bw.silver.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  
  # bandwidth information   
  print(" --- bandwith used for G  ---")
  print(" ---------------")
  print(summary(bw.g))
  tac  <- proc.time()
  duree.bw <-tac-tic
  save(bw.g, duree.bw, file = paste(Myroot,"Results/bw.g.",bwmethod,".",year,".RData", sep=""))
  
  
  ##  Np estimation of g 
  g.np <- npreg( bws = bw.g,
                 data = data.work,
                 gradients = TRUE,
                 residuals = TRUE )
  # summary(g.np)
  
  save(bw.f, bw.g, f.np, g.np,  file = paste(Myroot,"Results/np.",bwmethod,".",year,".RData", sep=""))
  
  # Step 3 -- Variance estimation 
  
  # tydat <- data.frame(e.f.2)
  # sigma2 <- npksum(txdat,
  #                  tydat, 
  #                   bws=bw.X$bw, 
  #                   bandwidth.divide=TRUE)$ksum
  
  # from Hardle 1990, page 100 and Kumbhakar and Tsionas (2009) p242 formula 42 we derive
  e2 <- e.f^2
  
  bw.sigma <- npregbw(ydat = e2, 
                      xdat = W,   
                      bandwidth.compute=FALSE, 
                      bws = bw.f$bw,
                      data = data.work)
  
  sigma2 <- npreg(ydat = e2, 
                  xdat = W,                      
                  bws = bw.sigma,
                  data = data.work)
  
  sigma2.fit <- fitted(sigma2)
  
  #### computing  theta's  ###
  
  # step 4 --  derivative of f and  g 
  # first retrive the names of the variables
  nom <- f.np$xnames
  idx.grain <- match("concvlr", nom) # <<<- correspond to  concvlr
  idx.irri <- match("hasfpirri", nom) # <<<- correspond to  hasfirri
  
  # gradient for grain 
  f1.np <- f.np$grad[,  idx.grain]  
  g1.np <- g.np$grad[,  idx.grain]  
  
  # gradient for irri
  f2.np <- f.np$grad[, idx.irri]  
  g2.np <- g.np$grad[, idx.irri] 
  

  # Step 4 -- computing theta1 
  w1 <- prixconcvlr
  p <-  prixlaiteriemoyen
  
  theta1 <-  -(f1.np - (w1/p))/g1.np
  
  #       --- computing theta2
  # ordre de grandeur pour le coût marginal de l'irrigation : environ 15 euros par hectare 
  # c'est une info moyenne pour la campagne 2010 pour la région du sud-ouest.
  
  w2 <- 15 
  theta2 <-  -(f2.np - (w2/p))/g2.np
  
  # computing theta 
  theta <- (theta1  + theta2) /2
  
  # Computing PROFIT  and sigma  (CB 27/11/2015)
  
  # Milk was divided by 1000 in the Stata file (EpurationOptilait.do)
  # concvlr is defined per cow !!! 
  
  Profit <-   p*laitproduit - w1*concvlr*eqvltotal - w2*hasfpirri 
  SigmaProf <- sd(Profit)
  AR <- - theta/SigmaProf

  # Computing RP    
  
  RP <- -0.5 * theta * SigmaProf
  RP.pc <- RP/Profit
  
  
  # Computing AR :::::::::::::! OLD VERSION ::::::::::::::::
  sigma <- sqrt(sigma2.fit)
  
  AR.old <-  - theta/sigma
  
  ####  RESULTS ############
  summary(theta)
  summary(AR)
  summary(AR.old)
  
  ###returned results ####
  
  #     risk.results <- list(sample.used =  data.work, 
  #                        theta = theta, 
  #                        AR = AR )  
  #     
  risk.results <- as.data.frame(cbind(ident, 
                                      annee,
                                      region, 
                                      f1.np,
                                      g1.np,
                                      f2.np,
                                      g2.np,
                                      theta1,
                                      theta2,
                                      theta,
                                      Profit,
                                      SigmaProf,
                                      AR,
                                      RP, 
                                      RP.pc,
                                      sigma, 
                                      AR.old)) 
  save( risk.results, file = paste(Myroot,"Results/Risknp.",bwmethod,".",year,".RData", sep=""))
  
  return(risk.results)
  
}

# Function inspired by NpAversion4 with pooled data
# 12/10/15 : test with speeding option in CV computation 
# 16/10/15 : Major change : Introducing "region" in  estimation
# 16/10/15 : Subsample option embedded 

NpPool2 <- function( Y.min, Y.max, 
                    Xvar = x,
                    N.sample = 0,
                    bwmethod = "CV",
                    bwcompute = "FALSE",
                    datafich = data.all,
                    ampli = 1,
                    ...
) 
{
  options(np.messages=FALSE)
  ################### TEST ##############
  #    year <-1996
  #      bwmethod <- "CV"
  #############################
  
  Myroot <- "D:/progs/Optilait/"
  
  # We take observations for the considered years 
  data.work <- subset(datafich,  annee >= Y.min & annee <= Y.max)
  N <- nrow(data.work)
  print(N)
  
  # we need to restrict the dataframe to complete observations
  data.work <- data.work[complete.cases(data.work),]
  
  #  Subssample of the data (for tests)
  if (N.sample > 0) { 
    # For testing  : subset of the data 
    set.seed(1234)
    S <-  N.sample   # Sample size of the training file
    ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
    data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
  }  
  
  attach(data.work)
  # NEW : All the variable are in here 
  W <- as.data.frame(data.work[,Xvar])      
  
  
  ##################   Estimation ###################  
  # New model without  "prixlaiteriemoyen" 
  ## Step 1 -  bandwidth for f 
  
  tic  <- proc.time()
  
  if(bwmethod == "CV"){
    if(bwcompute != "TRUE") {
      # 11/05/2015  we load the cross-validated bandwidths
      #load (paste("Results/CV/bw.pool.f.",Y.min,Y.max,".RData", sep=""))
      load(paste(Myroot,"Results/bw.pool.f.",bwmethod,".",Y.min,Y.max,".RData", sep=""))
      
    }
    if(bwcompute == "TRUE") {
      bw.f <- npregbw(ydat = data.work$laitproduit, 
                      xdat = W, 
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      ckertype = "epanechnikov",
                      ukertype = "liracine",
                      okertype = "liracine",
                      #  Options to change below 
                      nmulti = 5, 
                      tol = 0.1,
                      ftol = 0.1,
                      itmax = 20, 
                      data = data.work )
      
    }
  }
  
  if(bwmethod == "FIXED") { ### Will not work beacause not computed for "annee" !!
    # A déterminer 
    bw.fixed.f <- c(3.766250e+06, 4.150987e+05, 2.359682e+06, 4.393298e+06,
                    1.204398e+06, 8.197868e+07, 1.394551e+09, 1.174233e+00,
                    6.407930e+04, 1.398660e+05)
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = X,                       
                    regtype = "ll",
                    bws = bw.fixed.f,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  if(bwmethod == "SILVERMAN") {
    bw.silver.f <-  apply(W, 2, silverman.bw)  # function "Silverman.bw" applied here
    
    # A déterminer si on modifie la fenêtre de Silverman ... (ici x 2) 
    bw.silver.f <- ampli*bw.silver.f  ### <<<<<<<------- amplificateur  #####
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = W,                       
                    regtype = "ll",
                    bws = bw.silver.f,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  # bandwidth information   
  print(" --- bandwith used for F ---")
  print(" ---------------")
  print(summary(bw.f))
  
  tac  <- proc.time()
  duree.bw <-tac-tic
  
  save(bw.f, duree.bw,  file = paste(Myroot,"Results/bw.pool.f.",bwmethod,".",Y.min,Y.max,".RData", sep=""))
  
  # npplot(bw.f,  plot.errors.method="bootstrap")
  ##  Np estimation of f  (need to integrate to data.work)  
  
  f.np <- npreg( bws = bw.f,
                 gradients = TRUE,
                 residuals = TRUE)
  
  
  ## Step 2 -- Computing g
  # Residuals valeur absolue ou pas ????
  # e.f <- abs(f.np$resid)
  e.f <-f.np$resid
  
  # -- Computing bw for g  (same variables than for f )
  if(bwmethod == "CV") {
    if(bwcompute != "TRUE") {
      # 11/05/2015  we load the cross-validated bandwidths
      load (paste("Results/bw.pool.g.",bwmethod,".",Y.min,Y.max,".RData", sep=""))
    }
    if(bwcompute == "TRUE") {
      
      bw.g <- npregbw(ydat = e.f, 
                      xdat = W,  
                      regtype = "ll",
                      bwmethod = "cv.aic",
                      ckertype = "epanechnikov",
                      ukertype = "liracine",
                      okertype = "liracine",
                      #Options to change below !!!
                      nmulti = 5, 
                      tol = 0.1,
                      ftol = 0.1,
                      itmax = 20, 
                      data = data.work )
      
    }
  }
  
  if(bwmethod == "FIXED") {
    # A déterminer 
    bw.fixed.g <- c(1.522843e+07, 2.113754e+06 ,1.407951e+07, 5.513004e+06, 9.820912e+05 ,
                    1.740240e+08 , 2.020195e+08, 1.940852e-02 , 6.899923e+00, 7.041676e+00)
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = e.f, 
                    xdat = W,    
                    regtype = "ll",
                    bw <- bw.fixed.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  if(bwmethod == "SILVERMAN") {
    
    bw.silver.g <-  apply(W, 2, silverman.bw)
    # A déterminer si on modifie la fenêtre de Silverman .<<<---- ACHTUNG A MODIFIER ? #####
    bw.silver.g <- ampli*bw.silver.g 
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = e.f, 
                    xdat = W,         
                    regtype = "ll",
                    bws = bw.silver.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }
  
  
  # bandwidth information   
  print(" --- bandwith used for G  ---")
  print(" ---------------")
  print(summary(bw.g))
  tac  <- proc.time()
  duree.bw <-tac-tic
  save(bw.g, duree.bw, file = paste(Myroot,"Results/bw.pool.g.",bwmethod,".",Y.min,Y.max,".RData", sep=""))
  
  
  
  ##  Np estimation of f 
  g.np <- npreg( bws = bw.g,
                 data = data.work,
                 gradients = TRUE,
                 residuals = TRUE )
  # summary(g.np)
  
  save(bw.f, bw.g, f.np, g.np,  file = paste(Myroot,"Results/np.pool.",bwmethod,".",Y.min,Y.max,".RData", sep=""))
  
  # Step 3 -- Variance estimation 
  
  # tydat <- data.frame(e.f.2)
  # sigma2 <- npksum(txdat,
  #                  tydat, 
  #                   bws=bw.X$bw, 
  #                   bandwidth.divide=TRUE)$ksum
  
  # from Hardle 1990, page 100 and Kumbhakar and Tsionas (2009) p242 formula 42 we derive
  e2 <- e.f^2
  
  bw.sigma <- npregbw(ydat = e2, 
                      xdat = W,   
                      bandwidth.compute=FALSE, 
                      bws = bw.f$bw,
                      data = data.work)
  
  sigma2 <- npreg(ydat = e2, 
                  xdat = W,                      
                  bws = bw.sigma,
                  data = data.work)
  
  sigma2.fit <- fitted(sigma2)
  
  #### computing  theta's  ###
  
  # step 4 --  derivative of f and  g 
  
  f1.np <- f.np$grad[, 7]  # <<<- correspond to  concvlr
  g1.np <- g.np$grad[, 7]  # <<<- correspond to  concvlr
  
  f2.np <- f.np$grad[, 4]  # <<<- A MODIFIER !!!!!
  g2.np <- g.np$grad[, 4]  # <<<- A MODIFIER !!!!!
  
  # Step 4 -- computing theta1 
  w1 <- prixconcvlr
  p <-  prixlaiteriemoyen
  theta1 <-  -(f1.np - (w1/p))/g1.np
  
  #       --- computing theta2
  # ordre de grandeur pour le coût marginal de l'irrigation : environ 15 euros par hectare 
  # c'est une info moyenne pour la campagne 2010 pour la région du sud-ouest.
  
  w2 <- 15 
  theta2 <-  -(f2.np - (w2/p))/g2.np
  
  # computing theta 
  theta <- (theta1  + theta2) /2
  
  # Computing AR !! 
  sigma <- sqrt(sigma2.fit)
  
  AR <-  - theta/sigma
  
  ####  RESULTS ############
  summary(theta)
  summary(AR)
  
  ###returned results ####
  
  #     risk.results <- list(sample.used =  data.work, 
  #                        theta = theta, 
  #                        AR = AR )  
  #     
  risk.results <- as.data.frame(cbind(ident, 
                                      annee,
                                      region,
                                      theta,
                                      theta1,
                                      theta2,
                                      sigma, 
                                      AR)) 
  save(risk.results, file = paste(Myroot,"Results/Risknp.pool.",bwmethod,".",Y.min,Y.max,".RData", sep=""))
  
  return(risk.results)
  
}




## Fonction de calcul de la fenêtre de Silverman pour une variable
# From Specialreg/SpecialFunctions.R 

silverman.bw <- function(x){
  delta.k <- (3/(5*sqrt(5)))^(0.2)
  n.5 <- (length(x))^(-0.2)
  sigma.x <-min(sd(x, na.rm= TRUE), IQR(x, na.rm = TRUE, type = 7)/1.349)
  band <- 1.159* delta.k * sigma.x * n.5
  return(band)
}