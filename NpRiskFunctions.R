## Functons used to compute Risk aversion for each climatic year   
# 04/05/2015 : Version 0 recovering elements from NpModel.R


## Fonction de calcul de la fenêtre de Silverman pour une variable
# From Specialreg/SpecialFunctions.R 

silverman.bw <- function(x){
  delta.k <- (3/(5*sqrt(5)))^(0.2)
  n.5 <- (length(x))^(-0.2)
  sigma.x <-min(sd(x, na.rm= TRUE), IQR(x, na.rm = TRUE, type = 7)/1.349)
  band <- 1.159* delta.k * sigma.x * n.5
  return(band)
}

## Fonction de calcul ds l'AR , theta1, etc...
## 10/05/2015  : New version  : f(x, Z) , g(x)  x= factors,  z = weather factors...
##             :  Correction computing theta (denominator)
## 15/05/2015  : Adapted to parallel tratment
##             : Correction : computing g with abs(ei) on same regressors. 

## 03/06/2015  :Version 4 with computation of bandwidth alos for g

NpAversion4 <- function(year, 
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

  # We take observations for the considered year 
  data.work <- subset(datafich,  annee == year, 
                      select = c(annee, ident, prixconcvlr, prixlaiteriemoyen,
                                 laitproduit, sau, sfp, eqvltotal,
                                 hasfpirri, charpot, quotalait, concvlr,
                                 Tyear,ETPyear, DELTAyear) )
  N <- nrow(data.work)
  print(N)
  # we need to restrict the dataframe to complete observations
  
  data.work <- data.work[complete.cases(data.work),]
  
#       ############ For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
#             set.seed(1234)
#             S <-  100   # Sample size of the training file
#             ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
#             data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
#             
#       
# #       #############
  attach(data.work)
  X <- data.frame(data.work$sau, 
                  data.work$sfp,
                  data.work$eqvltotal,
                  data.work$hasfpirri, 
                  data.work$charpot,
                  data.work$quotalait,
                  data.work$concvlr) 
  Z <-  data.frame(data.work$Tyear,
                   data.work$ETPyear,
                   data.work$DELTAyear)
  W <- cbind(X,Z)
  
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
      # Residuals valeur absolue ou pas ????
      # e.f <- abs(f.np$resid)
  e.f <-f.np$resid

  
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
  
  
  ##  Np estimation of f 
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
                                      theta,
                                      theta1,
                                      theta2,
                                      sigma, 
                                      AR)) 
  save( risk.results, file = paste(Myroot,"Results/Risknp.",bwmethod,".",year,".RData", sep=""))

  return(risk.results)

}
# Function inspired by NpAversion4 with pooled data
# 12/10/15 : test with speeding option in CV computation 
# 16/10/15 : Major change : Introducing "region" in  estimation
# 16/10/15 : Subsample option embedded 

NpPool <- function( Y.min, Y.max, N.sample = 0,
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
  
  # We take observations for the considered year 
#   data.work <- subset(datafich,  annee >= Y.min & annee <= Y.max, 
#                       select = c(annee, ident, prixconcvlr, prixlaiteriemoyen,
#                                  laitproduit, sau, sfp, eqvltotal,
#                                  hasfpirri, charpot, quotalait, concvlr,region,
#                                  Tyear,ETPyear, DELTAyear) )
  data.work <- subset(datafich,  annee >= Y.min & annee <= Y.max, 
                      select = c(annee, ident, prixconcvlr, prixlaiteriemoyen,
                                 laitproduit, sau, sfp, eqvltotal,
                                 hasfpirri, charpot, quotalait, concvlr,region,
                                 Tprod,ETPprod, DELTAprod) )
  N <- nrow(data.work)
  print(N)
  # we need to restrict the dataframe to complete observations
  
  data.work <- data.work[complete.cases(data.work),]
  if (N.sample > 0) { 
# For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
      set.seed(1234)
      S <-  N.sample   # Sample size of the training file
      ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
      data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
              
  }  

  attach(data.work)
  X <- data.frame(data.work$sau, 
                  data.work$sfp,
                  data.work$eqvltotal,
                  data.work$hasfpirri, 
                  data.work$charpot,
                  data.work$quotalait,
                  data.work$concvlr,
                  ordered(data.work$annee)
                  ) 
  Z <-  data.frame(data.work$Tprod,
                   data.work$ETPprod,
                   data.work$DELTAprod)
#   Z <-  data.frame(data.work$Tyear,
#                    data.work$ETPyear,
#                    data.work$DELTAyear)
  U <- data.frame(data.work$region)
  W <- cbind(X,Z,U)
  
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




NpAversion3 <- function(year, 
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
  
  # We take observations for the considered year 
  data.work <- subset(datafich,  annee == year, 
                      select = c(annee, ident, prixconcvlr, prixlaiteriemoyen,
                                 laitproduit, sau, sfp, eqvltotal,
                                 hasfpirri, charpot, quotalait, concvlr,
                                 Tyear,ETPyear, DELTAyear) )
  N <- nrow(data.work)
  print(N)
  # we need to restrict the dataframe to complete observations
  
  data.work <- data.work[complete.cases(data.work),]
  
  #       ############ For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
  #             set.seed(1234)
  #             S <-  100   # Sample size of the training file
  #             ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
  #             data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
  #             
  #       
  # #       #############
  attach(data.work)
  X <- data.frame(data.work$sau, 
                  data.work$sfp,
                  data.work$eqvltotal,
                  data.work$hasfpirri, 
                  data.work$charpot,
                  data.work$quotalait,
                  data.work$concvlr) 
  Z <-  data.frame(data.work$Tyear,
                   data.work$ETPyear,
                   data.work$DELTAyear)
  W <- cbind(X,Z)
  
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
  
  # summary(f.np)
  
  ## Step 1-bis---  Li and Racine  significance test
  
  #   tic  <- proc.time()
  #   sigtest.model.np <- npsigtest(f.np, boot.num = 99 )
  #   tac  <- proc.time()
  #   duree.sig <-tac-tic
  #   summary( sigtest.model.np )
  #   save(bw.f, sigtest.model.np ,  file = "Results/sigtest.2006.RData")
  
  
  ## Step 2 -- Computing g
  # Residuals valeur absolue ou pas ????
  # e.f <- abs(f.np$resid)
  e.f <-f.np$resid
  
  
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
  save(bw.g, duree.bw, file = paste(Myroot,"Results/bw.g.",bwmethod,".",year,".RData", sep=""))
  
  
  ##  Np estimation of f 
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
                                      theta,
                                      theta1,
                                      theta2,
                                      sigma, 
                                      AR)) 
  save( risk.results, file = paste(Myroot,"Results/Risknp.",bwmethod,".",year,".RData", sep=""))
  
  return(risk.results)
  
}

##  Li, Maasoumi, and Racine (2009) test ###
TestDistLi <- function(eff1, eff2, nboot) {
  one <- data.frame(don=eff1)
  two <- data.frame(don=eff2)
  
  names(one) <- c("effF")  # variables must have the same name
  names(two) <- c("effF")  # variables must have the same name
  Test <- npdeneqtest(one, two, boot.num = nboot)
  list(Pval = Test$Tn.P , Val = Test$Tn)
}




###  Saved for future, but  f and ga should not have the same  regressors !!!

# NpAversion <- function(year, 
#                        bwmethod = "CV",
#                        bwcompute = "FALSE",
#                        ...
# ) 
# {
#   options(np.messages=FALSE)
#   ################### TEST ##############
#   # year <-1996
#   #   bwmethod <- "CV"
#   #############################
#   
#   # We take observations for the considered year 
#   data.work <- subset(data.all, annee == year, 
#                       select = c(annee, ident, prixconcvlr, prixlaiteriemoyen,
#                                  laitproduit, sau, sfp, eqvltotal,
#                                  hasfpirri, charpot, quotalait, concvlr,
#                                  Tyear,ETPyear, DELTAyear) )
#   N <- nrow(data.work)
#   # we need to restrict the dataframe to complete observations
#   
#   data.work <- data.work[complete.cases(data.work),]
#   
#   #     ############ For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
#   #           set.seed(1234)
#   #           S <-  100   # Sample size of the training file
#   #           ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
#   #           data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
#   #           
#   #     
#   #     #############
#   attach(data.work)
#   X <- data.frame(data.work$sau, 
#                   data.work$sfp,
#                   data.work$eqvltotal,
#                   data.work$hasfpirri, 
#                   data.work$charpot,
#                   data.work$quotalait,
#                   data.work$concvlr,
#                   data.work$Tyear,
#                   data.work$ETPyear,
#                   data.work$DELTAyear) 
#   
#   ##################   Estimation ###################  
#   # New model without  "prixlaiteriemoyen" 
#   ## Step 1 -  bandwidth for f 
#   
#   tic  <- proc.time()
#   
#   
#   if(bwmethod == "CV"){
#     if(bwcompute != "TRUE") {
#       # 11/05/2015  we load the cross-validated bandwidths
#       load (paste("Results/CV/bw.f.",year,".RData", sep=""))
#     }
#     if(bwcompute == "TRUE") {
#       bw.f <- npregbw(laitproduit~sau + sfp + eqvltotal + hasfpirri + charpot
#                       + quotalait  +  concvlr
#                       + Tyear + ETPyear + DELTAyear,  
#                       regtype = "ll",
#                       bwmethod = "cv.aic",
#                       ckertype = "epanechnikov",
#                       ukertype = "liracine",
#                       okertype = "liracine",
#                       #  Options to change below 
#                       nmulti = 5, 
#                       tol = 0.1,
#                       ftol = 0.1,
#                       itmax = 20, 
#                       data = data.work )
#     }
#   }
#   
#   if(bwmethod == "FIXED") {
#     # A déterminer 
#     bw.fixed.f <- c(3.766250e+06, 4.150987e+05, 2.359682e+06, 4.393298e+06,
#                     1.204398e+06, 8.197868e+07, 1.394551e+09, 1.174233e+00,
#                     6.407930e+04, 1.398660e+05)
#     
#     # Calcul de l'objet "bandwidth"
#     
#     bw.f <- npregbw(ydat = data.work$laitproduit, 
#                     xdat = X,                       
#                     regtype = "ll",
#                     bws = bw.fixed.f,
#                     bandwidth.compute=FALSE,  
#                     ckertype = "epanechnikov",
#                     ukertype = "liracine",
#                     okertype = "liracine",
#                     data = data.work )
#     
#   }
#   
#   if(bwmethod == "SILVERMAN") {
#     bw.silver.f <-  apply(X, 2, silverman.bw)
#     
#     # A déterminer si on modifie la fenêtre de Silverman ... (ici x 2) 
#     # bw.silver.f <- 2*bw.silver.f  ### <<<<<<<------- ACHTUNG A MODIFIER  #####
#     
#     # Calcul de l'objet "bandwidth"
#     
#     bw.f <- npregbw(ydat = data.work$laitproduit, 
#                     xdat = X,                       
#                     regtype = "ll",
#                     bws = bw.silver.f,
#                     bandwidth.compute=FALSE,  
#                     ckertype = "epanechnikov",
#                     ukertype = "liracine",
#                     okertype = "liracine",
#                     data = data.work )
#     
#   }
#   
#   # bandwidth information   
#   print(" --- bandwith used for F ---")
#   print(" ---------------")
#   print(summary(bw.f))
#   
#   tac  <- proc.time()
#   duree.bw <-tac-tic
#   
#   save(bw.f, duree.bw,  file = paste("Results/bw.f.",year,".RData", sep=""))
#   
#   # npplot(bw.f,  plot.errors.method="bootstrap")
#   ##  Np estimation of f  (need to integrate to data.work)  
#   
#   f.np <- npreg( bws = bw.f,
#                  gradients = TRUE,
#                  residuals = TRUE)
#   
#   # summary(f.np)
#   
#   # data.work$f.np.fit <- fitted(f.np)    #  peut fonctionner ici (pb de  nbre d'obs )
#   # attach(data.work)
#   ## Step 1-bis---  Li and Racine  significance test
#   
#   #   tic  <- proc.time()
#   #   sigtest.model.np <- npsigtest(f.np, boot.num = 99 )
#   #   tac  <- proc.time()
#   #   duree.sig <-tac-tic
#   #   summary( sigtest.model.np )
#   #   save(bw.f, sigtest.model.np ,  file = "Results/sigtest.2006.RData")
#   
#   
#   ## Step 2 -- Computing g
#   # Residuals
#   e.f <- f.np$resid
#   
#   # -- Computing bw for g  (same variable than for f )
#   if(bwmethod == "CV") {
#     if(bwcompute != "TRUE") {
#       # 11/05/2015  we load the cross-validated bandwidths
#       load (paste("Results/CV/bw.g.",year,".RData", sep=""))
#     }
#     if(bwcompute == "TRUE") {
#       bw.g <- npregbw(e.f ~ sau + sfp + eqvltotal + hasfpirri + charpot
#                       + quotalait  +  concvlr
#                       + Tyear + ETPyear + DELTAyear,  
#                       regtype = "ll",
#                       bwmethod = "cv.aic",
#                       ckertype = "epanechnikov",
#                       ukertype = "liracine",
#                       okertype = "liracine",
#                       nmulti = 5, 
#                       tol = 0.1,
#                       ftol = 0.1,
#                       itmax = 20, 
#                       data = data.work )
#     }
#   }
#   
#   if(bwmethod == "FIXED") {
#     # A déterminer 
#     bw.fixed.g <- c(1.522843e+07, 2.113754e+06 ,1.407951e+07, 5.513004e+06, 9.820912e+05 ,
#                     1.740240e+08 , 2.020195e+08, 1.940852e-02 , 6.899923e+00, 7.041676e+00)
#     
#     # Calcul de l'objet "bandwidth"
#     
#     bw.g <- npregbw(e.f~ sau + sfp + eqvltotal + hasfpirri + charpot
#                     + quotalait  +  concvlr
#                     + Tyear + ETPyear + DELTAyear,  
#                     regtype = "ll",
#                     bw <- bw.fixed.g,
#                     bandwidth.compute=FALSE,  
#                     ckertype = "epanechnikov",
#                     ukertype = "liracine",
#                     okertype = "liracine",
#                     data = data.work )
#     
#   }
#   
#   if(bwmethod == "SILVERMAN") {
#     
#     bw.silver.g <-  apply(X, 2, silverman.bw)
#     # A déterminer si on modifie la fenêtre de Silverman ... (ici x 2) 
#     bw.silver.g <- 4*bw.silver.g  ### <<<<<<<------- ACHTUNG A MODIFIER  #####
#     
#     # Calcul de l'objet "bandwidth"
#     
#     bw.g <- npregbw(ydat = data.work$laitproduit, 
#                     xdat = X,                       
#                     regtype = "ll",
#                     bws = bw.silver.g,
#                     bandwidth.compute=FALSE,  
#                     ckertype = "epanechnikov",
#                     ukertype = "liracine",
#                     okertype = "liracine",
#                     data = data.work )
#     
#   }
#   
#   
#   # bandwidth information   
#   print(" --- bandwith used for G  ---")
#   print(" ---------------")
#   print(summary(bw.g))
#   tac  <- proc.time()
#   duree.bw <-tac-tic
#   save(bw.g, duree.bw, file = paste("Results/bw.g.",year,".RData", sep=""))
#   
#   
#   ##  Np estimation of f 
#   g.np <- npreg( bws = bw.g,
#                  data = data.work,
#                  gradients = TRUE,
#                  residuals = TRUE )
#   # summary(g.np)
#   
#   save(bw.f, bw.g, f.np, g.np,  file = paste("Results/np.",bwmethod,"Y",year,".RData", sep=""))
#   
#   # Step 3 -- Variance estimation 
#   
#   txdat <- data.frame(sau,  sfp , eqvltotal , hasfpirri ,charpot,
#                       quotalait,  concvlr,  Tyear , ETPyear,  DELTAyear)
#   
#   # # We use Silverman for these bandwidth
#   # bw.X <- npudensbw(txdat ,
#   #                 bwmethod="normal-reference", 
#   #                 data= data.work) 
#   # 
#   
#   # tydat <- data.frame(e.f.2)
#   # sigma2 <- npksum(txdat,
#   #                  tydat, 
#   #                   bws=bw.X$bw, 
#   #                   bandwidth.divide=TRUE)$ksum
#   
#   # from Hardle 1990, page 100 and Kumbhakar and Tsionas (2009) p242 formula 42 we derive
#   sigma2 <- npreg(e.f^2~ sau + sfp + eqvltotal + hasfpirri + charpot
#                   + quotalait  +  concvlr
#                   + Tyear + ETPyear + DELTAyear,
#                   bws = bw.f$bw,
#                   data = data.work)
#   s2 <- fitted(sigma2)
#   
#   #### computing  theta's  ###
#   
#   # step 4 --  derivative of f and  g 
#   
#   f1.np <- f.np$grad[, 7]  # <<<- correspond to  concvlr
#   g1.np <- g.np$grad[, 7]  # <<<- correspond to  concvlr
#   
#   f2.np <- f.np$grad[, 7]  # <<<- A MODIFIER !!!!!
#   g2.np <- g.np$grad[, 7]  # <<<- A MODIFIER !!!!!
#   
#   # Step 4 -- computing theta1 
#   w1 <- prixconcvlr
#   p <-  prixlaiteriemoyen
#   theta1 <-  -(f1.np - (w1/p))/g1.np
#   
#   #       --- computing theta2
#   # ordre de grandeur pour le coût marginal de l'irrigation : environ 15 euros par hectare 
#   # c'est une info moyenne pour la campagne 2010 pour la région du sud-ouest.
#   
#   w2 <- 15 
#   theta2 <-  -(f2.np - (w2/p))/g2.np
#   
#   # computing tetha 
#   theta <- (theta1  + theta2) /2
#   
#   # Computing AR !! 
#   
#   AR <-  - theta/s2
#   
#   ####  RESULTS ############
#   summary(theta)
#   summary(AR)
#   
#   ###returned results ####
#   
#   #     risk.results <- list(sample.used =  data.work, 
#   #                        theta = theta, 
#   #                        AR = AR )  
#   #     
#   risk.results <- as.data.frame(cbind(ident, 
#                                       annee,
#                                       theta,
#                                       theta1,
#                                       theta2,
#                                       s2, 
#                                       AR))               
#   return(risk.results)
#   
# }
