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

NpAversion <- function(year, 
                       bwmethod = "CV",
                       ...
) 
{
  options(np.messages=FALSE)
  ###################TEST ##############
  # year <-1996
  #   bwmethod <- "CV"
  #############################
    
  # We take observations for the considered year 
    data.work <- subset(data.all, annee == year, 
                        select = c(annee, prixconcvlr, prixlaiteriemoyen,
                                   laitproduit, sau, sfp, eqvltotal,
                                   hasfpirri, charpot, quotalait, concvlr,
                                   Tyear,ETPyear, DELTAyear) )
    N <- nrow(data.work)
  # we need to restrict the dataframe to complete observations
  
  data.work <- data.work[complete.cases(data.work),]
  
#     ############ For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
#           set.seed(1234)
#           S <-  200   # Sample size of the training file
#           ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
#           data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
#           
#     
#     #############
    attach(data.work)
    X <- data.frame(data.work$sau, 
                    data.work$sfp,
                    data.work$eqvltotal,
                    data.work$hasfpirri, 
                    data.work$charpot,
                    data.work$quotalait,
                    data.work$concvlr,
                    data.work$Tyear,
                    data.work$ETPyear,
                    data.work$DELTAyear) 
    
    ##################   Estimation ###################  
    # New model without  "prixlaiteriemoyen" 
    ## Step 1 -  bandwidth for f 
    
    tic  <- proc.time()
    
    
    if(bwmethod == "CV"){
      bw.f <- npregbw(laitproduit~sau + sfp + eqvltotal + hasfpirri + charpot
                      + quotalait  +  concvlr
                      + Tyear + ETPyear + DELTAyear,  
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
    bw.silver.f <-  apply(X, 2, silverman.bw)
    
    bw.silver.f <- 2*bw.silver.f  ### <<<<<<<------- ACHTUNG A MODIFIER  #####
    
    # Calcul de l'objet "bandwidth"
    
    bw.f <- npregbw(ydat = data.work$laitproduit, 
                    xdat = X,                       
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

    save(bw.f, duree.bw,  file = paste("Results/bw.f.",year,".RData", sep=""))
    
    # npplot(bw.f,  plot.errors.method="bootstrap")
    ##  Np estimation of f  (need to integrate to data.work)  
   
    f.np <- npreg( bws = bw.f,
                   gradients = TRUE,
                   residuals = TRUE)
    
    # summary(f.np)
        
    # data.work$f.np.fit <- fitted(f.np)    #  peut fonctionner ici (pb de  nbre d'obs )
    # attach(data.work)
        ## Step 1-bis---  Li and Racine  significance test
    
    #   tic  <- proc.time()
    #   sigtest.model.np <- npsigtest(f.np, boot.num = 99 )
    #   tac  <- proc.time()
    #   duree.sig <-tac-tic
    #   summary( sigtest.model.np )
    #   save(bw.f, sigtest.model.np ,  file = "Results/sigtest.2006.RData")
    
    
    ## Step 2 -- Computing g
    # Residuals
    e.f <- f.np$resid
    
    # -- Computing bw for g  (same variable than for f )
  if(bwmethod == "CV"){
    bw.g <- npregbw(e.f ~ sau + sfp + eqvltotal + hasfpirri + charpot
                    + quotalait  +  concvlr
                    + Tyear + ETPyear + DELTAyear,  
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
  
  if(bwmethod == "FIXED") {
    # A déterminer 
    bw.fixed.g <- c(1.522843e+07, 2.113754e+06 ,1.407951e+07, 5.513004e+06, 9.820912e+05 ,
                    1.740240e+08 , 2.020195e+08, 1.940852e-02 , 6.899923e+00, 7.041676e+00)
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(e.f~ sau + sfp + eqvltotal + hasfpirri + charpot
                    + quotalait  +  concvlr
                    + Tyear + ETPyear + DELTAyear,  
                    regtype = "ll",
                    bw <- bw.fixed.g,
                    bandwidth.compute=FALSE,  
                    ckertype = "epanechnikov",
                    ukertype = "liracine",
                    okertype = "liracine",
                    data = data.work )
    
  }

  if(bwmethod == "SILVERMAN") {
    # A déterminer 
    bw.silver.g <-  apply(X, 2, silverman.bw)
    
    bw.silver.g <- 2*bw.silver.g  ### <<<<<<<------- ACHTUNG A MODIFIER  #####
    
    # Calcul de l'objet "bandwidth"
    
    bw.g <- npregbw(ydat = data.work$laitproduit, 
                    xdat = X,                       
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
  save(bw.g, duree.bw, file = paste("Results/bw.g.",year,".RData", sep=""))
  
 
    ##  Np estimation of f 
    g.np <- npreg( bws = bw.g,
                   data = data.work,
                   gradients = TRUE,
                   residuals = TRUE )
    # summary(g.np)
    
    save(bw.f, bw.g, f.np, g.np,  file = paste("Results/npfunctions.",year,".RData", sep=""))
    
    # Step 3 -- Variance estimation 
    
    txdat <- data.frame(sau,  sfp , eqvltotal , hasfpirri ,charpot,
                        quotalait,  concvlr,  Tyear , ETPyear,  DELTAyear)
    
    # # We use Silverman for these bandwidth
    # bw.X <- npudensbw(txdat ,
    #                 bwmethod="normal-reference", 
    #                 data= data.work) 
    # 
    
    # tydat <- data.frame(e.f.2)
    # sigma2 <- npksum(txdat,
    #                  tydat, 
    #                   bws=bw.X$bw, 
    #                   bandwidth.divide=TRUE)$ksum
    
    # from Hardle 1990, page 100 and Kumbhakar and Tsionas (2009) p242 formula 42 we derive
    sigma2 <- npreg(e.f^2~ sau + sfp + eqvltotal + hasfpirri + charpot
                    + quotalait  +  concvlr
                    + Tyear + ETPyear + DELTAyear,
                    bws = bw.f$bw,
                    data = data.work)
    s2 <- fitted(sigma2)
    
    #### computing  theta's  ###
    
    # step 4 --  derivative of f and  g 
    
    f1.np <- f.np$grad[, 7]  # <<<- correspond to  concvlr
    g1.np <- g.np$grad[, 7]  # <<<- correspond to  concvlr
    
    
    # Step 4 computing theta1
    w1 <- prixconcvlr
    p <-  prixlaiteriemoyen
    theta1 <-  f1.np - (w1/p)
    
    # computing tetha 
    theta <- theta1 # + theta 2) /2
    
    # Computing AR !! 
    
    AR <-  - theta/s2
    
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
                                        AR, 
                                        theta1))               
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


