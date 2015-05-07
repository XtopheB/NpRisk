# Estimation on Np regressions for the Bontemps-Couture paper
#  The impact of extreme climatic events on farmers' risk preferences (2015)
# 25/02/2015 


## Start 
rm(list=ls())
#setwd("C:/Chris/progs/Optilait")   
setwd("D:/progs/Optilait")

## libraries

library(np)
library(foreign)
library(AER)
library(Formula)
library(xtable)
library(texreg)
#library(VIM)
library(reporttools)
library(ggplot2)

options(np.messages=FALSE)

## DATA 
# File already merged with weather data
data.all <- read.dta("data/OptilaitMerge.dta")   # All years 

dim(data.all)

## Defining discrete variables in the dataframe 
# data.all$ident <- factor(data.all$ident)
# data.all$region <- factor(data.all$region)
# data.all$annee  <- ordered(data.all$annee)

# chosing the year
year <- 2003

data.work <- subset(data.all, annee ==year )
N <- nrow(data.work)

# ############ For testing  : subset of the data <<<<<<<<<<<<<<<<<<< to remove !!!!!
#       set.seed(1234)
#       S <-  111   # Sample size of the training file
#       ii <- sample(seq(1:nrow(data.work)),replace=FALSE) #  random values of iobs. index
#       data.work <- data.work[ii[1:S],]                # <<<<< ===== random sample from original file
#       
# 
# #############
attach(data.work)
## Just to compare, the linear model : 

lin.mod1 <- lm(formula = laitproduit~sau + sfp + eqvltotal + hasfpirri + charpot
               + quotalait  +  concvlr)

lin.mod2 <- lm(formula = laitproduit~sau + sfp + eqvltotal + hasfpirri + charpot
               + quotalait  +  concvlr
               + Tyear + ETPyear + DELTAyear)

lin.mod3 <- lm(formula = laitproduit~sau + sfp + eqvltotal + hasfpirri + charpot
                + quotalait  +  concvlr
               + Tyear + ETPyear + DELTAyear  + region)

screenreg(list(lin.mod1, lin.mod2, lin.mod3),
       stars = c(0.01, 0.05, 0.10), bold = 0.05, 
       custom.model.names = c("Model prod", "Model prod+weather", "Model prod+weather+region"),
       caption = paste("First linear regressions for year", year, ".")
)


##################   Estimation ###################  
# New model without  "prixlaiteriemoyen" 
## Step 1 -  bandwidth for f 

tic  <- proc.time()
bw.f <- npregbw(laitproduit~sau + sfp + eqvltotal + hasfpirri + charpot
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

tac  <- proc.time()
duree.bw <-tac-tic
save(bw.f, duree.bw,  file = "bw.f.RData")

summary(bw.f)
# npplot(bw.f,  plot.errors.method="bootstrap")
##  Np estimation of f  (nned to integrate to data.work)
f.np <- npreg( bws = bw.f,
                     data = data.work,
                     gradients = TRUE,
                     residuals = TRUE )

summary(f.np)
data.work$f.np.fit <- fitted(f.np)

# we neefd to restrict the datatframe to comuted values 
data.work <- data.work[complete.cases(data.work),]
attach(data.work)


## Step 1-bis---  Li and Racine  significance test

#   tic  <- proc.time()
#   sigtest.model.np <- npsigtest(f.np, boot.num = 99 )
#   tac  <- proc.time()
#   duree.sig <-tac-tic
#   summary( sigtest.model.np )
#   save(bw.f, sigtest.model.np ,  file = "results.2006.RData")


## Step 2 -- Computing g
# Residuals
e.f <- f.np$resid

# -- Computing bw for g  (same variable than for f )
bw.g <- npregbw(e.f~ sau + sfp + eqvltotal + hasfpirri + charpot
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

summary(bw.g)
##  Np estimation of f 
g.np <- npreg( bws = bw.g,
               data = data.work,
               gradients = TRUE,
               residuals = TRUE )
summary(g.np)

save(bw.f, bw.g, f.np, g.np,  file = "npfunctions.2003.RData")

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






#########################################

## Testing the function !!! 

## Start 
rm(list=ls())
#setwd("C:/Chris/progs/Optilait")   
setwd("D:/progs/Optilait")

## libraries

library(np)
library(foreign)
library(AER)
library(Formula)
library(xtable)
library(texreg)
#library(VIM)
library(reporttools)

options(np.messages=FALSE)

# loading the functions
source("Progs/NpRiskFunctions.R")

toto <- NpAversion(1996,  bwmethod = "SILVERMAN")

## DATA 
# File already merged with weather data
data.all <- read.dta("data/OptilaitMerge.dta")   # All years 

## Defining discrete variables in the dataframe 
# data.all$ident <- factor(data.all$ident)
# data.all$region <- factor(data.all$region)
# data.all$annee  <- ordered(data.all$annee)


dim(data.all)

# creating an empty dataframe 
risk.all<- read.csv(text="col1,col2")

yearlist<- c("1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006")

#yearlist<- c( "1999", "2001",  "2003",  "2005")

for (y in 1:length(yearlist)) {
  print(yearlist[y])
  risk.foo <- NpAversion(yearlist[y],  bwmethod = "SILVERMAN")
  risk.year  <- merge(data.all, risk.foo, by = c("ident","annee"))
  risk.all<- rbind.data.frame(risk.all, risk.year)
}


## Summary statistics 

r <- as.data.frame(risk.all)
Myvars<- with(risk.all,data.frame( "AR" = r$AR,
                            "Theta" = r$theta,
                            "Theta" = r$theta) 
            )

tableContinuous(vars = Myvars , group = r$annee,  stats = c("n",  "mean", "median","max"),
                cap = paste("Efficiencies on points."),  
                prec=3, longtable = FALSE
)


# Marking year 2001 with different color 
risk.all$color.2001 <- as.factor((risk.all$annee != 2001))

# First graph 
p <- ggplot(risk.all, aes(factor(annee), AR)) 
# Colors + no draw of outliers 
p <- p + geom_boxplot(aes(color=color.2001),outlier.colour = NA) 

# Compute lower and upper whisker limits
sts <- boxplot.stats(risk.all$AR)$stats  
p+ theme_bw() + coord_cartesian(ylim = c(min(sts)*0.95,max(sts)*1.05))


# NP Test of distribution change over time  !!



### Saving results  ##

save(risk.all,  file = paste("Results/Risk.all.SILVERx2"))


###  Outputs....
library(tables)
risk.all$years <- as.factor(risk.all$annee)
latex(
tabular( (years+1) ~ (n=1) + Format(digits=2)*(AR + theta1)*(median + sd), data=risk.all ), 
file = "Graphics/RiskSiverx2.tex"
)

X <- tabular( (years+1) ~ (n=1) + Format(digits=2)*(AR + theta1)*(median + sd), data=risk.all )
X
write.table.tabular(X,file = "Graphics/RiskSiverx2.tex")




