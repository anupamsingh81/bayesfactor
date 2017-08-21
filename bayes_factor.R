rm()

a = c(5)
b = c(10)
c = c(20)
d = c(12)

x =  data.frame(
a = c(20:21),
b = c(10:11),
c = c(25:26),
d = c(50:51)
)


prop.test(x0)

BFP = function(a,b,c,d){
  
  H = as.matrix(cbind(a,c),
  I = prop.test(H)
  
  k = I$p.value
  
  k
}

lapply(x$a,x$b,x$c,x$d,FUN = BFP)



# Dummy data set 
mydata=read.table(textConnection(" Marker Treatment Genotype1 Genotype2 Genotype3 
1        A         23        57        32 
1        B         43        59        12 
2        A         13        27        12 
2        B         23        29        22"),header=TRUE) 


# Chi-square 
res=lapply(split(mydata,mydata$Marker),function(x){ 
  res=chisq.test(x[,-c(1,2)]) # Deleting columns 1 and 2 of every list 
  res=c(res$statistic,res$p.value) 
  names(res)=c('statistic','pvalue') 
  res 
} 
) 

# The output 
res2=do.call(rbind,res) 
rownames(res2)=paste('marker_',1:length(res),sep="") 
res2 
  
# Bayes Factor
library(BayesFactor)

res3 =lapply(split(mydata,mydata$Marker),function(x){
  a = BayesFactor::contingencyTableBF(as.matrix(x[,-c(1,2)]),sampleType = "indepMulti", fixedMargin = "cols")  
  bf =  c(exp(a@bayesFactor$bf))
  bf
  } 
) 


res3$`2`@bayesFactor$code



t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  BF = BayesFactor::meta.ttestBF(t,n1,n2)
  bf =  c(exp(BF@bayesFactor$bf))
 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
 
}

m1 = c(10,20)
m2 = c(30,40)
s1 = c(10,12)
s2 = c(12,14)
n1 = c(16,26)
n2 = c(18,28) 

BT = t.test2(m1,m2,s1,s2,n1,n2)

X1 = data.frame(
t = c(2,3),
n1 = c(30,40),
n2 = c(40,60)
)

M = mapply(BayesFactor::meta.ttestBF,X1$t,X1$n1,X1$n2) # http://stackoverflow.com/questions/6827299/r-apply-function-with-multiple-parameters
W <- lapply( M , function(x) `@`( x , bayesFactor))[[1]]  # isolating @ from list of list, http://stackoverflow.com/questions/17971311/r-extracting-members-of-lists-and-sub-lists-s4-class
W1 = lapply(W,function(x)exp(x)) # isolating BF
X1$bf = unlist(W1) # unlisting and supplying bf as part of data frame

z1 =rnorm(20,60,10)
z2 = rnorm(40,40,20)
t.test(z1,z2)
ttestBF(z1,z2)
meta.ttestBF(4.08,20,40)

library(fmsb)
o = oddsratio(23,25,42,30)
selog =( log(o$conf.int[2]) - log(o$conf.int[1]))/3.92
logOR = log(o$estimate)
t = (log(1) - logOR)/selog

# imagining x i.e odds ratio as 1 , hence tstatistic is log1 - logOR/s.e

exp(ttest.tstat(t=1.119, n1=65, n2=55, rscale = 0.707)[['bf']]) # same calculation

bfdichotomous = function(a,b,c,d,r){
  o = oddsratio(a,b,c,d)
  p = o$p.value
  selog =( log(o$conf.int[2]) - log(o$conf.int[1]))/3.92
  logOR = log(o$estimate)
  t = (log(1) - logOR)/selog
  
  
 bfo = exp(ttest.tstat(t=t, n1= a+c, n2=b+d, rscale = r)[['bf']])
 
 g = c(o$estimate,o$conf.int,o$p.value,bfo)
 
  return(g)
}

bfdichotomous(32,33,216,225,0.7)

0.29/0.21
bfdichotomous(20,12,196,202,0.7)
20/216
12/214

oddsratio(4.4,12.7,95.6,87.3)
log(0.316)


oddsratio(32,33,216,225)

obtained =log(1.01)

sd = (log(1.7)-log(0.6))/3.92
exp(0.85)

Bfnormal(sd=0.26,obtained =0.0099,meanoftheory = 0,sdtheory = - 0.16 )
Bfchristie(sd=0.26,obtained =0.0099,meanoftheory = 0,sdtheory =  0.5 )

Bfnormal(sd=5,obtained =12,meanoftheory = 0,sdtheory = 11 ,tail = 1)
Bfchristie(sd=5,obtained =12,meanoftheory = 0,sdtheory = 11 ,tail = 1)



log(0.85)


s.e = sqrt(1/23 + 1/25 +1/42 + 1/30)
logOR

t = abs(logOR/s.e)

t1 = as.table(cbind(c(32,216),c(33,225)))
t1
tbf = contingencyTableBF(t1, sampleType = "indepMulti", fixedMargin = "cols")
tbf
1/0.0744
#Dienes calculator

#Code for Bt with normal likelihood (adapting the R code provided by Baguley and Kaye, 2010, for the Dienes 2008 calculator)

*************************************************************************
  # sd is standard error, also can calculate by se= t/meandifference;also sqrt(1/a +1/b+1/c+1/d), obtained=effect size, mean difference or log(OR),mean of theory is default set to zero, sd theory = log(OR) or mean difference of pilot study or study used to calculate power.)
  # Bfnormal with 1 tail is half normal.
  Bfnormal<-function(sd, obtained, meanoftheory=0, sdtheory,  dftheory = 1, tail=1)
  {
    area <- 0
    normarea <- 0
    theta <- meanoftheory - 10 * sdtheory
    incr <- sdtheory / 200
    for (A in -2000:2000){
      theta <- theta + incr
      tscore = (theta - meanoftheory)/sdtheory
      dist_theta <- dt(tscore, df=dftheory)
      
      if(identical(tail, 1)){
        if (theta <= 0){
          dist_theta <- 0
        } else {
          dist_theta <- dist_theta * 2
        }
      }
      height <- dist_theta * dnorm(obtained, theta, sd)
      area <- area + height * incr
      normarea <- normarea + dist_theta*incr
    }
    LikelihoodTheory <- area/normarea
    Likelihoodnull <- dnorm(obtained, 0, sd)
    BayesFactor <- LikelihoodTheory / Likelihoodnull 
    BayesFactor
  }

#Bft is t distribution with1 degree of freedom also called cauchy, if greater than 30 it is normal
Bft<-function(sd , obtained, dfdata , meanoftheory=0, sdtheory=1, dftheory = 1,tail=2)
{           
  area <- 0
  normarea <- 0
  
  theta <- meanoftheory - 10 * sdtheory
  incr <- sdtheory / 200
  for (A in -2000:2000){
    theta <- theta + incr
    dist_theta <- dnorm(theta, meanoftheory, sdtheory)
    dist_theta <- dt((theta-meanoftheory)/sdtheory, df=dftheory)
    if(identical(tail, 1)){
      if (theta <= 0){
        dist_theta <- 0
      } else {
        dist_theta <- dist_theta * 2
      }
    }
    height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta*incr
  }
  LikelihoodTheory <- area/normarea
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  BayesFactor <- LikelihoodTheory / Likelihoodnull 
  BayesFactor
}

# Usage
#http://www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/Bayes%20factor%20with%20t%20distribution.html

Bfnormal(sd = 1, obtained = 2,meanoftheory=0, sdtheory = 3) # half normal as default tail is zero
Bft(sd = 1, obtained = 2,dfdata=3,meanoftheory=0, sdtheory = 3)


Bfchristie(sd = 1, obtained = 2,meanoftheory=0, sdtheory = 3) # half normal, said to be more acurate as calculate by area under the curve. Uniform is set to false.
#Bfchristie
Bfchristie(sd = 1, obtained = 2,uniform=TRUE,lower=0 ,upper=10,meanoftheory=0) # uniform
Bfchristie<-function(sd, obtained, uniform = FALSE, lower=0, upper=1, meanoftheory=0, sdtheory=1, tails=2)
{
  #Version 2.0.1
  #modifications by John Christie
  # modification means that this does not exactly replicate Dienes.  What it does
  # is do the right thing instead.  :)  The current version is more accurate.
  #03/10/2011
  #Authors Danny Kaye & Thom Baguley
  #Version 1.0
  #19/10/2009
  # test data can be found starting at p100
  # notes on how to use these functions can be found at
  # http://danny-kaye.co.uk/Docs/Dienes_notes.pdf
  #raised from 2000 for better accuracy - the speed of the new code allows it
  slices <- 20000
  if(uniform){
    range <- upper - lower
    dist_theta <- 1 / range
    # incr <- range / slices
    # theta <- seq(lower + incr, by = incr, length.out = slices+1)
    # height <- dist_theta * dnorm(obtained, theta, sd)
    # area <- sum(height * incr)
    # the commented code above replicates the original result but at the
    # limit (slices <- 5e6) it's actually equivalent to the followingâ???¦
    area <- dist_theta * diff(pnorm(c(lower, upper), obtained, sd)) 
  }else{
    # the code below again doesnt' replicate the original code.
    # incrementing in a scalar loop was causing an accumulation of tiny fp errors
    # the lower end was incremented prior to starting (that's fixed above too)
    zlim <- 5
    incr <- sdtheory / (slices/(zlim*2))
    newLower <- meanoftheory - zlim * sdtheory
    theta <- seq(newLower, by = incr, length.out = slices+1)
    dist_theta <- dnorm(theta, meanoftheory, sdtheory)
    if (tails == 1){
      dist_theta <- dist_theta[theta > 0]	* 2
      theta <- theta[theta > 0]	
    }
    height <- dist_theta * dnorm(obtained, theta, sd)
    area <- sum(height * incr)
  }
  LikelihoodTheory <- area
  Likelihoodnull <- dnorm(obtained, 0, sd)
  BayesFactor <- LikelihoodTheory / Likelihoodnull
  return( list("LikelihoodTheory" = LikelihoodTheory, "Likelihoodnull" = Likelihoodnull, "BayesFactor" = BayesFactor) )
}
# Dienes paper
# case 2.5 The answer to the question should depend on the question
# reported in paper as BH(0,1.11) symbolic of mean of theory 0, sd theory 1.11 , reported 0.38
Bfnormal(sd = 0.25, obtained = 0.15,meanoftheory=0, sdtheory = 1.11)
Bfchristie(sd = 0.25, obtained = 0.15,meanoftheory=0, sdtheory = 1.11)

Bft(sd = 0.24, obtained = 0.15,dfdata=124,meanoftheory=1.11, sdtheory = 1.11)
Bfchristie(sd = 0.25, obtained = 0.15,uniform=TRUE,lower=0,upper =6) # reported in paper 0.09

# case 2.4 4 A high-powered significant result is not necessarily evidence for a theory

Bfnormal(sd = 0.14, obtained = -0.26,meanoftheory=0, sdtheory = 1.26) # effect in opposite direction

Bfnormal(sd = 0.14, obtained = 0.26,meanoftheory=0, sdtheory = 1.26) # quoted 1.18 ,effect in same direction

Bfchristie(sd = 0.14, obtained = 0.26,meanoftheory=0, sdtheory = 1.26)

Bfnormal(sd = 0.14, obtained = 0.28,meanoftheory=0, sdtheory = 1.26) # quoted 1.56

Bft(sd = 0.14, obtained = 0.26,dfdata=1,meanoftheory=0, sdtheory = 1.26) # trying t distribution with wider range 

Bfnormal(sd = 0.14, obtained = 0.28,meanoftheory=0, sdtheory = 0.40) #more precise theory, quoted 3.81

# case 2.3 A low-powered non-significant result is not necessarily insensitive

Bfnormal(sd = 4/1.15, obtained = -4,meanoftheory=0, sdtheory = 5) # quoted 0.31, effect in wrong direction hence favors null

Bfnormal(sd = 8/1.15, obtained = -4,meanoftheory=0, sdtheory = 5) # quoted 0.63, effcet in wrong direction , but twice as large standard error, hence evidence for null is poorer

# 2.2 A high powered non-significant result is not necessarily sensitive


Bfnormal(sd = 5.5/0.87, obtained = 5.5,meanoftheory=0, sdtheory = 13.3) #quoted 0.97
Bfchristie(sd = 5.5/0.87, obtained = 5.5,meanoftheory=0, sdtheory = 13.3) # christie in generalgives evidence more in favour of null with same parameters

#2.1 Often significance testing will provide adequate answers


Bfnormal(sd = 12/2.4, obtained = 12,meanoftheory=0, sdtheory = 11) # quoted 4.5


Bfchristie(sd = 12/2.4, obtained = 12,meanoftheory=0, sdtheory = 11) 

Bfnormal(sd = 0.14, obtained = -0.26,meanoftheory=0, sdtheory = 1.26) # quoted 0.04

Bfchristie(sd = 0.14, obtained = -0.26,meanoftheory=0, sdtheory = 1.26) # christie almost differs by a factor! strange!
# online calculator yields exactly 0.04

Bft(sd = 1, obtained = 2,dfdata=3,meanoftheory=0, sdtheory = 3) # mentioned as Bt(mean,se,df) of original study



