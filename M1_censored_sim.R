# rm(list=ls())

library(R2WinBUGS)
library(arm)
library(lme4)
library(MCMCpack)


####### Generate data #############
# number of lakes
n = 50
# intercept
a <- -0.24
# slope of the TP-CHL relationship
b <- 0.83  
# residual SD
sigma <- 0.40
# centered values of TP, covariate
tp <- seq(-0.9, 1.7, length=n) 

# linear predictor (deterministic part of the model)
mui <- a + b * tp

# Generate CHL values (stochastic, random part of the model)
# set.seed(1224) # so we all get the same random numbers
yi <- exp(rnorm(n, mean=mui, sd=sigma))

# Fit model using lm function
summary(lm(log(yi) ~ tp))


######## Censor data
# log-transform yi
log.yi <- log(yi)
# log-transform yi: "C" used for censored response
log.yiC <- log(yi)


# Proportion to censor
cenProp <- 0.3

### Censor data (data < 30th percentile is censored)
detect <- ifelse(log.yiC < quantile(log.yiC, cenProp),0,1)
# Reporting limit is same for all observations
rl <- as.numeric(rep(quantile(log.yiC,cenProp ),n))

# Create data frame
dat <- data.frame(yi,log.yi,log.yiC,tp,detect,rl)

# If an observation was not detected then set to NA
dat$log.yiC[dat$detect==0] <- NA

# Set rl to exp(1000) for detect == 1
dat$rl[dat$detect==1] <- 1000

head(dat)

####### Plot censorted data #################

symbol <- rep(16,length(dat$log.yi))
symbol[dat$log.yi < quantile(dat$log.yi, 0.3)] <- 21

res <- 6
name_figure <- "CensoredPlot.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1,1,1,1), oma=c(3,3,0.1,0.1) )

size.labels = 1
size.text = 1.3
x.label = 'Predictor'
y.label = 'Response'

plot(log.yi ~ tp,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(log.yi), max(log.yi)) )
axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0) )
axis(side=2,cex.axis=size.text, las=1, mgp=c(0,0.5,0),tck=-0.01)


points(tp, log.yi,pch=symbol,cex=1.7)

abline(lm(log.yi~tp), lwd=2)
abline(lm(log.yiC~tp, data=dat), lty=2, lwd=2)

mtext(x.label, line = 1, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 1.2, side = 2, cex = size.text, outer=T)

box()
par(def.par)
dev.off()


# Fit model using lm function to censored data
summary(lm(log.yiC ~ tp, data=dat))

# Reformat data for analysis using survival package


############## BUGS MODEL #######################


# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
model {

# Likelihood: 
for (i in 1:n){
   y[i] ~ dnorm(mu[i], tau)I(,upper[i])     # Distribution for random part
   mu[i] <- a + b * x[i]                    # Linear predictor
} # i


# Priors
a ~ dnorm(0, 0.001)
b ~ dnorm(0, 0.001)
sigma ~ dunif(0, 100)


# Derived quantities
tau <- pow(sigma,-2) # precision
sigma2 <- pow(sigma,2)



} # end model
",fill = TRUE)
sink()


# Load data
data <- list(y = dat$log.yiC, x = dat$tp, n = n, upper = dat$rl  )


# Initial values
yInit <- 0.5 * rl
yInit[dat$detect==1] <- NA
inits <- function (){
  list (a=rnorm(1), b=rnorm(1), sigma=runif(1),y=yInit )
}

# Parameters monitored
parameters <- c("a","b","sigma")


# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
nc <- 3


bugs.dir <- "C:/Program Files/WinBUGS14/"

start.time = Sys.time()         # Set timer 
# Call BUGS from R 

out <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb,debug = F, bugs.directory=bugs.dir)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time

# Summarize posteriors
print(out, dig = 3)



#################################################
######## GRAPH RESULTS #############

# Set new y to where nondetects are 0.5*DL
dat$ysub <- ifelse(is.na(dat$log.yiC), 0.5*dat$rl, dat$log.yiC)


### Create fake predictor variable
fake.x <- seq(min(tp), max(tp),length=50)

#### Obtain fitted line and 95% credible interval for fitted regression line
# Create container to hold estimates
est.line <- matrix(NA, ncol=length(fake.x), nrow=out$n.sims)
# dim(est.line)

# Start stopwatch
start.time <- proc.time()
# Predict response variable for all MCMC samples and for every value of fake.x
for(i in 1:out$n.sims){
  for(j in 1:length(fake.x) ){
    est.line[i, j] <- out$sims.list$a[i] + out$sims.list$b[i] * fake.x[j]
  }
}
# End stopwatch anc calculate time
end.time <- proc.time() - start.time
# end.time[3] # time in seconds
cat('This procedure took ',end.time[3]/60, ' minutes\n\n', sep='') 


### Obtain posterior mean fitted line and associated 
# 95% CRIs (upper and lower CRIs for predicted values in the matrix est.line
fitted.mean <-  apply(est.line, 2, mean )
CRIs <- apply(est.line, 2, quantile, probs=c(0.025,0.975) )

symbol <- rep(16,length(dat$log.yiC))
symbol[log.yiC < quantile(log.yiC, 0.3)] <- 21

res <- 6
name_figure <- "BayesFit.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 		# save default, for resetting...

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1,1,1,1), oma=c(3,3,0.1,0.1) )

size.labels = 1
size.text = 1.0
y.label = 'Response'
x.label = 'Predictor'

plot(dat$log.yi ~ dat$tp, type='n', axes=F, ylim=c(min(dat$log.yi), max(dat$log.yi)),
     xlab='',ylab='', xlim=c(min(dat$tp),max(dat$tp)) )
axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0) )
axis(side=2,cex.axis=size.text, las=1, mgp=c(0,0.5,0),tck=-0.01)

# Add credible region
i.for <- order(fake.x)
i.back <- order(fake.x, decreasing = TRUE )
x.polygon <- c( fake.x[i.for] , fake.x[i.back] )
y.polygon <- c( CRIs[1,][i.for] , CRIs[2,][i.back] )
polygon( x.polygon , y.polygon , col = "lightblue" , border = NA )

# Add fitted line
lines(fake.x,fitted.mean, lwd = 10, col="darkblue", lty = 1)

# Add data points
points(dat$tp, dat$log.yi,pch=symbol,cex=1.7) # Bayes estimate
# 0.5*DL points
points(dat$tp, dat$ysub, pch='*', cex=1.5)
abline(lm(yi~tp), lwd=2, col='white') # Truth
abline(lm(ysub~tp, data=dat), lty=2, lwd=2) # Ignore censored data

# Add axis labels
mtext(x.label, line = 1, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 0.8, side = 2, cex = size.text, outer=T)

legend(-0.75, 1.5, c('Bayes','0.5*DL','* = 0.5*DL',"O = censored data"), lty=c(1,2,NA,NA), lwd=c(5,2,NA,NA), 
       text.col=c('blue','black','black','gray') ,x.intersp = 0.5)

box()
par(def.par)
dev.off()
### END PLOT




