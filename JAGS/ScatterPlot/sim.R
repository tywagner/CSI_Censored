# rm(list=ls())

library(R2jags)
library(arm)
library(lme4)
library(MCMCpack)
library(colorspace)
library(RColorBrewer)

nsim <- 1
# Scenarios: (1) Account for censoring; (2) Substitute 0.5 * DL; (3) Substitute DL; (4) Substitute 0; (5) Drop (ignore) data
n.scenarios <- 5
# NUmber of parameters to compare (int, slope, sigma)
n.params <- 3
# Number of MCMC iterations saved
n.iter <- 6000

# Hold MCMC for intercept
outDiff <- array(NA, c(nsim, n.params, n.scenarios))



for(i in 1:nsim){

####### Generate data #############
# number of lakes
n = 10000
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

# # Fit model using lm function
# summary(lm(log(yi) ~ tp))


######## Censor data
# log-transform yi
log.yi <- log(yi)
# log-transform yi: "C" used for censored response
log.yiC <- log(yi)


# Proportion to censor
cenProp <- 0.30

### Censor data (data < xth percentile is censored)

# Create data frame
dat <- data.frame(yi,log.yi,log.yiC,tp)

dat$censorLimitVec <- as.numeric(rep(quantile(log.yiC,cenProp ),n))

dat$isCensored <- ( dat$log.yiC > dat$censorLimitVec ) # Must tell JAGS which observations are ABOVE Censoring limit

dat$log.yiC[!dat$isCensored] <- NA

# Set new y to where nondetects are 0.5*DL
dat$ysub0.5 <- ifelse(is.na(dat$log.yiC), 0.5*dat$censorLimit, dat$log.yiC)
# Set new y to where nondetects are DL
dat$ysubDL <- ifelse(is.na(dat$log.yiC), dat$censorLimit, dat$log.yiC)
# Set new y to where nondetects are 0
dat$ysub0 <- ifelse(is.na(dat$log.yiC), 0, dat$log.yiC)
# Create new data frame where we drop all censored data
dat.drop <- dat[!is.na(dat$log.yiC),]

# head(dat)
# head(dat.drop)

# Censored data model ###########################
############## BUGS MODEL #######################


# Define the model in the BUGS language and write a text file
sink("model.txt")
cat("
model {

# Likelihood: 
for (i in 1:n){
   isCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
   y[i] ~ dnorm(mu[i], tau)
   mu[i] <- a + b * x[i]                   
}


# Priors
a ~ dnorm(0, 0.001)
b ~ dnorm(0, 0.001)
sigma ~ dunif(0, 1)


# Derived quantities
tau <- pow(sigma,-2) # precision
sigma2 <- pow(sigma,2)



} # end model
",fill = TRUE)
sink()


# JAGS dinterval needs 0,1 not T,F for isCensored
# Load data
data <- list(y = dat$log.yiC, x = dat$tp, n = n, 
             censorLimitVec = dat$censorLimitVec,
             isCensored = as.numeric(dat$isCensored) )


# yInit = rep( NA , length(y) )
# yInit[isCensored] = censorLimitVec[isCensored]+1

# Initial values
yInit <- rep( NA , nrow(dat) )
yInit[!dat$isCensored] <- dat$censorLimitVec[!dat$isCensored]*0.5
inits <- function (){
  list (a=a, b=b, sigma=sigma,y=yInit )
}

# Parameters monitored
parameters <- c("a","b","sigma")


# MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 3


out <- jags(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb)
print(out)
################# END Censored data model

# Standard regression data model ###########################
############## BUGS MODEL #######################

# Define the model in the BUGS language and write a text file
sink("model2.txt")
cat("
    model {
    
    # Likelihood: 
    for (i in 1:n){
    y[i] ~ dnorm(mu[i], tau)     # Distribution for random part
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

#################################
# Load data
# Set to 0.5*DL
data <- list(y = dat$ysub0.5, x = dat$tp, n = n)


# Initial values
inits <- function (){
  list (a=rnorm(1), b=rnorm(1), sigma=runif(1) )
}

out0.5 <- jags(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb)

# out0.5$mean$a - a

############################################
# Set to DL
# Load data
data <- list(y = dat$ysubDL, x = dat$tp, n = n)
outDL <- jags(data = data, inits = inits, parameters.to.save = parameters, 
               model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
               n.burnin = nb)
############################################
# Set to zero
############################################
# Load data
data <- list(y = dat$ysub0, x = dat$tp, n = n)
out0 <- jags(data = data, inits = inits, parameters.to.save = parameters, 
              model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
              n.burnin = nb)
############################################
############################################
# Drop data
############################################
# Load data
data <- list(y = dat.drop$log.yi, x = dat.drop$tp, n = nrow(dat.drop) )
outDrop <- jags(data = data, inits = inits, parameters.to.save = parameters, 
             model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb)
############################################

# array(NA, c(nsim, n.params, n.scenarios))
# Scenarios: (1) Account for censoring; (2) Substitute 0.5 * DL; (3) Substitute DL; (4) Substitute 0; (5) Drop (ignore) data
# NUmber of parameters to compare (int, slope, sigma)
# Hold MCMC for intercept
outDiff[i,1,1] <- out$BUGSoutput$mean$a 
# Hold MCMC for slope
outDiff[i,2,1] <- out$BUGSoutput$mean$b 
# Hold MCMC for sigma
outDiff[i,3,1] <- out$BUGSoutput$mean$sigma 

outDiff[i,1,2] <- out0.5$BUGSoutput$mean$a
# Hold MCMC for slope
outDiff[i,2,2] <- out0.5$BUGSoutput$mean$b 
# Hold MCMC for sigma
outDiff[i,3,2] <- out0.5$BUGSoutput$mean$sigma 

outDiff[i,1,3] <- outDL$BUGSoutput$mean$a
# Hold MCMC for slope
outDiff[i,2,3] <- outDL$BUGSoutput$mean$b 
# Hold MCMC for sigma
outDiff[i,3,3] <- outDL$BUGSoutput$mean$sigma 

outDiff[i,1,4] <- out0$BUGSoutput$mean$a 
# Hold MCMC for slope
outDiff[i,2,4] <- out0$BUGSoutput$mean$b 
# Hold MCMC for sigma
outDiff[i,3,4] <- out0$BUGSoutput$mean$sigma 

outDiff[i,1,5] <- outDrop$BUGSoutput$mean$a 
# Hold MCMC for slope
outDiff[i,2,5] <- outDrop$BUGSoutput$mean$b 
# Hold MCMC for sigma
outDiff[i,3,5] <- outDrop$BUGSoutput$mean$sigma 

} # End for loop



#####################
# Percent bias
#####################
# library(hydroGOF)
# # Add back in parameter value
# outDiff2 <- array(NA, c(nsim, n.params, n.scenarios))
# # out
# outDiff2[,1,1] <- outDiff[,1,1]  + a
# # Hold MCMC for slope
# outDiff2[,2,1] <- outDiff[,2,1] + b
# # Hold MCMC for sigma
# outDiff2[,3,1] <- outDiff[,3,1] + sigma
# 
# # out0.5
# outDiff2[,1,2] <- outDiff[,1,2] + a
# # Hold MCMC for slope
# outDiff2[,2,2] <- outDiff[,2,2] + b
# # Hold MCMC for sigma
# outDiff2[,3,2] <- outDiff[,3,2] + sigma
# 
# # outDL
# outDiff2[,1,3] <- outDiff[,1,3] + a
# # Hold MCMC for slope
# outDiff2[,2,3] <- outDiff[,2,3] + b
# # Hold MCMC for sigma
# outDiff2[,3,3] <- outDiff[,3,3] + sigma
# 
# # out0
# outDiff2[,1,4] <- outDiff[,1,4] + a
# # Hold MCMC for slope
# outDiff2[,2,4] <- outDiff[,2,4] + b
# # Hold MCMC for sigma
# outDiff2[,3,4] <- outDiff[,3,4] + sigma
# 
# # outDrop
# outDiff2[,1,5] <- outDiff[,1,5] + a
# # Hold MCMC for slope
# outDiff2[,2,5] <- outDiff[,2,5] + b
# # Hold MCMC for sigma
# outDiff2[,3,5] <- outDiff[,3,5] + sigma
# 
# aa <- rep(a, nsim)
# bb <- rep(b, nsim)
# ss <- rep(sigma,nsim)
# 
# biasIntCen <- abs(pbias(outDiff2[,1,1], aa))
# biasSlopeCen <- pbias(outDiff2[,2,1], bb)
# biasSigCen <- pbias(outDiff2[,3,1], ss)
# 
# biasInt0.5 <- abs(pbias(outDiff2[,1,2], aa))
# biasSlope0.5 <- pbias(outDiff2[,2,2], bb)
# biasSig0.5 <- pbias(outDiff2[,3,2], ss)
# 
# biasIntDL <- abs(pbias(outDiff2[,1,3], aa))
# biasSlopeDL <- pbias(outDiff2[,2,3], bb)
# biasSigDL <- pbias(outDiff2[,3,3], ss)
# 
# biasInt0 <- abs(pbias(outDiff2[,1,4], aa))
# biasSlope0 <- pbias(outDiff2[,2,4], bb)
# biasSig0 <- pbias(outDiff2[,3,4], ss)
# 
# biasIntDrop <- abs(pbias(outDiff2[,1,5], aa))
# biasSlopeDrop <- pbias(outDiff2[,2,5], bb)
# biasSigDrop <- pbias(outDiff2[,3,5], ss)



# Output is scenarios in columns (censor, 0.5*DL, DL, 0, omit), parameters (int, slope, sd) in rows
mean.diff <- apply(outDiff, 2:3, mean)
quant.diff <- apply(outDiff, 2:3, quantile, c(0.025, 0.975))

# re-order for plotting
# censor, 0, 0.5*DL, DL, omit
mean.diff2 <- mean.diff[,c(1,4,2,3,5)]
quant.diff2 <- quant.diff[,,c(1,4,2,3,5)]

# saveRDS(mean.diff2, "mean5percent.rds")
# saveRDS(quant.diff2, "quant5percent.rds")

# #####################################################
# ########### PLOT ####################################
# #####################################################

# plotcol <- diverge_hcl(6, h = c(260, 0), c = 80, l = c(30, 90), power = 1.5)
# plotcol <- diverge_hsv(6)

plotcol <- brewer.pal(n = 6, name = "Dark2")

symbol <- rep(16,length(dat$log.yi))
symbol[log.yi < quantile(dat$log.yi, 0.3)] <- 21

res <- 6
name_figure <- "ScatterPlot.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 		# save default, for resetting...

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1,1,0,0), oma=c(2.0,2.0,0.1,0.1) )

size.labels = 1
size.text = 1.0
y.label = 'Response'
x.label = 'Predictor'

plot(log.yi ~ tp, type='n', axes=F, ylim=c(min(log.yi), max(log.yi)),
     xlab='',ylab='', xlim=c(min(tp),max(tp)), data=dat )
axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0) )
axis(side=2,cex.axis=size.text, las=1, mgp=c(0,0.5,0),tck=-0.01)

# Add data points
points(dat$tp, dat$log.yi,pch=symbol,cex=0.5) 
# # 0.5*DL points
# points(dat$tp, dat$ysub, pch='*', cex=1.5)


# Add fitted lines
# censor, 0, 0.5*DL, DL, omit
# Truth
abline(a,b, lwd = 5, col=plotcol[1], lty = 1)
# Censored
abline(mean.diff2[1,1],mean.diff2[2,1], lwd = 5, col=plotcol[2], lty = 1)
# Set to 0
abline(mean.diff2[1,2],mean.diff2[2,2], lwd = 5, col=plotcol[3], lty = 1)
# 0.5*DL
abline(mean.diff2[1,3],mean.diff2[2,3], lwd = 5, col=plotcol[4], lty = 1)
# DL
abline(mean.diff2[1,4],mean.diff2[2,4], lwd = 5, col=plotcol[5], lty = 1)
# Omit
abline(mean.diff2[1,5],mean.diff2[2,5], lwd = 5, col=plotcol[6], lty = 1)


# Add axis labels
mtext(x.label, line = 0.8, side = 1, cex = size.text, outer=T)
mtext(y.label, line = 0.9, side = 2, cex = size.text, outer=T)

legend(-0.95, 2.9, c('Truth','Censored','Set to 0', 'Set to 0.5*DL','Set to DL',"Omit"), lty=c(1,1,1,1,1,1), 
       lwd=c(5,5,5,5,5,5), col=c(plotcol[1],plotcol[2],plotcol[3],plotcol[4],plotcol[5],plotcol[6]),
       text.col=c('black','black','black','black','black','black') ,x.intersp = 0.5)

box()
par(def.par)
dev.off()
### END PLOT


