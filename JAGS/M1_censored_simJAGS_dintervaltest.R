# rm(list=ls())

library(R2jags)
library(arm)
library(lme4)
library(MCMCpack)

nsim <- 100
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

# # Fit model using lm function
# summary(lm(log(yi) ~ tp))


######## Censor data
# log-transform yi
log.yi <- log(yi)
# log-transform yi: "C" used for censored response
log.yiC <- log(yi)


# Proportion to censor
cenProp <- 0.3

### Censor data (data < xth percentile is censored)

# Create data frame
dat <- data.frame(yi,log.yi,log.yiC,tp)

dat$censorLimitVec <- as.numeric(rep(quantile(log.yiC,cenProp ),n))

dat$isCensored <- ( dat$log.yiC > dat$censorLimitVec ) # Must tell JAGS which observations are ABOVE Censoring limit

dat$log.yiC[dat$isCensored] <- NA

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
   y[i] ~ dnorm(mu, tau)
   # mu[i] <- a + b * x[i]                   
}


# Priors
mu ~ dnorm(0, 0.001)
# b ~ dnorm(0, 0.001)
sigma ~ dunif(0, 1)


# Derived quantities
tau <- pow(sigma,-2) # precision
sigma2 <- pow(sigma,2)



} # end model
",fill = TRUE)
sink()


# JAGS dinterval needs 0,1 not T,F for isCensored
# Load data
data <- list(y = dat$log.yiC,  n = n, 
             censorLimitVec = dat$censorLimitVec,
             isCensored = as.numeric(dat$isCensored) )


# yInit = rep( NA , length(y) )
# yInit[isCensored] = censorLimitVec[isCensored]+1

# Initial values
yInit <- rep( NA , nrow(dat) )
yInit[dat$isCensored] <- dat$censorLimitVec[dat$isCensored]*0.5
inits <- function (){
  list (mu=rnorm(1), sigma=runif(1),y=yInit )
}

# Parameters monitored
parameters <- c("mu","sigma")


# MCMC settings
ni <- 3000
nt <- 1
nb <- 1000
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
    y[i] ~ dnorm(mu, tau)     # Distribution for random part
    # mu[i] <- a + b * x[i]                    # Linear predictor
    } # i
    
    
    # Priors
    mu ~ dnorm(0, 0.001)
    # b ~ dnorm(0, 0.001)
    # sigma ~ dunif(0, 100)
    # 
    
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

out0.5 <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
            model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
            n.burnin = nb,debug = F, bugs.directory=bugs.dir)

# out0.5$mean$a - a

############################################
# Set to DL
# Load data
data <- list(y = dat$ysubDL, x = dat$tp, n = n)
outDL <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
               model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
               n.burnin = nb,debug = F, bugs.directory=bugs.dir)
############################################
# Set to zero
############################################
# Load data
data <- list(y = dat$ysub0, x = dat$tp, n = n)
out0 <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
              model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
              n.burnin = nb,debug = F, bugs.directory=bugs.dir)
############################################
############################################
# Drop data
############################################
# Load data
data <- list(y = dat.drop$log.yi, x = dat.drop$tp, n = nrow(dat.drop) )
outDrop <- bugs(data = data, inits = inits, parameters.to.save = parameters, 
             model.file = "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, 
             n.burnin = nb,debug = F, bugs.directory=bugs.dir)
############################################

# array(NA, c(nsim, n.params, n.scenarios))
# Scenarios: (1) Account for censoring; (2) Substitute 0.5 * DL; (3) Substitute DL; (4) Substitute 0; (5) Drop (ignore) data
# NUmber of parameters to compare (int, slope, sigma)
# Hold MCMC for intercept
outDiff[i,1,1] <- out$mean$a - a
# Hold MCMC for slope
outDiff[i,2,1] <- out$mean$b - b
# Hold MCMC for sigma
outDiff[i,3,1] <- out$mean$sigma - sigma

outDiff[i,1,2] <- out0.5$mean$a - a
# Hold MCMC for slope
outDiff[i,2,2] <- out0.5$mean$b - b
# Hold MCMC for sigma
outDiff[i,3,2] <- out0.5$mean$sigma - sigma

outDiff[i,1,3] <- outDL$mean$a - a
# Hold MCMC for slope
outDiff[i,2,3] <- outDL$mean$b - b
# Hold MCMC for sigma
outDiff[i,3,3] <- outDL$mean$sigma - sigma

outDiff[i,1,4] <- out0$mean$a - a
# Hold MCMC for slope
outDiff[i,2,4] <- out0$mean$b - b
# Hold MCMC for sigma
outDiff[i,3,4] <- out0$mean$sigma - sigma

outDiff[i,1,5] <- outDrop$mean$a - a
# Hold MCMC for slope
outDiff[i,2,5] <- outDrop$mean$b - b
# Hold MCMC for sigma
outDiff[i,3,5] <- outDrop$mean$sigma - sigma

} # End for loop


# array(NA, c(nsim, n.params, n.scenarios))
# mean(outDiff[,1,1])
# mean(outDiff[,2,1])
# mean(outDiff[,3,1])
# 
# mean(outDiff[,1,2])
# mean(outDiff[,2,2])
# mean(outDiff[,3,2])

#####################
# Percent bias
#####################
library(hydroGOF)
# Add back in parameter value
outDiff2 <- array(NA, c(nsim, n.params, n.scenarios))
# out
outDiff2[,1,1] <- outDiff[,1,1]  + a
# Hold MCMC for slope
outDiff2[,2,1] <- outDiff[,2,1] + b
# Hold MCMC for sigma
outDiff2[,3,1] <- outDiff[,3,1] + sigma

# out0.5
outDiff2[,1,2] <- outDiff[,1,2] + a
# Hold MCMC for slope
outDiff2[,2,2] <- outDiff[,2,2] + b
# Hold MCMC for sigma
outDiff2[,3,2] <- outDiff[,3,2] + sigma

# outDL
outDiff2[,1,3] <- outDiff[,1,3] + a
# Hold MCMC for slope
outDiff2[,2,3] <- outDiff[,2,3] + b
# Hold MCMC for sigma
outDiff2[,3,3] <- outDiff[,3,3] + sigma

# out0
outDiff2[,1,4] <- outDiff[,1,4] + a
# Hold MCMC for slope
outDiff2[,2,4] <- outDiff[,2,4] + b
# Hold MCMC for sigma
outDiff2[,3,4] <- outDiff[,3,4] + sigma

# outDrop
outDiff2[,1,5] <- outDiff[,1,5] + a
# Hold MCMC for slope
outDiff2[,2,5] <- outDiff[,2,5] + b
# Hold MCMC for sigma
outDiff2[,3,5] <- outDiff[,3,5] + sigma

aa <- rep(a, nsim)
bb <- rep(b, nsim)
ss <- rep(sigma,nsim)

biasIntCen <- abs(pbias(outDiff2[,1,1], aa))
biasSlopeCen <- pbias(outDiff2[,2,1], bb)
biasSigCen <- pbias(outDiff2[,3,1], ss)

biasInt0.5 <- abs(pbias(outDiff2[,1,2], aa))
biasSlope0.5 <- pbias(outDiff2[,2,2], bb)
biasSig0.5 <- pbias(outDiff2[,3,2], ss)

biasIntDL <- abs(pbias(outDiff2[,1,3], aa))
biasSlopeDL <- pbias(outDiff2[,2,3], bb)
biasSigDL <- pbias(outDiff2[,3,3], ss)

biasInt0 <- abs(pbias(outDiff2[,1,4], aa))
biasSlope0 <- pbias(outDiff2[,2,4], bb)
biasSig0 <- pbias(outDiff2[,3,4], ss)

biasIntDrop <- abs(pbias(outDiff2[,1,5], aa))
biasSlopeDrop <- pbias(outDiff2[,2,5], bb)
biasSigDrop <- pbias(outDiff2[,3,5], ss)



# Output is scenarios in columns (censor, 0.5*DL, DL, 0, omit), parameters (int, slope, sd) in rows
mean.diff <- apply(outDiff, 2:3, mean)
quant.diff <- apply(outDiff, 2:3, quantile, c(0.025, 0.975))

# re-order for plotting
# censor, 0, 0.5*DL, DL, omit
mean.diff2 <- mean.diff[,c(1,4,2,3,5)]
quant.diff2 <- quant.diff[,,c(1,4,2,3,5)]

saveRDS(mean.diff2, "mean30percent.rds")
saveRDS(quant.diff2, "quant30percent.rds")

# #####################################################
# ########### PLOT ####################################
# #####################################################
res <- 6
name_figure <- "SimDiffPlot30Percent.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:2), nrow = 1, ncol=2,byrow=TRUE), widths = c(0.7, 0.3))
layout.show(nf)
par(mar = c(1, 1, 1, 0) + 0.2,oma=c(2,5,1,0))
def.par <- par(no.readonly = TRUE)     #

size.labels = 1
size.text = 1

x.label <- 'Difference'
y.label <- 'Scenario'

nlake <- 5 # num
### axis label options
spc <- 0.06
# lab <- 1:nlake
# lab <- c("Censored \nmodel", "Set to \n0.5 * DL", "Set to \nDL", "Set to \n0", "Omit")
# re-order for plotting
# censor, 0, 0.5*DL, DL, omit
lab <- c("Censored \nmodel", "Set to \n0", "Set to \n0.5 * DL", "Set to \nDL", "Omit")
cex <- 0.5
adj <- 0
const <- 0.1
const2 <- 2*const
###

# quant.diff[,1,1]
# mean.diff

plot(c(-0.5, 0.5), c(1-const,nlake + const), 
     axes=F, xlab='',ylab='',type='n')
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01, cex=0.8) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:nlake),labels=lab,tck= -0.01, las=1, mgp=c(0,0.5,0), cex=0.8)

# text((par("usr")[1]-0.1) - spc,1:nlake,srt = 0, adj =adj,
#      labels = lab, xpd = TRUE,cex=0.8)
abline(v=0, lwd=1)

# 95% CIs for int
segments(x0=quant.diff2[,1,][1,], x1=quant.diff2[,1,][2,],
         y0=1:nlake, y1=1:nlake, col='black',lwd=1)
## mean diff for ints
points(mean.diff2[1,], 1:nlake, col='black',cex=1, pch=16)

# 95% CIs for slopes
segments(x0=quant.diff2[,2,][1,], x1=quant.diff2[,2,][2,],
         y0=1:nlake+const, y1=1:nlake+const, col='black',lwd=1, lty=2)
## mean diff for slopes
points(mean.diff2[2,], 1:nlake+const, col='blue',cex=1, pch=16)

# 95% CIs for sd
segments(x0=quant.diff2[,3,][1,], x1=quant.diff2[,3,][2,],
         y0=1:nlake+const2, y1=1:nlake+const2, col='black',lwd=1, lty=3)
## mean diff for sd
points(mean.diff2[3,], 1:nlake+const2, col='red',cex=1, pch=16)

# Add x- and y-axis lables
mtext(x.label, line = 0.4, side = 1, cex = size.text, outer=T, adj=0.35)
mtext(y.label, line = 3, side = 2, cex = size.text, outer=T)

# abline(h=0)
box()

### ADD legend
par(mar = c(0, 0, 1, 2))
plot(1:3, rnorm(3), pch = 1, lty = 1, ylim=c(-2,2), type = "n", axes = FALSE, ann = FALSE)
legend(1, 1, c("SD", "Slope", "Intercept"), pch = c(16,16,16), lty = c(3,2,1), cex=0.8, col=c('red','blue','black'),
       x.intersp = 0.5)
###

par(def.par)
dev.off()


