# rm(list=ls())

library(R2jags)
library(arm)
library(lme4)
library(MCMCpack)

# 5%
mean15 <- readRDS('mean5percent.rds')
quant15 <-readRDS('quant5percent.rds')

# 15%
mean30 <- readRDS('mean15percent.rds')
quant30 <-readRDS('quant15percent.rds')

# 30%
mean50 <- readRDS('mean30percent.rds')
quant50 <-readRDS('quant30percent.rds')



# #####################################################
# ########### PLOT ####################################
# #####################################################
res <- 6
name_figure <- "SimDiffPlotPanel.jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)     # save default, for resetting...

nf <- layout(matrix(c(1:3), nrow = 3, ncol=1,byrow=TRUE))
layout.show(nf)
par(mar = c(1, 1, 1, 0) + 0.2,oma=c(2.5,5.5,0.1,0), mai=c(0.1, 0, 0.1, 0))
def.par <- par(no.readonly = TRUE)     #

size.labels = 1
size.text = 1

x.label <- 'Difference (estimate - true value)'
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
const <- 0.15
const2 <- 2*const
###

## 15%
plot(c(-0.4, 0.4), c(1-const,nlake + const), 
     axes=F, xlab='',ylab='',type='n', ylim=c(c(1-const,nlake + 0.3)))
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01, cex=0.8, labels=FALSE) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:nlake),labels=lab,tck= -0.01, las=1, mgp=c(0,0.5,0), cex=0.8)

R# text((par("usr")[1]-0.1) - spc,1:nlake,srt = 0, adj =adj,
#      labels = lab, xpd = TRUE,cex=0.8)
abline(v=0, lwd=1)

# 95% CIs for int
segments(x0=quant15[,1,][1,], x1=quant15[,1,][2,],
         y0=1:nlake, y1=1:nlake, col='black',lwd=1)
## mean diff for ints
points(mean15[1,], 1:nlake, col='black',cex=1, pch=1)

# 95% CIs for slopes
segments(x0=quant15[,2,][1,], x1=quant15[,2,][2,],
         y0=1:nlake+const, y1=1:nlake+const, col='black',lwd=1, lty=1)
## mean diff for slopes
points(mean15[2,], 1:nlake+const, col='blue',cex=1, pch=2)

# 95% CIs for sd
segments(x0=quant15[,3,][1,], x1=quant15[,3,][2,],
         y0=1:nlake+const2, y1=1:nlake+const2, col='black',lwd=1, lty=1)
## mean diff for sd
points(mean15[3,], 1:nlake+const2, col='red',cex=1, pch=0)

# Add x- and y-axis lables
mtext(x.label, line = 1, side = 1, cex = size.text, outer=T, adj=0.5)
mtext(y.label, line = 4, side = 2, cex = size.text, outer=T)


legend(0.3, 5, c("SD", "Slope", "Intercept"), pch = c(0,2,1), lty = c(1,1,1), 
       cex=1.0, col=c('red','blue','black'),
       x.intersp = 0.6)

text(-0.4, 5, '(A)')

box()

## 30%
plot(c(-0.4, 0.4), c(1-const,nlake + const), 
     axes=F, xlab='',ylab='',type='n', ylim=c(c(1-const,nlake + 0.3)))
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01, cex=0.8, labels=FALSE) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:nlake),labels=lab,tck= -0.01, las=1, mgp=c(0,0.5,0), cex=0.8)

# text((par("usr")[1]-0.1) - spc,1:nlake,srt = 0, adj =adj,
#      labels = lab, xpd = TRUE,cex=0.8)
abline(v=0, lwd=1)

# 95% CIs for int
segments(x0=quant30[,1,][1,], x1=quant30[,1,][2,],
         y0=1:nlake, y1=1:nlake, col='black',lwd=1)
## mean diff for ints
points(mean30[1,], 1:nlake, col='black',cex=1, pch=1)

# 95% CIs for slopes
segments(x0=quant30[,2,][1,], x1=quant30[,2,][2,],
         y0=1:nlake+const, y1=1:nlake+const, col='black',lwd=1, lty=1)
## mean diff for slopes
points(mean30[2,], 1:nlake+const, col='blue',cex=1, pch=2)

# 95% CIs for sd
segments(x0=quant30[,3,][1,], x1=quant30[,3,][2,],
         y0=1:nlake+const2, y1=1:nlake+const2, col='black',lwd=1, lty=1)
## mean diff for sd
points(mean30[3,], 1:nlake+const2, col='red',cex=1, pch=0)

# Add x- and y-axis lables
# mtext(x.label, line = 0.4, side = 1, cex = size.text, outer=T, adj=0.35)
# mtext(y.label, line = 3, side = 2, cex = size.text, outer=T)


text(-0.4, 5, '(B)')

box()


## 50%
plot(c(-0.4, 0.4), c(1-const,nlake + const), 
     axes=F, xlab='',ylab='',type='n', ylim=c(c(1-const,nlake + 0.3)))
axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01, cex=0.8) #, at=xlab1, labels=round(xlab2,2)
#axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)
axis(side=2,at=c(1:nlake),labels=lab,tck= -0.01, las=1, mgp=c(0,0.5,0), cex=0.8)

# text((par("usr")[1]-0.1) - spc,1:nlake,srt = 0, adj =adj,
#      labels = lab, xpd = TRUE,cex=0.8)
abline(v=0, lwd=1)

# 95% CIs for int
segments(x0=quant50[,1,][1,], x1=quant50[,1,][2,],
         y0=1:nlake, y1=1:nlake, col='black',lwd=1)
## mean diff for ints
points(mean50[1,], 1:nlake, col='black',cex=1, pch=1)

# 95% CIs for slopes
segments(x0=quant50[,2,][1,], x1=quant50[,2,][2,],
         y0=1:nlake+const, y1=1:nlake+const, col='black',lwd=1, lty=1)
## mean diff for slopes
points(mean50[2,], 1:nlake+const, col='blue',cex=1, pch=2)

# 95% CIs for sd
segments(x0=quant50[,3,][1,], x1=quant50[,3,][2,],
         y0=1:nlake+const2, y1=1:nlake+const2, col='black',lwd=1, lty=1)
## mean diff for sd
points(mean50[3,], 1:nlake+const2, col='red',cex=1, pch=0)

# Add x- and y-axis lables
# mtext(x.label, line = 0.4, side = 1, cex = size.text, outer=T, adj=0.35)
# mtext(y.label, line = 3, side = 2, cex = size.text, outer=T)


text(-0.4, 5, '(C)')

box()


par(def.par)
dev.off()


