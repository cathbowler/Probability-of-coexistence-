source("BH.equalizing.stabilizing.R")
source("CB.equalising.stabilizing.distributions.R")

#install.packages("gridBase")
library(lattice)
library(gridBase)
library(grid) 
library(rstan)

# make coexistence line
xseq <- seq(from=0, to=.99, by=.01) #0.01
yseq <- 1/(1-xseq)
plot(xseq, yseq)
length(xseq)
#colors <- ifelse(coexists == 0, "grey", "black")
# resetting some things 
cols = c("#F8766D","#E7861B", "#95A900", "#39B600", "#00C19C",
         "#0088D8", "#00A5FF", "#7997FF", "#DC71FA",
         "brown", "black", "#FC61D5", "#FC717F","red", "blue")

pairs <- c("A. calendula & P. airoides", "A.calendula & P.canescens", "H.glutinosum & M.minima", 
           "H. glutinosum & P. airoides", "H. glutinosum & P.deibilis", "H. glutinosum & G. rosea", 
           "M.minima & P. airoides", "P.airoides & P. debilis", "P. airoides & P. canescens", 
           "P. airoides & T. cyanopetala", "P. airoides & G. rosea", "P. debilis & P. canescens", 
           "P.canescens & T. cyanopetala", "P. canescens & G. rosea")

do.transparent <- function(InCol, NewAlpha){
  rgb_vals <- col2rgb(InCol)[,1]
  OutCol <- rgb(red = rgb_vals[1], green = rgb_vals[2], blue = rgb_vals[3],
                alpha = NewAlpha, maxColorValue = 255)
  return(OutCol)
}

tgreen <- do.transparent("seagreen4", 100)

# make graph
pdf("Figures/EqualizingStabilizingNew.pdf")
plot(mean.corrected.stabilizing, log(mean.equalizing), cex=1.5, pch=1:14, col=cols, xlab=expression("Stabilizing Niche Differences "(1-rho)), 
     ylab=expression("Log Fitness Differences "(Kappa[j]/Kappa[i])), ylim=c(0,6))
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
#text("coexistence", x=0.8, y=1)
#text("comp. exclusion", x=0.5, y=3)
legend(-.15,7.2, legend=pairs, col=cols, ncol=4, bty="n", pch=1:14, border=FALSE, cex=0.7, xpd=T)

place.x = c(.28, .65, .28, .28, .79)
place.y=c(.37, .4, .23, .64, .57)
x <- list(corrected.stabilizing[,1], corrected.stabilizing[,5],  corrected.stabilizing[,8], corrected.stabilizing[,7], corrected.stabilizing[,9])
y = list(log(equalizing[,1]), log(equalizing[,5]), log(equalizing[,8]), log(equalizing[,7]), log(equalizing[,9]))
means.x <- c(mean.corrected.stabilizing[1], mean.corrected.stabilizing[5], mean.corrected.stabilizing[8], mean.corrected.stabilizing[7],  mean.corrected.stabilizing[9])
means.y <- c(log(mean.equalizing[1]), log(mean.equalizing[5]), log(mean.equalizing[8]), log(mean.equalizing[7]), log(mean.equalizing[9]))
change.cols <- c(cols[1], cols[5], cols[8], cols[7] , cols[9])
change.pch <- c(1,5,8,7,9)
## splitting the graphs ###
for (i in 1:5){
pushViewport(viewport(x=place.x[i], y=place.y[i] ,width=.15,height=.1))
grid.rect()
par(plt = gridPLT(), new=TRUE)
plot(x[[i]], y[[i]], pch=change.pch[i], col=alpha(change.cols[i], 0.2), xlab="", ylab="", xaxt="n",yaxt="n", ylim=c(0,6), cex=0.4)
points(means.x[i], means.y[i], pch=change.pch[i], col="black")
lines(xseq, log(yseq), col=tgreen, lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
#mtext(pairs[1], side=3)
axis(2, mgp=c(3, 0, 0), cex.axis=0.7, tick=F)
axis(1, mgp=c(3, 0, 0), cex.axis=0.7, tick=F)
popViewport()
}
dev.off()

# points falling within coexistence space 
stab.int <- corrected.stabilizing[,5][which(corrected.stabilizing[,5]>HPDinterval(mcmc(corrected.stabilizing[,5]))[1] &
                                              corrected.stabilizing[,5]<HPDinterval(mcmc(corrected.stabilizing[,5]))[2])]
eq.int <- equalizing[,14][which(equalizing[,5]>HPDinterval(mcmc(equalizing[,5]))[1] &
                                 equalizing[,5]<HPDinterval(mcmc(equalizing[,5]))[2])]


#yseq <- 1/(1-xseq)

sum(eq.int[which(eq.int>0)]<1/(1-stab.int))

sum(log(equalizing[,6][which(equalizing[,6]>0)])<log(yseq)/(1-corrected.stabilizing[,6]))

log(HPDinterval(mcmc(equalizing[,5])))

par(mar=c(4,4,2,2))
plot(corrected.stabilizing[,5], log(equalizing[,5]), pch=16, col="skyblue", xlab="Stabilizing ND", ylab="Log FD", xlim=c(0,1), ylim=c(0,12))
points(mean.corrected.stabilizing[5], log(mean.equalizing[5]), pch=19, col="black")
lines(xseq, log(yseq), col="grey", lwd=2)
abline(v=HPDinterval(mcmc(corrected.stabilizing[,5]))[2], h=log(HPDinterval(mcmc(equalizing[,5]))))
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)




      