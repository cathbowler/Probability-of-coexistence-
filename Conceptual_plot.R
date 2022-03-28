#############################
##### Conceptual figure #####
#############################
library(coda)
source("invasion.approach.ldgr.means.R")
source("invasions.approach.revisions.novar.R")

do.transparent <- function(InCol, NewAlpha){
  rgb_vals <- col2rgb(InCol)[,1]
  OutCol <- rgb(red = rgb_vals[1], green = rgb_vals[2], blue = rgb_vals[3],
                alpha = NewAlpha, maxColorValue = 255)
  return(OutCol)
}

tgreen <- do.transparent("seagreen4", 100)

xseq <- seq(from=0, to=.99, by=.01)
yseq <- 1/(1-xseq)


# CONCEPTUAL PLOT #####
pdf("Figures/conceptual_plot.pdf")
colour = c("turquoise", "palevioletred")
Species1 <- 0.8
Species2 <- 0.2
sp <- as.data.frame(rbind(Species1, Species2))
par(mfrow=c(2,2), mar=c(2,2,2,2), oma=c(2,2,1,1), cex.axis=1.5)
# Invasion 1
barplot(height=sp$V1, names.arg = c("Species 1", "Species 2"), ylab="", ylim=c(-1,1), col = colour,
        space = 0.5, border = FALSE, yaxt="n", cex.axis = 1.3)
abline(h=0, lty=2)
text(1.7, 0.9, "coexists", cex=1.5)
text(1.7, -0.3, "competitive \n exclusion", cex=1.5)
text(0.3, 0, "0", xpd=TRUE, cex=1.3)
mtext("Low density growth rate", side=2, line=2, cex=1.3)
mtext("(a)", side=3, adj=0, cex=1.3)
box()

# ES 1
plot(0.8, 1.4, pch=19, col="black", xlab="", ylab="", xlim=c(0,0.95), ylim=c(0.5,5), xaxt="n", yaxt="n", cex.axis=1.3)
mtext("Stabilizing Niche Differences", side=1, line=1, cex=1.3)
mtext("Fitness Differences", side=2, line=1, cex=1.3)
lines(xseq, log(yseq), col="grey", lwd=2)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
text(0.75, 0.8, "coexists", col="darkgreen", cex=1.5)
text(0.5, 3, "competitive \n exclusion", cex=1.5)
mtext("(b)", side=3, adj=0, cex=1.3)

# Invasion 2
plot(density(hygl.into.plde, na.rm = T), main="", xlab="", ylab = "", lwd=3, col="turquoise",
     xlim=c(-2,3.5), yaxt="n", xaxt="n", ylim=c(0,3), cex.axis=1.3)
lines(density(plde.into.hygl, na.rm = T), col="palevioletred", lwd=3)
abline(v=0, lty=2)
mtext("Low density growth rate", side=1, line=2, cex=1.3)
text(0, -0.3, "0", xpd = TRUE, cex=1.3)
points(mean.hygl.into.plde, y=0, pch=19, col="turquoise")
points(mean.plde.into.hygl+0.3, y=0, pch=19, col="palevioletred")
text(2,1.5, "coexists", cex=1.5)
text(-1, 1.5, "competitive \n exclusion", cex=.8, cex=1.5)
legend("topright", legend=c("Species 1", "Species 2"), cex=1.3, lty=1, lwd=3, col = c("turquoise", "palevioletred"), bty="n")
mtext("(c)", side=3, adj=0, cex=1.3)
mtext(side=2, line=2, "Probability density", cex=1.3)
#rect(-0,-0.5,6.5,2,col = rgb(0.5,0.5,0.5,1/4), border=NA)

# ES 2
x<-rnorm(0.8,0.08, n=20)
y<-rnorm(1.6,0.3, n=20)

# calculate HPD invtervals
x.int <- (as.numeric(HPDinterval(mcmc(x))[,2]) - as.numeric(HPDinterval(mcmc(x))[,1]))/2
y.int <- (as.numeric(HPDinterval(mcmc(y))[,2]) - as.numeric(HPDinterval(mcmc(y))[,1]))/2

plot(x, y, pch=19, col="black", xlab="", ylab="", xlim=c(0,0.95), ylim=c(0.5,5), xaxt="n", yaxt="n")
mtext("Stabilizing Niche Differences", side=1, line=1, cex=1.3)
mtext("Fitness Differences", side=2, line=1, cex=1.3)
arrows(mean(x), mean(y)+y.int, mean(x), mean(y)-y.int, length=0.05, angle=90, code=3)
arrows(mean(x)-x.int, mean(y), mean(x)+x.int, mean(y), length=0.05, angle=90, code=3)
polygon(c(xseq, rev(xseq)), c(log(yseq), rep(0, length(xseq))), col=tgreen)
text(0.7, 0.8, "coexists", col="darkgreen", cex=1.5)
text(0.5, 3, "competitive \n exclusion", cex=1.5)
mtext("(d)", side=3, adj=0, cex=1.3)

dev.off()

