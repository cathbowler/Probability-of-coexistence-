#########################################
# PLOTTING PAIRWISE LDGR DISTRIBUTIONS  #
#########################################
library(rstan)
#source("Invasion.approach.revisions.R")
source("Invasion.approach.ldgr.means.R")
source("Invasions.approach.revisions.novar.R")

species.names = c("dagl", "hygl", "plde", "poca", "trcy", "vero", "arca", "medi", "peai")
cols = c("#F8766D","#E7861B", "#95A900", "#39B600", "#00C19C",
         "#0088D8", "#00A5FF", "#7997FF", "#DC71FA",
         "brown", "black", "#FC61D5", "#FC717F","red", "blue")

xseq <- seq(from=0, to=.99, by=.01)
yseq <- 1/(1-xseq)

do.transparent <- function(InCol, NewAlpha){
        rgb_vals <- col2rgb(InCol)[,1]
        OutCol <- rgb(red = rgb_vals[1], green = rgb_vals[2], blue = rgb_vals[3],
                      alpha = NewAlpha, maxColorValue = 255)
        return(OutCol)
}

tgreen <- do.transparent("seagreen4", 100)

# MUTUAL INVASION PAIRS
# plot ####
pdf("Figures/Mutual_invasion_pairwise_notime.pdf")
par(mfrow=c(5,3), mar=c(2,1.5,2,1.5), oma=c(3,3,3,3))

#hygl and plde
means <- c(mean.hygl.into.plde, mean.plde.into.hygl)
plot(density(hygl.into.plde, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1,3.5), ylim=c(0,3.5), xlab="")
lines(density(plde.into.hygl, na.rm=T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,3)], pch=19)
mtext("case 4", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("P. debilis"))), col=c(cols[2], cols[3]),lwd=2, lty=1, bty="n")
mtext(side=2, line=3, adj=-15, "Density")
#auc <- ecdf(plde.into.hygl[which(plde.into.hygl<HPDinterval((as.mcmc(as.numeric(plde.into.hygl))))[,2])])
dat <- cbind(hygl.into.plde,plde.into.hygl)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)


# peai and plde
means <- c(mean.peai.into.plde, mean.plde.into.peai)
plot(density(peai.into.plde, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,4), ylim=c(0,3), xlab="")
lines(density(plde.into.peai, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,3)], pch=19)
mtext("case 6", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. debilis"))), col=c(cols[9], cols[3]),lwd=2, lty=1, bty="n")

dat <- cbind(peai.into.plde,plde.into.peai)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# hygl and peai
means <- c(mean.hygl.into.peai, mean.peai.into.hygl)
plot(density(hygl.into.peai, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1.5,3), ylim=c(0,2), xlab="")
lines(density(peai.into.hygl, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,9)], pch=19)
mtext("case 6", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("P. airoides"))), col=c(cols[2], cols[9]),lwd=2, lty=1, bty="n")

dat <- cbind(hygl.into.peai,peai.into.hygl)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)


# peai and poca
means <- c(mean.peai.into.poca, mean.poca.into.peai)
plot(density(peai.into.poca, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,7), ylim=c(0,1), xlab="")
lines(density(poca.into.peai, na.rm = T), col=cols[4], main ="", lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,4)], pch=19)
species.names = c("dagl", "hygl", "plde", "poca", "trcy", "vero", "arca", "medi", "peai")
mtext("case 2", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. canescens"))), col=c(cols[9], cols[4]),lwd=2, lty=1, bty="n")

dat <- cbind(peai.into.poca,poca.into.peai)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# arca and peai
means <- c(mean.arca.into.peai, mean.peai.into.arca)
plot(density(arca.into.peai, na.rm = T), col=cols[7], main ="", lwd=2,
     xlim=c(-1,5), ylim=c(0,5), xlab="")
lines(density(peai.into.arca, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(7,9)], pch=19)
mtext("case 5", side=3, adj=0)
legend("topright", legend=c(expression(italic("A. calendula")), expression(italic("P. airdoides"))), col=c(cols[7], cols[9]),lwd=2, lty=1, bty="n")

dat <- cbind(arca.into.peai,peai.into.arca)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# hygl and medi
means <- c(mean.hygl.into.medi, mean.medi.into.hygl)
plot(density(hygl.into.medi, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-0.5,3.5), ylim=c(0,8), xlab="")
lines(density(medi.into.hygl, na.rm = T), col=cols[8], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,8)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("M. minima"))), col=c(cols[2], cols[8]),lwd=2, lty=1, bty="n")

dat <- cbind(hygl.into.medi,medi.into.hygl)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# plde and poca
means <- c(mean.plde.into.poca, mean.poca.into.plde)
plot(density(plde.into.poca, na.rm = T), col=cols[3], main ="", lwd=2,
     xlim=c(-1,5), ylim=c(0,8), xlab="")
lines(density(poca.into.plde, na.rm = T), col=cols[4], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(3,4)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. debilis")), expression(italic("P. canescens"))), col=c(cols[3], cols[4]),lwd=2, lty=1, bty="n")

dat <- cbind(plde.into.poca,poca.into.plde)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# arca and poca
means <- c(mean.arca.into.poca, mean.poca.into.arca)
plot(density(arca.into.poca, na.rm = T), col=cols[7], main ="", lwd=2,
     xlim=c(-1,6), ylim=c(0,7), xlab="")
lines(density(poca.into.arca, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(7,3)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("A. calendula")), expression(italic("P. canescens"))), col=c(cols[7], cols[3]),lwd=2, lty=1, bty="n")

dat <- cbind(arca.into.poca,poca.into.arca)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

#hygl and vero
means <- c(mean.hygl.into.vero, mean.vero.into.hygl)
plot(density(hygl.into.vero, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1,3), ylim=c(0,6), xlab="")
lines(density(vero.into.hygl, na.rm = T), col=cols[6], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,6)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("G. rosea"))), col=c(cols[2], cols[6]),lwd=2, lty=1, bty="n")

dat <- cbind(hygl.into.vero,vero.into.hygl)
length(which(dat[,2]>0))/4500

coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)


# medi and peai 
means <- c(mean.medi.into.peai, mean.peai.into.medi)
plot(density(medi.into.peai, na.rm = T), col=cols[8], main ="", lwd=2,
     xlim=c(-0.5,5), ylim=c(0,8), xlab="")
lines(density(peai.into.medi, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(8,9)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("M. minima")), expression(italic("P. airoides"))), col=c(cols[8], cols[9]),lwd=2, lty=1, bty="n")

dat <- cbind(medi.into.peai,peai.into.medi)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# peai and trcy
means <- c(mean.peai.into.trcy, mean.trcy.into.peai)
plot(density(peai.into.trcy, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1.5,5.5), ylim=c(0,6), xlab="")
lines(density(trcy.into.peai, na.rm = T), col=cols[5], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,5)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("T. cyanopetala"))), col=c(cols[9], cols[5]),lwd=2, lty=1, bty="n")

dat <- cbind(peai.into.trcy,trcy.into.peai)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# peai and vero
means <- c(mean.peai.into.vero, mean.vero.into.peai)
plot(density(peai.into.vero, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,5.5), ylim=c(0,8), xlab="")
lines(density(vero.into.peai, na.rm = T), col=cols[6], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,6)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("G. rosea"))), col=c(cols[9], cols[6]),lwd=2, lty=1, bty="n")

dat <- cbind(peai.into.vero,vero.into.peai)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# poca and trcy 
means <- c(mean.poca.into.trcy, mean.trcy.into.poca)
plot(density(poca.into.trcy, na.rm = T), col=cols[4], main ="", lwd=2,
     xlim=c(-1.5,6), ylim=c(0,5.5), xlab="")
lines(density(trcy.into.poca, na.rm = T), col=cols[5], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(4,5)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. canescens")), expression(italic("T. cyanopetala"))), col=c(cols[4], cols[5]),lwd=2, lty=1, bty="n")

dat <- cbind(poca.into.trcy,trcy.into.poca)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)


# poca and vero 
means <- c(mean.poca.into.vero, mean.vero.into.poca)
plot(density(poca.into.vero, na.rm = T), col=cols[4], main ="", lwd=2,
     xlim=c(-1,6), ylim=c(0,6), xlab="")
lines(density(vero.into.poca, na.rm = T), col=cols[6], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(4,6)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. canescens")), expression(italic("G. rosea"))), col=c(cols[4], cols[6]),lwd=2, lty=1, bty="n")
mtext(side=1, line=3, adj=4, "Low Density Growth Rate")

dat <- cbind(poca.into.vero,vero.into.poca)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)


# summary barplot
par(mar=c(2,4,1.5,1.5))
case2 <- 1/14
case3 <- 9/14
case4 <- 1/14
case5 <- 1/14
case6 <- 2/14
cases <- cbind(case2, case3, case4, case5, case6)
a <- barplot(cases, ylab = "Fraction of sp. pairs", xaxt="n")
axis(1, at=a, tck=-0.01, labels=c("case \n 2", "case \n 3", "case \n 4", "case \n 5", "case\n  6"), cex.axis=0.9)

dev.off()



## Fraction of cases figure -----------------------------------------------------------------
pdf("Figures/Fraction_of_cases.pdf", width = 6, height = 6)
par(mfrow=c(2,1), mar=c(3,4,1,1), cex=1.3)
# Mean barplot (old way)
both.coex <- 1/14
one.comp.ex <- 12/14
neither.coex <- 1/14
old.way <- cbind(both.coex, one.comp.ex, neither.coex)
a <- barplot(old.way, ylab = "",  xaxt="n", ylim=c(0,1))
axis(1, at=a, tck=-0.01, labels=c("both \n coexist", "one comp. \n excluded", "neither \n invade"), cex.axis=1)
mtext(side=3, adj=0, "(a)  mean only", cex=1.3)
mtext(side = 2, line=3, "Fraction of sp. pairs", cex=1.3)

# New way 
both.coex <- 1/14
one.comp.ex <- 9/14
neither.coex <- 0
dist.dep <- 4/14
new.way <- cbind(both.coex, one.comp.ex, neither.coex, dist.dep)
b <- barplot(new.way, ylab = "", xaxt="n",  ylim=c(0,1))
axis(1, at=b, tck=-0.01, cex.axis=1, labels=c("both \n coexist", "one comp. \n excluded", "neither \n invade", "distribution \n dependent"))
mtext(side=3, adj=0, "(b) based on distributions", cex=1.3)
mtext(side = 2, line=3, "Fraction of sp. pairs", cex=1.3)

dev.off()


## Example cases #####
pdf("Figures/Example_of_cases.pdf", width = 10, height = 6)
par(mfrow=c(2,2), mar=c(3,3,1,1), cex=1.3)
# peai and poca
means <- c(mean.peai.into.hygl, mean.hygl.into.peai)
plot(density(peai.into.hygl, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-2,3), ylim=c(0,1.5), xlab="")
lines(density(hygl.into.peai, na.rm = T), col=cols[4], main ="", lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,4)], pch=19)
#species.names = c("dagl", "hygl", "plde", "poca", "trcy", "vero", "arca", "medi", "peai")
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("H. glutinosum"))), col=c(cols[9], cols[4]),lwd=2, lty=1, bty="n")
mtext(side=3, adj=0, "(a) case 2", cex=1.3)
mtext(side = 2, line=2, "Density", cex=1.3)
dat <- cbind(peai.into.hygl,hygl.into.peai)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

#hygl and plde
means <- c(mean.hygl.into.plde, mean.plde.into.hygl)
plot(density(hygl.into.plde, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1,3.5), ylim=c(0,3.5), xlab="", ylab="")
lines(density(plde.into.hygl, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,3)], pch=19)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("P. debilis"))), col=c(cols[2], cols[3]),lwd=2, lty=1, bty="n")
mtext(side=3, adj=0, "(b) case 4", cex=1.3)
dat <- cbind(hygl.into.plde,plde.into.hygl)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# arca and peai
means <- c(mean.arca.into.peai, mean.peai.into.arca)
plot(density(arca.into.peai, na.rm = T), col=cols[7], main ="", lwd=2,
     xlim=c(-1,5), ylim=c(0,5), xlab="")
lines(density(peai.into.arca, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(7,9)], pch=19)
mtext(side=3, adj=0, "(c) case 5", cex=1.3)
legend("topright", legend=c(expression(italic("A. calendula")), expression(italic("P. airdoides"))), col=c(cols[7], cols[9]),lwd=2, lty=1, bty="n")
mtext(side = 2, line=2, "Density", cex=1.3)
mtext("Low density growth rate", side=1, line=2, cex=1.3)
dat <- cbind(arca.into.peai,peai.into.arca)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

# peai and plde
means <- c(mean.peai.into.plde, mean.plde.into.peai)
plot(density(peai.into.plde, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,4), ylim=c(0,3), xlab="", ylab="")
lines(density(plde.into.peai, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,3)], pch=19)
mtext("Low density growth rate", side=1, line=2, cex=1.3)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. debilis"))), col=c(cols[9], cols[3]),lwd=2, lty=1, bty="n")
mtext(side=3, adj=0, "(d) case 6", cex=1.3)
dat <- cbind(peai.into.plde,plde.into.peai)
coex <- length(which(dat[,1]>0&dat[,2]>0))/4500
mtext(round(coex*100,2), side=3, adj=1, cex=0.8)

 dev.off()




# Cases figure ------------------------------------------------------------------

pdf("Figures/Cases_figure.pdf", width = 5, height = 3.5)

layout(matrix(c(1,1,1,
                2,3,4,
                5,6,7), 3, 3, byrow = TRUE), heights=c(.5,2,2))
par(mar=c(1,1,1,1))
plot.new()
legend("center", legend=c("Species 1", "Species 2"), ncol=2, col=c("black", "grey"),lwd=3, lty=1, bty="n", cex=1.3, xpd=T)
par(mar=c(3,2,1,1))
# case 1
plot(density(peai.into.poca-10, na.rm = T, bw = 0.8), col="black", lwd=3, yaxt="n", xaxt="n",
     xlim=c(-12,1), ylim=c(0,.4), main="CASE 1", cex=1.3)
lines(density(peai.into.poca-7, na.rm = T, bw = 0.8), col="grey", main ="", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)
mtext("Density", side=2)

# case 2
plot(density(peai.into.poca+5, na.rm = T, bw = 0.8), col="black", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-1,15), ylim=c(0,.4), main="CASE 2", cex=1.3)
lines(density(peai.into.poca+8, na.rm = T, bw = 0.8), col="grey", main ="", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)

# case 3
plot(density(peai.into.poca-9, na.rm = T, bw = 0.8), col="black", main ="CASE 3", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-12,12), ylim=c(0,.4))
lines(density(peai.into.poca+5, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)

# case 4
plot(density(peai.into.poca-2, na.rm = T, bw = 0.8), col="black", main ="CASE 4", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-5,12), ylim=c(0,.4))
lines(density(peai.into.poca+5, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
mtext("LDGR", side=1, line=2)
mtext("Density", side=2)
axis(side = 1, 0)

# case 5
plot(density(peai.into.poca-8, na.rm = T, bw = 0.8), col="black", main ="CASE 5", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-12,5), ylim=c(0,.4))
lines(density(peai.into.poca-2, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)
mtext("LDGR", side=1, line=2)

# case 6
plot(density(peai.into.poca-1, na.rm = T, bw = 0.8), col="black", main ="CASE 6", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-12,12), ylim=c(0,.4))
lines(density(peai.into.poca-2, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)
mtext("LDGR", side=1, line=2)

dev.off()


