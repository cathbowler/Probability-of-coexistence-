###############################################################
# PLOTTING PAIRWISE LDGR DISTRIBUTIONS from no lottery models #
###############################################################
library(coda)
library(rstan)
source("Invasion.approach.revisions.R")
source("Invasion.approach.ldgr.means.R")
#source("Invasions.approach.revisions.novar.R")

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
pdf("Figures/Nolot_Mutual_invasion_pairwise.pdf")
#pdf("Figures/Nolot_Mutual_invasion_pairwise_no_time_var.pdf")
par(mfrow=c(5,3), mar=c(2,1.5,2,1.5), oma=c(3,3,3,3))

#hygl and plde
means <- c(mean.hygl.into.plde, mean.plde.into.hygl)
plot(density(nolot.hygl.into.plde, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1,3.5), ylim=c(0,3.5), xlab="")
lines(density(nolot.plde.into.hygl, na.rm=T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,3)], pch=19)
mtext("case 4", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("P. debilis"))), col=c(cols[2], cols[3]),lwd=2, lty=1, bty="n")
mtext(side=2, line=3, adj=-15, "Density")

auc <- ecdf(nolot.plde.into.hygl[which(nolot.plde.into.hygl<HPDinterval((as.mcmc(as.numeric(nolot.plde.into.hygl))))[,2])])
1-auc(0) # gives percentage of distribution greater than zero and below upper bound of CI

# peai and plde
means <- c(mean.peai.into.plde, mean.plde.into.peai)
plot(density(nolot.peai.into.plde, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,4), ylim=c(0,3), xlab="")
lines(density(nolot.plde.into.peai, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,3)], pch=19)
mtext("case 6", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. debilis"))), col=c(cols[9], cols[3]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.peai.into.plde[which(nolot.peai.into.plde<HPDinterval((as.mcmc(as.numeric(nolot.peai.into.plde))))[,2])])
1-auc(0)
auc <- ecdf(nolot.plde.into.peai[which(nolot.plde.into.peai<HPDinterval((as.mcmc(as.numeric(nolot.plde.into.peai))))[,2])])
1-auc(0) 

# hygl and peai
means <- c(mean.hygl.into.peai, mean.peai.into.hygl)
plot(density(nolot.hygl.into.peai, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1.5,3), ylim=c(0,2), xlab="")
lines(density(nolot.peai.into.hygl, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,9)], pch=19)
mtext("case 6", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("P. airoides"))), col=c(cols[2], cols[9]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.hygl.into.peai[which(nolot.hygl.into.peai<HPDinterval((as.mcmc(as.numeric(nolot.hygl.into.peai))))[,2])])
1-auc(0)
auc <- ecdf(nolot.peai.into.hygl[which(nolot.peai.into.hygl<HPDinterval((as.mcmc(as.numeric(nolot.peai.into.hygl))))[,2])])
1-auc(0) 

# peai and poca
means <- c(mean.peai.into.poca, mean.poca.into.peai)
plot(density(nolot.peai.into.poca, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,7), ylim=c(0,1), xlab="")
lines(density(nolot.poca.into.peai, na.rm = T), col=cols[4], main ="", lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,4)], pch=19)
species.names = c("dagl", "hygl", "plde", "poca", "trcy", "vero", "arca", "medi", "peai")
mtext("case 2", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. canescens"))), col=c(cols[9], cols[4]),lwd=2, lty=1, bty="n")

# arca and peai
means <- c(mean.arca.into.peai, mean.peai.into.arca)
plot(density(nolot.arca.into.peai, na.rm = T), col=cols[7], main ="", lwd=2,
     xlim=c(-1,5), ylim=c(0,5), xlab="")
lines(density(nolot.peai.into.arca, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(7,9)], pch=19)
mtext("case 5", side=3, adj=0)
legend("topright", legend=c(expression(italic("A. calendula")), expression(italic("P. airdoides"))), col=c(cols[7], cols[9]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.arca.into.peai[which(nolot.arca.into.peai<HPDinterval((as.mcmc(as.numeric(nolot.arca.into.peai))))[,2])])
1-auc(0)
auc <- ecdf(nolot.peai.into.arca[which(nolot.peai.into.arca<HPDinterval((as.mcmc(as.numeric(nolot.peai.into.arca))))[,2])])
1-auc(0)

# hygl and medi
means <- c(mean.hygl.into.medi, mean.medi.into.hygl)
plot(density(nolot.hygl.into.medi, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-0.5,3.5), ylim=c(0,8), xlab="")
lines(density(nolot.medi.into.hygl, na.rm = T), col=cols[8], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,8)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("M. minima"))), col=c(cols[2], cols[8]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.hygl.into.medi[which(nolot.hygl.into.medi<HPDinterval((as.mcmc(as.numeric(nolot.hygl.into.medi))))[,2])])
1-auc(0)
auc <- ecdf(nolot.medi.into.hygl[which(nolot.medi.into.hygl<HPDinterval((as.mcmc(as.numeric(nolot.medi.into.hygl))))[,2])])
1-auc(0)

auc <- ecdf(nolot.hygl.into.medi)
1-auc(0)
auc <- ecdf(nolot.medi.into.hygl)
1-auc(0)


# plde and poca
means <- c(mean.plde.into.poca, mean.poca.into.plde)
plot(density(nolot.plde.into.poca, na.rm = T), col=cols[3], main ="", lwd=2,
     xlim=c(-1,5), ylim=c(0,8), xlab="")
lines(density(nolot.poca.into.plde, na.rm = T), col=cols[4], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(3,4)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. debilis")), expression(italic("P. canescens"))), col=c(cols[3], cols[4]),lwd=2, lty=1, bty="n")

# arca and poca
means <- c(mean.arca.into.poca, mean.poca.into.arca)
plot(density(nolot.arca.into.poca, na.rm = T), col=cols[7], main ="", lwd=2,
     xlim=c(-1,6), ylim=c(0,7), xlab="")
lines(density(nolot.poca.into.arca, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(7,3)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("A. calendula")), expression(italic("P. canescens"))), col=c(cols[7], cols[3]),lwd=2, lty=1, bty="n")

#hygl and vero
means <- c(mean.hygl.into.vero, mean.vero.into.hygl)
plot(density(nolot.hygl.into.vero, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1,3), ylim=c(0,6), xlab="")
lines(density(nolot.vero.into.hygl, na.rm = T), col=cols[6], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,6)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("G. rosea"))), col=c(cols[2], cols[6]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.hygl.into.vero[which(nolot.hygl.into.vero<HPDinterval((as.mcmc(as.numeric(nolot.hygl.into.vero))))[,2])])
1-auc(0)
auc <- ecdf(nolot.vero.into.hygl[which(nolot.vero.into.hygl<HPDinterval((as.mcmc(as.numeric(nolot.vero.into.hygl))))[,2])])
1-auc(0)

auc <- ecdf(nolot.vero.into.hygl)
1-auc(0)

# medi and peai 
means <- c(mean.medi.into.peai, mean.peai.into.medi)
plot(density(nolot.medi.into.peai, na.rm = T), col=cols[8], main ="", lwd=2,
     xlim=c(-0.5,5), ylim=c(0,8), xlab="")
lines(density(nolot.peai.into.medi, na.rm = T), col=cols[9], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(8,9)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("M. minima")), expression(italic("P. airoides"))), col=c(cols[8], cols[9]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.medi.into.peai[which(nolot.medi.into.peai<HPDinterval((as.mcmc(as.numeric(nolot.medi.into.peai))))[,2])])
1-auc(0)
auc <- ecdf(nolot.peai.into.medi[which(nolot.peai.into.medi<HPDinterval((as.mcmc(as.numeric(nolot.peai.into.medi))))[,2])])
1-auc(0)

# peai and trcy
means <- c(mean.peai.into.trcy, mean.trcy.into.peai)
plot(density(nolot.peai.into.trcy, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1.5,5.5), ylim=c(0,6), xlab="")
lines(density(nolot.trcy.into.peai, na.rm = T), col=cols[5], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,5)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("T. cyanopetala"))), col=c(cols[9], cols[5]),lwd=2, lty=1, bty="n")

# peai and vero
means <- c(mean.peai.into.vero, mean.vero.into.peai)
plot(density(nolot.peai.into.vero, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,5.5), ylim=c(0,8), xlab="")
lines(density(nolot.vero.into.peai, na.rm = T), col=cols[6], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,6)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("G. rosea"))), col=c(cols[9], cols[6]),lwd=2, lty=1, bty="n")

auc <- ecdf(nolot.peai.into.vero[which(nolot.peai.into.vero<HPDinterval((as.mcmc(as.numeric(nolot.peai.into.vero))))[,2])])
1-auc(0)
auc <- ecdf(nolot.vero.into.peai[which(nolot.vero.into.peai<HPDinterval((as.mcmc(as.numeric(nolot.vero.into.peai))))[,2])])
1-auc(0)

auc <- ecdf(nolot.peai.into.vero)
1-auc(0)
auc <- ecdf(nolot.vero.into.peai)
1-auc(0)

# poca and trcy 
means <- c(mean.poca.into.trcy, mean.trcy.into.poca)
plot(density(nolot.poca.into.trcy, na.rm = T), col=cols[4], main ="", lwd=2,
     xlim=c(-1.5,6), ylim=c(0,5.5), xlab="")
lines(density(nolot.trcy.into.poca, na.rm = T), col=cols[5], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(4,5)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. canescens")), expression(italic("T. cyanopetala"))), col=c(cols[4], cols[5]),lwd=2, lty=1, bty="n")

# poca and vero 
means <- c(mean.poca.into.vero, mean.vero.into.poca)
plot(density(nolot.poca.into.vero, na.rm = T), col=cols[4], main ="", lwd=2,
     xlim=c(-1,6), ylim=c(0,6), xlab="")
lines(density(nolot.vero.into.poca, na.rm = T), col=cols[6], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(4,6)], pch=19)
mtext("case 3", side=3, adj=0)
legend("topright", legend=c(expression(italic("P. canescens")), expression(italic("G. rosea"))), col=c(cols[4], cols[6]),lwd=2, lty=1, bty="n")
mtext(side=1, line=3, adj=4, "Low Density Growth Rate")

auc <- ecdf(nolot.vero.into.poca[which(nolot.vero.into.poca<HPDinterval((as.mcmc(as.numeric(nolot.vero.into.poca))))[,2])]) # change to >
1-auc(0)
auc <- ecdf(nolot.poca.into.vero[which(nolot.poca.into.vero<HPDinterval((as.mcmc(as.numeric(nolot.poca.into.vero))))[,2])])
1-auc(0)

auc <- ecdf(nolot.vero.into.poca)
1-auc(0)
auc <- ecdf(nolot.poca.into.vero)
1-auc(0)

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
one.comp.ex <- 10/14
neither.coex <- 0
dist.dep <- 3/14
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
means <- c(mean.peai.into.poca, mean.poca.into.peai)
plot(density(nolot.peai.into.poca, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,7), ylim=c(0,1), xlab="")
lines(density(nolot.poca.into.peai, na.rm = T), col=cols[4], main ="", lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,4)], pch=19)
#species.names = c("dagl", "hygl", "plde", "poca", "trcy", "vero", "arca", "medi", "peai")
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. canescens"))), col=c(cols[9], cols[4]),lwd=2, lty=1, bty="n")
mtext(side=3, adj=0, "(a) case 2", cex=1.3)
mtext(side = 2, line=2, "Density", cex=1.3)

# plde and poca
means <- c(mean.plde.into.poca, mean.poca.into.plde)
plot(density(nolot.plde.into.poca, na.rm = T), col=cols[3], main ="", lwd=2,
     xlim=c(-1,5), ylim=c(0,8), xlab="")
lines(density(nolot.poca.into.plde, na.rm = T), col=cols[4], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(3,4)], pch=19)
mtext(side=3, adj=0, "(b) case 3", cex=1.3)
legend("topright", legend=c(expression(italic("P. debilis")), expression(italic("P. canescens"))), col=c(cols[3], cols[4]),lwd=2, lty=1, bty="n")

#hygl and plde
means <- c(mean.hygl.into.plde, mean.plde.into.hygl)
plot(density(nolot.hygl.into.plde, na.rm = T), col=cols[2], main ="", lwd=2,
     xlim=c(-1,3.5), ylim=c(0,3.5), xlab="", ylab="")
lines(density(nolot.plde.into.hygl, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(2,3)], pch=19)
legend("topright", legend=c(expression(italic("H. glutinosum")), expression(italic("P. debilis"))), col=c(cols[2], cols[3]),lwd=2, lty=1, bty="n")
mtext("Low density growth rate", side=1, line=2, cex=1.3)
mtext(side=3, adj=0, "(c) case 4", cex=1.3)
mtext(side = 2, line=2, "Density", cex=1.3)

# peai and plde
means <- c(mean.peai.into.plde, mean.plde.into.peai)
plot(density(nolot.peai.into.plde, na.rm = T), col=cols[9], main ="", lwd=2,
     xlim=c(-1,4), ylim=c(0,3), xlab="", ylab="")
lines(density(nolot.plde.into.peai, na.rm = T), col=cols[3], lwd=2)
abline(v=0, lty=2)
points(x=means, y=rep(0,2),  col=cols[c(9,3)], pch=19)
mtext("Low density growth rate", side=1, line=2, cex=1.3)
legend("topright", legend=c(expression(italic("P. airoides")), expression(italic("P. debilis"))), col=c(cols[9], cols[3]),lwd=2, lty=1, bty="n")
mtext(side=3, adj=0, "(d) case 6", cex=1.3)

# case2 <- 1/14
# case3 <- 10/14
# case4 <- 1/14
# case6 <- 2/14
# cases <- cbind(case2, case3, case4, case6)
# a <- barplot(cases, ylab = "", xaxt="n", ylim=c(0,1))
# axis(1, at=a, tck=-0.01, labels=c("case \n 2", "case \n 3", "case \n 4", "case\n  6"), cex.axis=1)
# mtext(side=3, adj=0, "(c) distribution dependent", cex=1.3)
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
plot(density(nolot.peai.into.poca-10, na.rm = T, bw = 0.8), col="black", lwd=3, yaxt="n", xaxt="n",
     xlim=c(-12,1), ylim=c(0,.4), main="CASE 1", cex=1.3)
lines(density(nolot.peai.into.poca-7, na.rm = T, bw = 0.8), col="grey", main ="", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)
mtext("Density", side=2)

# case 2
plot(density(nolot.peai.into.poca+5, na.rm = T, bw = 0.8), col="black", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-1,15), ylim=c(0,.4), main="CASE 2", cex=1.3)
lines(density(nolot.peai.into.poca+8, na.rm = T, bw = 0.8), col="grey", main ="", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)

# case 3
plot(density(nolot.peai.into.poca-9, na.rm = T, bw = 0.8), col="black", main ="CASE 3", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-12,12), ylim=c(0,.4))
lines(density(nolot.peai.into.poca+5, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)

# case 4
plot(density(nolot.peai.into.poca-2, na.rm = T, bw = 0.8), col="black", main ="CASE 4", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-5,12), ylim=c(0,.4))
lines(density(nolot.peai.into.poca+5, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
mtext("LDGR", side=1, line=2)
mtext("Density", side=2)
axis(side = 1, 0)

# case 5
plot(density(nolot.peai.into.poca+5, na.rm = T, bw = 0.8), col="black", main ="CASE 5", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-5,12), ylim=c(0,.4))
lines(density(nolot.peai.into.poca-2, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)
mtext("LDGR", side=1, line=2)

# case 6
plot(density(nolot.peai.into.poca-1, na.rm = T, bw = 0.8), col="black", main ="CASE 6", lwd=3, yaxt="n",xaxt="n",
     xlim=c(-12,12), ylim=c(0,.4))
lines(density(nolot.peai.into.poca-2, na.rm = T, bw = 0.8), col="grey", lwd=3)
abline(v=0, lty=2, lwd=3)
axis(side = 1, 0)
mtext("LDGR", side=1, line=2)

dev.off()


