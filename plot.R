a <- read.table("homosapiens_feliscatus_ortho.txt", header=F)

lmmod = lm(a[,4]~a[,3])
pdf("correlationPlot.pdf")
plot(a[,3], a[,4], xlab="Neighborhood similarity", ylab="Alignment percentage identity")
abline(a=lmmod$coefficients[1], b=lmmod$coefficients[2], col='red', lwd=2)
dev.off()

summary(lmmod)
lmmod$coefficients
