a <- read.table("homosapiens_feliscatus_ortho.txt", header=F)
a[1,]
plot(a[,3], a[,4])
lmmod = lm(a[,4]~a[,3])
summary(lmmod)
lmmod$coefficients
abline(a=lmmod$coefficients[1], b=lmmod$coefficients[2], col='red', lwd=2)
