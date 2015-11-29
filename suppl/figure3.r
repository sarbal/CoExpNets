load("fig3_micr_data.Rdata")

# Panel A
plot(x=NULL, y=NULL, xlim=c(1,3), ylim=c(0.5,0.65), xlab="Sample size (log10)", ylab="Average GO AUROCs", sub="Microarray", axes=F)
axis(2)
axis(1, at=(1:3), labels=c(10,100,1000))

x = log10(as.numeric(final[filt,ii]))
y = as.numeric(final[filt,ik])
fi = x != -Inf  & x != Inf
points( x[fi], y[fi], pch=19 )
z = lm( y[fi]~x[fi])
abline( z, lty=1)
cor.temp = cor.test( x[fi], y[fi], method="s")

x = unlist(x.all)
y = unlist(y.all)

z = lm( y~x)
abline( z, lwd=4,lty=2,col=makeTransparent(1) )
cor.temp2 = cor.test( x, y, method="s")

legend( "bottomright", col=c(1, makeTransparent(1)), lty=c(1,2), lwd=c(1,4), legend=c(paste("All experiments, rho = ", round(cor.temp$estimate,3), "P = ",round(cor.temp$p.value,5)),
                                                                                            paste("Downsampled, rho = ", round(cor.temp2$estimate,3), "P = ",round(cor.temp2$p.value,5))), bty="n")

load("fig3_rna_data.Rdata")

# Panel B
plot(x=NULL, y=NULL, xlim=c(1,3), ylim=c(0.5,0.65), xlab="Sample size (log10)", ylab="Average GO AUROCs", sub="RNA-seq", axes=F)
axis(2)
axis(1, at=(1:3), labels=c(10,100,1000))

x = log10(as.numeric(final[filt,ii]))
y = as.numeric(final[filt,ik])
fi = x != -Inf  & x != Inf
points( x[fi], y[fi], pch=19 )
z = lm( y[fi]~x[fi])
abline( z, lty=1)
cor.temp = cor.test( x[fi], y[fi], method="s")


x = unlist(x.all)
y = unlist(y.all)

z = lm( y~x)
abline( z, lwd=4,lty=2,col=makeTransparent(1) )
cor.temp2 = cor.test( x, y, method="s")

legend( "bottomright", col=c(1, makeTransparent(1)), lty=c(1,2), lwd=c(1,4), legend=c(paste("All experiments, rho = ", round(cor.temp$estimate,3), "P = ",round(cor.temp$p.value,5)),
                                                                                            paste("Downsampled, rho = ", round(cor.temp2$estimate,3), "P = ",round(cor.temp2$p.value,5))), bty="n")


