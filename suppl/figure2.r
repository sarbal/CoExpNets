load("fig2_data.Rdata")
rocs.goexp.means  = colMeans(rocs.go)
rocs.goterm.means = rowMeans(rocs.go)
rocs.goexp.sd = apply(rocs.go, 2, sd)
rocs.goterm.sd = apply(rocs.go,1, sd)

rocs.go.test = rocs.go > agg.go[,2]
frac.goexp = colSums(rocs.go.test)/dim(rocs.go)[1]
frac.goterm = rowSums(rocs.go.test)/dim(rocs.go)[2]

h1a = get_density(hist( rocs.goexp.means, plot=F))
h2a = get_density(hist( rocs.goterm.means, plot=F))
h1 = density( rocs.goexp.means)
h2 = density( rocs.goterm.means)


means.filt = (rocs.goterm.means > agg.go[,2])
#sum(means.filt)
#colMeans( cbind(rocs.goterm.means[means.filt], agg.go[means.filt,2]) )
#apply( cbind(rocs.goterm.means[means.filt], agg.go[means.filt,2]), 2 ,sd )
max.gos = apply(rocs.go, 1, max)

# Panel A
rtest =  wilcox.test( rocs.goexp.means, rocs.goterm.means)
pval = paste("P = ",format(rtest$p.value))
h1 = get_density( hist(rocs.goexp.means,breaks=10,plot=F))
h2 = get_density( hist(rocs.goterm.means,breaks=20,plot=F))
plot(h2, xlab="GO AUROCs", ylab="Distribution density", main="", type="l", lwd=4, xlim=c(0.45,0.95), ylim=c(0,15), bty="n")
lines(h1, lwd=2, col=makeTransparent(1))
abline(v= mean(rocs.goexp.means))
abline(v= mean(rocs.goterm.means))
legend("topright", legend=c("Experiments","GO terms"), col=c(1,makeTransparent(1)), lwd=c(4,2), bty="n")
legend("right", legend=pval, pch=19)



# Panel B
plot( sort(max.gos), (length(max.gos)-1:length(max.gos))/length(max.gos), type="l", lwd=2, bty="n", xlab="GO AUROCs", ylab="" )
abline(v=min(max.gos) )
text(min(max.gos)+0.05, 0, round(min(max.gos),3)  )


