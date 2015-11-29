load("fig2_data.Rdata")


rocs.go.test = rocs.go > agg.go[,2]
frac.goexp = colSums(rocs.go.test)/dim(rocs.go)[1]
frac.goterm = rowSums(rocs.go.test)/dim(rocs.go)[2]

h1 = get_counts(hist( frac.goexp, plot=F))
h2 = get_counts(hist( frac.goterm, plot=F))

# Panel A
plot(h2, xlab="Fraction of experiments outperforming aggregate", ylab="Distribution counts of GO terms", main="", bty="n", type="l", lwd=4, xlim=c(0,1))
par(new=T)
plot( frac.goterm,  rocs.goterm.means, col="lightgrey", xlab="", ylab="", pch=19, axes=F, bty="n")
axis(4)
mtext("Average GO AUROC across experiments", side=4, padj=3)


# Panel B
plot(h1, main="", xlab="Fraction of experiments outperforming aggregate", ylab="Distribution counts of experiments",type="l", bty="n" )
par(new=T)
plot(frac.goexp,  log10(stats.ids),xlim=range(h1[,1]), pch=19, xlab="", ylab="", axes=F, bty="n")
axis(1)
axis(4, at=(0:3), labels=c(1,10,100,1000))
mtext("Experiment sample size log10", side=4, padj=3)

