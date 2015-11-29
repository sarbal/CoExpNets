load("fig4_data.Rdata")

# Panel A
method= methods[2]
mode = modes[2]
type = types[1]

part = 1:29
filt = final.rna.agg[,5]==type & final.rna.agg[,4] == method & final.rna.agg[,1] == mode
b.all  = (boxplot( as.numeric( final.rna.agg[filt,6]) ~as.numeric(final.rna.agg[filt,3]),  plot=F ))
ma.all = (tapply( as.numeric( final.rna.agg[filt,6]),as.numeric(final.rna.agg[filt,3]), mean))[part]
sd.all = (tapply( as.numeric( final.rna.agg[filt,6]),as.numeric(final.rna.agg[filt,3]),sd))[part]
se.all = (sd.all/ sqrt(b.all$n[part]))
b = as.numeric(b.all$names)[part]

b=b

filtm = final.micr.agg[,5]==type & final.micr.agg[,1] == method & final.micr.agg[,4] == mode
b.allm  = boxplot( as.numeric( final.micr.agg[filtm,6]) ~as.numeric(final.micr.agg[filtm,2]),  plot=F )
mm.all = tapply( as.numeric( final.micr.agg[filtm,6]),as.numeric(final.micr.agg[filtm,2]), mean)
sdm.all = tapply( as.numeric( final.micr.agg[filtm,6]),as.numeric(final.micr.agg[filtm,2]),sd)
sem.all = sdm.all/ sqrt(b.allm$n)
bm = as.numeric(b.allm$names)



plot(b, ma.all,lwd=3,lty=1,type="l", xlim = range(part+1),  ylab="Average AUROCs", xlab="Number of networks in aggregate", bty="n")
segments(b, ma.all-se.all,b, ma.all+se.all, lwd=1, col=makeTransparent(1) )
lines(bm, mm.all, lwd=3, col="grey" )
segments(bm, mm.all-sem.all,bm, mm.all+sem.all, lwd=1, col=makeTransparent("grey") )
legend("bottomright", c("RNA-seq", "Microarray"), lwd=3,col=c(1,"grey"), bty="n" )



# Panel B
m.agg  = round(mean(agg.roc[,2], na.rm=T),2)
m.micr =  round(mean(agg.roc.micr[,2], na.rm=T),2)
m.sub.agg  = round(mean(agg.sub.roc[,2], na.rm=T),2)

h.roc.agg = get_density( hist(agg.roc[,2], breaks=50) )
h.roc.micr = get_density( hist(agg.roc.micr[,2], breaks=50) )
h.roc.sub.agg = get_density( hist(agg.sub.roc[,2], breaks=50) )


rtest =  wilcox.test( agg.sub.roc[,2], agg.roc.micr[,2])

pval = paste("P = ",format(rtest$p.value))



plot( h.roc.micr, xlim=c(0.35,1), col="grey",type="l", lwd=3, xlab="GO AUROCs", ylab="Density", bty="n")
lines(h.roc.sub.agg, lwd=3 )
abline(v=m.sub.agg,lty=2)
abline(v=m.micr, col="grey",lty=2)
legend("topright", c("RNA-seq aggregate", "Microarray aggregate"), col=c(1,"grey"), lwd=3, bty="n")
legend("right", pval, pch=19, bty="n")


# inset
plot( agg.sub.roc[,2], agg.roc.micr[,2], xlab="GO AUROCs in RNA-seq aggregate", ylab="GO AUROCs in microarray aggregate", bty="n", pch=19)
cor.temp = cor.test( agg.roc[,2], agg.roc.micr[,2], method="s")
legend( "bottomright", legend=paste("rho = ", round(cor.temp$estimate,2), "P = ",round(cor.temp$p.value,5)),bty="n")






