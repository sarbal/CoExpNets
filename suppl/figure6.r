load("fig6_data.Rdata")

h1 = get_density(hist( rank(node_degree.sub.rna.rank)[filt.common.poorly]/n,  plot=F))
h2 = get_density(hist( rank(node_degree.sub.micr.rank)[filt.common.poorly]/n,plot=F))
 ymax=max(h1[,2], h2[,2])

# Panel A

# Panel B

# Panel C
