load("~/shared/bin/run_GBA.Rdata")
load("~/shared/GO.mouse.Rdata")
load("~/shared/GO.human.Rdata")
load("~/shared/coexp/agg.rank.20A.coexp.spearman.Rdata") # sample human network



# Get properties of experiment
N = dim(exprs)[1]
S = dim(exprs)[2] 
genes = rownames(exprs)
samples = colnames(exprs)
avg.exprs = rowMeans(exprs, na.rm=T)
sd.exprs = apply( exprs, 1, sd)


# Generate co-expression network from expression data
network = make_network(exprs)

# Check the node degrees of the network
nd = node_degree(network)
hist(nd)

# Compare node degrees to average expression levels
plot( nd, avg.exprs )

# Run GBA cross validation to get AUROC scores
GBA = run_GBA(network, GO.labels)



