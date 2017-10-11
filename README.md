# CoExpNets
How-to generate and analyze co-expression aggregate networks


## Basics 
1. Get expression experiment 
``` 
exprs_file = "GSE12946_expression_FPKM.parsed"
exprs = read.table(file=exprs_file, header=T)
```
2. Rank standardize by sample 
``` 
exprs = rank(exprs, na.last="keep")
exprs  exprs/max(exprs, na.rm=T)
```

3. Calculate network edges using spearman correlations 
``` 
load("bin/run_GBA.Rdata")
load("data/GO.human.Rdata")
network = make_network(exprs) 
```

4. Rank standardize network
```
n = dim(network)[1]
network = matrix(rank(network, na.last="keep",ties.method="average"), nrow=n, ncol=n)
```
5. Options: 

  a. Take top 0.1% of connections 
    ```
    network = threshold_network_top_genes(network, 0.001)
    ```
  b. Repeat on other networks and aggregate 
   see ```agg_network.r```
  
  c. As is

6. Test for network functionality 
``` 
GBA = run_GBA(network, GO.labels)
``` 

