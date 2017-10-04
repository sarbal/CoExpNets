
source("helper_functions.r") 
load("run_GBA.Rdata")
load("GO.human.Rdata")

# Minimum number of smaples needed in expression experiment 
MIN_SAMP = 10 

# List of expression files. Each expression file should have the column names as the sample ids and the first column the gene ids   
exprs_ids = read.table(file="expression_file_names.txt") 

# Genes you want to aggregate over 
agg.genes = unlist(read.table(file="gene_list.txt"))

# Initialize aggregate and other variables  
agg = diag(length(agg.genes)) 
exprs_used = "" 
roc = list() 
i = 1 

# Aggregate across 
for (exprs_file in exprs_ids){
          exprs = read.table(file=exprs_file, header=T)
          # number of samples in experiment 
          n = dim(exprs)[2]
          
          # Only use experiments with more than MIN_SAMP
          if(n <= MIN_SAMP ) { next }
          exprs_used = append(exprs_file,exprs_used)  
          # Get names and labels from file 
          samples = colnames(exprs[,2:n])
          genes = exprs[,1]
	
          m = match( genes, agg.genes)
          f1 = !is.na(m)
          f2 = m[f1]
          	  
          net = build_coexp_network(data.matrix, genes[f1])
	  sub = agg[f2,f2]
	  sub = net + sub 
	  agg[f2,f2] = sub 
	
	  # Rerank agg 
	  agg.r = matrix(rank(agg, na.last="keep",ties.method="average"), nrow=dim(agg)[1], ncol=dim(agg)[2])
	  rownames(agg.r) =genes 
	  colnames(agg.r) =genes 
	
	  # Run GBA 
	  roc[[i]] = run_GBA(sub, GO.labels)
	  i = i + 1 
}

agg.r = matrix(rank(agg , na.last="keep",ties.method="average"), nrow=dim(agg)[1], ncol=dim(agg)[2])
rownames(agg.r) =genes 
colnames(agg.r) =genes 

# Run GBA 
roc_final = run_GBA(agg.r, GO.labels)

# Save binaries 
save(roc, roc_final, file="GBA.Rdata") 
save(sub, agg, genes, exprs_used, file="aggregate.Rdata") 


