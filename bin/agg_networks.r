
source("helper_functions.r") 
load("~/shared/bin/run_GBA.Rdata")
load("~/shared/GO.human.Rdata")

# Minimum number of smaples needed in expression experiment 
MIN_SAMP = 10 

# List of expression files. Each expression file should have the column names as the sample ids and the first column the gene ids   
exprs_ids = read.table(file="expression_file_names.txt") 

# Genes you want to aggregate over 
agg.genes = unlist(read.table(file="gene_list.txt"))

# Initialize aggregate 
agg.r = diag(length(agg.genes)) 

# Aggregate across 
for (exprs_file in exprs_ids){
          exprs = read.table(file=exprs_file, header=T)
          # number of samples in experiment 
          n = dim(exprs)[2]
          
          # Only use experiments with more than MIN_SAMP
          if(n <= MIN_SAMP ) { next }
          
          # Get names and labels from file 
          samples = colnames(exprs[,2:n])
          genes = exprs[,1]
	
          m = match( genes, agg.genes)
          f1 = !is.na(m)
          f2 = m[f1]
          	  
          net = build_coexp_network(data.matrix, genes[f1])
	  sub = agg.r[f2,f2]
	  sub = net + sub 
	  agg.r[f2,f2] = sub 
	  # Rerank agg.r 
	  sub = matrix(rank(agg.r, na.last="keep",ties.method="average"), nrow=dim(agg.r)[1], ncol=dim(agg.r)[2])
	 
	  # Run GBA 
	  roc = run_GBA()
}


