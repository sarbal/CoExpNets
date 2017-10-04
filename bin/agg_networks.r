


MIN_SAMP = 10 
## List of expression files. Each expression file should have the column names as the sample ids and the first column the gene ids   
exprs_ids = read.table(file="expression_file_names.txt") 

    for (exprs_file in exprs_ids){
          exprs = read.table(file=exprs_file, header=T)
          # number of samples in experiment 
          n = dim(exprs)[2]
          
          # Only use experiments with more than MIN_SAMP
          if(n <= MIN_SAMP ) { next }
          
          # Get names and labels from file 
          samples = colnames(exprs[,2:n])
          genes = exprs[,1]
          m = match( attr$ensemblID, genes)
          f1 = !is.na(m)
          f2 = m[f1]
          f3 = !is.na(attr$entrezID[f1])
          data.matrix= exprs[f2,2:n][f3,]
          genes = attr$entrezID[f1][f3]
          filtercols = colSums(data.matrix) != 0
          data.matrix = data.matrix[,filtercols]        
          tmp <- aggregate(data.matrix, list(genes), median)
          data.matrix <- as.matrix(tmp[, -1])
          rownames(data.matrix) <- tmp[, 1]
          colnames(data.matrix) <- samples[filtercols] 
          
          net = make_coexp_network(indir, id, data.matrix, n, rownames(data.matrix), "")


	}

