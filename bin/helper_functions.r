library(GEOquery, quietly = TRUE)
library(arrayQualityMetrics, quietly = TRUE)
library(impute, quietly = TRUE)
library(limma, quietly = TRUE)
library(gplots)
require(Biobase)
library(parmigene)
library(WGCNA)

build_coexp_network <- function(data.all, genes.subset){

        ### data.all should contain: info (sample information), genes (entrez gene IDs), data.matrix ( gene expressions )
        info  = data.all[[1]]
        genes = data.all[[2]]
        data  = data.all[[3]]

        # If gene subset exists, filter expression table on subset
        if( ! missing(genes.subset) ){
                # Make filters
                m = match(genes.subset[,1], genes)
                f.g = !is.na(m)
                f.d = m[f.g]

                data = data[f.d,]
                genes = genes[f.d]

        }

        # Get correlations and rank then standardize
        network = cor( t(data), method="s")

        # Fill in missing genes
        if( ! missing(genes.subset) ){
                n = dim(genes.subset)[1]
                network.temp = matrix( NA, ncol = n, nrow= n )
                rownames(network.temp) = genes.subset[,1]
                colnames(network.temp) = genes.subset[,1]

                genes = genes.subset[,1] 
                
                network.temp[f.g,f.g] = network
                network = network.temp
                diag(network)  = 1

        }

        # Number of genes
        n = dim(network)[1]

        network.rank = matrix(rank(network, na.last="keep",ties.method="average"), nrow=n, ncol=n)
        network.rank = network.rank/max(network.rank, na.rm=T)

        rownames(network.rank) = genes
        colnames(network.rank) = genes

        network.rank[is.na(network.rank)] = 1
        diag(network.rank) = 1

        return(network.rank)
}

get_expression_matrix_from_GEO <- function(gseid){

        # Get the family_soft.gz file for the GSEid
	gseSOFT <- getGEO(GEO=gseid, GSEMatrix=F)

        # Get the properties of the microarray samples
        descrip = Columns(GSMList(gseSOFT)[[1]])$Description[2]
        names = names(GSMList(gseSOFT))
        platforms = lapply(GSMList(gseSOFT),function(x) {Meta(x)$platform})
        gplid = unlist(unique(platforms))[1]  # take the first platform, sometimes mutliple platforms, this will fail s

        # Get more sample info
	pD1 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$title)
        pD2 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$description)
        pD3 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$source_name_ch1)
        pD4 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$type)
        pD5 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$characteristics_ch1)

        # Parse the sample info
        pD1 <- get_charcs(pD1,names)
        pD2 <- get_charcs(pD2,names)
	pD3 <- get_charcs(pD3,names)
	pD4 <- get_charcs(pD4,names)
	pD5 <- get_charcs(pD5,names)

	info <- cbind(names,pD1,pD2,pD3,pD4,pD5)

        # Get the probeset and corresponding ORFs
        ### This is not consistent across platforms (ie the gplids), may cause issues
        probesets <- Table(GPLList(gseSOFT)[[gplid]])$ID
	ORF <- Table(GPLList(gseSOFT)[[gplid]])$ENTREZ_GENE_ID

        # For each sample, get the expression levels for the probes
        data.matrix <- do.call('cbind',lapply(GSMList(gseSOFT),function(x)
                                      {tab <- Table(x)
                                       mymatch <- match(probesets,tab$ID_REF)
                                       return(tab$VALUE[mymatch])
                                     }))

        # Sometimes the data is not a numeric variable, convert
	data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

        # Label rows and columns
	rownames(data.matrix) <- probesets
	colnames(data.matrix) <- names(GSMList(gseSOFT))

	# Use probes with gene ids, filter out probes with multiple ORFs and probes with no ORFs
	use_probe <- which( is.na(ORF) == F & match(ORF, "", nomatch = 0) == 0 & regexpr("///", ORF) == -1 )
	data.matrix <- data.matrix[use_probe, ]
	rownames(data.matrix) <- ORF[use_probe]

	# Check to see if data is log2 transformed
	Med <- median(data.matrix, na.rm = T)
	if (Med > 16) data.matrix <- log2(data.matrix)


	# Normalize between arrays so that the intensities or log-ratios have similar distributions
	na.length <- length(which(is.na(data.matrix) == T))
	if (na.length > 0) data.matrix <- impute.knn(data.matrix)$data
	data.matrix <- normalizeBetweenArrays(data.matrix)

	# Aggregate multiple genes, using median expression value
	tmp <- aggregate(data.matrix, list(rownames(data.matrix)), median)
	data.matrix <- as.matrix(tmp[, -1])
	genes = tmp[,1]
        rownames(data.matrix) <- genes

	# Clean up
        rm(tmp,na.length,use_probe,ORF)


        ### Write to text files
        # write( t(info), ncol=dim(final)[2], file=paste(GEOout,"info",sep="."), sep="\t")
        # write.table(data.matrix, file=paste(GEOout,"data",sep="."))


	return( list(info, genes, data.matrix) )
}

get_expression_matrix_from_file <- function(gseid){

        dir = "/home/sballouz/data/RNAseq/gemma_rnaseq_data/data_new/"
        # Get the FPKM file for the GSEid
	data.matrix <- read.table( paste(dir,gseid, "_expression_FPKM.parsed", sep="") )

        genes = data.matrix[,1] 
        data.matrix = data.matrix[,-1]

        # Get the meta data from GEO
        gseSOFT <- getGEO(GEO=gseid, GSEMatrix=F)

        # Get the properties of the RNAseq samples
        descrip = Columns(GSMList(gseSOFT)[[1]])$Description[2]
        names = names(GSMList(gseSOFT))
        platforms = lapply(GSMList(gseSOFT),function(x) {Meta(x)$platform})
        gplid = unlist(unique(platforms))[1]  # take the first platform, sometimes mutliple platforms, this will fail s

        # Get more sample info
	pD1 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$title)[[1]]
        pD2 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$description)[[1]]
        pD3 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$source_name_ch1)[[1]]
        pD4 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$type)[[1]]
        pD5 <- sapply(names, function(gsm)GSMList(gseSOFT)[[gsm]]@header$characteristics_ch1)[[1]]

        # Parse the sample info
        pD1 <- get_charcs(pD1,names)
        pD2 <- get_charcs(pD2,names)
	pD3 <- get_charcs(pD3,names)
	pD4 <- get_charcs(pD4,names)
	pD5 <- get_charcs(pD5,names)

	info <- cbind(names,pD1,pD2,pD3,pD4,pD5)


	return( list(info, genes, data.matrix) )
}

get_expression_matrix_from_gemma <- function(gseid){

        # Get the text file from gemma (via url)
        filtered = "true"
        #filtered = "false"
        url = paste("http://chibi.ubc.ca/Gemma/rest/experimentData/findExpressionDataByEeName/", gseid, ",", filtered, sep="")

	data.matrix <- read.table( url, sep="\t", header=T, quote="")

        write.table(data.matrix, file=paste(gseid,".data", sep=""))

        # Entrez gene IDs column 6
        genes = data.matrix[,6]
        use_genes <- which( is.na(genes) == F & match(genes, "", nomatch = 0) == 0 & regexpr("\\|", genes) == -1 )
        genes = genes[use_genes]

        # Remove extra ID columns (1:6) and missing genes
        S = dim(data.matrix)[2]
        data.matrix = data.matrix[use_genes,7:S]

        # Remove and aggregate duplicate genes
        tmp <- aggregate(data.matrix, list(genes), median)
	data.matrix <- as.matrix(tmp[, -1])
	genes = tmp[,1]
        rownames(data.matrix) <- genes

        # Sample info in column names
        info = colnames(data.matrix)

        return( list(info, genes, data.matrix) )
}

save_expression_matrix_from_gemma <- function(gseid){

        # Get the text file from gemma (via url)
        filtered = "true"
        #filtered = "false"
        url = paste("http://chibi.ubc.ca/Gemma/rest/experimentData/findExpressionDataByEeName/", gseid, ",", filtered, sep="")

	data.matrix <- read.table( url, sep="\t", header=T, quote="")

        write.table(data.matrix, file=paste(gseid,".data", sep=""), sep="\t", quote=F, row.names=F)
}

qc_expression_matrix_from_gemma <- function(gseid){

	data.matrix <- read.table(paste(gseid,".data", sep=""), quote="", header=T, sep="\t")
        genes = data.matrix[,6]
        use_genes <- which( is.na(genes) == F & match(genes, "", nomatch = 0) == 0 & regexpr("\\|", genes) == -1 )
        genes = genes[use_genes]

        # Remove extra ID columns (1:6) and missing genes
        S = dim(data.matrix)[2]
        data.matrix = data.matrix[use_genes,7:S]

        # Remove and aggregate duplicate genes
        tmp <- aggregate(data.matrix, list(genes), median)
	data.matrix <- as.matrix(tmp[, -1])
	genes = tmp[,1]
        rownames(data.matrix) <- genes

        # Sample info in column names
        info = colnames(data.matrix)


        sample.cor = cor(data.matrix, method="s")
        png(paste("./imgs/",gseid,".sample.cor.png", sep="") )
        heatmap.2(sample.cor, density="none", trace="none", main=gseid)
        dev.off()
        

        return( list(info, genes, data.matrix) )

}

# Transparent colors
makeTransparent<-function(someColor, alpha=100)
{
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
	blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}





get_charcs <- function(pD,names){

       	if( class(pD) == "matrix")
	       if( dim(pD)[1] == length(names)) return(pD)
	       else return(t(pD))
        if( class(pD) == "character") return(pD)
	if( class(pD) == "list"){
	        max = 0
	       	for (i in names){
		        t <- length(pD[[i]])
			if( t > max ) max = t
		}

		nrows = length(names)
		ncols = max


		char <- matrix("-", ncol=ncols, nrow=nrows)
		rownames(char) = names

		for (i in names){
			if( length(pD[[i]]) < max ){
				missing = max - length(pD[[i]])
			} else {
				missing = 0
			}
			char[i,] <- c(unlist(pD[[i]]), rep("-", missing))
		}
	}

	return(char)
}


make_network <- function(exprs){

        n = dim(exprs)[1]
        s = dim(exprs)[2]
        genes = rownames(exprs)
        samples = colnames(exprs)

        network = cor( t(exprs), method="s")

        network.rank = matrix(rank(network, na.last="keep",ties.method="average"), nrow=n, ncol=n)
        network.rank = network.rank/max(network.rank, na.rm=T)

        rownames(network.rank) = genes
        colnames(network.rank) = genes

        network.rank[is.na(network.rank)] = 0
        network = network.rank 
        rm(network.rank)
        return(network)

}

node_degree <- function(network){
        return( colSums(network, na.rm=T) )
}


threshold_network_top_genes <- function(network, p){
	# threshold final network to top 0.1 genes
	diag(network) <- 0
	n = dim(network)[1]
	x = n-1
        top = ceiling(x*p)

	ord = order(network, decreasing=T)
	ord2 = apply(network, 2, order, decreasing=T)
	ranks <-  apply( network, 2, rank,na.last="keep",ties.method="first")
	x = which(ranks < (n-top))
	network[x]=0

	a = network != 0
	b =  t(network) != 0
	c = a+b
	d = c==0
	c[d] = 1

	network = (network + t(network))/c
	return(network)
}


make_coexp_network_thresh <- function(indir, id, exprs, n, gene.ids, method,p){

        network = cor( t(exprs), method=method)
        diag(network) = 0

        network[is.na(network)] = 0
        network = abs(network)
        temp = sort(network, method="q")
        i = temp[length(temp)*p]
        network[network < i] = 0
        network[network>0] = 1

        rownames(network) = gene.ids
        colnames(network) = gene.ids

        save(network,exprs, file=paste( indir,"/coexp/thresh/",id,".", method, ".Rdata",sep=""))
}



make_coexp_network_MI <- function(data.matrix,type,out){
        # Run Mutual Information
        print (system.time(mi <- knnmi.all(data.matrix)))

	# Run ARACNE-a
        print(system.time(grn.a.a <- aracne.a(mi, 0.05)))

        # Run ARACNE-m
        print(system.time(grn.a.m <- aracne.m(mi,tau=0.15)))

        # Run CLR
        print(system.time(grn.clr <- clr(mi)))

        # Run MRNET
        print(system.time(grn.mrnet <- mrnet(mi)))

        save(mi, grn.a.m, grn.a.a, grn.clr, grn.mrnet, file=paste(out,type, "Rdata", sep="."))
}

make_coexp_network_WGCNA <- function(data.matrix,type,out){

        names = rownames(data.matrix)
        #
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(t(data.matrix), powerVector = powers, verbose = 5)
        softPower = sft$powerEstimate

        if( is.na(softPower)){ softPower = 6 }
        adjacency = adjacency( t(data.matrix), power=softPower)

        test = is.na(adjacency)
        adjacency[test] = 0

        network <- TOMsimilarity(adjacency)
        rownames(network) <- names
        colnames(network) <- names
        

        save(network, file=paste(out,type,softPower, "Rdata", sep="."))
}


