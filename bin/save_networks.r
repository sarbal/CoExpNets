source("helper_functions.r")

# Example for GEO data, given GSE id
genes.subset = read.table("entrezGeneID_Human.txt")
gseid = "GSE16615"

# Get expression data from GEO
data.all = get_expression_matrix_from_GEO(gseid)
network.rank = build_coexp_network(data.all, genes.subset)

# Get expression data from RNAseq data file
data.all = get_expression_matrix_from_file(gseid)
network.rank = build_coexp_network(data.all, genes.subset)

# Get expression data gemma
data.all = get_expression_matrix_from_gemma(gseid)
network.rank = build_coexp_network(data.all, genes.subset)


