#' Bipartite Edge Reconstruction from Expression data
#'
#' This function generates a complete bipartite network from gene expression data and sequence motif data 
#' 
#' @param motif A motif dataset, a data.frame, matrix or exprSet containing 3 columns. Each row describes an motif associated with a transcription factor (column 1) a gene (column 2) and a score (column 3) for the motif.
#' @param expr An expression dataset, as a genes (rows) by samples (columns) data.frame
#' @param verbose logical to indicate printing of output for algorithm progress.
#' @param cpp logical use C++ for maximum speed, set to false if unable to run.
#' @keywords keywords
#' @export
#' @return An object of class "bere"
#' @examples
#' data(yeast)

regpredict <- function(motif.data, 
                      expr.data,
                      verbose=T,
                      randomize="none",
                      cpp=T){
  if(verbose)
    print('Initializing and validating')
  # Create vectors for TF names and Gene names from Motif dataset
  tf.names   <- sort(unique(motif.data[,1]))
  num.TFs    <- length(tf.names)
  if (is.null(expr.data)){
    stop("Error: Expression data null")
  } else {
    # Use the motif data AND the expr data (if provided) for the gene list
    gene.names <- sort(intersect(motif.data[,2],rownames(expr.data)))
    num.genes  <- length(gene.names)
    
    # Filter out the expr genes without motif data
    expr.data <- expr.data[rownames(expr.data) %in% gene.names,]
    
    # Keep everything sorted alphabetically
    expr.data      <- expr.data[order(rownames(expr.data)),]
    num.conditions <- ncol(expr.data);
    if (randomize=='within.gene'){
      expr.data <- t(apply(expr.data, 1, sample))
      if(verbose)
        print("Randomizing by reordering each gene's expression")
    } else if (randomize=='by.genes'){
      rownames(expr.data) <- sample(rownames(expr.data))
      expr.data           <- expr.data[order(rownames(expr.data)),]
      if(verbose)
        print("Randomizing by reordering each gene labels")
    }
  }
  
  # Bad data checking
  if (num.genes==0){
    stop("Error validating data.  No matched genes.\n  Please ensure that gene names in expression file match gene names in motif file.")
  }
  
  strt<-Sys.time()
  if(num.conditions==0) {
    stop("Error: Number of samples = 0")
    gene.coreg <- diag(num.genes)
  } else if(num.conditions<3) {
    stop('Not enough expression conditions detected to calculate correlation.')
  } else {
    if(verbose)
      print('Verified adequate samples, calculating correlation matrix')
    if(cpp){
      # C++ implementation
      gene.coreg <- rcpp_ccorr(t(apply(expr.data, 1, function(x)(x-mean(x))/(sd(x)))))
      
    } else {
      # Standard r correlation calculation
      gene.coreg <- cor(t(expr.data), method="pearson", use="pairwise.complete.obs")      
    }
  }
  
  print(Sys.time()-strt)
  
  if(verbose)
    print('More data cleaning')
  # Convert 3 column format to matrix format
  colnames(motif.data) <- c('TF','GENE','value')
  regulatory.network <- tidyr::spread(motif.data, GENE, value, fill=0)
  rownames(regulatory.network) <- regulatory.network[,1]
  # sort the TFs (rows), and remove redundant first column
  regulatory.network <- regulatory.network[order(rownames(regulatory.network)),-1]
  # sort the genes (columns)
  regulatory.network <- as.matrix(regulatory.network[,order(colnames(regulatory.network))])
  
  # Filter out any motifs that are not in expr dataset (if given)
  if (!is.null(expr.data)){
    regulatory.network <- regulatory.network[,colnames(regulatory.network) %in% gene.names]
  }
  
  # store initial motif network (alphabetized for rows and columns)
#   starting.motifs <- regulatory.network
  

  if(verbose)
    print('Main calculation')
  ########################################

strt<-Sys.time()
  correlation.dif <- sweep(regulatory.network,1,rowSums(regulatory.network),`/`)%*%gene.coreg-sweep(1-regulatory.network,1,rowSums(1-regulatory.network),`/`)%*%gene.coreg
  result <- sweep(correlation.dif, 2, apply(correlation.dif, 2, sd),'/')
#   regulatory.network <- ifelse(res>quantile(res,1-mean(regulatory.network)),1,0)

print(Sys.time()-strt)
  ########################################
  
  return(result)
}
