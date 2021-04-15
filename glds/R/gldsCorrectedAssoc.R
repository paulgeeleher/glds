#' Calculate assocations GLDS corrected assocations between drug sensitivity data and some predictive marker (e.g. somatic mutations)
#' 
#' Calculate assocations GLDS corrected assocations between drug sensitivity data and some predictive marker (e.g. somatic mutations)
#' 
#' @param drugMat Matrix of drug sensitivity data, with drug names as column names and samples (e.g. cell lines) as row names.
#' @param drugRelationshipList A list, which for each drug contains a vector of 1's and 0's, indicating whether the drugs are related. I.e. if both drugs were an inhibitor of ERBB2, that position would contain a 1. This order of the drugs should be the same as the order of the drug sensitivty matrix, i.e. drugMat.
#' @param markerMat A matrix containing the data for which you are looking for an association with drug sensitivity, i.e. a matrix of somatic mutation status. This matrix has samples names (i.e. cell lines) as columns and marker names (e.g. gene names) as row names.
#' @param numCorDrugsExclude When calculating GLDS, one of the steps involves searching for drugs with a high correlation with the drug on interest and excluding them from the calculation of GLDS (so as to not remove drug specific signal). As a rule of thumb, about 25% of the size of the dataset may be suitable.
#' @param minMuts For each marker, i.e. row of markerMat, what is the minimum number of non-zero enteries required so that we calculate a P-value. E.g. how many somatic mutations need to be present.
#' @param additionalCovariateMatrix A matrix containing additional covariates to be fit in the drug biomarker association models. This could be, for example, tissue of origin or cancer type. Columns should be sample names (i.e. cell lines) and rows should be the row names.
#' 
#' @export
gldsCorrectedAssoc <- function(drugMat, drugRelationshipList, markerMat, numCorDrugsExclude=100, minMuts=5, additionalCovariateMatrix=NULL)
{
  results_gldsPs <- list()
  results_gldsBetas <- list()
  results_naivePs <- list()
  results_naiveBetas <- list()  
  numDrugs <- ncol(drugMat)
  pairCor <- cor(drugMat, method="spearman")
  comNames <- colnames(markerMat)[colnames(markerMat) %in% rownames(drugMat)] # cell lines for which we have both mutation and drug data....
  
  if(!is.null(additionalCovariateMatrix)) # if additional co variate matrix was provided, then also subset to those samples
  {
    comNames <- comNames[comNames %in% rownames(additionalCovariateMatrix)]
  }
  
  pb <- txtProgressBar(min = 0, max = ncol(drugMat), style = 3) # create progress bar
  for(i in 1:ncol(drugMat))
  {
    # Calculate 10 PCs on non-related sets of drugs....
    negControlDrugs <- colnames(drugMat)[!as.logical(drugRelationshipList[[i]])]
    pairwiseCorNear <- names(sort(abs(pairCor[, colnames(drugMat)[i]]))[(numDrugs-numCorDrugsExclude):numDrugs]) # NB also remove very correlated drugs. Number defined by "numCorDrugsExclude".
    negControlDrugs <- setdiff(negControlDrugs, pairwiseCorNear) # remove very highly correlated drugs from "negative controls"
    controlPCsAll <- prcomp(drugMat[, negControlDrugs])$x
    controlPCsAllCom <- controlPCsAll[comNames, ]
    
    # Calculate the P-values and beta values for each marker for this drug, controlling for GLDS and not controlling for GLDS
    results_gldsPs[[i]] <- numeric()
    results_gldsPs[[i]] <- rep(NA, nrow(markerMat))
    results_gldsBetas[[i]] <- numeric()
    results_gldsBetas[[i]] <- rep(NA, nrow(markerMat))
    results_naivePs[[i]] <- numeric()
    results_naivePs[[i]] <- rep(NA, nrow(markerMat))
    results_naiveBetas[[i]] <- numeric()
    results_naiveBetas[[i]] <- rep(NA, nrow(markerMat))
    for(j in 1:nrow(markerMat))
    {
      if(sum(markerMat[j, comNames]) > minMuts)
      {
	if(is.null(additionalCovariateMatrix)) # if no additional covariate have been provided....
	{
	  theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames])))
	  results_naivePs[[i]][j] <- theCoefs[2, 4]
	  results_naiveBetas[[i]][j] <- theCoefs[2, 1]
	  
	  theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+controlPCsAllCom[,1:10])))
	  results_gldsPs[[i]][j] <- theCoefs[2, 4]
	  results_gldsBetas[[i]][j] <- theCoefs[2, 1]
	}
	else # if there are other covariates, include them in the models.
	{
	  theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+additionalCovariateMatrix[comNames,])))
	  results_naivePs[[i]][j] <- theCoefs[2, 4]
	  results_naiveBetas[[i]][j] <- theCoefs[2, 1]
	  
	  theCoefs <- coef(summary(lm(drugMat[comNames, i]~markerMat[j, comNames]+controlPCsAllCom[,1:10]+additionalCovariateMatrix[comNames,])))
	  results_gldsPs[[i]][j] <- theCoefs[2, 4]
	  results_gldsBetas[[i]][j] <- theCoefs[2, 1]
	}
      }
    }
    # cat(paste(i, " ", sep=""))
    names(results_gldsPs[[i]]) <- rownames(markerMat)
    names(results_naivePs[[i]]) <- rownames(markerMat)
    names(results_gldsBetas[[i]]) <- rownames(markerMat)
    names(results_naiveBetas[[i]]) <- rownames(markerMat)
    
    setTxtProgressBar(pb, i) # update progress bar
    
  }
  close(pb)
  names(results_gldsPs) <- colnames(drugMat)
  names(results_naivePs) <- colnames(drugMat)
  names(results_gldsBetas) <- colnames(drugMat)
  names(results_naiveBetas) <- colnames(drugMat)
  
  outList <- list(pGlds=results_gldsPs, betaGlds=results_gldsBetas, pNaive=results_naivePs, betaNaive=results_naiveBetas)
  return(outList)  
}
