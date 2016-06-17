#' Impute missing values in a matrix using an expectation maximization approach based on lasso regression models.
#' 
#' Impute missing values in a matrix using an expectation maximization approach based on lasso regression models.
#' 
#' @import glmnet
#'
#' @param allIc50 A matrix of drug sensitivity data with missing values, 
#' @param nPerms The number of iterations that the EM-algorithm should run.
#' 
#' @export
completeMatrix <- function(allIc50, nPerms=50)
{
  cat("\nNumber of iterations:")
  #' To initialize the algorithm, all missing values are first imputed to the median.
  numCellLinesNotScreened <- apply(allIc50, 2, function(r)return(sum(is.na(r))))
  numDrugsNotScreened <- apply(allIc50, 1, function(r)return(sum(is.na(r))))
  indexNotScreened <- apply(allIc50, 2, function(r)return(which(is.na(r))))
  drugMed <- apply(allIc50, 2, function(col)return(median(na.omit(as.numeric(col)))))
  hundIc50sImpute <- data.matrix(allIc50)
  for(i in 1:ncol(allIc50))
  {
    datImpute <- as.numeric(allIc50[,i])
    datImpute[is.na(datImpute)] <- drugMed[i]
    hundIc50sImpute[,i] <- datImpute
  }

  #' Sort the matrix by least number of missing values to most...
  hundIc50sImputeSort <- hundIc50sImpute[, order(numCellLinesNotScreened[colnames(hundIc50sImpute)])]

  #' Calculate PCs of this matrix.
  pcs <- prcomp(hundIc50sImputeSort)$x

  #' Now fit a lasso model for all other drugs, for all samples other than those NOT screened by the drug we are predicting for. We will iterate 50 times.
  imputeSortList <- list()
  imputeSortList[[1]] <- hundIc50sImputeSort
  medianDistance <- numeric()
  for(j in 1:nPerms)
  {
    pb <- txtProgressBar(min = 0, max = ncol(hundIc50sImputeSort), style = 3) # create progress bar
    for(i in 1:ncol(hundIc50sImputeSort))
    {
      setTxtProgressBar(pb, i) # update progress bar
      pcs <- prcomp(hundIc50sImputeSort)$x # calcualte the PCs of the current matrix
      indexNotScreenedThisDrug <- indexNotScreened[[colnames(hundIc50sImputeSort)[i]]] # index of "imputed" cell lines for this drug
      
      if(length(indexNotScreenedThisDrug) > 0)
      {
	# make the training (design matrix) for the non-imputed cell lines for this drug
	pcMat <- pcs[-indexNotScreenedThisDrug, ] 
	dmRaw <-  model.matrix(~pcMat)

	# calcualte lambda using CV
	FitRaw <- cv.glmnet(dmRaw, hundIc50sImputeSort[-indexNotScreenedThisDrug, i], nfolds=20)
	getLambdasRaw <- FitRaw$lambda.min # 0.1132716
	
	# fit the model and extract the co-efficients
	theModRaw <- glmnet(dmRaw, hundIc50sImputeSort[-indexNotScreenedThisDrug, i], lambda=getLambdasRaw)
	coef(theModRaw)[,1][coef(theModRaw)[,1] != 0] # 
	betas <- coef(theModRaw)[,1][coef(theModRaw)[,1] != 0][-1]
	intercept <- coef(theModRaw)[,1][coef(theModRaw)[,1] != 0][1]
	names(betas) <- substring(names(betas), 6)

	# use the model to update the estimates for the imputed values for this drug.
	if(length(indexNotScreenedThisDrug) == 1)
	{
	  predictions <- intercept + apply((t(pcs[indexNotScreenedThisDrug, names(betas)]) * betas), 1, sum); # predict new values
	}
	else
	{
	  predictions <- intercept + apply(t(t(pcs[indexNotScreenedThisDrug, names(betas)]) * betas), 1, sum);
	}
	
	hundIc50sImputeSort[indexNotScreenedThisDrug, i] <- predictions  # update matrix
      }

    }
    close(pb)
    imputeSortList[[j + 1]] <- hundIc50sImputeSort
    medianDistance[j] <- mean(imputeSortList[[j]] - imputeSortList[[j + 1]])
    cat(paste("\nIteration: ", j, "\n"), "")
    plot(medianDistance, main=j, ylab="Median Distance from previous imputed matrix")
  }
  cat("\nDone\n")
  return(hundIc50sImputeSort[,colnames(allIc50)])
}