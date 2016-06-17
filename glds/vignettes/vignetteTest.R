# Prototype the code for the vignette here.


# Load the package
library(glds)
library(glmnet)
# source("../R/gldsCorrectedAssoc.R")
# source("../R/completeMatrixFunction.R")

# Load the CTRP drug data....
# load(file="/home/pgeeleher/Dropbox/mdrPaper/glds_package/glds/data/aucMat.RData") # lVecs, markerMat, aucMat, fullAucMat
data(aucMat)

# Impute the CTRP drug data. This data is also a loaded above. Uncomment this line to re-run this algorithm. Note that the algorithm contains a stochastic (i.e. random) component, therefore the results will not be identical when the algorithm is re-run.
# fullAucMat <- completeMatrix(aucMat, nPerms=3)

# to save time, lets just look at a subset of the mutations, these will include the top gene-drug associations.
mutationSubset <- c("BRAF", "TP53", "NRAS", "PIK3CA", "KRAS", "FGF4", "LASP1", "FNTA", "SSX1", "IL2", "SDHC")

# Now run the GLDS corrected P-values for this subset.
assocsOut <- gldsCorrectedAssoc(fullAucMat, lVecs, markerMat[mutationSubset,], numCorDrugsExclude=100, minMuts=5, additionalCovariateMatrix=NULL)

# assocsOut will be a list with 4 elements "pGlds", "betaGlds", "pNaive", "betaNaive", which contains the p-values and beta values (i.e. effect size) for each gene-drug combination for the uncorrected (naive) and GLDS corrected approaches.
# lets arrange these results in a table ordered by pGLDS
theOrd <- names(sort(unlist(assocsOut$pGlds)))
dfOut <- data.frame(pGlds=unlist(assocsOut$pGlds)[theOrd], betaGlds=unlist(assocsOut$betaGlds)[theOrd], pUncorr=unlist(assocsOut$pNaive)[theOrd], betaUncorr=unlist(assocsOut$betaNaive)[theOrd])

# print the top 10 associations, ordered by GLDS p-value
print(dfOut[1:10,])


