#' The print and summary method for parsimonious class objects. 
#' 
#' @param x An object of class parsimonious.  
#' @param ... unused. 
#' @method print parsimonious
#' @export 
print.parsimonious <- function(x,...){	

	cat('\n		------------------------------		\n')
	cat('Parsimonious time varying parameter models.\n')
	cat('Estimator: ',x$estimator,'.\n',sep = '')
	cat('BIC selected penalty parameter: ',x$lambda,'.\n')
	cat('Residual variance: ',var(residuals(x)),'.\n\n')
	
	cat('Variable selection:\nTotal number of non-zero parameters: ',sum(x$rawcoef!=0),'.\n')
	cat('Total number of regressors: ',length(x$rawcoef),'\n')
	cat('Non zero parameters per variable: \n')
	cat('\t\t',colnames(coefficients(x)),'\n')
	cat('\t\t',colSums(matrix(x$rawcoef,ncol=x$r,byrow=TRUE)!=0))
	cat('\n\n')	
	
	invisible(x)
}