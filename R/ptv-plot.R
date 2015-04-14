#' A function plotting the estimated parameter path in a parsimonious class object.
#' 
#' @param x An object of class parsimonious.
#'  @param ... unused.
#'  
#' @examples
#' \dontrun{
#' data(economics, package="ggplot2")
#' psr <- tail(economics$psavert,-1) # personal saving rate.
#' lur <- head(log(economics$unemploy/economics$pop),-1) # lagged log unemployment rate.
#' ptvmod <- ptvfit(y=psr,x=lur)
#' plot(ptvmod)
#' }
#'  
#' @rdname plot
#' @method plot parsimonious
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export 
plot.parsimonious <- function(x,...){	
	# checking class
	if(attr(x, "class")!='parsimonious') stop('x must be an object of class parsimonious.')
	# Initialization	
	r  <- x$r
	vn <- colnames(coefficients(x))
	if(is.null(vn))vn <- paste0('Variable ',1:r)	
	
	mcb <- melt(coefficients(x))
	colnames(mcb)[1:2] <- c('Time','Variable')
	mcb$Estimator <- x$estimator
	
	# checking if there's a post Lasso OLS estimator 
	if(!is.null(x$postcoef)){
		post <- matrix(x$postcoef,ncol=x$r)
		colnames(post) <- vn
		mpst <- melt(post)
		mpst$Estimator <- 'post Lasso OLS'
		colnames(mpst) <- colnames(mcb)
		mcb <- rbind(mcb,mpst)
	}
	
	# plotting
	gpar <- ggplot(mcb,aes(x=Time,y=value,colour=Estimator,shape=Estimator)) + geom_line() + geom_point() + facet_wrap(~Variable) 
	gpar <- gpar + theme_bw() + theme(legend.box = "horizontal",legend.position="bottom") 
	
	# print and return
	print(gpar)
	invisible(gpar)
}
