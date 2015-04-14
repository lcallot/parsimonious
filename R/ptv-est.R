#' A function to estimate a model with parsimoniously time varying parameters.
#' 
#' @name ptvfit
#' 
#' @param formula an object of class formula. An intercept should always be included in the model, do not explicitely remove the intercept from the formula. 
#' @param data an optional data frame.
#' @param weights an optional vector of weights (of length r(T-r+1)) used to compute an adaptive lasso. If a weight is equal to zero the variable is automatically excluded from the model.  
#' @param alpha the mixing parameter of elastic net. Default 1 (Lasso), set to 0 for ridge regression.
#' @param lambda user-defined lambda sequence.
#' @param peninit boolean: should the initial value of the parameters be penalized or only the increments in the random walk?
#' @param postest boolean: should the post Lasso OLS estimation be performed? 
#' @param ptvcste boolean: should the intercept be parsimoniously time varying?
#' @param pmax: maximum number of non-zero parameters to estimate.
#' @param ... unused.
#' 
#' @return An object of the parsimonious class.
#' 
#' @examples
#' \dontrun{
#' data(economics, package="ggplot2")
#' psr <- tail(economics$psavert,-1) # personal saving rate.
#' lur <- head(log(economics$unemploy/economics$pop),-1) # lagged log unemployment rate.
#' pce <- head(log(economics$pce/economics$pop),-1) # lagged log consumption expenditures per capita.
#' ptvmod <- ptvfit(psr~lur+pce)
#' }
#' 
#' @importFrom glmnet glmnet
#' @export
ptvfit <- function(formula, data, wada=NULL, alpha=1, lambda=NULL, peninit=TRUE, postest=TRUE, ptvcste = FALSE, pmax = NULL,...){
	
	# checking optional param
	if(!is.numeric(alpha))stop('alpha must be between 0 and 1')
	# if (!is.null(stable) && !is.numeric(stable)) stop('stable must be a numeric vector.')
	
	# From lm with some cleaning
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "na.action"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- quote(stats::model.frame)
	mf <- eval(mf, parent.frame())	
	mt <- attr(mf, "terms")
	
	intercept <- as.logical(attr(mt, 'intercept'))
	# getting the model components
	y <- model.response(mf, "numeric")
	x <- model.matrix(mt, mf, contrasts)
	if (intercept) x <- x[,-1]
	
	# A few useful variables
	nT <- length(y)
	r <- length(x)/nT
	x <- matrix(x,nrow=nT,ncol=r)
	
	#names
	if(!is.null(colnames(x)))	vn <- colnames(x) # works only if the rhs variables have individual names
	else {
		if(r>1)vn <- paste('x',1:ncol(x),sep='.') # works if they don't
		if(r==1)vn <- 'x'
	}
	colnames(x) <- vn
	if(ptvcste) vn <- c('cste',vn)

	# Checking the weights for the adaptive Lasso
	if (!is.null(wada) && !is.numeric(wada)) 
	{
		stop("'wada' must be a numeric vector.")
	}
		
	#Time varying intercept?
	if(ptvcste) {
		#intercept <- FALSE
		r <- r+1
		x <- cbind(1,x)
	}
	# Creating the W matrix
	w <- matrix(1,ncol=nT,nrow=nT)
	w[upper.tri(w)] <- 0
	if(r>1) w <- w[,-c((nT-r+2):nT)] # 1 too many column was removed.
	W <- kronecker(w,diag(r))
	
	# row-diagonalization of X
	dX <- matrix(0,nrow=nT,ncol=nT*r)
	for(n in 1:nT) dX[n,((n-1)*r+1):(n*r)] <- x[n,] 
	XDW <- dX%*%W

	if(ptvcste) XDW[,1]<-0
		
	if(is.null(pmax))pmax <- ncol(XDW)
	
	if(!is.null(wada))	{
		if(!((length(wada)==nT*r)|(length(wada)==(nT-r+1)*r)))stop('wada must have length r*nT or r*(nT-r+1).')
		wvec <- matrix(wada,ncol=1,byrow=TRUE)		
		
		if((r>1)&(length(wvec)==r*nT)) wvec<- head(wvec,-r*(r-1))
		
		# exclude infinite weights
		exclude <- which(wvec==Inf)		
		# removind the Infs
		wvec[which(wvec==Inf)] <- 0
		if (!peninit) wvec[1:r] <- 0
		fb <- glmnet(XDW,y,family='gaussian',intercept=intercept,standardize=FALSE,lambda=lambda,alpha=alpha,penalty.factor=wvec,exclude=exclude,pmax=pmax)
	}
	else 
	{
		penfac <- rep(1, ncol(XDW))
		if (!peninit) penfac[1:r] <- 0
		fb <- glmnet(XDW,y,family='gaussian',intercept=intercept,standardize=FALSE,lambda=lambda,alpha=alpha,penalty.factor=penfac,pmax=pmax)
	}

	
	# Getting fit and residuals
	pr.fb <- predict(fb,XDW)	
	res	<-matrix(rep(y,ncol(pr.fb)),ncol=ncol(pr.fb))-pr.fb
	
	# Computing BIC
	RSS	<-colSums(res^2)
	if(alpha>0)bic <- log(RSS/nT) + log(nT)/nT * fb$df# * log(log(ncol(XDW)))
	# BIC with ridge df
	if(alpha==0){rdf <- .ridge.df(XDW,fb$lambda); bic <- log(RSS/nT) + log(nT)*rdf/nT}
	
	# Selecting BIC minimizing beta
	bv <- fb$beta[,which.min(bic)]
	# adding the last obs and possibly setting up the inital value of the ptv intercept. 
	bvec <- c(bv,rep(0,r*(r-1))) 
	if(ptvcste) bvec[1] <- fb$a0[which.min(bic)]
	beta <- matrix(bvec,ncol=r,byrow=TRUE)
	cb <- apply(beta,2,cumsum)
	cb <- matrix(cb,ncol=r)
	colnames(cb) <- vn
	
	# Creating the ptv object 
	ptv <- list('y'=y,'x'=x,'lambda'=fb$lambda[which.min(bic)],'intercept'=fb$a0[which.min(bic)],'coefficients'=cb,'rawcoef'=bvec,'residuals'=res[,which.min(bic)],'wada'=wada,'alpha'=alpha,'r'=r,'nT'=nT,'estimator'=ifelse(is.null(wada),'Lasso','Adaptive Lasso'),'call'=cl,'post'=FALSE, 'glmnetobj'=fb)
	
	# Post Lasso OLS + stderr
	if(sum(bv!=0)>0 & postest){
		ptv$post <- TRUE
		post <- lm.fit(y=y,x=cbind(1,(XDW)[,bv!=0]))
		cpst <- matrix(0,nrow=length(bv),ncol=1)
		cpst[bv!=0,] <- coefficients(post)[-1]
		cpst <- rbind(cpst,matrix(0,nrow=r*(r-1),ncol=1))
		cpst <- matrix(cpst,ncol=r,byrow=TRUE)
		cpst <- matrix(apply(cpst,2,cumsum),ncol=r)
		cstpst <- coefficients(post)[1]
		sepost <- try(sum(residuals(post)^2) *diag(solve(t((XDW)[,bv!=0])%*%((XDW)[,bv!=0])))/(nT-sum(bv!=0)))
		ptv<- c(ptv,list('postcoef'=cpst,'postse'=sepost,'postintercept'=cstpst,'postres'=residuals(post)))
	}
	
	class(ptv) <- 'parsimonious'
	return(ptv)
}

